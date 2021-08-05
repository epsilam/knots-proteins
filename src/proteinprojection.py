import prody as pd
import pylab as pl
from matplotlib import collections  as mc
import numpy as np
from scipy.linalg import expm

class ProteinProjection:
                                         # angle1 and angle2 determine the 2d
                                         # plane onto which the protein is
                                         # projected. they are along the z-axis
                                         # and y-axis, respectively.
    def __init__(self, proteinName, planeNormal):
        self.planeNormal = planeNormal / np.linalg.norm(planeNormal)

        print("Attempting to find or download PDB file...")
        self.p = pd.parsePDB(proteinName)
        print("Parsing PDB file...")
        self.backbone = self.p.select('protein') # removes irrelevant structures
                                            # (such as water molecules,
                                            # carbohydrates, and ligands) and
                                            # leaves only the polymers


        self.resNums = np.unique(self.backbone.getResnums())

        self.resCoords = {} # dictionary of residue coordinates as averages
                            # of atom coordinates in the respective residue

        self.projectedResCoords = {} # a dictionary containing the residue
                                     # number together with the 2d coordinates
                                     # of the residue projected onto the given
                                     # plane.

        self.backboneCrossings = {}      # a dictionary containing crossings of
                                         # bonds on the backbone. the key
                                         # (i,i+1,j,j+1) is included if the bond
                                         # between i and i+1 crosses the bond
                                         # between j and j+1. the value is a
                                         # pair (position, isPositive)
                                         # where position is the (x,y) position
                                         # of the crossing, and isPositive
                                         # is a boolean equal to True if the
                                         # crossing is positive and False if
                                         # negative.

        self.proteinStructureBonds = set() # a set of pairs of residues which
                                           # are bonded via a non-peptide bond.

        self.structureWithBackboneCrossings = {} # crossings between a bond on
                                         # backbone (peptide) and a secondary/
                                         # tertiary structure bond.
                                         # the element (i1,i2,j1,j2) is included
                                         # if the bond between i1 and i2 crosses
                                         # the bond between j1 and j2

                                         # crossings between two secondary/
                                         # tertiary structure bonds. same format
                                         # as above.
        self.structureCrossings = {}

                          # given two pairs (resi1,resi2) and (resj1,resj2)
                          # of res nums, determine if the bond between resi1
                          # and resi2 crosses the bond between resj1 and resj2
                          # in a given projection using some
                          # simple linear algebra.
        self.projectResCoords()
        self.findBackbonePrimaryCrossings()
        self.findProteinStructureBonds()
        self.findStructureCrossings()
        self.findStructureWithBackboneCrossings()

    def projectResCoords(self):
        """Project all residue coordinates onto the plane with normal specified
        at init. Produces 2D coordinates by rotating the plane to the x,y-plane
        and taking the first two coordinates of each point."""
        t = np.transpose

        n = self.planeNormal
        z = np.array([0,0,1])           # z-axis
        r = t(np.matrix(np.cross(n,z))) # rotation axis between n and z
        a = np.arccos(np.dot(n,z))      # rotation angle between n and z

        n = t(np.matrix(n))
        z = t(np.matrix(z))
        cpr = z * t(n) - n * t(z)       # cross product matrix of n and z
        M = np.cos(a) * np.eye(3) + np.sin(a) * cpr + (1-np.cos(a)) * r * t(r)
        for i in self.resNums:
            res = self.backbone['A',i]
            self.resCoords[i] = np.average(res.getCoords(),axis=0)
            rc = self.resCoords[i]

            # 3D coordinates of residue projected onto plane
            proj = rc - n.A1 * np.dot(rc, n.A1) / np.dot(n.A1, n.A1)

            # coordinates after rotating the plane to be perpendicular to z-axis
            rotCoords = M * t(np.matrix(proj))
            self.projectedResCoords[i] = rotCoords.A1[0:2]


    def _checkBondCrossing(self, resi1, resi2, resj1, resj2, crossingDict):
        Ri1 = self.projectedResCoords[resi1]
        Ri2 = self.projectedResCoords[resi2]
        Rj1 = self.projectedResCoords[resj1]
        Rj2 = self.projectedResCoords[resj2]
        E = Rj1 - Ri1
        R = Ri2 - Ri1
        S = Rj2 - Rj1
        crossRS = np.cross(R,S)
        t = np.cross(E,S) / crossRS
        u = np.cross(E,R) / crossRS
        n = self.planeNormal
        if crossRS != 0:
            if 0 <= t <= 1 and 0 <= u <= 1:
                Ri1_3D = self.resCoords[resi1]
                Ri2_3D = self.resCoords[resi2]
                Rj1_3D = self.resCoords[resj1]
                Rj2_3D = self.resCoords[resj2]
                x1 = Ri1_3D + t * (Ri2_3D - Ri1_3D) # 3D location of projected
                x2 = Rj1_3D + u * (Rj2_3D - Rj1_3D) # crossing along each strand
                # project x1 and x2 onto the plane normal. whichever projection
                # has the greater norm will be the over strand
                x1p = np.dot(x1, n) / np.dot(n, n)
                x2p = np.dot(x2, n) / np.dot(n, n)
                if x1p >= x2p:
                    isPositive = (crossRS > 0)
                    resPairs = (resi1, resi2, resj1, resj2)
                else:
                    isPositive = not (crossRS > 0)
                    resPairs = (resj1, resj2, resi1, resi2)
                crossingDict[resPairs] = (Ri1 + t * R, isPositive)
        else:
            raise Exception("Bonds appear tangent in projection.\
            Use different projection.")

                          # for each pair of two successive residues (i.e., two
                          # residues with a peptide bond), check if there is a
                          # crossing with another
    def findBackbonePrimaryCrossings(self):
        for idx in range(len(self.resNums)-1):
            resi1 = self.resNums[idx]
            resi2 = self.resNums[idx+1]
            for jdx in range(idx+2,len(self.resNums)-1):
                resj1 = self.resNums[jdx]
                resj2 = self.resNums[jdx+1]
                self._checkBondCrossing(resi1,resi2,resj1,resj2,
                                        self.backboneCrossings)

                          # find all bonds corresponding to the secondary and
                          # tertiary structure of the protein
                          # (excludes peptide bonds)
    def findProteinStructureBonds(self):
        maxRadius = 2 * 1.221
        self.p.inferBonds(max_bond = maxRadius)
        for bondCandidate in self.backbone.iterBonds():
            bondAtomPair = bondCandidate.getAtoms()
            a1 = bondAtomPair[0]
            a2 = bondAtomPair[1]
            r1 = a1.getResnum()
            r2 = a2.getResnum()
            if r1 != r2 and r1 != r2 + 1 and r1 != r2 - 1:
                if r2 < r1: # to ensure no duplicates, force r1 < r2
                    temp = r1; r1 = r2; r2 = temp
                self.proteinStructureBonds.add((r1,r2))


                            # find crossings between secondary/tertiary
                            # structure bonds.
    def findStructureCrossings(self):
        bonds = list(self.proteinStructureBonds)
        for i in range(len(bonds)):
            for j in range(i+1,len(bonds)):
                resi1 = bonds[i][0]
                resi2 = bonds[i][1]
                resj1 = bonds[j][0]
                resj2 = bonds[j][1]
                if resi1 in (resj1,resj2) or resi2 in (resj1,resj2):
                    continue
                self._checkBondCrossing(resi1,resi2,resj1,resj2,
                                        self.structureCrossings)

                            # find crossings between peptide bonds and
                            # secondary/tertiary structure bonds.
    def findStructureWithBackboneCrossings(self):
        bonds = list(self.proteinStructureBonds)
        for i in range(len(bonds)):
            for j in range(len(self.resNums)-1):
                resi1 = bonds[i][0]
                resi2 = bonds[i][1]
                resj1 = self.resNums[j]
                resj2 = self.resNums[j+1]
                if resi1 in (resj1,resj2) or resi2 in (resj1,resj2):
                    continue
                self._checkBondCrossing(resi1,resi2,resj1,resj2,
                                        self.structureWithBackboneCrossings)


    def plotCrossings(self):
        fig, ax = pl.subplots()

        rc = []                      # plot backbone
        for key,value in self.projectedResCoords.items():
            rc.append(tuple(value))
        lines = []
        for i in range(len(rc)-1):
            lines.append((rc[i],rc[i+1]))
        lc = mc.LineCollection(lines, linewidths=2, colors=[(0.1,0.7,0.1,1)])
        ax.add_collection(lc)

        lines = []                  # plot secondary/tertiary structure bonds
        for bondedResPair in self.proteinStructureBonds:
            rc1 = tuple(self.projectedResCoords[bondedResPair[0]])
            rc2 = tuple(self.projectedResCoords[bondedResPair[1]])
            lines.append((rc1,rc2))
        lc = mc.LineCollection(lines, linewidths=2, colors=[(0,0,0,0.1)])
        ax.add_collection(lc)

        crossingDict = {}
        crossingDict.update(self.backboneCrossings)
        crossingDict.update(self.structureCrossings)
        crossingDict.update(self.structureWithBackboneCrossings)
        crossCoords = []                # plot crossings
        for key,value in crossingDict.items():
            crossCoords.append(value[0])
        xCrossCoords, yCrossCoords = zip(*crossCoords)
        ax.scatter(xCrossCoords,yCrossCoords, c=[(1,0,0,1)], s=30, zorder=3)

        ax.autoscale()
        ax.margins(0.1)
        fig.show()


