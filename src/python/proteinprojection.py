import prody as pd
import pylab as pl
from matplotlib import collections  as mc
import igraph as ig
import numpy as np
from scipy.linalg import expm

#

class ProteinProjection:
                                         # angle1 and angle2 determine the 2d
                                         # plane onto which the protein is
                                         # projected. they are along the z-axis
                                         # and y-axis, respectively.
    def __init__(self, proteinName: str, angle1: float, angle2: float):

        self.p = pd.parsePDB(proteinName)
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
        for i in self.resNums:
            res = self.backbone['A',i]
            self.resCoords[i] = np.average(res.getCoords(),axis=0)
                                     # we rotate the residue coordinates by
                                     # the two angles and then project onto
                                     # the xy-plane.

            M1 = expm(np.cross(np.eye(3),  # rotation matrix along z-axis
                            np.array([0,0,1])*angle1))
            M2 = expm(np.cross(np.eye(3),  # rotation matrix along y-axis
                            np.array([0,1,0])*angle2))
            rot = np.dot(M1,np.dot(M2,self.resCoords[i]))
            self.projectedResCoords[i] = rot[0:2] # take the x,y coords

                          # given two pairs (resi1,resi2) and (resj1,resj2)
                          # of res nums, determine if the bond between resi1
                          # and resi2 crosses the bond between resj1 and resj2
                          # in a given projection using some
                          # simple linear algebra.


        self.backboneCrossings = {}      # a dictionary containing crossings of
                                         # bonds on the backbone. the key
                                         # (i,i+1,j,j+1) is included if the bond
                                         # between i and i+1 crosses the bond
                                         # between j and j+1.


        self.proteinStructureBonds = set() # a set of pairs of residues which
                                           # are bonded via a non-peptide bond.

                                         # crossings between a bond on the
                                         # backbone (peptide) and a secondary/
                                         # tertiary structure bond.
                                         # the element (i1,i2,j1,j2) is included
                                         # if the bond between i1 and i2 crosses
                                         # the bond between j1 and j2
        self.structureWithBackboneCrossings = {}

                                         # crossings between two secondary/
                                         # tertiary structure bonds. same format
                                         # as above.
        self.structureCrossings = {}

    def _checkBondCrossing(self, resi1, resi2, resj1, resj2, crossingDict):
        A = self.projectedResCoords[resi1]
        B = self.projectedResCoords[resi2]
        C = self.projectedResCoords[resj1]
        D = self.projectedResCoords[resj2]
        E = C - A
        R = B - A
        S = D - C
        crossRS = np.cross(R,S)
        t = np.cross(E,S) / crossRS
        u = np.cross(E,R) / crossRS
        if crossRS != 0:
            if 0 <= t <= 1 and 0 <= u <= 1:
                crossingDict[(resi1,resi2,resj1,resj2)] = A + t * R

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
        self.p.inferBonds(max_bond=2*1.221)
        for bond in self.backbone.iterBonds():
            bondAtomPair = bond.getAtoms()
            r1 = bondAtomPair[0].getResnum()
            r2 = bondAtomPair[1].getResnum()
            if r1 != r2 and r1 != r2 + 1 and r1 != r2 - 1:
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
        lc = mc.LineCollection(lines, linewidths=2, colors=[(0.2,1,0,1)])
        ax.add_collection(lc)

        lines = []                  # plot secondary/tertiary structure bonds
        for bondedResPair in self.proteinStructureBonds:
            rc1 = tuple(self.projectedResCoords[bondedResPair[0]])
            rc2 = tuple(self.projectedResCoords[bondedResPair[1]])
            lines.append((rc1,rc2))
        lc = mc.LineCollection(lines, linewidths=2, colors=[(0,0,0,0.3)])
        ax.add_collection(lc)

        crossingDict = {}
        crossingDict.update(self.backboneCrossings)
        crossingDict.update(self.structureCrossings)
        crossingDict.update(self.structureWithBackboneCrossings)
        crossCoords = []                # plot crossings
        for key,value in crossingDict.items():
            crossCoords.append(value)
        xCrossCoords, yCrossCoords = zip(*crossCoords)
        ax.scatter(xCrossCoords,yCrossCoords, c=[(1,0,0,1)], s=30, zorder=3)

        ax.autoscale()
        ax.margins(0.1)
        fig.show()


