from graph import *
from scipy.linalg import norm
from collections import OrderedDict
import proteinprojection as proproj

class ProteinGraph(Graph):
    def __init__(self, proteinName: str, angle1: float, angle2: float):
        Graph.__init__(self)
        self.ProProj = proproj.ProteinProjection(proteinName, angle1, angle2)
        self.initialVerts = {} # dictionary to keep track of S-vertices
                               # representing actual residues.

    def initResVerts(self):
        """Create initial S-vertices corresponding to each residue
        in the chain."""
        resNums = self.ProProj.resNums
        v = self._addSvert(resNums[0], None, None) # initial vertex
        self.initialVerts[resNums[0]] = v
        for i in range(1,len(resNums)):
            v = self._addSvert(resNums[i], v, None)
            self.initialVerts[resNums[i]] = v
        self._addSvert(resNums[-1], v, None) # terminal vertex

    def resolveHverts(self, res):
        """At the residue res, create all necessary H-vertices."""
        HBondedResidues = set()
        for (res1,res2) in self.ProProj.proteinStructureBonds:
            if res == res1:
                HBondedResidues.add(res2)
            elif res == res2:
                HBondedResidues.add(res1)
        HBondedResidues = list(HBondedResidues)
        vert = self.initialVerts[res]
        vert1 = vert
        for res2 in HBondedResidues:
            vert2 = self.initialVerts[res2]
            if vert < vert2: # prevent H-vertices being added twice
                vert1 = self._addHvert(vert1, vert2)

    def resolveXvertsBB(self, res):
        """Resolve all crossings along the bond between the residue res and the
        next residue on the backbone."""
        res2 = res + 1 # next residue on backbone

        # First, we create an ordered dictionary of crossings, ordered by
        # distance to the residue res.
        crossings = {} # format:
                       # (resi1,resi2,resj1,resj2)
                       # : (isOver, isOtherStrandHbond, dist, isPositive)
        for key,value in self.ProProj.backboneCrossings.items():
            isPositive = value[1]
            dist = norm(self.ProProj.projectedResCoords[res] - value[0])
            if (res,res2) == key[0:2]:
                crossings[key] = (True, False, dist, isPositive)
            if (res,res2) == key[2:4]:
                crossings[key] = (False, False, dist, isPositive)
        for key,value in self.ProProj.structureWithBackboneCrossings.items():
            isPositive = value[1]
            dist = norm(self.ProProj.projectedResCoords[res] - value[0])
            if (res,res2) == key[0:2]:
                crossings[key] = (True, True, dist, isPositive)
            if (res,res2) == key[2:4]:
                crossings[key] = (False, True, dist, isPositive)
        crossings = OrderedDict(sorted(crossings.items(),key=lambda x:x[1][2]))

        vert2 = self.initialVert[res2]
        vert1 = self.Sverts[vert2][1]
        for resPairs in list(crossings)[:-1]:
            isOver = crossings[resPairs][0]
            isOtherStrandHbond = crossings[resPairs][1]
            isPositive = crossings[resPairs][3]
            vert1 = self._addXvertsAlongBB(resPairs,
                                           vert1, None,
                                           isOver,
                                           isOtherStrandHbond,
                                           isPositive)
        resPairs = crossings[-1] # last crossing
        isOver = crossings[resPairs][0]
        isOtherStrandHbond = crossigns[resPairs][1]
        isPositive = crossings[resPairs][3]
        vert1 = self._addXvertsAlongBB(resPairs,
                                       vert1, vert2,
                                       isOver,
                                       isOtherStrandHbond,
                                       isPositive)
        self.Sverts[vert2][1] = vert1

    def resolveXvertsHbond(self, res1, res2):
        """Resolve all crossings on the H-bond between residues res1 and res2.
        """
        crossings = {} # format:
                       # (resi1,resi2,resj1,resj2)
                       # : (isOver, isOtherStrandHbond, dist, isPositive)
        for key,value in self.ProProj.structureWithBackboneCrossings.items():
            isPositive = value[1]
            dist = norm(self.ProProj.projectedResCoords[res] - value[0])
            if (res1,res2) == key[0:2]:
                crossings[key] = (True, False, dist, isPositive)
            if (res1,res2) == key[2:4]:
                crossings[key] = (False, False, dist, isPositive)
        for key,value in self.ProProj.structureCrossings.items():
            isPositive = value[1]
            dist = norm(self.ProProj.projectedResCoords[res] - value[0])
            if (res1,res2) == key[0:2]:
                crossings[key] = (True, True, dist, isPositive)
            if (res1,res2) == key[2:4]:
                crossings[key] = (False, True, dist, isPositive)
        crossings = OrderedDict(sorted(crossings.items(),key=lambda x:x[1][2]))

        # find the relevant H-vertex and the corresponding
        # two S-vertices that match the H-bond between the two given residues.
        vert = self.initialVerts[res1]
        nxt = self.Sverts[vert][2]
        try:
            while True:
                svert = self.Hverts[nxt][3]
                if self.Sverts[svert][0] == res2: # check that the H-vertex is
                                                  # along a strand which goes
                                                  # through res2
                    break
                nxt = self.Hverts[nxt][2]
        except KeyError:
            raise Exception("Error in resolving X-vertices along H-bond.\
             No H-bond found between these two residues")
        Hvert = nxt
        # we use the convention that the strand travelling away from the
        # H-vertex is labelled 1, and the one going toward it is labelled 2.
        Svert1 = self.Hverts[Hvert][4]
        Svert2 = self.Hverts[Hvert][3]


        for resPairs in list(crossings)[:-1]:
            isOver = crossings[resPairs][0]
            isOtherStrandHbond = crossings[resPairs][1]
            isPositive = crossings[resPairs][3]
            (v1,v2) = self._addXvertsAlongHbond(resPairs,
                                                v1,None,
                                                None,v2,
                                                isOver,
                                                isOtherStrandHbond,
                                                isPositive)
        resPairs = crossings[-1]
        isOver = crossings[resPairs][0]
        isOtherStrandHbond = crossigns[resPairs][1]
        isPositive = crossings[resPairs][3]
        (v1,v2) = self._addXvertsAlongHbond(resPairs,
                                            v1,Svert1,
                                            Svert2,v2,
                                            isOver,
                                            isOtherStrandHbond,
                                            isPositive)
