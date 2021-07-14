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
                       # : (isOver, isOtherStrandHbond, dist)
        for key,value in self.ProProj.backboneCrossings.items():
            dist = norm(self.ProProj.projectedResCoords[res] - value[0])
            if (res,res2) == key[0:2]:
                crossings[key] = (True, False, dist)
            if (res,res2) == key[2:4]:
                crossings[key] = (False, False, dist)
        for key,value in self.ProProj.structureWithBackboneCrossings.items():
            dist = norm(self.ProProj.projectedResCoords[res] - value[0])
            if (res,res2) == key[0:2]:
                crossings[key] = (True, True, dist)
            if (res,res2) == key[2:4]:
                crossings[key] = (False, True, dist)
        crossings = OrderedDict(sorted(crossings.items(),key=lambda x:x[1][2]))

        vert2 = self.initialVert[res2]
        vert1 = self.Sverts[vert2][1]
        for resPairs in list(crossings)[:-1]:
            isOver = crossings[resPairs][0]
            isOtherStrandHbond = crossigns[resPairs][1]
            vert1 = self._addXvertsAlongBB(resPairs, vert1, None,
                                           isOver, isOtherStrandHbond)
        resPairs = crossings[-1] # last crossing
        isOver = crossings[resPairs][0]
        isOtherStrandHbond = crossigns[resPairs][1]
        vert1 = self._addXvertsAlongBB(resPairs, vert1, vert2,
                                           isOver, isOtherStrandHbond)
        self.Sverts[vert2][1] = vert1

    def resolveXvertsHbond(self, res1, res2):
        """Resolve all crossings on the H-bond between residues res1 and res2.
        """
        crossings = {} # format:
                       # (resi1,resi2,resj1,resj2)
                       # : (isOver, isOtherStrandHbond, dist)
        for key,value in self.ProProj.structureWithBackboneCrossings.items():
            dist = norm(self.ProProj.projectedResCoords[res] - value)
            if (res1,res2) == key[0:2]:
                crossings[key] = (True, False, dist)
            if (res1,res2) == key[2:4]:
                crossings[key] = (False, False, dist)
        for key,value in self.ProProj.structureCrossings.items():
            dist = norm(self.ProProj.projectedResCoords[res] - value)
            if (res1,res2) == key[0:2]:
                crossings[key] = (True, True, dist)
            if (res1,res2) == key[2:4]:
                crossings[key] = (False, True, dist)
        crossings = OrderedDict(sorted(crossings.items(),key=lambda x:x[1][2]))

        # now, we need to find the relevant H-vertex and the corresponding
        # two S-vertices that match the H-bond between the two given residues.
