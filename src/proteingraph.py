from src.graph import *
from scipy.linalg import norm
from collections import OrderedDict
import src.proteinprojection as proproj

class ProteinGraph(Graph):
    def __init__(self, proteinName, planeNormal):
        Graph.__init__(self)
        self.ProProj = proproj.ProteinProjection(proteinName, planeNormal)
        self.initialVerts = {} # dictionary to keep track of S-vertices
                               # representing actual residues.
        print("Creating graph structure...")
        self.resolve()
        print("Number of atoms:          " + str(self.ProProj.p.numAtoms()))
        print("Number of residues:       " + str(self.ProProj.p.numResidues()))
        print("Number of H-contacts:     " + str(len(self.Hverts)))
        print("Number of real crossings: " \
              + str(  len(self.ProProj.backboneCrossings) \
                    + len(self.ProProj.structureWithBackboneCrossings) \
                    + len(self.ProProj.structureCrossings)))
        print("Number of vertices:       " + str(  len(self.Hverts) \
                                                 + len(self.Xverts)))

    def resolve(self):
        """Create the entire graph corresponding to the given protein,
        resolving all H-bonds and crossings."""
        self.initResVerts()

        resNums = self.ProProj.resNums
        for res in resNums[:-1]:
            self.resolveHverts(res)

        for res in resNums[:-1]:
            self.resolveXvertsBB(res)

        for (res1,res2) in self.ProProj.proteinStructureBonds:
            self.resolveXvertsHbond(res1,res2)

        for vert in list(self.Sverts)[1:]: # we prune all S-vertices except the
            self._pruneSvert(vert)         # initial vertex.

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

        vert2 = self.initialVerts[res2]
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
        if len(crossings) != 0: # deal with last crossing, if a crossing exists
            resPairs = list(crossings)[-1]
            isOver = crossings[resPairs][0]
            isOtherStrandHbond = crossings[resPairs][1]
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
            dist = norm(self.ProProj.projectedResCoords[res1] - value[0])
            if (res1,res2) == key[0:2]:
                crossings[key] = (True, False, dist, isPositive)
            if (res1,res2) == key[2:4]:
                crossings[key] = (False, False, dist, isPositive)
        for key,value in self.ProProj.structureCrossings.items():
            isPositive = value[1]
            dist = norm(self.ProProj.projectedResCoords[res1] - value[0])
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
                svert = nxt + 1 # the corresponding S-verts are created right
                                # after the H-verex, so their ID's are created
                                # one after another.
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
        Svert2 = Hvert + 1
        Svert1 = Hvert + 2
        v1 = Hvert
        v2 = Hvert
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
        if len(crossings) != 0: # deal with last crossing, if a crossing exists
            resPairs = list(crossings)[-1] # last crossing
            isOver = crossings[resPairs][0]
            isOtherStrandHbond = crossings[resPairs][1]
            isPositive = crossings[resPairs][3]
            (v1,v2) = self._addXvertsAlongHbond(resPairs,
                                                v1,Svert1,
                                                Svert2,v2,
                                                isOver,
                                                isOtherStrandHbond,
                                                isPositive)
        self.Sverts[Svert1][1] = v1
        self.Sverts[Svert2][2] = v2
