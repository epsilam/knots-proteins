import proteinprojection as pp
from scipy.linalg import norm
import numpy as np


class ProteinGraph():
    def __init__(self, proteinName: str, angle1: float, angle2: float):

        #
        self.PP = pp.ProteinProjection(proteinName,angle1,angle2)
        self.PP.findBackbonePrimaryCrossings()
        self.PP.findProteinStructureBonds()
        self.PP.findStructureCrossings()
        self.PP.findStructureWithBackboneCrossings()

        # dicts for (S)imple vertiecs, (H)-vertices, and (C)rossing vertices.
        # formats of the dictionary elements are given below.
        # n denotes the unique index for that vertex, res is the residue number
        # associated with that vertex, and prev and next are the unique indices
        # of the vertices connected to the vertex.
        self.Sverts = {}  # n : [res, prev, next]
        self.Hverts = {}  # n : [res, prev1, next1, prev2, next2]
        self.Cverts = {}  # n : [res, prevOver, nextOver, prevUnder, nextUnder]

        # add simple vertices. the first vertex has no previous vertex, and
        # the last vertex has no next vertex.
        numVerts = len(self.PP.resNums)
        self.Sverts[0] = [self.PP.resNums[0], None, 1]
        for i in range(1,numVerts-1):
            self.Sverts[i] = [self.PP.resNums[i], i-1, i+1]
        lastVert = numVerts - 1
        self.Sverts[lastVert] = [self.PP.resNums[lastVert], lastVert - 1, None]

        # add H-vertices and connect appropriately.
        for i in range(numVerts):
            self._resolveHVertices(i)

    def _findHBondedResidues(self, resnum: int):
        """Return the set of residues bonded to resnum via an H-bond."""
        residues = set()
        for (res1,res2) in self.PP.proteinStructureBonds:
            if resnum == res1:
                residues.add(res2)
            elif resnum == res2:
                residues.add(res1)
        return list(residues)

    def _resolveHVertices(self, n):
        """Create necessary H-vertices at each given simple vertex n based on
        H-bonds."""

        nResNum = self.Sverts[n][0] # residue number corresp. to vertex n
        HbondedResidues = self._findHBondedResidues(nResNum)
        if len(HbondedResidues) == 0:
            return
        numVerts = len(self.Sverts) + len(self.Hverts)
        firstHvert = numVerts
        lastHvert = numVerts + len(HbondedResidues) - 1
        if len(HbondedResidues) == 1:
            self.Hverts[firstHvert] = [nResNum,
                                       n, self.Sverts[n][2],
                                       None, None]
        else: # here, we are guaranteed to have at least two H-vertices
            self.Hverts[firstHvert] = [nResNum,
                                       n, firstHvert + 1,
                                       None, None]
            self.Hverts[lastHvert]  = [nResNum,
                                       lastHvert - 1, self.Sverts[n][2],
                                       None, None]
            # for all H-vertices except the first and last
            for i in range(firstHvert+1, lastHvert):
                self.Hverts[i] = [nResNum, i-1, i+1, None, None]

        # properly connect the vertex n and the next vertex on the protein chain
        nNext = self.Sverts[n][2]
        if nNext != None: # then n is not the last vertex in the chain
            self.Sverts[nNext][1] = lastHvert
        self.Sverts[n][2] = firstHvert # next vertex of n is firstHvert


        numVerts = len(self.Sverts) + len(self.Hverts)
        # add second strands for each of the H-vertices.
        # this requires adding 2*len(HbondedResidues) many simple vertices
        for i in range(len(HbondedResidues)):
            nBonded = HbondedResidues[i] # the residue which n is H-bonded to
            nBondedIdx = list(self.PP.resNums).index(nBonded)
            HvertIdx = firstHvert + i

            # create two additional vertices next to the bonded residue and
            # add edges accordingly
            v1 = numVerts + 2 * i
            v2 = numVerts + 2 * i + 1
            self.Sverts[v1] = [nBonded,
                               nBondedIdx,
                               HvertIdx]
            self.Sverts[v2] = [nBonded,
                               HvertIdx,
                               self.Sverts[nBondedIdx][2]]
            self.Hverts[HvertIdx][3] = v1
            self.Hverts[HvertIdx][4] = v2
