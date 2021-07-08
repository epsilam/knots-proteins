from graph import *
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
