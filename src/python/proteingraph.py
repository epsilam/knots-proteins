import proteinprojection as pp

class ProteinGraph:
    def __init__(self, proteinName: str, angle1: float, angle2: float):

        PP = ProteinProjection(proteinName,angle1,angle2)
        PP.findBackbonePrimaryCrossings()
        PP.findProteinStructureBonds()
        PP.findStructureCrossings()
        PP.findStructureWithBackboneCrossings()

        self.vertices = {}     # vertices from backbone, which are connected
                               # to precisely two other vertices, crossings,
                               # or H-vertices, or two of a mix of the above.
                               # The initial and terminal residues may have
                               # None as one of their connections.
                               # A crossing in a value in self.vertices is
                               # identified via its index in self.crossings,
                               # and is denoted as a tuple: (i,'X'), and
                               # H-vertices are denoted (i,'H')
        self.Hvertices = []
        self.crossings = []

        self.vertices[PP.resNums[0]] = (None,PP.resNums[1]) # initial vertex
        for i in range(1,len(PP.resNums)-1):
            self.vertices[PP.resNums[i]] = (PP.resNums[i-1],PP.resNums[i+1])
        self.vertices[PP.resNums[-1]] = (PP.resNums[-2],None)

    def evalStructureBond(self, pair):
        i = pair[0]
        j = pair[1]

    def evalBackboneCrossings(self):
        for key,value in PP.backboneCrossings:
            resi1 = key[0]
            resi2 = key[1]
            resj1 = key[2]
            resj2 = key[3]

            crossingIndex = len(self.crossings)
            self.crossings.append((resi1,resi2,resj1,resj2))
            self.vertices[resi1][1] = (crossingIndex,'X')
            self.vertices[resi2][0] = (crossingIndex,'X')
            self.vertices[resj1][1] = (crossingIndex,'X')
            self.vertices[resj2][0] = (crossingIndex,'X')


    def exportTangle(self): # exports crossings and H-vertices into a format
        pass                # readable by Mathematica
