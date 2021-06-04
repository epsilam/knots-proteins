import proteinprojection as pp
from scipy.linalg import norm
from collections import OrderedDict

class ProteinGraph:
    def __init__(self, proteinName: str, angle1: float, angle2: float):

        PP = ProteinProjection(proteinName,angle1,angle2)
        PP.findBackbonePrimaryCrossings()
        PP.findProteinStructureBonds()
        PP.findStructureCrossings()
        PP.findStructureWithBackboneCrossings()

        self.Bvertices = {}   # vertices corresp. to backbone residues
        self.Hvertices = {}   # H-vertices
        self.BBcrossings = {} # crossings between backbone bonds
        self.BScrossings = {} # crossings between backbone and structure bonds
        self.SScrossings = {} # crossings between structure bonds

                              # the above dictionaries have key,value pairs of
                              # the format
                              # ("Bv",res) : (v1,v2)
                              # ("Hv",res,n) : (v1,v2,v3,v4)
                              # ("BB",(resi1,resi2,resj1,resj2)) : (v1,v2,v3,v4)
                              # ("BS",(resi1,resi2,resj1,resj2)) : (v1,v2,v3,v4)
                              # ("SS",(resi1,resi2,resj1,resj2)) : (v1,v2,v3,v4)

                               # note: BB denotes a crossing where the backbone
                               # crosses itself, BS denotes a crossing between
                               # the backbone and a structure bond,
                               # and SS denotes a crossing between two structure
                               # bonds.
                               # vertices from backbone, which are connected
                               # to precisely two other vertices, crossings,
                               # or H-vertices, or two of a mix of the above.
                               # The initial and terminal residues may have
                               # None as one of their connections.
                               # A crossing in a value in self.vertices is
                               # identified via its index in self.crossings,
                               # and is denoted as a tuple: (i,'X'), and
                               # H-vertices are denoted (i,'H')

        self.vertices[PP.resNums[0]] = (None,PP.resNums[1]) # initial vertex
        for i in range(1,len(PP.resNums)-1):
            self.vertices[PP.resNums[i]] = (PP.resNums[i-1],PP.resNums[i+1])
        self.vertices[PP.resNums[-1]] = (PP.resNums[-2],None)

    def evalStructureBond(self, pair):
        i = pair[0]
        j = pair[1]


                                    # find all crossings intersecting a
                                    # backbone (peptide) bond between
                                    # resi1 and resi2, and return a sorted
                                    # dictionary containing distances to
                                    # the first residue resi1.
    def _findBBandBScrossingsOnBackboneBond(resi1,resi2):
        cr = {} # dictionary whose keys are crossing labels of crossings
                # intersecting the peptide bond, and whose values are
                # distances of crossings to the first res in the crossing.
        for key,value in PP.backboneCrossings.items():
            if resi1 == key[0] or resi2 == key[2]:
                d = norm(PP.projectedRescoords[resi1] - value)
                cr[("BB",key)] = d
        for key,value in PP.structureWithBackboneCrossings.items():
            if (resi1,resi2) == key[0:2] or (resi1,resi2) == key[2:4]:
                d = norm(PP.structureWithBackboneCrossings[resi1] - value)
                cr[("BS",key)] = d
        cr = OrderedDict(sorted(cr.items(), key = lambda x:x[1])) # sort by dist
        return cr

                                    # traveling along the backbone, add vertices
                                    # to dicts for all types of crossings along
                                    # each eptide bond & connect appropriately.
    def setCrossingVertices(self):
        resi1 = PP.resNums[0] # dealing with the startpoint residue of the chain
        resi2 = PP.resNums[1]
        cr = self._findBBandBScrossingsOnBackboneBond(resi1,resi2)
        if not cr: # if there are no crossings on this bond
            self.Bvertices[("Bv",resi1)] = (None, ("Bv",resi2))
        else:
            crlist = list(cr) # allow us to index crossings

            self.Bvertices[("Bv",resi1)] = (None,crlist[0]) # deal w/ 1st cross
            X = crlist[0]
            if (resi1,resi2) == X[0][0:2]: isOver = True
            if (resi1,resi2) == X[0][2:4]: isOver = False
            if X[0] == "BB":
                self.BBcrossings[X] = (None, None, None, None)
                if isOver:
                    self.BBcrossings[X][0] = ("Bv",resi1)
                    self.BBcrossings[X][1] = crlist[i+1]
                else:
                    self.BBcrossings[X][2] = ("Bv",resi1)
                    self.BBcrossings[X][3] = crlist[i+1]

            if X[0] == "BS":
                self.BScrossings[X] = (None, None, None, None)
                if isOver:
                    self.BScrossings[X][0] = ("Bv",resi1)
                    self.BScrossings[X][1] = crlist[i+1]
                else:
                    self.BScrossings[X][2] = ("Bv",resi1)
                    self.BScrossings[X][3] = crlist[i+1]


            for i in range(1,len(crlist)-1):
                X = crlist[i]
                                    # check if the current backbone bond goes
                                    # over or under this crossing.
                if (resi1,resi2) == X[0][0:2]: isOver = True
                if (resi1,resi2) == X[0][2:4]: isOver = False
                if X[0] == "BB":
                    self.BBcrossings[X] = (None, None, None, None)
                    if isOver:
                        self.BBcrossings[X][0] = crlist[i-1]
                        self.BBcrossings[X][1] = crlist[i+1]
                    else:
                        self.BBcrossings[X][2] = crlist[i-1]
                        self.BBcrossings[X][3] = crlist[i+1]

                if X[0] == "BS":
                    self.BScrossings[X] = (None, None, None, None)
                    if isOver:
                        self.BScrossings[X][0] = crlist[i-1]
                        self.BScrossings[X][1] = crlist[i+1]
                    else:
                        self.BScrossings[X][2] = crlist[i-1]
                        self.BScrossings[X][3] = crlist[i+1]



        for i in range(1,len(PP.resNums)-1):
            resi1 = PP.resNums[i]
            resi2 = PP.resNums[i+1]
            cr = self.findBBandBScrossingsOnBackboneBond(resi1,resi2)





    def exportTangle(self): # exports crossings and H-vertices into a format
        pass                # readable by Mathematica
