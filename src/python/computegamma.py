import prody as pd
import pylab as pl
import igraph as ig

class ProteinProjection:
    def __init__(self, proteinName: str, normal):
        self.normal = normal # normal is a 3d vector, the normal to the plane
                             # to which we intend to
                             # project the 3d protein structure

        self.p = parsePDB(proteinName)
        backbone = p.select('protein') # removes unnecessary water molecules,
                                       # leaving only the backbone residues

        resNums = backbone.getResnums()

        self.resCoords = {} # dictionary of residue coordinates as averages
                            # of atom coordinates in the respective residue

        self.projectedResidues = {} # a dictionary containing the residue number
                                    # together with the 2d coordinates of the
                                    # residue projected onto the given plane

                                    # for now, we just project onto the xy plane
                                    # for simplicity
        for i in range(resNums[0],resNums[-1]+1):
            projectedResidues[i] =

    def

    def exportTangle():
        pass
