from src.proteingraph import *

class ProteinGamma(ProteinGraph):
    def __init__(self, proteinName: str, angle1: float, angle2: float):
        ProteinGraph.__init__(self, proteinName, angle1, angle2)

        self.strandSegments = {}
        self.merges = {}

    def __nextVertex(self, prvStrand, prv, vert):
        """Returns a pair (nxt, s) where nxt is the ID of the next vertex from
        vert along the strand connected to prev, and s is the strand (1 or 2)
        along which we travel through vert. prvStrand is the strand (1 or 2)
        travelled through prv."""
        if vert in self.Hverts:
            (prv1,nxt1,prv2,nxt2) = self.Hverts[vert][1:5]
            if prv1 == prv2:
                strandInfo = self.Hverts[vert][6]
        elif vert in self.Xverts:
            (prv1,nxt1,prv2,nxt2) = self.Xverts[vert][1:5]
            if prv1 == prv2:
                strandInfo = self.Xverts[vert][6]
        else:
            raise Exception("Error in obtaining next vertex: argument is not\
            a vertex")

        # if the previous vertex along both strands is the same, we use the
        # information passed by prvStrand.
        if prv1 == prv2:
            si0 = int(strandInfo[0])
            si2 = int(strandInfo[2])
            if prvStrand == si0:
                if si2 == 1:
                    return (nxt1, 1)
                elif si2 == 2:
                    return (nxt2, 2)
            else:
                if si2 == 1:
                    return (nxt2, 2)
                elif si2 == 2:
                    return (nxt1, 1)

        elif prv == prv1:
            return (nxt1, 1)
        elif prv == prv2:
            return (nxt2, 2)
        else:
            raise Exception("Error in obtaining next vertex")


    def _vertexSequence(self):
        """Returns a generator that yields all H- and X-vertices in the graph
        in order from the initial to the terminal vertex, together with the
        strand of that vertex along which the generator travels. Return format
        is (v,s) with v the vertex and s the strand."""

        ## check that the graph has been adequately pruned.
        #if len(self.Sverts) != 1 or self.Sverts[0][1] is not None:
        #    raise Exception("Error: unexpected vertices. Graph must be pruned.")

        prvStrand = None
        prv = 0
        vert = self.Sverts[0][2]
        (nxt, strand) = self.__nextVertex(prvStrand, prv, vert)
        yield (vert, (None, nxt))

        numVertices = len(self.Hverts) + len(self.Xverts)
        #for i in range(2 * numVertices - 1):
        while nxt is not None:
            prvStrand = strand
            prv = vert
            vert = nxt
            (nxt, strand) = self.__nextVertex(prvStrand, prv, vert)
            yield (vert, (prv, nxt))

    def _createMerge(self):
        pass
