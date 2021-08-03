from src.proteingraph import *
import sys
from shutil import copyfile
from datetime import datetime

class ProteinGamma(ProteinGraph):
    def __init__(self, proteinName, planeNormal):
        ProteinGraph.__init__(self, proteinName, planeNormal)

        self._labels = {} # dictionary containing labels for strands and
                         # H-contacts. Converts from the vertex IDs of the
                         # Graph class to labels used in Mathematica.
        self._labelCounter = 0

    def __strandLabel(self, strand):
        """Return ID of strand. If strand ID has not been encountered, a new ID
        is created."""
        if strand in self._labels:
            return self._labels[strand]
        ID = self._labelCounter
        self._labelCounter = ID + 1
        self._labels[strand] = ID
        return ID

    def __HcontactLabel(self, hvert):
        if hvert in self._labels:
            return self._labels[hvert]
        ID = self._labelCounter
        self._labelCounter = ID + 1
        self._labels[hvert] = ID
        return ID


    def __nextVertex(self, prvStrandNum, prv, vert):
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
        # information passed by prvStrandNum.
        if prv1 == prv2:
            si0 = int(strandInfo[0])
            si2 = int(strandInfo[2])
            if prvStrandNum == si0:
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


    def vertexSequence(self):
        """Returns a generator that yields all H- and X-vertices in the graph
        in order from the initial to the terminal vertex, together with the
        strand of that vertex along which the generator travels."""

        # Output format: [vType, vTypeAnnot, strand1, strand2, mergeStrand]
        # where vType is 'H' or 'X',
        # vTypeAnnot is either the label of the H-contact or '+' or '-'
        # strand1,strand2 are the two strands of the vertex
        # mergeStrand is the strand which we merge to the main strand

        prvStrandNum = None
        prv = 0
        vert = self.Sverts[0][2]
        (nxt, strandNum) = self.__nextVertex(prvStrandNum, prv, vert)

        if vert in self.Hverts:
            (prv1,nxt1,prv2,nxt2) = self.Hverts[vert][1:5]
            vType = ("H",self.__HcontactLabel(vert))
        elif vert in self.Xverts:
            (prv1,nxt1,prv2,nxt2) = self.Xverts[vert][1:5]
            if self.Xverts[vert][5] == True:
                vType = ("X","+")
            else:
                vType = ("X","-")

        sID1 = self.__strandLabel((prv1,nxt1))
        sID2 = self.__strandLabel((prv2,nxt2))
        vertexReturn = (*vType, sID1, sID2)
        ret = (*vType, sID1, sID2)
        yield [*ret,ret[1+strandNum]]

        numVertices = len(self.Hverts) + len(self.Xverts)
        #for i in range(2 * numVertices - 1):
        while nxt is not None:
            prvStrandNum = strandNum
            prv = vert
            vert = nxt
            (nxt, strandNum) = self.__nextVertex(prvStrandNum, prv, vert)

            if vert in self.Hverts:
                (prv1,nxt1,prv2,nxt2) = self.Hverts[vert][1:5]
                vType = ("H",self.__HcontactLabel(vert))
            elif vert in self.Xverts:
                (prv1,nxt1,prv2,nxt2) = self.Xverts[vert][1:5]
                if self.Xverts[vert][5] == True:
                    vType = ("X","+")
                else:
                    vType = ("X","-")

            if strandNum == 1:
                strand = (prv1, nxt1)
            elif strandNum == 2:
                strand = (prv2, nxt2)
            else:
                raise Exception("Error")

            sID1 = self.__strandLabel((prv1,nxt1))
            sID2 = self.__strandLabel((prv2,nxt2))
            ret = (*vType, sID1, sID2)
            yield [*ret,ret[1+strandNum]]

    def formatVertexMerge(self, vertexmergelist, seenStrandPairs):
        """Takes the output of vertexSequence and prints formatted Mathematica
        code to stdout."""
        vType = vertexmergelist[0]
        vTypeAnnot = vertexmergelist[1]
        strand1 = vertexmergelist[2]
        strand2 = vertexmergelist[3]
        mergeStrand = vertexmergelist[4]

        if (strand1,strand2) not in seenStrandPairs:
            seenStrandPairs.add((strand1,strand2))
            if vType == 'H':
                output = " // Subscript[mH," \
                    + str(strand1) + "," \
                    + str(strand2) + "," \
                    + str(mergeStrand) + "]"
                return output

            if vType == 'X':
                if vTypeAnnot == '+':
                    output = " // Subscript[mXp," \
                        + str(strand1) + "," \
                        + str(strand2) + "," \
                        + str(mergeStrand) + "]"
                    return output
                elif vTypeAnnot == '-':
                    output = " // Subscript[mXn," \
                        + str(strand1) + "," \
                        + str(strand2) + "," \
                        + str(mergeStrand) + "]"
        else: # case that the vertex has not been seen
            output = " // Subscript[m, s, " + str(mergeStrand) + " -> s]"
        return output

    def printGamma(self):
        seenStrandPairs = set() # to avoid adding vertices twice, we keep track
                                # of previously seen pairs of strands.

        outputFileName = "gammaOutput-" + str(datetime.now()) + ".m"
        copyfile("gamma.m", outputFileName)
        with open(outputFileName, 'a') as file:
            file.write("(* ::Input:: *)\n")
            file.write("(*Subscript[i, s]")
            gen = self.vertexSequence()
            for v in gen:
                file.write(self.formatVertexMerge(v, seenStrandPairs))
            file.write("*)\n")
