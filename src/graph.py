import numpy as np

class Graph():
    def __init__(self):
        self._vertexIDTracker = 0

        self.Sverts = {}  # n : [res, prev, next]
        self.Hverts = {}  # n : [res, prev1, next1, prev2, next2]
        self.Xverts = {}  # n : [(resi1, resi2, resj1, resj2),
                          #      prevOver, nextOver,
                          #      prevUnder, nextUnder,
                          #      isPositive]

        self._seenXings = {} # Seen but unresolved crossings. As we
                             # process all the crossings on a bond between
                             # two vertices at once, we only process one of
                             # the two strands, and add the crossing to the
                             # "seen crossings" dict. When the crossing is
                             # encountered for the second time, it is
                             # resolved and removed from the set.
                             # Element format:
                             # (resi1, resi2, resj1, resj2) : vList
                             # where vList is a list of X-vertices (1, 2, or 4)
                             # associated with the crossing.

    def __newVertexID(self):
        """Return a new unique ID to keep track of vertices. To be used only
        when creating new vertices."""
        ID = self._vertexIDTracker
        self._vertexIDTracker = ID + 1
        return ID

    def _addSvert(self, res, prv, nxt):
        """Add and appropriately connect a new S-vertex. Arguments prv and nxt
        can be None in the case of an initial or terminal vertex.
        prv : an S-vertex or None.
        nxt : an S-vertex or None.
        Returns ID of new S-vertex."""
        vertexID = self.__newVertexID()
        self.Sverts[vertexID] = [res, prv, nxt]
        if prv is not None:
            self.Sverts[prv][2] = vertexID
        if nxt is not None:
            self.Sverts[nxt][1] = vertexID
        return vertexID

    def _addHvert(self, vert1, vert2):
        """Add new H-vertex to the graph representing the H-bond between the
        residues corresponding to vertices vert1 and vert2, automatically
        creating two new S-vertices in the process to preserve the graph
        structure, and connecting these vertices accordingly.
        vert1 : an S-vertex or an H-vertex.
        vert2 : an S-vertex.
        Returns ID of new H-vertex."""
        vertexID = self.__newVertexID()
        prv1 = vert1
        try: # check if vert1 is an S-vertex or an H-vertex
            nxt1 = self.Sverts[vert1][2]
            self.Sverts[vert1][2] = vertexID
            res1 = self.Sverts[vert1][0]
        except KeyError:
            nxt1 = self.Hverts[vert1][2]
            self.Hverts[vert1][2] = vertexID
            res1 = self.Hverts[vert1][0]
        self.Sverts[nxt1][1] = vertexID # note: nxt1 must be an S-vertex
        prv2 = self.__newVertexID()
        nxt2 = self.__newVertexID()
        res2 = self.Sverts[vert2][0]
        vert2nxt = self.Sverts[vert2][2]
        self.Sverts[vert2nxt][1] = nxt2
        self.Hverts[vertexID] = [res1, prv1, nxt1, prv2, nxt2]
        self.Sverts[vert2][2] = prv2
        self.Sverts[prv2] = [res2, vert2, vertexID]
        self.Sverts[nxt2] = [res2, vertexID, vert2nxt]
        return vertexID

    def __resolvePrvVertOnAddXvertAlongBB(self, prv, xvert):
        """Since X-vertices are created in a chain, we modify the previous
        vertex according to whether it is an S-, H-, or X-vertex. Only to be
        used in _addXvertsAlongBB().
        prv : any vertex.
        xvert : an X-vertex."""
        try: # case that prv is an S-vertex
            self.Sverts[prv][2] = xvert
        except KeyError:
            try: # case that prv is an H-vertex
                self.Hverts[prv][2] = xvert
            except KeyError: # case that prv is an X-vertex. check if over.
                if self.Xverts[prv][1] is not None and \
                   self.Xverts[prv][2] is None:
                    self.Xverts[prv][2] = xvert
                elif self.Xverts[prv][3] is not None and \
                     self.Xverts[prv][4] is None:
                    self.Xverts[prv][4] = xvert

    def _addXvertsAlongBB(self,
                          resPairs,
                          prv, nxt,
                          isOver,
                          isOtherStrandHbond,
                          isPositive):
        """If the crossing has not been seen before, add 1 or 2 new X-vertices
        along a strand on the backbone depending on if the strand crosses the
        backbone or an H-bond. If the crossing has been seen before, connect it
        appropriately. Also connect the previous vertex prv, but not nxt.
        resPairs : a 4-tuple of the form (resi1, resi2, resj1, resj2).
        prv    : either an S-vertex, H-vertex, or X-vertex.
        nxt    : any vertex or None.
        isOver : bool denoting if the strand between prv, nxt is over or under.
        isOtherStrandHbond : bool
        Returns ID of furthest of the two new X-vertex created on the strand."""

        if resPairs not in self._seenXings:
            if isOtherStrandHbond:
                X1 = self.__newVertexID()
                X2 = self.__newVertexID()
                self._seenXings[resPairs] = [X1, X2]
                if isOver:
                    self.Xverts[X1] = [resPairs,prv,X2,None,None,False]
                    self.Xverts[X2] = [resPairs,X1,nxt,None,None,True]
                else:
                    self.Xverts[X1] = [resPairs,None,None,prv,X2,True]
                    self.Xverts[X2] = [resPairs,None,None,X1,nxt,False]
                self.__resolvePrvVertOnAddXvertAlongBB(prv, X1)
                return X2
            else:
                X1 = self.__newVertexID()
                self._seenXings[resPairs] = [X1]
                if isOver:
                    self.Xverts[X1] = [resPairs,prv,nxt,None,None,isPositive]
                else:
                    self.Xverts[X1] = [resPairs,None,None,prv,nxt,isPositive]
                self.__resolvePrvVertOnAddXvertAlongBB(prv, X1)
                return X1
        else: # case that crossing has been seen before.
            # here, we only need to consider travelling along the backbone,
            # because crossings along H-bonds are resolved after crossings along
            # the backbone. thus, there will be only element in the seenXings
            # entry.
            if len(self._seenXings[resPairs]) != 1: # DEBUG
                raise Exception("error in seen crossings")
            X1 = self._seenXings[resPairs][0]
            if isOver:
                self.Xverts[X1][1] = prv
                self.Xverts[X1][2] = nxt
            else:
                self.Xverts[X1][3] = prv
                self.Xverts[X1][4] = nxt

            del self._seenXings[resPairs]
            self.__resolvePrvVertOnAddXvertAlongBB(prv, X1)
            return X1

    def __resolvePrvVertOnAddXvertAlongHbond(self, prv1, nxt2, xvert1, xvert2):
        """Same as __resolvePrvVertOnAddXvertAlongBB(), but for along Hbonds.
        Since the two strands go in opposite directions, we use nxt2 in place of
        prv2. We have two cases: either prv1 and nxt2 are both the same
        H-vertex, or they are two distinct X-vertices.
        """
        try:
            self.Hverts[prv1][4] = xvert1
            self.Hverts[prv1][3] = xvert2
            if prv1 != nxt2:
                raise Exception("unexpected H-vertex mismatch")
        except KeyError:
            # check if over. since strands will be either both over or both
            # under, we only check one vertex.
            if self.Xverts[prv1][1] is not None and \
               self.Xverts[prv1][2] is None: # on the over strand
                self.Xverts[prv1][2] = xvert1
                self.Xverts[nxt2][1] = xvert2
            elif self.Xverts[prv1][3] is not None and\
                 self.Xverts[prv1][4] is None: # on the under strand
                self.Xverts[prv1][4] = xvert1
                self.Xverts[nxt2][3] = xvert2


    def _addXvertsAlongHbond(self,
                             resPairs,
                             prv1, nxt1, prv2, nxt2,
                             isOver,
                             isOtherStrandHbond,
                             isPositive):
        """Same as _addXvertsAlongBB, but goes along H-bonds, and returns ID of
        two X-vertices. Creates either 2 or 4 vertices, depending on whether
        the H-bond crosses a peptide bond or another H-bond.
        Recall that the two strands along an H-bond travel in opposite
        directions. By convention, we take the 1 strand to be travelling away
        from the H-vertex, and the 2 strand to be going toward the H-vertex.
        """
        if not isOtherStrandHbond:
            # here, the other strand is a peptide bond, so it will have been
            # seen before.
            [X1, X2] = self._seenXings[resPairs]


            if isOver:
                strandPrv = 1; strandNxt = 2
                if not isPositive:
                    Xtemp = X1; X1 = X2; X2 = Xtemp
            else:
                strandPrv = 3; strandNxt = 4
                if isPositive:
                    Xtemp = X1; X1 = X2; X2 = Xtemp

            self.Xverts[X1][strandPrv] = prv1
            self.Xverts[X1][strandNxt] = nxt1
            self.Xverts[X2][strandPrv] = prv2
            self.Xverts[X2][strandNxt] = nxt2

            self.__resolvePrvVertOnAddXvertAlongHbond(prv1, nxt2, X1, X2)
            return (X1, X2)

        else: # crossing of two H-bonds.
            if resPairs not in self._seenXings:
                # create 4 vertices.
                # by convention, we will take X-verticess 1 and 3 to be on the
                # strand going away from the H-vertex, and 2 and 4 going toward
                # the H-vertex.
                X1 = self.__newVertexID()
                X2 = self.__newVertexID()
                X3 = self.__newVertexID()
                X4 = self.__newVertexID()
                self._seenXings[resPairs] = [X1, X2, X3, X4]
                if isOver:
                    self.Xverts[X1] = [resPairs,prv1,X3,None,None,False]
                    self.Xverts[X3] = [resPairs,X1,nxt1,None,None,True]
                    self.Xverts[X2] = [resPairs,X4,nxt2,None,None,True]
                    self.Xverts[X4] = [resPairs,prv2,X2,None,None,False]
                else:
                    self.Xverts[X1] = [resPairs,None,None,prv1,X3,True]
                    self.Xverts[X3] = [resPairs,None,None,X1,nxt1,False]
                    self.Xverts[X2] = [resPairs,None,None,X4,nxt2,False]
                    self.Xverts[X4] = [resPairs,None,None,prv2,X2,True]
                self.__resolvePrvVertOnAddXvertAlongHbond(prv1, nxt2, X1, X2)
                return (X3, X4)
            else: # case that we have encountered this crossing before.
                  # notice that the ordering of X1,X2,X3,X4 changes (we rotate
                  # the square by 90 degrees).
                [X1, X2, X3, X4] = self._seenXings[resPairs]

                if isOver:
                    strandPrv = 1; strandNxt = 2
                    if not isPositive:
                        Xtemp = X1; X1 = X4; X4 = Xtemp
                        Xtemp = X2; X2 = X3; X3 = Xtemp

                else:
                    strandPrv = 3; strandNxt = 4
                    if isPositive:
                        Xtemp = X1; X1 = X4; X4 = Xtemp
                        Xtemp = X2; X2 = X3; X3 = Xtemp

                self.Xverts[X2][strandPrv] = prv1
                self.Xverts[X2][strandNxt] = X1

                self.Xverts[X1][strandPrv] = X2
                self.Xverts[X1][strandNxt] = nxt1

                self.Xverts[X3][strandPrv] = prv2
                self.Xverts[X3][strandNxt] = X4

                self.Xverts[X4][strandPrv] = X3
                self.Xverts[X4][strandNxt] = nxt2

                self.__resolvePrvVertOnAddXvertAlongHbond(prv1,nxt2,X2,X4)
                del self._seenXings[resPairs]
                return (X1, X3)

    def _pruneSvert(self, vert):
        """Prune the S-vertex vert from the graph, if it exists, and
        appropriately connect its neighbours."""
        if vert not in self.Sverts:
            raise Exception("Error in pruning: given vertex is not an S-vertex")

        prv = self.Sverts[vert][1]
        nxt = self.Sverts[vert][2]

        del self.Sverts[vert]

        if prv is not None:
            if prv in self.Sverts:
                self.Sverts[prv][2] = nxt
            elif prv in self.Hverts:
                if vert == self.Hverts[prv][2]: # vert is on the backbone strand
                    self.Hverts[prv][2] = nxt
                elif vert == self.Hverts[prv][4]: # vert is on the H-bond strand
                    self.Hverts[prv][4] = nxt
                else: # vert is on neither strand
                    raise Exception("Error in pruning: H-vertex error")
            elif prv in self.Xverts:
                if vert == self.Xverts[prv][2]:
                    self.Xverts[prv][2] = nxt
                elif vert == self.Xverts[prv][4]:
                    self.Xverts[prv][4] = nxt
                else:
                    raise Exception("Error in pruning: X-vertex error")
            else:
                raise Exception("Error in pruning: prv not found")

        if nxt is not None:
            if nxt in self.Sverts:
                self.Sverts[nxt][1] = prv
            elif nxt in self.Hverts:
                if vert == self.Hverts[nxt][1]:
                    self.Hverts[nxt][1] = prv
                elif vert == self.Hverts[nxt][3]:
                    self.Hverts[nxt][3] = prv
                else:
                    raise Exception("Error in pruning: H-vertex error")
            elif nxt in self.Xverts:
                if vert == self.Xverts[nxt][1]:
                    self.Xverts[nxt][1] = prv
                elif vert == self.Xverts[nxt][3]:
                    self.Xverts[nxt][3] = prv
                else:
                    raise Exception("Error in pruning: X-vertex error")
            else:
                raise Exception("Error in pruning: nxt not found")
