import numpy as np

class Graph():
    def __init__(self):
        self._vertexIDTracker = 0

        self.Sverts = {}  # n : [res, prev, next]
        self.Hverts = {}  # n : [res, prev1, next1, prev2, next2]
        self.Xverts = {}  # n : [(resi1, resi2, resj1, resj2),
                          #      prevOver, nextOver,
                          #      prevUnder, nextUnder]

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
        prv2 = self.__newVertexID()
        nxt2 = self.__newVertexID()
        res2 = self.Sverts[vert2][0]
        vert2nxt = self.Sverts[vert2][2]
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
                          isOtherStrandHbond):
        """If the crossing has not been seen before, add 1 or 2 new X-vertices
        along a strand on the backbone depending on if the strand crosses the
        backbone or an H-bond. If the crossing has been seen before, connect it
        appropriately. Also connect the previous vertex prv, but not nxt.
        resPairs : a 4-tuple of the form (resi1, resi2, resj1, resj2).
        prv    : either an S-vertex, H-vertex, or X-vertex.
        nxt    : an S-vertex or X-vertex.
        isOver : bool denoting if the strand between prv, nxt is over or under.
        isOtherStrandHbond : bool
        Returns ID of furthest of the two new X-vertex created on the strand."""

        if resPairs not in self._seenXings:
            self._seenXings.add(resPairs)
            if isOtherStrandHbond:
                X1 = self.__newVertexID()
                X2 = self.__newVertexID()
                self._seenXings[resPairs] = [X1, X2]
                if isOver:
                    self.Xverts[X1] = [resPairs, prv, X2, None, None]
                    self.Xverts[X2] = [resPairs, X1, nxt, None, None]
                else:
                    self.Xverts[X1] = [resPairs, None, None, prv, X2]
                    self.Xverts[X2] = [resPairs, None, None, X1, nxt]
                self.__resolvePrvVertOnAddXvertAlongBB(prv, X2)
                return X2
            else:
                X1 = self.__newVertexID()
                self._seenXings[resPairs] = [X1]
                if isOver:
                    self.Xverts[X1] = [resPairs, prv, nxt, None, None]
                else:
                    self.Xverts[X1] = [resPairs, None, None, prv, nxt]
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
