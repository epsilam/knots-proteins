import numpy as np

class Graph():
    def __init__(self):
        self._vertexIDTracker = 0

        self.Sverts = {}  # n : [res, prev, next]
        self.Hverts = {}  # n : [res, prev1, next1, prev2, next2]
        self.Cverts = {}  # n : [(resi1, resi2, resj1, resj2),
                          #      prevOver, nextOver,
                          #      prevUnder, nextUnder]
    def __newVertexID(self):
        """Return a new unique ID to keep track of vertices. To be used only
        when creating new vertices."""
        ID = self._vertexIDTracker
        self._vertexIDTracker = ID + 1
        return ID


    def _addSvert(self, res, prv, nxt):
        """Add and appropriately connect a new S-vertex. Arguments prv and nxt
        must be S-vertices, or None (in case of initial or terminal vertex)."""
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
        vert1 is either an S-vertex or an H-vertex.
        vert2 is an S-vertex."""
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

    def _addCvert(self, ):
