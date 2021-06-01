import prody as pd
import pylab as pl
import igraph as ig
import numpy as np
from scipy.linalg import expm

class ProteinProjection:
                                         # angle1 and angle2 determine the 2d
                                         # plane onto which the protein is
                                         # projected. they are along the z-axis
                                         # and y-axis, respectively.
    def __init__(self, proteinName: str, angle1: float, angle2: float):

        self.p = pd.parsePDB(proteinName)
        backbone = self.p.select('protein') # removes unnecessary H2O molecules,
                                            # leaving only the backbone residues

        self.backboneCrossings = {}      # a dictionary containing crossings of
                                         # bonds on the backbone. the element
                                         # (i,j) is included if the bond
                                         # between i and i+1 crosses the bond
                                         # between j and j+1.

        self.resNums = np.unique(backbone.getResnums())

        self.resCoords = {} # dictionary of residue coordinates as averages
                            # of atom coordinates in the respective residue

        self.projectedResCoords = {} # a dictionary containing the residue
                                     # number together with the 2d coordinates
                                     # of the residue projected onto the given
                                     # plane.
        for i in self.resNums:
            res = backbone['A',i]
            self.resCoords[i] = np.average(res.getCoords(),axis=0)
                                     # we rotate the residue coordinates by
                                     # the two angles and then project onto
                                     # the xy-plane.

            M1 = expm(np.cross(np.eye(3),  # rotation matrix along z-axis
                            np.array([0,0,1])*angle1))
            M2 = expm(np.cross(np.eye(3),  # rotation matrix along y-axis
                            np.array([0,1,0])*angle2))
            rot = np.dot(M1,np.dot(M2,self.resCoords[i]))
            self.projectedResCoords[i] = rot[0:2] # take the x,y coords

                          # given two pairs (resi1,resi2) and (resj1,resj2)
                          # of res nums, determine if the bond between resi1
                          # and resi2 crosses the bond between resj1 and resj2
                          # in a given projection using some
                          # simple linear algebra.
    def checkBondCrossing(self, resi1, resi2, resj1, resj2):
        A = projectedResCoords[resi1]
        B = projectedResCoords[resi2]
        C = projectedResCoords[resj1]
        D = projectedResCoords[resj2]
        E = C - A
        R = B - A
        S = D - C
        crossRS = np.cross(R,S)
        t = np.cross(E,S) / crossRS
        u = np.cross(E,R) / crossRS
        if crossRS != 0:
            if 0 <= t <= 1 and 0 <= u <= 1:
                backboneCrossings[(resi1,resj1)] = A + t * R

                          # for each pair of two successive residues (i.e., two
                          # residues with a peptide bond), check if there is a
                          # crossing with another
    def findBackboneCrossings(self):
        for idx in range(len(self.resNums))-1:
            resi1 = resNums[idx]
            resi2 = resNums[idx+1]
            for jdx in range(idx+1,len(self.resNums)-1):
                resj1 = resNums[jdx]
                resj2 = resNums[jdx+1]
                checkBondCrossing(resi1,resi2,resj1,resj2)



    def exportTangle(self): # exports crossings and H-vertices into a format
                        # readable by Mathematica
        pass
