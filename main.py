from src.gamma import *

PROTEIN_NAME = '1R02'
PROJECTION_PLANE_NORMAL = np.array([1,-1,1])
PLOT_GRAPH = False

P = ProteinGamma(PROTEIN_NAME, PROJECTION_PLANE_NORMAL)
if PLOT_GRAPH:
    P.ProProj.plotCrossings()
P.printGamma()
