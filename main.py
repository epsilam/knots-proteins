from src.gamma import *

P = ProteinGamma('1R02', np.array([1,-1,1]))
P.ProProj.plotCrossings()
P.printGamma()
