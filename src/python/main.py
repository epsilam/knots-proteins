from computegamma import *

#PP = ProteinProjection('1L2Y',0,np.pi/3)
PP = ProteinProjection('7LT8',0,0)
PP.findBackbonePrimaryCrossings()
PP.findProteinStructureBonds()
PP.plotCrossings()
