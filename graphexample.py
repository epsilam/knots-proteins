from src.proteingraph import *

P = ProteinGraph('1L2Y',0,np.pi/3)

print("\n LIST OF H-VERTICES:")
print(P.Hverts)
print("\n LIST OF X-VERTICES:")
print(P.Xverts)
print()

P.printVertex(21)
P.printVertex(72)
