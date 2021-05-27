from prody import *
from pylab import *

p = parsePDB('7LT8')

# Bonds need to be set. You can set them automatically with inferBonds,
# which infers bonds based on distances between atoms.
p.inferBonds()

# You can select specific parts of a protein (the actual protein, the water
# molecules, etc). Below we select only the protein, so water molecules are not
# shown.
showProtein(p.select('protein'))
show()

# We can access atoms in a protein with the index operator, e.g., p[0].
# After bonds have been set, we can iterate through all the bonds of a
# particular atom with iterBonded().
# We can also find the number of atoms an atom is bonded with via numBonds()
print(p[1].numBonds())
for atom in p[1].iterBonded():
    print(atom)


#to iterate over residues, use
for res in p.iterResidues():
    pass
