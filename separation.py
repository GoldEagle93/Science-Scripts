from prody import *
import matplotlib
from numpy import ones, zeros
import sys

# generate PDB

fetchPDB(sys.argv[1])

atoms = parsePDB(sys.argv[1])
sel1 = atoms.select('protein')
sel2 = atoms.select('nucleic')
sel3 = atoms.select('protein or nucleic')

moveAtoms(sel1, by=ones(3) * 30)

writePDB(sys.argv[1] + '-sep.pdb', sel3)

# Generate sequence.dat

sel4 = sel2.select("name C1'")
characters = "{'}"
sequence = ""

for i in sel4:
     x = str(set(i.getSequence()))
     for j in characters:
         x = x.replace(j, "")
     sequence = sequence + x.upper()

seqLength = int(len(sequence)/2)

# f = open('sequence.dat', "w")
# f.write(sequence[:-seqLength])
# f.close()

file='sequence.dat' 
with open(file, 'w') as filetowrite:
     filetowrite.write(sequence[:-seqLength])