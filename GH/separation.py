from prody import *
import matplotlib
from numpy import ones, zeros
import sys

pdb = str(sys.argv[1])

# generate PDB
def separate(pdbID):

     """
     Generates a new PDB file from the code provided as pdbID through sys.argv[1] with the following modifications:
     - Writes only protein and nucelic acid atoms, ommiting water molecules, ions, etc.
     - Separates the protein and DNA using Prody by a certain distance specified in moveAtoms.

     Also Generates a pdbID-sequnce.dat file including the sequence of the DNA the first DNA strand.
     """

     atoms = parsePDB("./" + pdbID + "-ref.pdb")
     sel1 = atoms.select('protein')
     sel2 = atoms.select('nucleic or resname DA5 or resname DT5 or resname DC5 or resname DG5 or resname DA3 or resname DT3 or resname DC3 or resname DG3')
     sel3 = atoms.select('protein or nucleic or resname DA5 or resname DT5 or resname DC5 or resname DG5 or resname DA3 or resname DT3 or resname DC3 or resname DG3')

     moveAtoms(sel1, by=ones(3) * 50)

     writePDB(pdbID + '-sep.pdb', sel3)

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

     print(sequence[:-seqLength])

     file=pdbID + '-seq.dat' 
     with open(file, 'w') as filetowrite:
          filetowrite.write(sequence[:-seqLength])

separate(pdb)