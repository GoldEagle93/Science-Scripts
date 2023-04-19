import mdtraj as md
import sys

"""
Determines residue number of first and last protein and nucleic acid residue. Then generates a contact file based on distance of each residue i with i+6.
"""

pdb = md.load_pdb('ref.pdb')
topology = pdb.topology
CBs = []
O5s = []
for i in topology.atoms:
    if i.name == 'CA':
        CBs.append(i.index)
    elif i.name == "O5'":
        O5s.append(i.index)

pairs = []        
        
for i in CBs:
    for j in O5s:
        dist = md.compute_distances(pdb,[[i, j]])[0]
        if dist < 1.0:
            pairs.append(str(topology.atom(i).residue.index+1) + ' ' + topology.atom(i).name + ' ' + str(topology.atom(j).residue.index+1) + ' ' + topology.atom(j).name + ' 10')

with open('protein-DNA-contacts.dat' %(sys.argv[1]), 'w') as out:
    for i in pairs:
        out.write(i+'\n\n')
out.close()

print(len(pairs))
print(pairs[0])