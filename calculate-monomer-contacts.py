import mdtraj as md
import sys

"""
Determines residue number of first and last protein and nucleic acid residue. Then generates a contact file based on distance of each residue i with i+6.
"""

pair=[]
pdb = md.load_pdb(sys.argv[1])
topology = pdb.topology
CAs = []
for i in topology.atoms:
    if i.name == 'CA':
        CAs.append(i.index)


pairs = []        

for i in CAs:
    ii = topology.atom(i).residue.index
    for j in CAs:
        ji = topology.atom(j).residue.index
        dist = md.compute_distances(pdb,[[i, j]])[0]
        if dist < 0.8 and ((abs(ii - ji)) > 6) and \
        ((ii > 366 and ji > 366) or \
        (ii < 366 and ji < 366) or \
        (ii > 246 and ii < 487 and ji > 246 and ji < 487)):
            pairs.append(str(ii+1) + ' ' + topology.atom(i).name + ' ' + str(ji+1) + ' ' + topology.atom(j).name + ' ' + str(dist[0]))


with open('protein-contacts.dat', 'w') as out:
    for i in pairs:
        out.write(i+'\n\n')
out.close()

print(len(pairs))
print(pairs[0])
