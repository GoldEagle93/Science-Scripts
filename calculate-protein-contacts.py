import mdtraj as md

"""
Generates a contact file based on distance of each residue i with i+6.
"""

pair=[]
pdb = md.load_pdb('ref.pdb')
topology = pdb.topology
CAs = []
for i in topology.atoms:
    if i.name == 'CA':
        CAs.append(i.index)


pairs = []        
        
for i in CAs:
    ii = topology.atom(i).residue.index
    for j in CAs:
        jj = topology.atom(j).residue.index
        dist = md.compute_distances(pdb,[[i, j]])[0]
        if dist < 1.0 and ((abs(ii - jj)) > 6):
            pairs.append(str(ii+1) + ' ' + topology.atom(i).name + ' ' + str(jj+1) + ' ' + topology.atom(j).name + ' ' + str("%.2f" % (dist[0]*10)))

with open('protein-contacts.dat', 'w') as out:
    for i in pairs:
        out.write(i+'\n\n')
out.close()

print('Total pairs: ' + str(len(pairs)))
print('Last pairs: ' + str(pairs[-1]))
