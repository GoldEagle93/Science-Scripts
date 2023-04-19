import mdtraj as md
import sys

"""
Determines residue number of first and last protein and nucleic acid residue. Then generates a contact file based on distance of each residue i with i+6.
"""

pair=[]
pdb = md.load_pdb(sys.argv[1] + '-sep.pdb')
aminoacids = ['A', 'C', 'D', 'E', 'F',
              'G', 'H', 'I', 'K', 'L',
              'M', 'N', 'O', 'P', 'Q',
              'R', 'S', 'T', 'U', 'V',
              'W', 'Y',
              'ALA', 'ARG', 'ASN', 'ASP', 'CYS',
              'GLU', 'GLN', 'GLY', 'HIS', 'ILE',
              'LEU', 'LYS', 'MET', 'PHE', 'PRO',
              'SER', 'THR', 'TRP', 'TYR', 'VAL',
              'SEC', 'PYL',]
nucleobases = ['DA', 'DT', 'DC', 'DG', 'DU',
               'DA5', 'DT5', 'DC5', 'DG5', 'DU5',
               'DA3', 'DT3', 'DC3', 'DG3', 'DU3']

nucleic_start = 0
nucleic_end = 0
protein_start = 0
protein_end = 0
output = open('%s-contacts.dat' %(sys.argv[1]), 'w')

for line in open(sys.argv[1] + '-sep.pdb'):
    list = line.split()
    if list[0] == 'ATOM' and list[3] in nucleobases and nucleic_start == 0:
        nucleic_start = int(list[4])
    elif list[0] == 'ATOM' and list[3] in nucleobases and int(list[4]) >= nucleic_end:
        nucleic_end = int(list[4])
    elif list[0] == 'ATOM' and list[3] in aminoacids and protein_start == 0:
        protein_start = int(list[4])
    elif list[0] == 'ATOM' and list[3] in aminoacids and int(list[4]) >= protein_end:
        protein_end = int(list[4])

# for i in range(protein_start,protein_end):
#     print(i)

for i in range(protein_start-1,protein_end):
    for j in range(i+6,protein_end):
        pair.append([i,j])
#a= md.compute_contacts(pdb,contacts=pair,scheme='CA')
z=md.compute_contacts(pdb,contacts=pair,scheme='ca')
dist=z[0]

indices = dist < 0.8
(a,b) =  indices.shape
native_contacts = z[1][indices[0],:]
native_distances = dist[0,indices[0]]
# print(native_contacts)
for i,d in zip(native_contacts,native_distances):
   a,b = i
   output.write(str(a+1) + ' CA ' + str(b+1) + ' CA ' + str(d*10) + '\n\n')
#    print(a+1,'CA',b+1,'CA',d*10)
#    print()