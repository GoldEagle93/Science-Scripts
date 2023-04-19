with open('ref-pdb4amber.pdb', 'r') as f:
    file = f.readlines()
f.close()

ref_atoms = ['N', 'C', 'CA', "O5'", "C4'", "C1'"]

# with open('reference.pdb', 'w') as o:
#     for line in file:
#         if line.split()[0] == 'ATOM' and line.split()[2] in ref_atoms:
#             o.write(line)
# o.close()

# DNA_resnames = ['DT', 'DC', 'DG', 'DA', 'DT5', 'DC5', 'DG5', 'DA5']
# removed
# with open('reference.pdb', 'w') as o:
#     for i in range(len(file)):
#         if file[i].startswith('TER') and file[i+1].split()[3] in DNA_resnames:
#             fiveprime = file[i+1].split()[5]
#             for line in file[i:i+50]:
#                 if line.split()[5] == fiveprime and 
#         if line.split()[0] == 'ATOM' and line.split()[2] in ref_atoms:
#             o.write(line)
# o.close()


# with open('strip.cpptraj', 'w') as o:
#     o.write(
# """parm topol.top
# trajin topol.crd
# strip !@N,C,CA,O5',C4',C1' parmout reference.top
# go

# """
#     )
# o.close()

for i in ['0', '1', '2', '3', '4']:
    with open('rmsd%s.cpptraj'%i, 'w') as o:
        o.write(
"""parm topol.top [system]
parm reference.top [reftop]
reference reference.pdb parm [reftop] [ref]
trajin trajectory.0%s.dcd parm [system]

align @O5',C4',C1' ref [ref]
rms fit @O5',C4',C1'
rms traj%s ref [ref] @N,C,CA nofit out rmsd%s.dat

go

"""%(i, i, i)
        )
    o.close()

# :45,52,55,56,59,3,53,4,124,7,123,125,126,8,133,136,137,10,134,140,11,135,138,141,18,132,19,57,60,25,51,26@N,C,CA