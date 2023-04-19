import subprocess

def bash(command):
    subprocess.run(command, shell=True, executable='/bin/bash')

# bash('deactivate')
# bash('conda deactivate')
# bash('ls -l')
bash('amber')
# bash('module restore Amber20')



# bashCommand = "deactivate; conda deactivate; module purge; module restore Amber20"
# process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
# output, error = process.communicate()
# subprocess.call('amber')
# subprocess.run(bashCommand.split(),shell=True, stdout=subprocess.PIPE, universal_newlines=True)


# # subprocess.call('amber')
# subprocess.run('pdb4amber -i ref.pdb -o ref-amber.pdb')

# with open('ref-amber.pdb', 'r') as f:
#     file = f.readlines()
# f.close()

# ref_atoms = ['N', 'C', 'CA', "O5'", "C4'", "C1'"]

# with open('reference.pdb', 'w') as o:
#     for line in file:
#         if line.split()[0] == 'ATOM' and line.split()[2] in ref_atoms:
#             o.write(line)
# o.close()

# with open('strip.cpptraj', 'w') as o:
#     o.write(
# """
# cpptraj topol.top<<EOF
# trajin topol.crd
# strip !@N,C,CA,O5',C4',C1' parmout reference.top
# go
# EOF

# """
#     )
# o.close()

# subprocess.run('cpptraj -i strip.cpptraj')

# for i in ['0', '1', '2', '3', '4']:
#     with open('rmsd%s.cpptraj'%i, 'w') as o:
#         o.write(
# """
# cpptraj<<EOF
# parm topol.top [system]
# parm reference.top [reftop]
# reference reference.pdb parm [reftop] [ref]
# trajin trajectory.0%s.dcd [system]

# rms fit @O5',C4',C1'
# rms traj1 @N,C,CA ref ref nofit out rmsd%s.dat

# go
# EOF
# """%(i, i)
#         )
#     o.close()

