import subprocess
import mdtraj as md
import sys
import numpy as np

"""
Determines residue number of first and last protein and nucleic acid residue. Then generates a contact file based on distance of each residue i with i+6.
"""

with open('sequence.dat', 'r') as inputfile:
    sequence = inputfile.readlines()[0]
    SL = 2*len(sequence)
inputfile.close()

SL = int(2*len(sequence))
midpoint = int(SL/4)
left = list(range(1, midpoint+1))
right = list(range(midpoint+1, (int(SL/2))+1))
left = left + list(range(SL-len(left)+1, SL+1))
right = right + list(range(int(SL/2)+1, SL-int(len(left)/2)+1))


pdb = md.load_pdb('ref.pdb')
topology = pdb.topology
CAs = []
N1s = []
threshold = 12 #in A
for i in topology.atoms:
    if i.name == 'CA':
        CAs.append(i.index)
    elif i.name == "N1" and ('A' in i.residue.name or 'G' in i.residue.name):
        N1s.append(i.index)

# poslist = [i for i in pdb.xyz]
# reference = md.formats.PDBTrajectoryFile('ref-mdtraj.pdb', mode='w', standard_names=False)
# reference.write(poslist[0], topology)

projector = []
right_pairs = []
left_pairs = []
mask = []
c = 1
last_res = topology.atom(CAs[-1]).residue.index + 1

for i in N1s:
    n = 0
    collection = []
    for j in CAs:
        dist = md.compute_distances(pdb,[[i, j]])[0]
        if dist < threshold*.1:
            # collection.append(str(topology.atom(i).residue.index+1) + ' ' + topology.atom(i).name + ' ' + str(topology.atom(j).residue.index+1) + ' ' + topology.atom(j).name + ' 10')
            collection.append(str(last_res+c) + ' S0 ' + str(topology.atom(j).residue.index+1) + ' ' + topology.atom(j).name + ' ' + str(threshold))
    # print(len(collection))
    if len(collection)>3 and topology.atom(i).residue.index+1 in left:
        left_pairs.append(collection)
        for z in collection:
            if z.split()[2] not in mask:
                mask.append(z.split()[2])
        if str(topology.atom(i).residue.index+1) not in mask:
            mask.append(str(topology.atom(i).residue.index+1))
        if n == 0:
            projector.append(str(topology.atom(i).residue.index+1) + ' N1 ' + str(last_res+c) + ' S0')
            n += 1
        c += 1
    elif len(collection)>3 and topology.atom(i).residue.index+1 in right:
        right_pairs.append(collection)
        for z in collection:
            if z.split()[2] not in mask:
                mask.append(z.split()[2])
        if str(topology.atom(i).residue.index+1) not in mask:
            mask.append(str(topology.atom(i).residue.index+1))
        if n == 0:
            projector.append(str(topology.atom(i).residue.index+1) + ' N1 ' + str(last_res+c) + ' S0')
            n += 1
        c += 1

with open('dummy-contacts-left.dat', 'w') as out:
    for i in left_pairs:
        for j in i:
            out.write(j+'\n')
        out.write('\n')
out.close()

with open('dummy-contacts-right.dat', 'w') as out:
    for i in right_pairs:
        for j in i:
            out.write(j+'\n')
        out.write('\n')
out.close()

with open('projector.dat', 'w') as out:
    for i in projector:
        out.write(i+'\n')
    out.write('\n')
out.close()

print('groups on the left: %d'%len(left_pairs))
print('groups on the right: %d'%len(right_pairs))

with open('cluster-contacts.sh', 'w') as out:
    out.write(
"""#!/bin/bash
#SBATCH --job-name=cls      # Job name
#SBATCH --mail-type=END,FAIL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=reza.esmaeeli@ufl.edu    # Where to send mail
#SBATCH --ntasks=1                  # Number of MPI ranks
#SBATCH --cpus-per-task=1            # Number of cores per MPI rank 
#SBATCH --nodes=1                    # Number of nodes
#SBATCH --ntasks-per-node=1         # How many tasks on each node
#SBATCH --ntasks-per-socket=1        # How many tasks on each CPU or socket
#SBATCH --distribution=cyclic:cyclic # Distribute tasks cyclically on nodes and sockets
#SBATCH --mem-per-cpu=4000mb          # Memory per processor
#SBATCH --time=12:00:00              # Time limit hrs:min:sec
#SBATCH --output=equil%%j.log     # Standard output and error log
pwd; hostname; date

deactivate
conda deactivate
module purge
ml cuda/11.0.207 gcc/9.3.0 openmpi/4.0.4 mkl/2020.0.166 meld/0.4.19
source /home/alberto.perezant/Source/amber/amber.sh

for a in 0 1 2 3 4 
do
    b=`perl -e 'printf("%%02i", $ARGV[0]);' $a`
    extract_trajectory extract_traj_dcd --replica $a trajectory.$b.dcd
done

cat<<EFO>get_top.py
#! /usr/bin/env python
import pickle as cPickle
x = cPickle.load(open('Data/system.dat','rb'))
f = open('topol.top', 'w')
f.write(x.top_string)
g = open('topol.crd', 'w')
g.write(x._mdcrd_string)
EFO
python get_top.py

mkdir contact-clusters
cd contact-clusters
mkdir DNA PROT DNA-PROT

cpptraj ../topol.top<<EOF
trajin ../trajectory.00.dcd  10000 20000
trajin ../trajectory.01.dcd  10000 20000
trajin ../trajectory.02.dcd 10000 20000
trajin ../trajectory.03.dcd 10000 20000
trajin ../trajectory.04.dcd 10000 20000
rms PredSSE first @P,OP1,OP2 out trajrmsd.dat
cluster hieragglo epsilon 1.5 linkage rms @P,OP1,OP2 nofit  sieve 10 summary DNA/summary singlerepout representative repout DNA/unique repfmt pdb clusterout DNA/clusttraj clusterfmt netcdf avgout DNA/avg avgfmt pdb out DNA/frame_vs_cluster.txt
cluster hieragglo epsilon 1.5 linkage rms @CA,N,O,C nofit  sieve 10 summary PROT/summary singlerepout representative repout PROT/unique repfmt pdb clusterout PROT/clusttraj clusterfmt netcdf avgout PROT/avg avgfmt pdb out PROT/frame_vs_cluster.txt
cluster hieragglo epsilon 1.5 linkage rms :%s@CA,N,O,C,P,OP1,OP2 nofit  sieve 10 summary DNA-PROT/summary singlerepout representative repout DNA-PROT/unique repfmt pdb clusterout DNA-PROT/clusttraj clusterfmt netcdf  avgout DNA-PROT/avg avgfmt pdb out DNA-PROT/frame_vs_cluster.txt
go
EOF

"""%','.join(mask)
    )
    out.write('\n')
out.close()

for i in ['00', '01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29']:
    with open('rmsd%s.cpptraj'%i, 'w') as o:
        o.write(
"""parm topol.top [system]
parm reference.prmtop [reftop]
reference reference.inpcrd parm [reftop] [ref]
trajin trajectory.%s.dcd parm [system]

align @O5',C4',C1' ref [ref]
rms fit @O5',C4',C1'
rms traj%s ref [ref] :%s@N,C,CA nofit out rmsd%s.dat

go

"""%(i, i, ','.join(mask), i)
        )
    o.close()


# def bash(command):
#     subprocess.run(command, shell=True)

