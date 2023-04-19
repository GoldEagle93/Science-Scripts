#!/bin/bash
#SBATCH --job-name=equil0      # Job name
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
#SBATCH --output=equil%j.log     # Standard output and error log
pwd; hostname; date

deactivate
conda deactivate
module purge
ml cuda/11.0.207 gcc/9.3.0 openmpi/4.0.4 mkl/2020.0.166 meld/0.4.19
source /home/alberto.perezant/Source/amber/amber.sh

#for a in 0 1 2 3 4 
#do
#    b=`perl -e 'printf("%02i", $ARGV[0]);' $a`
#    extract_trajectory extract_traj_dcd --replica $a trajectory.$b.dcd
#done



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

mkdir Cpptraj_linkage_sieve_eps_1.5
cd Cpptraj_linkage_sieve_eps_1.5
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
cluster hieragglo epsilon 1.5 linkage rms @CA,N,O,C,P,OP1,OP2 nofit  sieve 10 summary DNA-PROT/summary singlerepout representative repout DNA-PROT/unique repfmt pdb clusterout DNA-PROT/clusttraj clusterfmt netcdf  avgout DNA-PROT/avg avgfmt pdb out DNA-PROT/frame_vs_cluster.txt
go
EOF

