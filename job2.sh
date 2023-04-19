#!/bin/bash
#SBATCH --job-name=script                 # Job name
#SBATCH --mail-type=FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=reza.esmaeeli@ufl.edu   # Where to send mail
#SBATCH --ntasks=1                          # Number of MPI ranks
#SBATCH --cpus-per-task=1                   # Number of cores per MPI rank 
#SBATCH --nodes=1                           # Number of nodes
#SBATCH --ntasks-per-node=1                 # How many tasks on each node
#SBATCH --ntasks-per-socket=1               # How many tasks on each CPU or socket
#SBATCH --distribution=cyclic:cyclic        # Distribute tasks cyclically on nodes and sockets
#SBATCH --mem-per-cpu=4000mb                # Memory per processor
#SBATCH --time=7-00:00:00                   # Time limit hrs:min:sec
#SBATCH --output=script2.log                # Standard output and error log
#SBATCH --qos=alberto.perezant

pwd; hostname; date

deactivate
conda deactivate
module purge
# module load conda
# conda activate RF2NA

module load python


python funnel2.py


