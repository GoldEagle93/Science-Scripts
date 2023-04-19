#!/bin/bash
#SBATCH --job-name=migration      # Job name
#SBATCH --mail-type=BEGIN,END,FAIL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=reza.esmaeeli@ufl.edu    # Where to send mail
#SBATCH --ntasks=1                  # Number of MPI ranks
#SBATCH --cpus-per-task=1            # Number of cores per MPI rank 
#SBATCH --nodes=1                    # Number of nodes
#SBATCH --ntasks-per-node=1         # How many tasks on each node
#SBATCH --ntasks-per-socket=1        # How many tasks on each CPU or socket
#SBATCH --distribution=cyclic:cyclic # Distribute tasks cyclically on nodes and sockets
#SBATCH --mem-per-cpu=4000mb          # Memory per processor
#SBATCH --time=2-00:00:00              # Time limit hrs:min:sec
#SBATCH --output=sbatch.log     # Standard output and error log
#SBATCH --qos=alberto.perezant-b

cd /blue/alberto.perezant/reza/TF/1a74-op/

date

cp -r ./* /orange/alberto.perezant/Reza/TF/1a74-optimization/

date
