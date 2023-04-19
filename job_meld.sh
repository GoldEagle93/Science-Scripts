#!/bin/bash
#SBATCH --job-name=meldTest      # Job name
#SBATCH --mail-user=reza.esmaeeli@ufl.edu
#SBATCH --mail-type=BEGIN,END,FAIL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --ntasks=30      
#SBATCH --gpus-per-task=1
#SBATCH --cpus-per-gpu=1
#SBATCH --mem-per-cpu=2000mb            
#SBATCH --partition=gpu
#SBATCH --distribution=cyclic:cyclic
#SBATCH --reservation=perez
#SBATCH --mem-per-cpu=2000mb          # Memory per processor
#SBATCH --time=100:00:00              # Time limit hrs:min:sec
#SBATCH --output=equil%j.log     # Standard output and error log



source ~/.load_OpenMM

export OPENMM_CUDA_COMPILER=''

if [ -e remd.log ]; then             #If there is a remd.log we are conitnuing a killed simulation
    prepare_restart --prepare-run  #so we need to prepare_restart
      fi

srun --mpi=pmix_v1  launch_remd --debug
