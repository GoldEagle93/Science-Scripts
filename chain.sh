#!/bin/bash

script=$1
njob=$2

echo I\'ll submit the script $script will restart $njob times
echo If this is not what you want you have 10s to kill this process
echo The positional arguments are: PBS_script name number_of_restarts
#sleep 10


PID=$(sbatch $script|awk '{print $NF}')
echo $PID
for i in `seq 1 $njob`; do
        PID=$(sbatch --dependency=afterany:$PID  $script| awk '{print $NF}')
        echo $PID
    done
