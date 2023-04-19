#!/bin/bash

pdblink="https://files.rcsb.org/download/$1.pdb"
wget $pdblink

python pdb-stripper.py $1
