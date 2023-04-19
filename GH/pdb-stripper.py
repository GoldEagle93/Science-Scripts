#!/usr/bin/env python

import fileinput
import sys
import requests

pdb = str(sys.argv[1])

def strip(pdbID):

    """
    Fetches PDB file from RCSB and writes only lines with ATOM, TER and MASTER records.
    """

    pdbfile = requests.get('https://files.rcsb.org/download/' + pdbID + '.pdb').text
    newfile = open('%s-stripped.pdb' %(pdbID), 'w')

    for line in pdbfile.splitlines():
        # print(line)
        if line.startswith('ATOM') or line.startswith('TER') or line.startswith('MASTER') or line.startswith('END'):
            
            newfile.write(line + '\n')
    newfile.close()

strip(pdb)