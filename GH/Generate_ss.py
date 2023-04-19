#! /usr/bin/env python
import fileinput
import sys
import requests

pdb = str(sys.argv[1])

def generate_ss(pdbID):

    """
    Generates a sequence file from the PDB file generated by previous scripts which is needed by PSIPRED for secondary structure prediction.
    """

    fasta = requests.get('https://www.rcsb.org/fasta/entry/' + pdbID).text
    index = 0
    aminoacids = ['A', 'C', 'D', 'E', 'F',
                'G', 'H', 'I', 'K', 'L',
                'M', 'N', 'O', 'P', 'Q',
                'R', 'S', 'T', 'U', 'V',
                'W', 'Y',
                'ALA', 'ARG', 'ASN', 'ASP', 'CYS',
                'GLU', 'GLN', 'GLY', 'HIS', 'ILE',
                'LEU', 'LYS', 'MET', 'PHE', 'PRO',
                'SER', 'THR', 'TRP', 'TYR', 'VAL',
                'SEC', 'PYL', 'HIE']
    aminoacidsVSdna = ['D', 'E', 'F',
                'H', 'I', 'K', 'L',
                'M', 'N', 'O', 'P', 'Q',
                'R', 'S', 'U', 'V',
                'W', 'Y']
    aa3to1letter = {'ALA' : 'A', 'ARG' : 'R', 'ASN' : 'N', 'ASP' : 'D', 'CYS' : 'C',
                    'GLU' : 'E', 'GLN' : 'Q', 'GLY' : 'G', 'HIS' : 'H', 'ILE' : 'I',
                    'LEU' : 'L', 'LYS' : 'K', 'MET' : 'M', 'PHE' : 'F', 'PRO' : 'P',
                    'SER' : 'S', 'THR' : 'T', 'TRP' : 'W', 'TYR' : 'Y', 'VAL' : 'V',
                    'SEC' : 'U', 'PYL' : 'O', 'HIE' : 'H'}

    sequence = ""
    output = open('%s-fasta.fa' %(sys.argv[1]), 'w')
    output.write('>' + sys.argv[1] + '\n')

    for line in open(sys.argv[1] + '-sep.pdb'):
        list = line.split()
        if list[0] == 'ATOM' and list[3] in aminoacids and list[4] != index:
            sequence += aa3to1letter[list[3]]
            index = list[4]

    output.write(sequence)
    # print(sequence)

    # for line in fasta.splitlines():
    #     # print(line)
    #     if '>' not in line and line.isupper() and any(aa in line for aa in aminoacidsVSdna):
    #         output.write(line)

    output.close()

generate_ss(pdb)