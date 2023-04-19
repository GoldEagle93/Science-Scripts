#! /usr/bin/env python
import fileinput
import sys
import requests
import colorama
from colorama import Fore, Style

pdb = str(sys.argv[1])

def standardize(protein_sequence, horiz):

    """
    Generates a special format for the secondary structure based on sequence and PSIPRED output in which random coils are representated by ".", helices are represented by "H" and strands are represented by "E".
    """

    with open(protein_sequence, 'r') as inputfile:
        inputseq = inputfile.readlines()
    inputfile.close()
    with open(horiz, 'r') as inputfile:
        inputss = inputfile.readlines()
    inputfile.close()
    
    conf = ""
    pred = ""
    aa = ""
    standardSS = []
    finalSS = ""

    for line in inputss:
        list = line.split()
        if len(list) >=1:
            if list[0] == 'Conf:':
                conf += list[1]
            elif list[0] == 'Pred:':
                pred += list[1]
            elif list[0] == 'AA:':
                aa += list[1]

    for i in range(len(conf)):
        if int(conf[i]) <= 3:
            standardSS.append('.')
        elif int(conf[i]) > 3 and pred[i] == 'C':
            standardSS.append('.')
        elif int(conf[i]) > 3 and pred[i] == 'H':
            standardSS.append('H')
        elif int(conf[i]) > 3 and pred[i] == 'E':
            standardSS.append('E')

    print(Fore.BLUE + "WARNING: SS of last 10 residues are: %s" %(''.join(standardSS[-10:])) + Fore.BLACK)

    for i in range(len(standardSS)-4):
        if standardSS[i] == '.':
            if standardSS[i+1] == '.' or standardSS[i+2] == '.' or standardSS[i+3] == '.' or standardSS[i+4] == '.':
                standardSS[i+1] = '.'

    for i in range(len(str(inputseq[0]))):
        finalSS += '..'

    finalSS += ''.join(standardSS)

    output = open('ss.dat', 'w')
    output.write(finalSS)

standardize(sys.argv[1], sys.argv[2])