import os
import sys

def determine_complex_sequence(pdb):
    # seed = [i for i in os.listdir('TEMPLATES') if '_' not in i and i.endswith('-BDNA.pdb')][0]
    lines = [i for i in open(pdb,'r').readlines()]
    lines = [i for i in lines if i.startswith('ATOM')]
    complex_sequence = []
    for line in range(len(lines)):
        split = lines[line].split()
        if split[4].isdigit():
            if lines[line-1].split()[4] != lines[line].split()[4]:
                if split[3] == 'HIE':
                    complex_sequence.append('HIS')
                else:
                    complex_sequence.append(split[3])
        elif not split[4].isdigit() and split[5].isdigit():
            if lines[line-1].split()[5] != lines[line].split()[5]:
                if split[3] == 'HIE':
                    complex_sequence.append('HIS')
                else:
                    complex_sequence.append(split[3])
    # print(complex_sequence)
    return complex_sequence
# seq = determine_complex_sequence()

def three2one(sequence):
    NAs = 'DA DT DC DG DA5 DT5 DC5 DG5 DA3 DT3 DC3 DG3'.split(' ')
    AAs = 'ALA ARG ASP ASN CYS GLU GLN GLY HIS HIE HID HIP ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL SEL'.split(' ')
    one_NAs = {'DA':'A',  'DT':'T',  'DC':'C',  'DG':'G',
               'DA5':'A', 'DT5':'T', 'DC5':'C', 'DG5':'G',
               'DA3':'A', 'DT3':'T', 'DC3':'C', 'DG3':'G',
               'A':'A',   'U':'U',   'C':'C',   'G':'G',
               'A5':'A',  'U5':'U',  'C5':'C',  'G5':'G',
               'A3':'A',  'U3':'U',  'C3':'C',  'G3':'G'}

    one_AAs = {'ALA':'A', 'ARG':'R', 'ASP':'D', 'ASN':'N',
               'CYS':'C', 'GLU':'E', 'GLN':'Q', 'GLY':'G',
               'HIS':'H', 'HIE':'H', 'HID':'H', 'HIP':'H',
               'ILE':'I', 'LEU':'L', 'LYS':'K', 'MET':'M',
               'PHE':'F', 'PRO':'P', 'SER':'S', 'THR':'T',
               'TRP':'W', 'TYR':'Y', 'VAL':'V', 'SEL':'Z'}
    one_letter_sequence = []

    for residue in sequence:
        if residue in NAs:
            one_letter_sequence.append(one_NAs[residue].ljust(1))
        elif residue in AAs:
            one_letter_sequence.append(one_AAs[residue].ljust(1))
            # one_letter_sequence.append(one_AAs[residue]+'--')

    if len(sequence) == len(one_letter_sequence):
        return ''.join(one_letter_sequence)
    else:
        print('input/output sequence length mismatch')
        # print(','.join(sequence))
        # print(','.join(one_letter_sequence))
        print(len(sequence))
        print(len(one_letter_sequence))

forward  = three2one(determine_complex_sequence(sys.argv[1]))[:20]
with open('sequence.dat', 'w') as of:
    of.write(forward)
of.close()

