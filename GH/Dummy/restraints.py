import numpy as np
from meld.remd import ladder, adaptor, leader
from meld import comm, vault
from meld import system
from meld import parse
from meld.system.restraints import LinearRamp,ConstantRamp

def non_interacting(g,s,i,j,name_i,name_j,scaler=None):
    g.append(s.restraints.create_restraint('distance', scaler,ramp=LinearRamp(0,100,0,1), r1=0.0, r2=3.0, r3=100.0, r4=101.0, k=250.0,
            atom_1_res_index=i, atom_1_name=name_i, atom_2_res_index=j, atom_2_name=name_j))

def exclude_restraint(g,s,i,j,name_i,name_j,scaler=None):
    g.append(s.restraints.create_restraint('distance', scaler,ramp=LinearRamp(0,100,0,1), r1=4.0, r2=5.0, r3=7.0, r4=8.0, k=250.0,
            atom_1_res_index=i, atom_1_name=name_i, atom_2_res_index=j, atom_2_name=name_j))

def append_restraint(g,s,i,j,name_i,name_j,scaler=None):
    try:
        g.append(s.restraints.create_restraint('distance', scaler,ramp=LinearRamp(0,100,0,1), r1=0.0, r2=0.0, r3=0.8, r4=1.0, k=250.0,
            atom_1_res_index=i, atom_1_name=name_i, atom_2_res_index=j, atom_2_name=name_j))
    except:
        pass


# Restraints on Protein CA
def make_cartesian_collections(s, scaler, residues, atoms = ['CA'], delta=0.35, k=250.):
    cart = []
    atoms_restraint = atoms
    sequence = [(i,j) for i,j in zip(s.residue_numbers,s.residue_names)]
    sequence = sorted(set(sequence))
    sequence = dict(sequence)
    #Residues are 1 based
    #index of atoms are 1 base
    for i in residues:
        # print i
        if sequence[i] in 'ACE':
            atoms_restraint = ['CH3','C','O']
        elif sequence[i] in 'NME':
            atoms_restraint = ['CH3','N']
        else:
            atoms_restraint = atoms
        for b in atoms_restraint:
            try:
                #print(i,b)
                atom_index = s.index_of_atom(i,b) - 1
                x,y,z = s.coordinates[atom_index]/10.
                rest = s.restraints.create_restraint('cartesian',scaler, res_index=i, atom_name=b,
                    x=x, y=y, z=z, delta=delta,force_const=k)
                cart.append(rest)
            except:
                pass
    return cart

def create_contacts(s,scaler=None,group1=[],group2=[],group1sel=None,group2sel=None,group_active=1,collection_active=2):
    #Group1 will be any heavy atom in protein 
    #Group2 will be backbone atoms in peptide
    #Should have two contacts working out of all possibilities.
    scaler = scaler if scaler else s.restraints.create_scaler('nonlinear', alpha_min=0.4, alpha_max=1.0, factor=4.0)
    contact_restraints = []
    #Heavy atoms for each residue
    atoms = {"ALA":['N','C','O','CA','CB'],
             "VAL":['N','C','O','CA','CB','CG1','CG2'],
             "LEU":['N','C','O','CA','CB','CG','CD1','CD2'],
             "ILE":['N','C','O','CA','CB','CG1','CG2','CD1'],
             "PHE":['N','C','O','CA','CB','CG','CD1','CE1','CZ','CE2','CD2'],
             "TRP":['N','C','O','CA','CB','CG','CD1','NE1','CE2','CZ2','CH2','CZ3','CE3','CD2'],
             "MET":['N','C','O','CA','CB','CG','SD','CE'],
             "PRO":['N','C','O','CD','CG','CB','CA'],
             "ASP":['N','C','O','CA','CB','CG','OD1','OD2'],
             "GLU":['N','C','O','CA','CB','CG','CD','OE1','OE2'],
             "LYS":['N','C','O','CA','CB','CG','CD','CE','NZ'],
             "ARG":['N','C','O','CA','CB','CG','CD','NE','CZ','NH1','NH2'],
             "HIS":['N','C','O','CA','CB','CG','ND1','CE1','NE2','CD2'],
             "HID":['N','C','O','CA','CB','CG','ND1','CE1','NE2','CD2'],
             "HIE":['N','C','O','CA','CB','CG','ND1','CE1','NE2','CD2'],
             "HIP":['N','C','O','CA','CB','CG','ND1','CE1','NE2','CD2'],
             "GLY":['N','C','O','CA'],
             "SER":['N','C','O','CA','CB','OG'],
             "THR":['N','C','O','CA','CB','CG2','OG1'],
             "CYS":['N','C','O','CA','CB','SG'],
             "CYX":['N','C','O','CA','CB','SG'],
             "TYR":['N','C','O','CA','CB','CG','CD1','CE1','CZ','OH','CE2','CD2'],
             "ASN":['N','C','O','CA','CB','CG','OD1','ND2'],
             "GLN":['N','C','O','CA','CB','CG','CD','OE1','NE2'],
             "0EH":['N','C','O','CA'],
             "MK8":['N','C','O','CA'],
             "2JH":['N','C','O','CA','CB','CG','CD1','CD2','CE'],
             "BCK":['N','C','O','CA'],
             "PO2":['P','OP1','OP2'],
             "CB":['CB'],
             "DA":["C1'","C2","C2'","C3'","C4","C4'","C5","C5'","C6","C8","DA3","N1","N3","N6","N7","N9","O3'","O4'","O5'","OP1","OP2","P"],
             "DC":["C1'","C2","C2'","C3'","C4","C4'","C5","C5'","C6","N1","N3","N4","O2","O3'","O4'","O5'","OP1","OP2","P"],
             "DG":["C1'","C2","C2'","C3'","C4","C4'","C5","C5'","C6","C8","N1","N2","N3","N7","N9","O3'","O4'","O5'","O6","OP1","OP2","P"],
             "DT":["C1'","C2","C2'","C3'","C4","C4'","C5","C5'","C6","C7","N1","N3","O2","O3'","O4","O4'","O5'","OP1","OP2","P"]}


    #s.residue_numbers and s.residue_names have as many instances as atoms
    #we first use zip to create tuples and then set to create unique (disordered) lists of tuples
    #the sorted organizes them. Finally the dict creates an instance (one based) of residue to residue name
    sequence = [(i,j) for i,j in zip(s.residue_numbers,s.residue_names)]
    sequence = sorted(set(sequence))
    sequence = dict(sequence)

    for index_i in group1:
        for index_j in group2:
            res_i = sequence[index_i]
            res_j = sequence[index_j]
            if group1sel is None:
                atoms_i = atoms[res_i]
            else:
                atoms_i = atoms[group1sel]
                if group1sel == 'CB' and res_i == 'GLY':
                    atoms_j = ['CA']
            if group2sel is None:
                atoms_j = atoms[res_j]
            else:
                atoms_j = atoms[group2sel]
                if group2sel == 'CB' and res_j == 'GLY':
                    atoms_j = ['CA']

            local_contact = []
            for a_i in atoms_i:
                for a_j in atoms_j:
                    #print(index_i,a_i,index_j,a_j)
                    local_contact.append(s.restraints.create_restraint('distance', scaler, LinearRamp(0,100,0,1),r1=0.0, r2=0.0, r3=0.60, r4=0.70, k=250.0,
                atom_1_res_index=index_i, atom_1_name=a_i, atom_2_res_index=index_j, atom_2_name=a_j))

            contact_restraints.append(s.restraints.create_restraint_group(local_contact,int(group_active)))

            #print 'Contact:',index_i,index_j,res_i,res_j

    all_rest = len(contact_restraints)
    active = int( collection_active )
    #print active,all_rest
    #print contact_restraints
    s.restraints.add_selectively_active_collection(contact_restraints, active)


def cleft_contacts(s,scaler=None,groups=[],collection_active=2):
    #groups  will be tuples of pairs of lists and number of interactions to satisfy
    #will look for contacts between matchups in two lists. Activate on each list a number of interactions
    scaler = scaler if scaler else s.restraints.create_scaler('nonlinear', alpha_min=0.4, alpha_max=1.0, factor=4.0)
    contact_restraints = []
    #Heavy atoms for each residue
    atoms = {"ALA":['CB'],
             "VAL":['CB','CG1','CG2'],
             "LEU":['CB','CG','CD1','CD2'],
             "ILE":['CB','CG1','CG2','CD1'],
             "PHE":['CB','CG','CD1','CE1','CZ','CE2','CD2'],
             "TRP":['CB','CG','CD1','NE1','CE2','CZ2','CH2','CZ3','CE3','CD2'],
             "MET":['CB','CG','SD','CE'],
             "PRO":['CD','CG','CB'],
             "ASP":['CB','CG','OD1','OD2'],
             "GLU":['CB','CG','CD','OE1','OE2'],
             "LYS":['CB','CG','CD','CE','NZ'],
             "ARG":['CB','CG','CD','NE','CZ','NH1','NH2'],
             "HIS":['CB','CG','ND1','CE1','NE2','CD2'],
             "HID":['CB','CG','ND1','CE1','NE2','CD2'],
             "HIE":['CB','CG','ND1','CE1','NE2','CD2'],
             "HIP":['CB','CG','ND1','CE1','NE2','CD2'],
             "SER":['CB','OG'],
             "THR":['CB','CG2','OG1'],
             "CYS":['CB','SG'],
             "CYX":['CB','SG'],
             "TYR":['CB','CG','CD1','CE1','CZ','OH','CE2','CD2'],
             "ASN":['CB','CG','OD1','ND2'],
             "GLN":['CB','CG','CD','OE1','NE2'],
             "2JH":['CB','CG','CD1','CD2','CE']}

    #s.residue_numbers and s.residue_names have as many instances as atoms
    #we first use zip to create tuples and then set to create unique (disordered) lists of tuples
    #the sorted organizes them. Finally the dict creates an instance (one based) of residue to residue name
    sequence = [(i,j) for i,j in zip(s.residue_numbers,s.residue_names)]
    sequence = sorted(set(sequence))
    sequence = dict(sequence)

    for (group1,group2,group_active) in groups:
        local_contact = []
        for index_i in group1:
            for index_j in group2:
                res_i = sequence[index_i]
                res_j = sequence[index_j]
                atoms_i = atoms[res_i]
                atoms_j = atoms[res_j]

                for a_i in atoms_i:
                    for a_j in atoms_j:
                        #print a_i,a_j,group1,group2,group_active,res_i,res_j
                        local_contact.append(s.restraints.create_restraint('distance', scaler, LinearRamp(0,100,0,1),r1=0.0, r2=0.0, r3=0.65, r4=0.80, k=250.0,
                    atom_1_res_index=index_i, atom_1_name=a_i, atom_2_res_index=index_j, atom_2_name=a_j))

        #In this case we cannot enforce one restraint per pair of residues, so we are outside the loop
        contact_restraints.append(s.restraints.create_restraint_group(local_contact,int(group_active)))


    all_rest = len(contact_restraints)
    active = int( collection_active )
    #print active,all_rest
    #print contact_restraints
    s.restraints.add_selectively_active_collection(contact_restraints, active)


def create_specific_contacts(s,scaler=None,groups=[],offset=2,mapping=None,active=None):
    '''
    Groups is a list of lists. Each lists represents a different peptide residue interacting with different 
    atoms in the protein this will be a group. Each list starts with the number of restraints to satisfy in 
    each group. Then is a series of tuples, containing the residue number in the protein, the atom name and then the residue number in the peptide and the atom name in the peptide. We need a mapping to make sure numbering is accurate. Finally, the offset tells us how much to add to the peptide numbers to represent also the second peptide.
    '''
   
    scaler = scaler if scaler else s.restraints.create_scaler('nonlinear', alpha_min=0.4, alpha_max=1.0, factor=4.0)

    contact_restraints = []

    for g in groups:
        active = g[0]
        local_contact = []
        for r in g[1:]:
            #print(r)
            res_i,atoms_i,res_j,atoms_j = r
            res_i = mapping[res_i]
            res_j = mapping[res_j]

            local_contact.append(s.restraints.create_restraint('distance', scaler, LinearRamp(0,100,0,1),r1=0.0, r2=0.0, r3=0.5, r4=0.70, k=250.0,
                    atom_1_res_index=res_i, atom_1_name=atoms_i, atom_2_res_index=res_j, atom_2_name=atoms_j))
            #and the other peptide
            local_contact.append(s.restraints.create_restraint('distance', scaler, LinearRamp(0,100,0,1),r1=0.0, r2=0.0, r3=0.5, r4=0.70, k=250.0,
                    atom_1_res_index=res_i, atom_1_name=atoms_i, atom_2_res_index=int(res_j)+offset, atom_2_name=atoms_j))

        #In this case we cannot enforce one restraint per pair of residues, so we are outside the loop
        contact_restraints.append(s.restraints.create_restraint_group(local_contact,int(active)))


    s.restraints.add_selectively_active_collection(contact_restraints, active)


def keep_fixed_distance(filename, s, scaler):
    '''
    Used to keep a specific contact between residue 1 atom1 residue2 atom2 at a distance that is +- 1 the original distance
    '''
    dists = []
    rest_group = []
    lines = open(filename).read().splitlines()
    lines = [line.strip() for line in lines]
    for line in lines:
        if not line:
            dists.append(s.restraints.create_restraint_group(rest_group, 1))
            rest_group = []
        else:
            cols = line.split()
            i = int(cols[0])
            name_i = cols[1]
            j = int(cols[2])
            name_j = cols[3]
            dist = float(cols[4]) / 10.

            rest = s.restraints.create_restraint('distance', scaler,LinearRamp(0,100,0,1),
                                                 r1=dist-0.2, r2=dist-0.1, r3=dist+0.1, r4=dist+0.2, k=250,
                                                 atom_1_res_index=i, atom_2_res_index=j,
                                                 atom_1_name=name_i, atom_2_name=name_j)
            rest_group.append(rest)
    return dists

def contacts(filename, s, scaler):
    '''
    Used to keep a specific contact between residue 1 atom1 residue2 atom2 at a distance that is +- 1 the original distance
    '''
    dists = []
    rest_group = []
    lines = open(filename).read().splitlines()
    lines = [line.strip() for line in lines]
    for line in lines:
        if not line:
            dists.append(s.restraints.create_restraint_group(rest_group, 1))
            rest_group = []
        else:
            cols = line.split()
            i = int(cols[0])
            name_i = cols[1]
            j = int(cols[2])
            name_j = cols[3]
            dist = float(cols[4]) / 10.

            rest = s.restraints.create_restraint('distance', scaler,LinearRamp(0,100,0,1),
                                                 r1=0., r2=0., r3=dist, r4=dist+0.1, k=250,
                                                 atom_1_res_index=i, atom_2_res_index=j,
                                                 atom_1_name=name_i, atom_2_name=name_j)
            rest_group.append(rest)
    return dists

def make_hbond_restraint_file(sequence,offset,name='hbondsDNA.dat'):
    offset = int(offset)
    tot = len(sequence)
    print(sequence,tot)

    AT = '''{0} N1 {1} N3 3.0
    
{0} N6 {1} O4 3.0
    
'''
    
    CG = '''{0} O2 {1} N2 3.0
    
{0} N3 {1} N1 3.0
    
{0} N4 {1} O6 3.0
    
'''
    
    GC = '''{0} N2 {1} O2 3.0
    
{0} N1 {1} N3 3.0
    
{0} O6 {1} N4 3.0
    
'''
    
    TA = '''{0} N3 {1} N1 3.0
    
{0} O4 {1} N6 3.0
    
'''
    
    values = {'A':AT,'C':CG,'G':GC,'T':TA}
    
    with open(name,'w') as fo:
        for i,j in enumerate(sequence):
            ni = i + 1
            ncomplementary = tot*2 - ni +1
            ni = ni + offset
            ncomplementary += offset
            fo.write(values[j].format(ni,ncomplementary))


def readSeq(fname):
    '''
    Reads a Levitt like sequence file with the following format
    SEQ A >CCTTGGCTGACGTCAGCCAAG
    SEQ B >CTTGGCTGACGTCAGCCAAGG
    We want the first strand

    Or read a sequence file CCCCTACACT
    '''
    fi = open(fname,'r').readlines()[0].rstrip()
    if 'SEQ' in fi:
        return fi.split(">")[1]
    else:
        return fi


def get_dist_restraints(filename, s, scaler):
    dists = []
    rest_group = []
    lines = open(filename).read().splitlines()
    lines = [line.strip() for line in lines]
    for line in lines:
        if not line:
            dists.append(s.restraints.create_restraint_group(rest_group, 1))
            rest_group = []
        else:
            cols = line.split()
            i = int(cols[0])
            name_i = cols[1]
            j = int(cols[2])
            name_j = cols[3]
            dist = float(cols[4]) / 10.

            rest = s.restraints.create_restraint('distance', scaler,LinearRamp(0,100,0,1),
                                                 r1=0.0, r2=0.0, r3=dist, r4=dist+0.2, k=250,
                                                 atom_1_res_index=i, atom_2_res_index=j,
                                                 atom_1_name=name_i, atom_2_name=name_j)
            rest_group.append(rest)
    return dists
