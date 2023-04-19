import os
import glob
import dummy
import numpy as np
from meld.remd import ladder, adaptor, leader
import meld.system.montecarlo as mc
from meld import system
from meld.system import patchers
from meld import comm, vault
from meld import parse
from meld import remd
from meld.system import param_sampling
from openmm import unit as u
import mdtraj as md

N_REPLICAS = 30
N_STEPS = 20000
BLOCK_SIZE = 50

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

def make_cartesian_collections(s, scaler, ramp, residues, atoms, delta, k):
    cart = []
    atoms_restraint = atoms
    # sequence = [(i,j) for i,j in zip(s.residue_numbers,s.residue_names)]
    # sequence = sorted(set(sequence))
    # sequence = dict(sequence)
    # #Residues are 1 based
    # #index of atoms are 1 base
    for i in residues:
    #     # print i
    #     if sequence[i] in 'ACE':
    #         atoms_restraint = ['CH3','C','O']
    #     elif sequence[i] in 'NME':
    #         atoms_restraint = ['CH3','N']
    #     else:
    #         atoms_restraint = atoms
        for b in atoms_restraint:
            # print('i in residues is: ', i)
            # print('sequence at i is: ', seq[i])
            # print('sequence at i-1 is: ', seq[i-1])
            atom_index = s.index.atom(i-1, b, expected_resname=seq[i-1])
            # print(s.coordinates[atom_index])
            x,y,z = s.coordinates[atom_index]
            rest = s.restraints.create_restraint('cartesian',scaler, ramp, atom_index=atom_index,
                x=x*u.nanometer, y=y*u.nanometer, z=z*u.nanometer, delta=delta*u.nanometer,force_const=k*u.kilojoule_per_mole/u.nanometer **2)
            cart.append(rest)
    # print(x, y, z)
    # print(rest)
    # print(cart)
    return cart

def define_dummy_cartesian(filename, s, scaler, ramp, seq, p1, p2, p3, p4, trust_in_group=2):
    dists = []
    rest_group = []
    lines = open(filename).read().splitlines()
    lines = [line.strip() for line in lines]
    for line in lines:
        if not line:
            dists.append(s.restraints.create_restraint_group(rest_group, trust_in_group))
            rest_group = []
        else:
            cols = line.split()
            i = int(cols[0])-1
            name_i = cols[1]
            j = int(cols[2])-1
            name_j = cols[3]

            rest = s.restraints.create_restraint('distance', scaler, ramp,
                                                 r1=p1, r2=p2, r3=p3, r4=p4,
                                                 k=250*u.kilojoule_per_mole/u.nanometer **2,
                                                 atom1=s.index.atom(i,name_i, expected_resname=seq[i]),
                                                 atom2=s.index.atom(j,name_j, expected_resname='SDM'))
            rest_group.append(rest)
    return dists

def define_dummy_contacts(filename, s, scaler, ramp, seq, p1, p2, p3, p4, trust_in_group=2):
    dists = []
    rest_group = []
    lines = open(filename).read().splitlines()
    lines = [line.strip() for line in lines]
    for line in lines:
        if not line:
            dists.append(s.restraints.create_restraint_group(rest_group, trust_in_group))
            rest_group = []
        else:
            cols = line.split()
            i = int(cols[0])-1
            name_i = cols[1]
            j = int(cols[2])-1
            name_j = cols[3]

            rest = s.restraints.create_restraint('distance', scaler, ramp,
                                                 r1=p1, r2=p2, r3=p3, r4=p4,
                                                 k=250*u.kilojoule_per_mole/u.nanometer **2,
                                                 atom1=s.index.atom(i,name_i, expected_resname='SDM'),
                                                 atom2=s.index.atom(j,name_j, expected_resname=seq[j]))
            rest_group.append(rest)
    return dists

def keep_fixed_distance(filename, s, scaler, ramp, seq):
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
            i = int(cols[0])-1
            name_i = cols[1]
            j = int(cols[2])-1
            name_j = cols[3]
            dist = float(cols[-1])*.1

            rest = s.restraints.create_restraint('distance', scaler, ramp,
                                                 r1=(dist-0.2)*u.nanometer, r2=(dist-0.1)*u.nanometer, r3=(dist+0.1)*u.nanometer, r4=(dist+0.2)*u.nanometer,
                                                 k=250*u.kilojoule_per_mole/u.nanometer **2,
                                                 atom1=s.index.atom(i,name_i, expected_resname=seq[i]),
                                                 atom2=s.index.atom(j,name_j, expected_resname=seq[j]))
            rest_group.append(rest)
    return dists


def get_dist_restraints_hydrophobe(filename, s, scaler, ramp, seq, p1, p2, p3, p4, trust_in_group=2):
    dists = []
    rest_group = []
    lines = open(filename).read().splitlines()
    lines = [line.strip() for line in lines]
    for line in lines:
        if not line:
            dists.append(s.restraints.create_restraint_group(rest_group, trust_in_group))
            rest_group = []
        else:
            cols = line.split()
            i = int(cols[0])-1
            name_i = cols[1]
            j = int(cols[2])-1
            name_j = cols[3]

            rest = s.restraints.create_restraint('distance', scaler, ramp,
                                                 r1=p1, r2=p2, r3=p3, r4=p4,
                                                 k=250*u.kilojoule_per_mole/u.nanometer **2,
                                                 atom1=s.index.atom(i,name_i, expected_resname=seq[i]),
                                                 atom2=s.index.atom(j,name_j, expected_resname=seq[j]))
            rest_group.append(rest)
    return dists

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


with open('projector.dat', 'r') as inputfile:
    projection = inputfile.readlines()
inputfile.close()
projection = [i for i in projection if i != '\n']
num_of_dum = len(projection)

# Determine the coordinates on which to put the dummy particles
def project(map_list, s):
    coords = []
    for contact in map_list:
        split = contact.split()
        resid = int(split[0])-1
        resname = split[1]
        atom = s.index.atom(resid, resname, expected_resname=seq[resid])
        # atom = system.indexing.AtomIndex(int(i.split()[0])) - 1
        print(atom)
        x,y,z = s.coordinates[atom]
        atom_coords = [x*10,y*10,z*10]
        print(atom_coords)
        # print(s.coordinates[-1])
        coords.append(atom_coords)    
    return coords

# Determine the sequence of the system from PDB seed, in three-letter code
def determine_complex_sequence():
    seed = [i for i in os.listdir('TEMPLATES') if '_' not in i and i.endswith('-BDNA.pdb')][0]
    lines = [i for i in open('TEMPLATES/'+seed,'r').readlines()]
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
seq = determine_complex_sequence()

def gen_state(s, index):
    state = s.get_state_template()
    state.alpha = index / (N_REPLICAS - 1.0)
    return state

def setup_system():
    # load the sequence
    templates = glob.glob('TEMPLATES/*.pdb')
    # sequence = parse.get_sequence_from_AA1(filename='sequence.dat')
    # n_res = len(sequence.split())

    # build the system
    # p = system.subsystem.SubSystemFromSequence(sequence)
    # b = system.builder.SystemBuilder()
    p = system.subsystem.SubSystemFromPdbFile(templates[0])
    b = system.builder.SystemBuilder(forcefield="ff14sbside")
    #rdc_patcher = patchers.RdcAlignmentPatcher(n_tensors=1)
    #s = b.build_system([p], patchers=[rdc_patcher])
    s = b.build_system([p])
    coords = project(projection, s)
    dummy_patcher = dummy.DummyPatcher(n_dummies=num_of_dum,coord_dummies=coords)
    #coords = project(projection, s)
    s = b.build_system([p], patchers=[dummy_patcher])
    n_res = s.residue_numbers[-1]

    s.temperature_scaler = system.temperature.GeometricTemperatureScaler(0, 0.3, 300.*u.kelvin, 550.*u.kelvin)

    ramp = s.restraints.create_scaler('nonlinear_ramp', start_time=1, end_time=200,
                                        start_weight=1e-3, end_weight=1, factor=4.0)

    ss_scaler = s.restraints.create_scaler('constant')
    ss_rests = parse.get_secondary_structure_restraints(filename='ss.dat', system=s, scaler=ss_scaler,
            ramp=ramp, torsion_force_constant=0.01*u.kilojoule_per_mole/u.degree **2, distance_force_constant=2.5*u.kilojoule_per_mole/u.nanometer **2, quadratic_cut=2.0*u.nanometer)
    n_ss_keep = int(len(ss_rests) * 0.85)
    s.restraints.add_selectively_active_collection(ss_rests, n_ss_keep)

    const_scaler = s.restraints.create_scaler('constant')
    for i in projection:
        # print('i in projection is: ', int(i.split()[0]),int(i.split()[0])+1)
        rest = make_cartesian_collections(s, const_scaler, ramp, range(int(i.split()[0]),int(i.split()[0])+1),atoms=["N1"], delta=0.6, k=250.)
        s.restraints.add_as_always_active_list(rest)
    #
    # Dummy Contacts
    #
    dummy_cartesian = define_dummy_cartesian('projector.dat', s, const_scaler, ramp, seq,  0*u.nanometer, 0.1*u.nanometer, 0.2*u.nanometer, 0.3*u.nanometer, trust_in_group=len(projection))
    s.restraints.add_selectively_active_collection(dummy_cartesian,int(len(dummy_cartesian)))



    # create_hydrophobes(s,group_1=subset1,group_2=subset1,CO=False)

    dist = keep_fixed_distance('protein-contacts.dat', s, const_scaler, ramp, seq)
    s.restraints.add_selectively_active_collection(dist,int(len(dist)))

    sequenceDNA = readSeq('sequence.dat')
    make_hbond_restraint_file(sequenceDNA,0)
    dist = keep_fixed_distance('hbondsDNA.dat', s, const_scaler, ramp, seq)
    s.restraints.add_selectively_active_collection(dist,int(len(dist)))

    dummy_scaler = s.restraints.create_scaler('nonlinear',alpha_min=0.0,alpha_max=1.0, factor=4.0, strength_at_alpha_min=0.5, strength_at_alpha_max=1.0)

    pos1 = s.restraints.create_scaler('linear_positioner',alpha_min=0.0, alpha_max=1.0, pos_min=0.0*u.nanometer, pos_max=6.7*u.nanometer)
    pos2 = s.restraints.create_scaler('linear_positioner',alpha_min=0.0, alpha_max=1.0, pos_min=0.0*u.nanometer, pos_max=6.8*u.nanometer)
    pos3 = s.restraints.create_scaler('linear_positioner',alpha_min=0.0, alpha_max=1.0, pos_min=1.2*u.nanometer, pos_max=7.2*u.nanometer)
    pos4 = s.restraints.create_scaler('linear_positioner',alpha_min=0.0, alpha_max=1.0, pos_min=1.3*u.nanometer, pos_max=7.3*u.nanometer)

    dummy_contacts_left = define_dummy_contacts('dummy-contacts-left.dat', s, dummy_scaler, ramp, seq, pos1, pos2, pos3, pos4, trust_in_group=2)
    dummy_contacts_right = define_dummy_contacts('dummy-contacts-right.dat', s, dummy_scaler, ramp, seq, pos1, pos2, pos3, pos4, trust_in_group=2)


    #creates parameter for sampling for hydrophobic contacts
    # dists = get_dist_restraints_hydrophobe('hydrophobe.dat', s, scaler, ramp, seq)
    prior_left = param_sampling.ScaledExponentialDiscretePrior(u0=1.0, temperature_scaler=s.temperature_scaler, scaler=dummy_scaler)
    sampler_left = param_sampling.DiscreteSampler(int(1), int(1.00 * len(dummy_contacts_left)), 1)
    param_left = s.param_sampler.add_discrete_parameter("left", int(0.3 * num_of_dum), prior_left, sampler_left)
    s.restraints.add_selectively_active_collection(dummy_contacts_left, param_left)

    prior_right = param_sampling.ScaledExponentialDiscretePrior(u0=1.0, temperature_scaler=s.temperature_scaler, scaler=dummy_scaler)
    sampler_right = param_sampling.DiscreteSampler(int(1), int(1.00 * len(dummy_contacts_right)), 1)
    param_right = s.param_sampler.add_discrete_parameter("right", int(0.3 * num_of_dum), prior_right, sampler_right)
    s.restraints.add_selectively_active_collection(dummy_contacts_right, param_right)

    print('N1s:')
    print(s.coordinates[17])
    print(s.coordinates[80])
    print(s.coordinates[115])
    print(s.coordinates[459])
    print('dummies:')
    print(s.coordinates[-4:])
    # create the options
    options = system.options.RunOptions()
    options.implicit_solvent_model = 'gbNeck2'
    options.use_big_timestep = True
    options.cutoff = 1.8*u.nanometers
    options.remove_com = False
    options.use_amap = False
    options.amap_beta_bias = 1.0
    options.timesteps = 14286
    options.minimize_steps = 20000
    # options.min_mc = sched
    options.param_mcmc_steps=200

    # create a store
    store = vault.DataStore(s.get_state_template(),N_REPLICAS, s.get_pdb_writer(), block_size=BLOCK_SIZE)
    store.initialize(mode='w')
    store.save_system(s)
    store.save_run_options(options)

    # create and store the remd_runner
    l = ladder.NearestNeighborLadder(n_trials=48 * 48)
    policy_1 = adaptor.AdaptationPolicy(2.0, 50, 50)
    a = adaptor.EqualAcceptanceAdaptor(n_replicas=N_REPLICAS, adaptation_policy=policy_1, min_acc_prob=0.02)

    remd_runner = remd.leader.LeaderReplicaExchangeRunner(N_REPLICAS, max_steps=N_STEPS,
                                                            ladder=l, adaptor=a)
    store.save_remd_runner(remd_runner)

    # create and store the communicator
    c = comm.MPICommunicator(s.n_atoms, N_REPLICAS, timeout=60000)
    store.save_communicator(c)

    # create and save the initial states
    states = [gen_state(s, i) for i in range(N_REPLICAS)]
    store.save_states(states, 0)

    # save data_store
    store.save_data_store()

    return s.n_atoms


setup_system()
