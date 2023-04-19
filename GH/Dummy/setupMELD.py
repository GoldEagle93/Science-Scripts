#!/usr/bin/env python
# encoding: utf-8

import numpy as np
from meld.remd import ladder, adaptor, leader
from meld import comm, vault
from meld import system
from meld import parse
from meld.system.restraints import LinearRamp,ConstantRamp
import glob
from restraints import *
import dummy

# number of replicas
N_REPLICAS = 30
# number of steps (units of exchange period)
N_STEPS = 20000
# controles frequence of output
BLOCK_SIZE = 100

def gen_state_templates(index, templates,coor):
    n_templates = len(templates)
    # print index,n_templates,index%n_templates
    a = system.ProteinMoleculeFromPdbFile(templates[index%n_templates])
    b = system.SystemBuilder(forcefield="ff14sbside")
    dummy_patcher = dummy.DummyPatcher(n_dummies=1,coord_dummies=coor)
    c = b.build_system_from_molecules([a],patchers=[dummy_patcher])
    pos = c._coordinates
    vel = np.zeros_like(pos)
    alpha = index / (N_REPLICAS - 1.0)
    energy = 0
    return system.SystemState(pos, vel, alpha, energy,[999,999,999] )



# MAIN routine
def setup_system():
    # create the system starting from coordinates in template.pdb
    templates = glob.glob('TEMPLATES/*.pdb')
    p = system.ProteinMoleculeFromPdbFile(templates[0])
    b = system.SystemBuilder(forcefield="ff14sbside")
    # load non-standard AA force field params, bonds
    s = b.build_system_from_molecules([p])

    #Amber indexing of residues
    DNA1_at = s.index_of_atom(10,'N1') - 1
    # DNA2_at = s.index_of_atom(50,'N1') - 1
    Prot_at = s.index_of_atom(130,'CA') - 1
    DNA1_x,DNA1_y,DNA1_z = s.coordinates[DNA1_at]
    # DNA2_x,DNA2_y,DNA2_z = s.coordinates[DNA2_at]
    Prot_x,Prot_y,Prot_z = s.coordinates[Prot_at]

    #the z coordinate is close on all cases
    z = np.average([DNA1_z,DNA1_z])
    # print("z_average:", z)
    #Solve linear system of equations, make distance 200 angs
    from scipy.optimize import fsolve

    def equations(p):
        x, y = p
        return ( (x-DNA1_x)**2 + (y-DNA1_y)**2 + (z-DNA1_z)**2 - 0**2, (x-DNA1_x)**2 + (y-DNA1_y)**2 + (z-DNA1_z)**2 - 0**2)

    x, y =  fsolve(equations, (5, 5))

    print('DNA Coords: ')
    print(s.coordinates[DNA1_at])
    print('Dummy Coords: ')
    print("x: ",x)
    print("y: ",y)
    print("z: ",z)

    dummy_patcher = dummy.DummyPatcher(n_dummies=1,coord_dummies=[x,y,z])
    s = b.build_system_from_molecules([p], patchers=[dummy_patcher])

 

 
    # Create temperature ladder
    s.temperature_scaler = system.GeometricTemperatureScaler(0.0, 0.3, 300., 450.)

    # Keep protein dimer conformation fairly constant
    dist_scaler = s.restraints.create_scaler('nonlinear', alpha_min=0.4, alpha_max=1.0, factor=4.0)
    dist_scaler2 = s.restraints.create_scaler('nonlinear', alpha_min=0.3, alpha_max=0.9, factor=4.0)
    dist_scaler2 = s.restraints.create_scaler('nonlinear', alpha_min=0.4, alpha_max=1.0, factor=4.0)
    const_scaler = s.restraints.create_scaler('constant')
    dist = keep_fixed_distance('protein-contacts.dat',s,scaler=const_scaler)
    s.restraints.add_selectively_active_collection(dist,int(len(dist)))

    # Keep DNA hbonds
    #Read sequence file
    sequenceDNA = readSeq('sequence.dat')
    # sequenceDNA2 = readSeq('sequence2.dat')
    #Generate hbondsDNA.dat
    make_hbond_restraint_file(sequenceDNA,0)
    dist = keep_fixed_distance('hbondsDNA.dat',s,scaler=const_scaler)
    s.restraints.add_selectively_active_collection(dist,int(len(dist)))
    # make_hbond_restraint_file(sequenceDNA2,40,name='hbondsDNA2.dat')
    # dist = keep_fixed_distance('hbondsDNA2.dat',s,scaler=const_scaler)
    # s.restraints.add_selectively_active_collection(dist,int(len(dist)))

    # Keep DNA close to starting conformation
    rest = make_cartesian_collections(s, const_scaler, range(1,41),atoms=["C1'","C2","C2'","C3'","C4","C4'","C5","C5'","C6","C7","C8","DA3","N1","N2","N3","N4","N6","N7","N9","O2","O3'","O4","O4'","O5'","O6","OP1","OP2","P"])
    #These are the common atoms to all DNA bases including ends:
    #C1' C2 C2' C3' C4 C4' C5 C5' C6 N1 N3 O3' O4' O5'
    s.restraints.add_as_always_active_list(rest)

    # Create Contacts between protein and DNA
    dist_scaler = s.restraints.create_scaler('nonlinear', alpha_min=0.3, alpha_max=0.5, factor=4.0)
    dom1 = get_dist_restraints('DNA-contacts.dat',s,scaler= dist_scaler)
    s.restraints.add_selectively_active_collection(dom1,12)


    #constraint the dummy particle
    # rest = make_cartesian_collections(s, const_scaler, [365],atoms=["S0"])
    # s.restraints.add_as_always_active_list(rest)

    # Make sure protein is far from DNA at higher replicas for homogeneous sampling
    scaler3 = s.restraints.create_scaler('nonlinear',alpha_min=0.2,alpha_max=1.0, factor=4.0, strength_at_alpha_min=0.2, strength_at_alpha_max=1.0)

    conf_rest = []
    group1 = []
    group2 = []
    group1.append( (365,"S0") )
    for j in range(41,364):
        group2.append( (j,"CA") )
    positioner = s.restraints.create_scaler('linear_positioner',alpha_min=0.2, alpha_max=1.0, pos_min=1.5, pos_max=5.)
    conf_rest.append(s.restraints.create_restraint('com', scaler3,ramp=LinearRamp(0,100,0,1), 
                                                       force_const=2.0,group1=group1,group2=group2,
                                                       distance =positioner,weights1=None, weights2=None, dims='xyz'))
    s.restraints.add_as_always_active_list(conf_rest)





    #
    # Secondary Structure
    #
    ss_scaler = s.restraints.create_scaler('constant')
    ss_rests = parse.get_secondary_structure_restraints(filename='ss.dat', system=s,ramp=LinearRamp(0,100,0,1), scaler=ss_scaler,
            torsion_force_constant=2.5, distance_force_constant=2.5)
    n_ss_keep = int(len(ss_rests) * 0.96) 
    s.restraints.add_selectively_active_collection(ss_rests, n_ss_keep)



    # create the options
    options = system.RunOptions()
    options.implicit_solvent_model = 'gbNeck2'
    options.remove_com = False
    options.use_big_timestep = False # MD timestep (3.3 fs)
    options.use_bigger_timestep = True # MD timestep (4.0 fs)
    options.cutoff = 1.8 # cutoff in nm
    options.soluteDielectric = 1.
    #options.implicitSolventSaltConc = None

    options.use_amap = False # correction to FF12SB
    options.amap_beta_bias = 1.0
    options.timesteps = 11111 # number of MD steps per exchange
    options.minimize_steps = 20000 # init minimization steps

    # create a store
    store = vault.DataStore(s.n_atoms, N_REPLICAS, s.get_pdb_writer(), block_size=BLOCK_SIZE)
    store.initialize(mode='w')
    store.save_system(s)
    store.save_run_options(options)

    # create and store the remd_runner, sets up replica exchange details
    l = ladder.NearestNeighborLadder(n_trials=48)
    policy = adaptor.AdaptationPolicy(2.0, 50, 50)
    a = adaptor.EqualAcceptanceAdaptor(n_replicas=N_REPLICAS, adaptation_policy=policy)
    remd_runner = leader.LeaderReplicaExchangeRunner(N_REPLICAS, max_steps=N_STEPS, ladder=l, adaptor=a)
    store.save_remd_runner(remd_runner)

    # create and store the communicator
    c = comm.MPICommunicator(s.n_atoms, N_REPLICAS)
    store.save_communicator(c)

    # create and save the initial states
    # create and save the initial states, initialize each replica with a different template
    states = [gen_state_templates(i,templates,[x,y,z]) for i in range(N_REPLICAS)]
    store.save_states(states, 0)

    # save data_store
    store.save_data_store()

    return s.n_atoms

# RUN THE SETUP
setup_system()
