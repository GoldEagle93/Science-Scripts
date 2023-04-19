#!/usr/bin/env python
# encoding: utf-8

import sys
import numpy as np
from meld.remd import ladder, adaptor, master_runner
from meld import comm, vault
from meld import system
from meld import parse
from meld.system.restraints import LinearRamp,ConstantRamp
import glob
from restraints import *

# number of replicas
N_REPLICAS = 30
# number of steps (units of exchange period)
N_STEPS = 10000
# controles frequence of output
BLOCK_SIZE = 100

def gen_state_templates(index, templates):
    n_templates = len(templates)
    # print index,n_templates,index%n_templates
    a = system.ProteinMoleculeFromPdbFile(templates[index%n_templates])
    b = system.SystemBuilder(forcefield="ff14sbside")
    c = b.build_system_from_molecules([a])
    pos = c._coordinates
    vel = np.zeros_like(pos)
    alpha = index / (N_REPLICAS - 1.0)
    energy = 0
    return system.SystemState(pos, vel, alpha, energy,[999,999,999] )


# MAIN routine
def setup_system():
    # create the system starting from coordinates in template.pdb
    templates = glob.glob('%s-sep.pdb' %(sys.argv[1]))
    p = system.ProteinMoleculeFromPdbFile(templates[0])
    b = system.SystemBuilder(forcefield="ff14sbside")
    # load non-standard AA force field params, bonds
    s = b.build_system_from_molecules([p])
 
    # Create temperature ladder
    s.temperature_scaler = system.GeometricTemperatureScaler(0.0, 0.5, 300., 500.)

    # Keep protein dimer conformation fairly constant
    dist_scaler = s.restraints.create_scaler('nonlinear', alpha_min=0.4, alpha_max=1.0, factor=4.0)
    
    const_scaler = s.restraints.create_scaler('constant')
    dist = keep_fixed_distance('%s-contacts.dat' %(sys.argv[1]),s,scaler=const_scaler)
    s.restraints.add_selectively_active_collection(dist,int(len(dist)))

    # Keep DNA hbonds
    #Read sequence file
    sequenceDNA = readSeq('%s-seq.dat' %(sys.argv[1]))
    #Generate hbondsDNA.dat
    make_hbond_restraint_file(sequenceDNA,0)
    dist_scaler3 = s.restraints.create_scaler('nonlinear', alpha_min=0.9, alpha_max=1.0, factor=4.0)
    dist = keep_fixed_distance('hbondsDNA.dat',s,scaler=const_scaler)
    s.restraints.add_selectively_active_collection(dist,int(len(dist)))

    # Keep DNA close to starting conformation
    rest = make_cartesian_collections(s, const_scaler, range(1,43),atoms=["C1'","C2","C2'","C3'","C4","C4'","C5","C5'","C6","C7","C8","DA3","N1","N2","N3","N4","N6","N7","N9","O2","O3'","O4","O4'","O5'","O6","OP1","OP2","P"])
    # rest = make_cartesian_collections(s, const_scaler, range(1,16),atoms=["C1'", "C2", "C2'", "C3'", "C4", "C4'", "C5", "C5'", "C6", "N1", "N3", "O3'", "O4'"])
    #These are the common atoms to all DNA bases including ends:
    #C1' C2 C2' C3' C4 C4' C5 C5' C6 N1 N3 O3' O4' O5'
    s.restraints.add_as_always_active_list(rest)

    # Create Contacts between protein and DNA
    dom1 = get_dist_restraints('%s-DNA-contacts.dat' %(sys.argv[1]),s,scaler= dist_scaler)
    s.restraints.add_selectively_active_collection(dom1,int(len(dom1)))

    # Find Glycines and Restrain peptide within reasonable distance from DNA
    names  = np.array(s.atom_names)
    resid = np.array(s.residue_numbers)    
    # resnames = np.array(s.residue_names)
    select = names == 'CB'
    non_gly = resid[select]

    # scaler3 = s.restraints.create_scaler('nonlinear',alpha_min=0.7,alpha_max=1.0, factor=4.0, strength_at_alpha_min=1.0, strength_at_alpha_max=0.5)

    # conf_rest = []
    # group1 = []
    # group2 = []
    # for i in range(2,21):
    #     group1.append( (i,"O5'") )
    # for i in range(22,41):
    #     group1.append( (i,"O5'") )
    # for j in non_gly:
    #     group2.append( (j,"CB") )
    # positioner = s.restraints.create_scaler('linear_positioner',alpha_min=0.7, alpha_max=1.0, pos_min=10., pos_max=15.) 
    # conf_rest.append(s.restraints.create_restraint('com', scaler3,ramp=LinearRamp(0,100,0,1), 
    #                                                    force_const=75.0,group1=group1,group2=group2,
    #                                                    distance =positioner,weights1=None, weights2=None, dims='xyz'))
    # s.restraints.add_as_always_active_list(conf_rest)

    dist_scaler2 = s.restraints.create_scaler('nonlinear', alpha_min=0.7, alpha_max=1.0, factor=4.0)
    res_groups = get_distance_rests('%s-res_groups.dat' %(sys.argv[1]),s,scaler= dist_scaler2)
    s.restraints.add_selectively_active_collection(res_groups,int(len(res_groups)-10))


    #
    # Secondary Structure
    #
    ss_scaler = s.restraints.create_scaler('constant')
    ss_rests = parse.get_secondary_structure_restraints(filename='%s-ss.dat' %(sys.argv[1]), system=s,ramp=LinearRamp(0,100,0,1), scaler=ss_scaler,
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
    remd_runner = master_runner.MasterReplicaExchangeRunner(N_REPLICAS, max_steps=N_STEPS, ladder=l, adaptor=a)
    store.save_remd_runner(remd_runner)

    # create and store the communicator
    c = comm.MPICommunicator(s.n_atoms, N_REPLICAS)
    store.save_communicator(c)

    # create and save the initial states
    # create and save the initial states, initialize each replica with a different template
    states = [gen_state_templates(i,templates) for i in range(N_REPLICAS)]
    store.save_states(states, 0)

    # save data_store
    store.save_data_store()

    return s.n_atoms

# RUN THE SETUP
setup_system()
