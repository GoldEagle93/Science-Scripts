import os
import mdtraj as md
import datetime
import numpy as np
import matplotlib.pyplot as plt

trajname = 'md.nc'
topname = 'topol.prmtop'
stride = 50
cutoff = 2.0
frames = range(0, 1000)
trajectory = md.load(trajname, top=topname)
t1 = datetime.datetime.now()
topology = trajectory.topology
waters = topology.select("name O and (resname =~ 'HOH*')")
reference = topology.select("resname F86 and name O2")
waters_resids = []
for i in waters:
    waters_resids.append(topology.atom(i).residue.index)
print(len(topology.select("name O and (resname =~ 'HOH*')")) == len(topology.select("name O and water")))
print('There are '+str(len(waters)+' water molecules')

#Make a list of water molecules that are nearby

inshell = []
visiting_waters = []

for water in range(len(waters)):
    time = 0
    watertimes = []
    for frame in frames:
        if calculate_distance(waters[water], reference, frame) < cutoff and calculate_distance(waters[water], reference, frame-1) < cutoff:
            time += stride
        elif time != 0:
            watertimes.append(time)
            time = 0
    if time != 0:
        watertimes.append(time)
        visiting_waters.append(waters[water])
        time = 0
        inshell.append(watertimes)

#See if any water moleucles resided for the whole time
constant_waters = []
for i in range(len(inshell)):
    if inshell[i] > [(max(frames)-min(frames))*stride]:
        constant_waters.append(visiting_waters[i])
        
constant_waters_resids = []
for i in constant_waters:
    constant_waters_resids.append(topology.atom(i).residue.index)
