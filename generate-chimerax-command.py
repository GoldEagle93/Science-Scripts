import os
import sys
import pyperclip

source = '/Users/HQ/orange/Reza/TF/NAR/' #1.14-two-sided-dummy-contacts/'
# systems = sorted([i for i in os.listdir(source) if '.' not in i])
# systems = ['1a74', '1azp', '1by4', '1bgb',  '2b0d', '1cdw', '1dh3', '1jj4', '1r4o', '1r4r', '1zme', '2dgc', '2r1j', '3cro', '1ysa', '2dgc']
systems = [sys.argv[1]]
bold = ['#fcff5d', '#7dfc00', '#0ec434', '#228c68', '#8ad8e8', '#235b54', '#29bdab', '#3998f5', '#37294f', '#277da7', '#3750db', '#f22020', '#991919', '#ffcba5', '#e68f66', '#c56133', '#96341c', '#632819', '#ffc413', '#f47a22', '#2f2aa0', '#b732cc', '#772b9d', '#f07cab', '#d30b94', '#edeff3', '#c3a5b4', '#946aa2', '#5d4c86', '#201923']
uf = ['#0021A5', '#F2A900', '#FA4616', '#22884C']

os.chdir(source)
for system in systems:
    ref = open(system+'/14-two-sided-dummy-contacts/'+'ref.pdb', 'r').readlines()
    temp = open(system+'/14-two-sided-dummy-contacts/TEMPLATES/'+system+'-BDNA.pdb', 'r').readlines()
    rmsd = open(system+'/14-two-sided-dummy-contacts/rmsd00.cpptraj', 'r').readlines()
    mask = ','.join([str(z) for z in (sorted([int(j) for j in [i for i in rmsd if 'rms traj00 ref [ref] :' in i][0].split(':')[1].split('@')[0].split(',')]))])
    chimera_command = "select #1; color sel #FA4616; select #2; color sel #0021A5; align #2:%s@N,C,CA,O5',C4',C1' toAtoms #1:%s@N,C,CA,O5',C4',C1'; select clear"%(mask, mask)
    pyperclip.copy(chimera_command)
        

