import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import os
from os.path import exists
# matplotlib.use('Agg')
import glob

source = '/Users/HQ/orange/Reza/TF/NAR/' #14-two-sided-dummy-contacts/'
source = '/orange/alberto.perezant/Reza/TF/NAR'
systems = sorted([i for i in os.listdir(source) if '.' not in i])
systems = ['1a74', '1azp', '1by4', '1bgb',  '2b0d', '1cdw', '1dh3', '1jj4', '1r4o', '1r4r', '1zme', '2r1j', '3cro', '1ysa', '2dgc']
bold = ['#fcff5d', '#7dfc00', '#0ec434', '#228c68', '#8ad8e8', '#235b54', '#29bdab', '#3998f5', '#37294f', '#277da7', '#3750db', '#f22020', '#991919', '#ffcba5', '#e68f66', '#c56133', '#96341c', '#632819', '#ffc413', '#f47a22', '#2f2aa0', '#b732cc', '#772b9d', '#f07cab', '#d30b94', '#edeff3', '#c3a5b4', '#946aa2', '#5d4c86', '#201923']
uf = ['#0021A5', '#F2A900', '#FA4616', '#22884C']


os.chdir(source)

print('starting')

fig, axs = plt.subplots(3, 5, sharex='col',sharey='row',figsize=(24,18))
for i in range(15):
    row = int(i/5)
    column = i % 5
    if exists('%s/14-two-sided-dummy-contacts/contact-clusters-all-traj-pd/DNA-PROT/summary'%(systems[i])) and exists('%s/14-two-sided-dummy-contacts/contact-clusters-all-traj/DNA-PROT/summary'%(systems[i])) and exists('%s/14-two-sided-dummy-contacts/contact-clusters-all-traj-li/DNA/summary'%(systems[i])):
        print('Now doing %s\n'%(systems[i]))
        # rmsd    = pd.read_csv('%s/14-two-sided-dummy-contacts/contact-clusters-all-traj/DNA-PROT/summary'%(systems[i]),    delim_whitespace=True)
        # rmsd_pd = pd.read_csv('%s/14-two-sided-dummy-contacts/contact-clusters-all-traj-pd/DNA-PROT/summary'%(systems[i]), delim_whitespace=True)
        # rmsd_li = pd.read_csv('%s/14-two-sided-dummy-contacts/contact-clusters-all-traj-li/DNA/summary'%(systems[i]), delim_whitespace=True)
        clusters    = int(open('%s/14-two-sided-dummy-contacts/contact-clusters-all-traj/DNA-PROT/summary'%(systems[i]),   'r').readlines()[-1].split()[0])
        clusters_pd = int(open('%s/14-two-sided-dummy-contacts/contact-clusters-all-traj-pd/DNA-PROT/summary'%(systems[i]),'r').readlines()[-1].split()[0])
        clusters_li = int(open('%s/14-two-sided-dummy-contacts/contact-clusters-all-traj-li/DNA/summary'%(systems[i]),     'r').readlines()[-1].split()[0])
        cluster_number = clusters + 1
        
        cluster = np.loadtxt('%s/14-two-sided-dummy-contacts/contact-clusters-all-traj/DNA-PROT/frame_vs_cluster.txt'%(systems[i]), skiprows=1)[:,1]

        rmsd = np.loadtxt('%s/14-two-sided-dummy-contacts/contact-clusters-all-traj/trajrmsd.dat'%(systems[i]), skiprows=1)[:,1]
        replica_index = np.zeros(len(rmsd))
        #  summary_pd[-1].split()[0])
        # summary_pd = open(f+g+'Funnel_linkage_sieve_eps_new_1.5/summary','r').readlines()

        number_replicas=30
        chunks = int(float(len(rmsd))/number_replicas)
        for i in range(number_replicas):
            replica_index[chunks*i:] = i 
        # np.savetxt('%s/14-two-sided-dummy-contacts/contact-clusters-all-traj/chunks.dat'%(systems[i]), replica_index)

        cluster_avg_rmsd = np.zeros(cluster_number)
        cluster_avg_replica = np.zeros(cluster_number)
        pop = np.zeros(cluster_number)
        color = np.chararray(cluster_number,itemsize=5)
# plt.figure()
        for i in range(cluster_number):
            index_cluster = cluster == i
            cluster_avg_rmsd[i] = np.average(rmsd[index_cluster])
            cluster_avg_replica[i] = np.average(replica_index[index_cluster])
            pop[i] = index_cluster.sum()
    #print(i, cluster_avg_rmsd[i], cluster_avg_replica[i], pop[i])
            col = 'red'
            if cluster_avg_replica[i] < 15.:
                col = 'blue'
                if cluster_avg_rmsd[i] < 5.:
                    col = 'green'
            axs[row,column].scatter(cluster_avg_rmsd[i],cluster_avg_replica[i],s=(pop[i])/100.,c=col,alpha=1)
#             axs[j,k].set_xticks([0,5,10,15,20,25])
#             #axs[j,k].set_xticklabels([0,5,10,15,20,25],fontsize=12)
#             axs[j,k].set_ylim([-1,33])
#             axs[j,k].set_xlim([0,30])
#             axs[j,k].set_yticks([0,10,20,30])
#             #axs[j,k].set_title(f[:-1],x=0.5,y=0.9,fontsize=12)
# axs[0,0].set_yticklabels([0,10,20,30],fontsize=25)
# axs[1,0].set_yticklabels([0,10,20,30],fontsize=25)
# axs[2,0].set_yticklabels([0,10,20,30],fontsize=25)
# axs[2,0].set_xticklabels([0,5,10,15,20,25],fontsize=25)
# axs[2,1].set_xticklabels([0,5,10,15,20,25],fontsize=25)
# axs[2,2].set_xticklabels([0,5,10,15,20,25],fontsize=25)
# axs[2,1].set_xlabel('RMSD (Å)', fontsize=25)
# axs[1,0].set_ylabel('Replica Index',fontsize=25)
# #plt.show()
# axs[0,0].set_title('CSP',fontsize=25)
# axs[0,1].set_title('CSP+TALOS',fontsize=25)
# axs[0,2].set_title('CSP+TALOS+NOE',fontsize=25)
# axs[0,2].yaxis.set_label_position("right")
# axs[1,2].yaxis.set_label_position("right")
# axs[2,2].yaxis.set_label_position("right")
# axs[0,2].set_ylabel('BRD3-TP',fontsize=25)
# axs[1,2].set_ylabel('BRD3-NSD3',fontsize=25)
# axs[2,2].set_ylabel('BRD3-JMJD6',fontsize=25)
print('saving figure')
plt.savefig('funnel2.png', bbox_inches='tight', dpi=300)
# plt.savefig('funnel_plot_new2.pdf')
# plt.show()

