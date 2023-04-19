#! /usr/bin/env python

import numpy
import matplotlib 
matplotlib.use('Agg')
from matplotlib import pyplot
import glob
import os
from scipy import stats

pwd = os.getcwd()
#print(pwd)
directories = pwd.split("/")
def make_plot():
    fig, axs = pyplot.subplots(3, 3, sharex='col',sharey='row',figsize=(24,18))
    pyplot.subplots_adjust(hspace=0.001,wspace=0.001)
    files=['peptide_TP/','peptide_NSD3/','peptide_JMJD6/']
    files2=['csp/','csp_talos/','csp_talos_noe/']
    for f in files:
        j=files.index(f)
        for g in files2:
        
            k = files2.index(g)
            with open(f+g+'Funnel_linkage_sieve_eps_new_1.5/summary','r') as fi:
                clusters = int(fi.readlines()[-1].split()[0])
            cluster_number = clusters + 1

            rmsd = numpy.loadtxt(f+g+'Funnel_linkage_sieve_eps_new_1.5/trajrmsd_2.out', skiprows=1)[:,1]
            rmsd[rmsd>30] = 30
            cluster = numpy.loadtxt(f+g+'Funnel_linkage_sieve_eps_new_1.5/frame_vs_cluster.txt', skiprows=1)[:,1]
            replica_index = numpy.zeros(len(rmsd))
            #replicas = glob.glob(f+'trajectory.*.dcd')
            #number_replicas = len(replicas)
            number_replicas=30
            chunks = int(float(len(rmsd))/number_replicas)
            for i in range(number_replicas):
                replica_index[chunks*i:] = i 

            numpy.savetxt(f+g+'Funnel_linkage_sieve_eps_new_1.5/chunks.txt',replica_index)

#there are cluster_number clusters 0..3350
            cluster_avg_rmsd = numpy.zeros(cluster_number)
            cluster_avg_replica = numpy.zeros(cluster_number)
            pop = numpy.zeros(cluster_number)
            color = numpy.chararray(cluster_number,itemsize=5)
   # pyplot.figure()
            for i in range(cluster_number):
                index_cluster = cluster == i
                cluster_avg_rmsd[i] = numpy.average(rmsd[index_cluster])
                cluster_avg_replica[i] = numpy.average(replica_index[index_cluster])
                pop[i] = index_cluster.sum()
        #print(i, cluster_avg_rmsd[i], cluster_avg_replica[i], pop[i])
                col = 'red'
                if cluster_avg_replica[i] < 15.:
                    col = 'blue'
                    if cluster_avg_rmsd[i] < 5.:
                        col = 'green'
                axs[j,k].scatter(cluster_avg_rmsd[i],cluster_avg_replica[i],s=(pop[i])/100.,c=col,alpha=1)
                axs[j,k].set_xticks([0,5,10,15,20,25])
                #axs[j,k].set_xticklabels([0,5,10,15,20,25],fontsize=12)
                axs[j,k].set_ylim([-1,33])
                axs[j,k].set_xlim([0,30])
                axs[j,k].set_yticks([0,10,20,30])
                #axs[j,k].set_title(f[:-1],x=0.5,y=0.9,fontsize=12)
    axs[0,0].set_yticklabels([0,10,20,30],fontsize=25)
    axs[1,0].set_yticklabels([0,10,20,30],fontsize=25)
    axs[2,0].set_yticklabels([0,10,20,30],fontsize=25)
    axs[2,0].set_xticklabels([0,5,10,15,20,25],fontsize=25)
    axs[2,1].set_xticklabels([0,5,10,15,20,25],fontsize=25)
    axs[2,2].set_xticklabels([0,5,10,15,20,25],fontsize=25)
    axs[2,1].set_xlabel('RMSD (Å)', fontsize=25)
    axs[1,0].set_ylabel('Replica Index',fontsize=25)
    #pyplot.show()
    axs[0,0].set_title('CSP',fontsize=25)
    axs[0,1].set_title('CSP+TALOS',fontsize=25)
    axs[0,2].set_title('CSP+TALOS+NOE',fontsize=25)
    axs[0,2].yaxis.set_label_position("right")
    axs[1,2].yaxis.set_label_position("right")
    axs[2,2].yaxis.set_label_position("right")
    axs[0,2].set_ylabel('BRD3-TP',fontsize=25)
    axs[1,2].set_ylabel('BRD3-NSD3',fontsize=25)
    axs[2,2].set_ylabel('BRD3-JMJD6',fontsize=25)
    pyplot.savefig('funnel_plot_new2.png')
    pyplot.savefig('funnel_plot_new2.pdf')

make_plot()









