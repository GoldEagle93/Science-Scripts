### The required libraries and packages ###
import networkx as nx
import requests
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
# from matplotlib import cm
from nxviz.plots import CircosPlot
import random

seq=open('sequence.dat','r')
seq=seq.readlines()
for i in seq:
    L = len(i)
print(L)

def plot_type(wt,output,output2):
    G=nx.Graph(name='Protein Interaction Graph')
    for a in range(1,L):
        if a==1:
            G.add_node(a,group=1)
        if 1<a<int(L/4):
            G.add_node(a,group=int(L/8))
        if int(L/4)<=a<int(L/2):
            G.add_node(a,group=int(3*L/8))
        if int(L/2)<=a<int(3*L/4):
            G.add_node(a,group=int(5*L/8))
        if int(3*L/4)<=a<L-1:
            G.add_node(a,group=int(7*L/8))
        if a==L-1:
            G.add_node(a,group=L-1)

    good=[]
    f=open('good.dat','r').readlines()
    for i in f:
        i=i.strip()
        i=i.split()
        good.append(i)
    
    #    G.add_edge(int(i[0]),int(i[2]),color='r',weight=0,edgeprops={'alpha':1})

    f=open('4_satisfied.dat','r').readlines()
    for i in f:
        i=i.strip()
        i=i.split()
        if i in good:
            G.add_edge(int(i[0]),int(i[2]),color='r',weight=wt,edgeprops={'alpha':1})
        if i not in good:
            G.add_edge(int(i[0]),int(i[2]),color='pink',weight=wt,edgeprops={'alpha':1})
    


#c = CircosPlot(G,edge_color='color',edge_width='weight',alpha=1,node_label_layout="rotation").draw()
    c = CircosPlot(G,node_grouping='group',node_order='group',group_label_position="middle",edge_color='color',edge_width='weight',alpha=1,node_label_layout="rotation",fontsize=18).draw()
#c = CircosPlot(G,node_grouping='group',node_order='group',edge_color='color',edge_width='weight')
#annotate.circos_group(G, group_by="group")

    plt.savefig(output)
    plt.savefig(output2)
#plot_type(2,'ambiguous_new.pdf','ambiguous_new.png')
#plot_type(0,'true_new.pdf','true_new.png')
plot_type(5,'satisfied_new.pdf','satisfied_new.png')
