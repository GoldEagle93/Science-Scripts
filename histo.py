#! /usr/bin/env python

import matplotlib
from matplotlib import pyplot

import numpy

data = numpy.loadtxt('/orange/alberto.perezant/Reza/TF/NAR/1azp/12-dummy-contacts-70%-trust-1/rmsd.dat',skiprows=1)
data = data[:,1:]
bold = ['#fcff5d', '#7dfc00', '#0ec434', '#228c68', '#8ad8e8', '#235b54', '#29bdab', '#3998f5', '#37294f', '#277da7', '#3750db', '#f22020', '#991919', '#ffcba5', '#e68f66', '#c56133', '#96341c', '#632819', '#ffc413', '#f47a22', '#2f2aa0', '#b732cc', '#772b9d', '#f07cab', '#d30b94', '#edeff3', '#c3a5b4', '#946aa2', '#5d4c86', '#201923']
fig,ax = pyplot.subplots(nrows=6,ncols=5,sharex=True,sharey=True)
for i in range(30):
    remainder = i % 5
    row = int(i/5)
    ax[row,remainder].hist(data[:,i],color=bold[i],range=(0,20))
    #pyplot.violinplot(data[:,i],positions=i,color=bold[i])
pyplot.show()
