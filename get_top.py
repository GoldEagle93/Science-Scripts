#! /usr/bin/env python
import pickle as cPickle
x = cPickle.load(open('Data/system.dat','rb'))
f = open('topol.top', 'w')
f.write(x.top_string)
g = open('topol.crd', 'w')
g.write(x._mdcrd_string)
