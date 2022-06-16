#IMPORTING PACKAGES:
import os
import random
import numpy as np
import matplotlib.pyplot as plt
from scipy import pi
import h5py
import sys
from classes import universe,particle
import time



#PROGRAM VARIABLES:
OUTDIR = 'output/'
ICSDIR = 'ics/'
FIGDIR = 'figures/'
U_T = 365*24*3600




#EXERCICES 1 & 2:

    #Arrays for diferent Nparts, methods and plot styles:
n = [1,10,100,1000]
sys = ['sun','load10','load100','load1000']
method = ['Euler','RK4']
c = ['blue','tab:orange']
ticks = [1,10,100,1000]

    #Computing different universes for each case:
#for j in range(len(method)):
#    t = np.array([])
#    nt = np.array([])
#    for i in range(len(sys)):
#        ti = time.time()
#        u = universe(method[j],sys[i],'no')
#        u.whole(18.6e8*U_T)
#        u.plot_trace()
#        tf = time.time()
#        t = np.append(t, (tf-ti))
#        nt = np.append(nt, n[i])
#    np.savetxt(OUTDIR+'timevsN_'+method[j]+'.txt',np.transpose([t,nt]),header='Time(s) Nparts')

    #Computing time vs Npart relation:
data = np.array([[],[]])
plt.figure(dpi=150)
for j in range(len(method)):
    data = np.loadtxt(OUTDIR+'timevsN_'+method[j]+'.txt', skiprows=1)
    plt.plot(data[:,1],data[:,0],color=c[j],ls='-',lw=2,label=method[j])
    plt.plot(data[:,1],data[:,0],color=c[j],marker='.',ms=5)
plt.xticks(ticks)
plt.xlabel('Nparts')
plt.ylabel('t (s)')
plt.legend()
plt.savefig(FIGDIR+'timevsN.png')
plt.show()



#EXERCISE 3:


