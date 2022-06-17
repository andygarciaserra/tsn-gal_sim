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




#EXERCICE 3:

    #Arrays for diferent Nparts, methods and plot styles:
n = [10,100]
sys = ['load10','load100']
c = ['blue','tab:orange']
ticks = [1,10,100]
j=1

#    #Computing different universes for each case:
#t = np.array([])
#nt = np.array([])
#for i in range(len(sys)):
#    ti = time.time()
#    u = universe('RK4',sys[i],'yes')
#    u.whole(18.6e8*U_T)
#    u.plot_trace()
#    tf = time.time()
#    t = np.append(t, (tf-ti))
#    nt = np.append(nt, n[i])
#np.savetxt(OUTDIR+'timevsN_'+'RK4'+'_yesint.txt',np.transpose([t,nt]),header='Time(s)\tNparts')

    #Plotting computing time vs Npart (no interaction between stars):
data = np.array([[],[]])
plt.figure(dpi=150)
data = np.loadtxt(OUTDIR+'timevsN_'+'RK4'+'_yesint.txt', skiprows=1)
plt.plot(data[:,1],data[:,0],color=c[j],ls='-',lw=2,label='RK4')
plt.plot(data[:,1],data[:,0],color=c[j],marker='.',ms=5)
plt.xticks(ticks)
plt.xlabel('Nparts')
plt.ylabel('t (s)')
plt.legend()
plt.savefig(FIGDIR+'timevsN_yesint.png')
plt.show()
