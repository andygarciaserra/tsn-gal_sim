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
n = [1,10,100,1000]
sys = ['sun','load10','load100','load1000']
method = ['Euler','RK4']
c = ['blue','tab:orange']
ticks = [1,10,100,1000]

u = universe('RK4','load10','no')
u.whole(18.6e8*U_T)
u.plot_trace()
