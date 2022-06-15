#IMPORTING PACKAGES:
import os
import random
import numpy as np
import matplotlib.pyplot as plt
from scipy import pi
import h5py
import sys
from classes import universe,particle
import random

U_T = 365*24*3600


#CREATING A TEST UNIVERSE OF N=1000:
sys = ['sun','load10','load100','load1000']
method = ['Euler','RK4']

for i in range(len(sys)):
    for j in range(len(method)):
        u = universe(method[j],sys[i])
        u.whole(18.6e8*U_T)
        u.plot_trace()
