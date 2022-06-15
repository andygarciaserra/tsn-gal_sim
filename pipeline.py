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

u = universe('Euler','load10')
u.whole(18.6e8*U_T)

#Plotting traces:
frames = u.frames
u.plot_trace()
