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

u = universe('RK4','sun')
u.whole(18e8*U_T)

#Plotting traces:
v = u.frames
for i in range(len(v)):
    u.plot(i)
