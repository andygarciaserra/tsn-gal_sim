#IMPORTING PACKAGES:
import os
import random
import numpy as np
import matplotlib.pyplot as plt
from scipy import pi
import h5py
import sys
import random
from matplotlib.offsetbox import TextArea, DrawingArea, OffsetImage, AnnotationBbox
import matplotlib.pyplot as plt
import matplotlib.image as mpimg



#General variables:
LOAD_DIR = './'

#Natural to cgs convertion units:
U_T = 3600*24*365 #yrs to s
U_VEL = 1e5 #[km/s] to [cm/s]
U_DIST = 3.085678e21 #[kpc] to [cm]
U_MASS = 1.989e33 #[Msol] to [g]
G = 6.67e-8 #cgs units

#NFW parameters:
RHOo = 5932371.*U_MASS*(U_DIST**(-3))
Rs = 20.*U_DIST

TIMESTEP = 9.8e8*U_T

# Particle class:
class particle:
    m = 0
    xpos = 0
    ypos = 0
    zpos = 0
    vx = 0
    vy = 0
    vz = 0

    def __init__(self,mass,x,y,z,vx,vy,vz):
        self.m = mass
        self.xpos = x
        self.ypos = y
        self.zpos = z
        self.vx = vx
        self.vy = vy
        self.vz = vz

    def pos(self):
        return print('[x,y,z]= '+str(np.array([self.x,self.y,self.z])))
    
    def vel(self):
        return print('[vx,vy,vz]= '+str(np.array([self.vx,self.vy,self.vz]))) 



# Universe class:
class universe:
    part = np.array([])
    frames = np.array([])
    
    tbin = 0.001*(18e8*365*24*3600)      # Time bin for orbit integration.
    method = ''     # Integration method ('Euler' for Euler and 'RK4' for Runge-Kutta)
    sys = ''        # System type ('sun' for sin-SagA / 'loadN' to load ics from loadN.txt)

    def __init__(self,method,sys):          #Definition example: u = Universe(100, 10 , 10 ,  'RK4' ,'fig','binary')  
        self.sys = sys
        self.method = method

        # Creating the actual Universe:
        if self.sys=='sun':
            self.add(particle(1*U_MASS, 8*U_DIST,0,0, 0,127*U_VEL,0))
        if self.sys=='load10':
            data = self.load(10)
        if self.sys=='load100':
            data = self.load(100)
        if self.sys=='load1000':
            data = self.load(1000)
    
    def load(self,N):
        data = np.loadtxt(LOAD_DIR+'disk'+str(N)+'.txt')
        for i in range(len(data)):
            self.part = np.array([])
            self.add(particle(data[i,0],data[i,1],data[i,2],data[i,3],data[i,4],data[i,5],data[i,6]))
        return data

    def show(self):
        return self.part
    
    def add(self,newpart):
        self.part = np.append(self.part,newpart)
    
    def plot(self):
        figdir='figures/'
            
            #getting data
        x = [i.xpos/U_DIST for i in self.part]
        y = [i.ypos/U_DIST for i in self.part]
            #plotting
        plt.plot(0,0,'tab:orange',marker='o',ms=15,alpha=0.5)
        plt.plot(x,y,'ko',ms=3)
            #formatting plot
        plt.xlim((-20,20))
        plt.ylim((-20,20)) 
        plt.grid()
            #animating and showing
        plt.show(block=False)
        plt.pause(.0000001)
        plt.clf()

    def nextframe(self):
        if (self.method=='Euler'):
            for i in range(len(self.part)):
                mass = self.part[i].m
                newpartx = self.part[i].xpos + self.part[i].vx * self.tbin
                newparty = self.part[i].ypos + self.part[i].vy * self.tbin
                newpartz = self.part[i].zpos + self.part[i].vz * self.tbin
                r = np.sqrt(newpartx**2 + newparty**2 + newpartz**2)
                mass_r = 4*np.pi*RHOo*(Rs**3)*(np.log((Rs+r)/Rs)-(r/(Rs+r)))
                newpartvx = self.part[i].vx - (G*mass_r/(r**3)) * newpartx * self.tbin
                newpartvy = self.part[i].vy - (G*mass_r/(r**3)) * newparty * self.tbin
                newpartvz = self.part[i].vz - (G*mass_r/(r**3)) * newpartz * self.tbin
                self.part[i] = particle(mass,newpartx,newparty,newpartz,newpartvx,newpartvy,newpartvz)
            return self.part
                
        if (self.method=='RK4'):
            for i in range(len(self.part)):
                mass = self.part[i].m
                [newpartx,newparty,newpartz,newpartvx,newpartvy,newpartvz] = rk4(self.part.x,self.part.y,self.part.z,self.part.vx,self.part.vy,self.part.vz)
                self.part[i] = particle(mass,newpartx,newparty,newpartz,newpartvx,newpartvy,newpartvz)
            return self.part
    
    def rk4(x,y,z,vx,vy,vz):
        def f(x,y,z):
            r = np.sqrt(x**2 + y**2 + z**2)
            mass_r = 4*np.pi*RHOo*(Rs**3)*(np.log((Rs+r)/Rs)-(r/(Rs+r)))
            f = [-(G*mass_r/(r**3))*x,-(G*mass_r/(r**3))*y,-(G*mass_r/(r**3))*z]
            return f
        pos = [x,y,z]
        v = [vx,vy,vz]
        h = self.tbin
        i=0
        while i<3:
            f1=f(pos[0],pos[1],pos[2])
            k1[i] = h*f1[i]
            f2=f(pos[0],pos[1],pos[2])
            k2[i] = h*f2[i]
            f3=f(pos[0],pos[1],pos[2])
            k3[i] = 
            f4=f(pos[0],pos[1],pos[2])
            k4[i] = 
            i+=1


        return newpartx,newparty,newpartz,newpartvx,newpartvy,newpartvz 
   
    def whole(self,t):
        totalt = t
        i = 0
        self.show()
        while i < totalt:
            self.frames = np.append(self.frames,self.part)
            self.part = self.nextframe()
            self.plot()
            i+=self.tbin
