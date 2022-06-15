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
    
    tbin = 0.001*(18.6e8*U_T)      # Time bin for orbit integration.
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
    
    def plot(self,n):       #Plots the n-element inside the frames vector of particles through time
        figdir='figures/'
            
            #getting data
        x = [i.xpos/U_DIST for i in self.frames[n]]
        y = [i.ypos/U_DIST for i in self.frames[n]]
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
                parti = self.part[i]
                [newpartx,newparty,newpartz,newpartvx,newpartvy,newpartvz] = self.rk4(parti.xpos,parti.ypos,parti.zpos,parti.vx,parti.vy,parti.vz)
                self.part[i] = particle(parti.m,newpartx,newparty,newpartz,newpartvx,newpartvy,newpartvz)
            return self.part
   
    def mass_r(self,x,y,z):
        r = np.sqrt(x**2 + y**2 + z**2)
        mass = 4*np.pi*RHOo*(Rs**3)*(np.log((Rs+r)/Rs)-(r/(Rs+r)))
        return mass,r

    def acc(self,x,y,z,vx,vy,vz):
        m,r = self.mass_r(x,y,z)
        acc = [vx,vy,vz,-(G*m/(r**3))*x,-(G*m/(r**3))*y,-(G*m/(r**3))*z]
        return acc
    
    def rk4(self,x,y,z,vx,vy,vz):
        h = self.tbin    
        
        k1=self.acc(x,y,z,vx,vy,vz)
        k2=self.acc(x+k1[0]*h/2, y+k1[1]*h/2, z+k1[2]*h/2, vx+k1[3]*h/2, vy+k1[4]*h/2, vz+k1[5]*h/2)
        k3=self.acc(x+k2[0]*h/2, y+k2[1]*h/2, z+k2[2]*h/2, vx+k2[3]*h/2, vy+k2[4]*h/2, vz+k2[5]*h/2)
        k4=self.acc(x+k3[0], y+k3[1], z+k3[2], vx+k3[3], vy+k3[4], vz+k3[5])

        newx = x + (h/6.)*(k1[0] + 2*k2[0] + 2*k3[0] + k4[0])
        newy = y + (h/6.)*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1])
        newz = z + (h/6.)*(k1[2] + 2*k2[2] + 2*k3[2] + k4[2])
        newvx = vx + (h/6.)*(k1[3] + 2*k2[3] + 2*k3[3] + k4[3])
        newvy = vy + (h/6.)*(k1[4] + 2*k2[4] + 2*k3[4] + k4[4])
        newvz = vz + (h/6.)*(k1[5] + 2*k2[5] + 2*k3[5] + k4[5])

        return newx,newy,newz,newvx,newvy,newvz 
   
    def whole(self,t):
        totalt = t
        i = 0
        while i < totalt:
            self.frames = np.append(self.frames,self.part)
            self.part = self.nextframe()
            i+=self.tbin
