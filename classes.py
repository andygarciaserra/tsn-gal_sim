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
from mpl_toolkits.mplot3d import Axes3D


#General variables:
LOAD_DIR = 'ics/'
SAVEDIR = 'figures/'

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
    oldpart = np.array([])
    acc_int = np.array([])      # Acceleration tensor.
    tbin = 0.001*(18.6e8*U_T)   # Time bin for orbit integration.
    method = ''                 # Integration method ('Euler' for Euler and 'RK4' for Runge-Kutta)
    sys = ''                    # System type ('sun' for sun-SagA / 'loadN' to load ics from loadN.txt)
    inter = ''                  # Boolean for grav interaction between stars. ('yes'/'no')
    softlength = 1*U_DIST       # System softening length
    
    def __init__(self,method,sys,inter):          #Definition example: u = Universe('Euler','load10','yes')  
        self.sys = sys
        self.method = method
        self.inter = inter
        
        # Creating the actual Universe:
        if self.sys=='sun':
            self.add(particle(1*U_MASS, 8*U_DIST,0,0, 0,127*U_VEL,0))
        if self.sys=='load10':
            self.load(10)
        if self.sys=='load100':
            self.load(100)
        if self.sys=='load1000':
            self.load(1000)  

    def load(self,N):
        data = np.loadtxt(LOAD_DIR+'disk'+str(N)+'.txt')
        for i in range(len(data)):
            self.add(particle(data[i,0]*U_MASS,data[i,1]*U_DIST,data[i,2]*U_DIST,data[i,3]*U_DIST,data[i,4]*U_VEL,data[i,5]*U_VEL,data[i,6]*U_VEL))

    def show(self):
        print(self.part)
    
    def add(self,newpart):
        self.part = np.append(self.part,newpart)
    
    def plot_trace(self):       #Plots the n-element inside the frames vector of particles through time
        
        fig = plt.figure(figsize=(10,5))
        ax1 = fig.add_subplot(121, projection='3d')
        ax2 = fig.add_subplot(122)
        
        for i in range(len(self.part)):
            x = [self.frames[j,i].xpos/U_DIST for j in range(len(self.frames[:,i]))]
            y = [self.frames[j,i].ypos/U_DIST for j in range(len(self.frames[:,i]))]
            z = [self.frames[j,i].zpos/U_DIST for j in range(len(self.frames[:,i]))]
                #plotting
            ax1.scatter3D(0,0,0,c='grey',marker='o',alpha=0.5)
            ax1.plot3D(x,y,z,lw=0.5)
            ax2.plot(x,y,ms=0.5)

        #formatting plot
        ax1.set_xlabel('x (kpc)')
        ax1.set_ylabel('y (kpc)')
        ax1.set_zlabel('z (kpc)')
        ax2.set_xlabel('x (kpc)')
        ax2.set_ylabel('y (kpc)')
        ax2.yaxis.set_label_position('right')
        ax2.yaxis.tick_right()
        fig.tight_layout(pad=3.0)
        fig.subplots_adjust(hspace=4)
        plt.savefig(SAVEDIR+str(self.sys)+'_'+str(self.method)+'_'+str(self.inter)+'int'+'.png',dpi=200)
        #plt.show()

    def nextframe(self):
        if (self.method=='Euler'):
            for i in range(len(self.oldpart)):
                mass = self.oldpart[i].m
                newpartx = self.oldpart[i].xpos + self.oldpart[i].vx * self.tbin
                newparty = self.oldpart[i].ypos + self.oldpart[i].vy * self.tbin
                newpartz = self.oldpart[i].zpos + self.oldpart[i].vz * self.tbin
                mass_r,r = mass_r(newpartx,newparty,newpartz)
                newpartvx = self.oldpart[i].vx - (G*mass_r/(r**3)) * newpartx * self.tbin
                newpartvy = self.oldpart[i].vy - (G*mass_r/(r**3)) * newparty * self.tbin
                newpartvz = self.oldpart[i].vz - (G*mass_r/(r**3)) * newpartz * self.tbin
                self.add(particle(mass,newpartx,newparty,newpartz,newpartvx,newpartvy,newpartvz))
                
        if (self.method=='RK4'):
            for i in range(len(self.oldpart)):
                [newpartx,newparty,newpartz,newpartvx,newpartvy,newpartvz] = self.rk4(i)
                self.add(particle(self.oldpart[i].m,newpartx,newparty,newpartz,newpartvx,newpartvy,newpartvz))
   
    def mass_r(self,x,y,z):
        r = np.sqrt(x**2 + y**2 + z**2)
        mass = 4*np.pi*RHOo*(Rs**3)*(np.log((Rs+r)/Rs)-(r/(Rs+r)))
        return mass,r

    def acc(self,x,y,z,vx,vy,vz,i):
        if self.inter=='no':
            m,r = self.mass_r(x,y,z)
            acc = [vx,vy,vz,-(G*m/(r**3))*x,-(G*m/(r**3))*y,-(G*m/(r**3))*z]
            return acc
        if self.inter=='yes':
            m,r = self.mass_r(x,y,z)
            accx = 0.
            accy = 0.
            accz = 0.
            for j in range(len(self.oldpart)):
                if j!=i:
                    dx = self.oldpart[j].xpos-x 
                    dy = self.oldpart[j].ypos-y
                    dz = self.oldpart[j].zpos-z
                    intr = np.sqrt(dx**2+dy**2+dz**2+self.softlength**2)
                    accx += (self.oldpart[j].m)*(G*dx*(intr**(-3)))
                    accy += (self.oldpart[j].m)*(G*dy*(intr**(-3)))
                    accz += (self.oldpart[j].m)*(G*dz*(intr**(-3)))
            acc = [vx,vy,vz,accx-(G*m/(r**3))*x,accy-(G*m/(r**3))*y,accz-(G*m/(r**3))*z]
            return acc
            

    def rk4(self,i):
        h = self.tbin
        x = self.oldpart[i].xpos ;  y = self.oldpart[i].ypos ;  z = self.oldpart[i].zpos
        vx = self.oldpart[i].vx  ;  vy = self.oldpart[i].vy  ;  vz = self.oldpart[i].vz
        k1=self.acc(x,y,z,vx,vy,vz,i)
        k2=self.acc(x+k1[0]*h/2, y+k1[1]*h/2, z+k1[2]*h/2, vx+k1[3]*h/2, vy+k1[4]*h/2, vz+k1[5]*h/2,i)
        k3=self.acc(x+k2[0]*h/2, y+k2[1]*h/2, z+k2[2]*h/2, vx+k2[3]*h/2, vy+k2[4]*h/2, vz+k2[5]*h/2,i)
        k4=self.acc(x+k3[0]*h, y+k3[1]*h, z+k3[2]*h, vx+k3[3]*h, vy+k3[4]*h, vz+k3[5]*h,i)

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
        self.oldpart = self.part
        self.frames = np.append(self.frames,self.oldpart)
        self.part = np.array([])
        self.nextframe()
        i+=self.tbin
        while i < totalt:
            self.oldpart = self.part
            self.frames = np.vstack((self.frames,self.oldpart))
            self.part = np.array([])
            self.nextframe()
            i+=self.tbin
