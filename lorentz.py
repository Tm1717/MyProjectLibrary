# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 23:34:17 2022

@author: tanzi
"""

 


""" Section 1:  Start by importing relevant libraries
#---------------------------------------------------------------------------"""
import numpy as np
import matplotlib.pyplot as plt
from  matplotlib import animation as animation
import time as time

""" Section 2:  Define functions for the script
#---------------------------------------------------------------------------"""



def B_bottle(M, r):
    #M is Dipole Moment (3d array)
    #R vector will give 3 coordinates xyz.
    
    r[2] -= 10
    rMag = np.dot(r,r)**(0.5)
    B1top = 3.0 * r * np.dot(M,r) / (0.00000001+rMag**5)
    B2bot = -1.0 * M / (0.00000001+rMag**3)
    
    Btop = mu0/(4.0*np.pi) * (B1top + B2bot)
    
    r[2] += 20
    rMag = np.dot(r,r)**(0.5)
    B1top = 3.0 * r * np.dot(M,r) / (0.0000001+rMag**5)
    B2bot = -1.0 * M / (0.00000001+rMag**3)
    Bbot = mu0/(4.0*np.pi) * (B1top + B2bot)
    
    r[2] -= 10
    
    return Btop + Bbot


""" Section 3a:  Main body - Define Constants
#---------------------------------------------------------------------------"""
# Define physical constants
m_p = 1.67E-27       # mass of proton: kg
qe = 1.602E-19       # charge of proton: C
mu0 = 4*np.pi*1e-7   # permeability of free space

# alpha particle 
mass = 4.0*m_p
q = 2.0*qe
Q_M = q/mass

#magnetic dipole moment
M = 10000.0 * np.array([0.0, 0.0, 1.0])

""" Section 3b:  Main body - Set up Euler's Method
#---------------------------------------------------------------------------"""
# trajectory setup
dt = 1.0E-5                       # time step
endTime = 2.0                   # run from zero for this many seconds
T = np.arange(0.0, endTime, dt)   # array of times for trajectory

# Set up arrays for position and velocity
rp = np.zeros( (T.size,3) )
vp = np.zeros( (T.size,3) )

v0 = 100
rp[0,:] = np.array([2.0, 1.0, 0.0])
vp[0,:] = np.array([v0, 0.0, v0])

clockStart = time.time()
# Euler time steps
for i in range(0, len(T)-1):
    # Acceleration Ap due to dipole field
    accelB = Q_M*np.cross( vp[i,:] , B_bottle(M,rp[i,:]) )
    
    # update position and velocity according to Euler method
    vtemp1 = vp[i,:] + accelB*dt
    rtemp1 = rp[i,:] + vp[i,:]*dt
    
    # Modified Euler:  WTFysics is going on here?
    accelE = Q_M*np.cross( vp[i,:] , B_bottle(M,rtemp1) )
    vtemp2 = vp[i,:] + accelE*dt
    vp[i+1,:] = vp[i,:] + 0.5*(accelB+accelE)*dt
    rp[i+1,:] = rp[i,:] + 0.5*(vtemp1+vtemp2)*dt

print("duration of for loop: %5f" % (time.time()-clockStart) )

""" Section 4 plotting
#---------------------------------------------------------------------------"""
fig = plt.figure(figsize=(20,10))
ax1 = fig.add_subplot(231)
ax2 = fig.add_subplot(232, projection='3d')
ax3 = fig.add_subplot(233)
ax4 = fig.add_subplot(234)
ax5 = fig.add_subplot(235)
ax6 = fig.add_subplot(236)
ax2.view_init(azim=-110.0, elev=10.0)
ax2.dist = 10

# Set up plot aesthetics
ax1.set_xlabel('x position (m)')
ax1.set_ylabel('z position (m)')
ax1.set_title('trajectory of alpha particle \n in a magnetic field')

ax1.plot(rp[:,0],rp[:,2],'b')

# axis labels for 3D plot
ax2.set_xlabel('X position',fontsize=14)
ax2.set_ylabel('Y position',fontsize=14)
ax2.set_zlabel('Z position',fontsize=14)
ax2.plot(rp[:,0], rp[:,1], rp[:,2], 'r')

ax3.set_title('X, Y, Z coordinates vs time')
ax3.plot(T,rp[:,0],'r', label = 'X-coordinate')
ax3.plot(T,rp[:,1],'b', label = 'Y-coordinate')
ax3.plot(T,rp[:,2],'g', label = 'Z-coordinate')
ax3.legend()


ax4.plot(T,rp[:,0],'r')
ax4.set_xlabel('Time (s)')
ax4.set_ylabel('X-position (m)')
ax4.set_title('X-coordinate vs time graph')

ax5.plot(T,rp[:,1],'b')
ax5.set_xlabel('Time (s)')
ax5.set_ylabel('Y-position (m)')
ax5.set_title('Y-coordinate vs time graph')

ax6.plot(T,rp[:,2],'g')
ax6.set_xlabel('Time (s)')
ax6.set_ylabel('Z-position (m)')
ax6.set_title('Z-coordinate vs time graph')

# This will overlay B-field lines on graph:  Used when switching from
# uniform field to dipole field
if True:
    # Set up the meshgrid for 2D plotting
    x = np.arange(-20,20,0.1)
    z = np.arange(-20,20,0.1)
    X,Z = np.meshgrid(x,z)
    ilen, klen = np.shape(X)
    
    # 3D array of magnetic field values - 3rd dimension for x,y,z components of the
    # magnetic field vector
    Bfield = np.zeros((ilen,klen,3))
    
    for i in range(ilen):
        for k in range(klen):
            R = np.array( [X[i,k], 0.0, Z[i,k]] )
            Bfield[i,k,:] = B_bottle(M, R)
    
    ax1.streamplot(X, Z, Bfield[:,:,0], Bfield[:,:,2], density=0.5, color='k')
    
    # streamline plot for magnetic field lines
plt.savefig('bottle.png', dpi=300)    


