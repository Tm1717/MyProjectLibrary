# -*- coding: utf-8 -*-
"""
Created on Tue Oct 26 20:25:40 2021

@author: tanzi
"""

import matplotlib.pyplot as plt
import numpy as np
import random as rnd
from matplotlib import animation as am
from scipy.optimize import curve_fit

def energyFlow(eSolid,numExch):
    maxN = eSolid.size-1
    for i in range(numExch):
        gain = rnd.randint(0,maxN) 
        lose = rnd.randint(0,maxN) 
        while eSolid[lose] == 0:
            lose = rnd.randint(0,maxN)
        eSolid[gain] += 1
        eSolid[lose] -= 1

    return eSolid

def evolve1(*args):
    energyFlow(solid1,1000)
    imgPlot1.set_data( np.reshape(solid1, (L,L) ) )
    return [imgPlot1]    

def evolve2(*args):
    energyFlow(solid2,1000)
    imgPlot2.set_data( np.reshape(solid2, (L,B) ) )
    return [imgPlot2]    

def exponential(x, A, b): ##FIXME: 
    return A * np.exp(-b * x)

L = 50 
B = 25
N1 = L*L 
N2 = L*B
qavg1 = 15 
qavg2 = 30

solid1 = qavg1*np.ones(N1)
solid2 = qavg2*np.ones(N2)
solid3 = np.append(solid1, solid2)

 
    
animate = True

if animate:
    fig1 = plt.figure(figsize=(10,10))          
    
    img1 = np.reshape(solid1, (L,L))
    
    
    
    
    imgPlot1 = plt.imshow(img1, interpolation='none', vmin=0, vmax=50, cmap='coolwarm')
    plt.colorbar(imgPlot1)
    
    
    
    anim1 = am.FuncAnimation(fig1, evolve1, \
    						interval=10, frames=500, repeat=False, blit=True)
        
        
    fig2 = plt.figure(figsize=(10,10))
    img2 = np.reshape(solid2, (L,B))
    
    imgPlot2 = plt.imshow(img2, interpolation='none', vmin=0, vmax=50, cmap='coolwarm')
    plt.colorbar(imgPlot2)
    
    anim2 = am.FuncAnimation(fig2, evolve2, \
    						interval=10, frames=500, repeat=False, blit=True)    
    

else:
    solid1 = energyFlow(solid1, 500*solid1.size)
    
    energy1, counts1 = np.unique(solid1, return_counts=True)
    
    solid2 = energyFlow(solid2, 500*solid2.size)
    
    
    
    energy2, counts2 = np.unique(solid2, return_counts=True)
    
    
    
    solid3 = energyFlow(solid3, 1000*solid3.size) 
    
    energy3, counts3 = np.unique(solid3, return_counts=True)
    
    
    fig3 = plt.figure(figsize=(20,20))
    ax1 = fig3.add_subplot(321)
    ax2 = fig3.add_subplot(322)
    ax3 = fig3.add_subplot(323)
    ax4 = fig3.add_subplot(324)
    ax5 = fig3.add_subplot(326)
    ax6 = fig3.add_subplot(325)
    
    popt1, pcov1 = curve_fit(exponential, energy1, counts1)
    popt2, pcov2 = curve_fit(exponential, energy2, counts2)
    popt3, pcov3 = curve_fit(exponential, energy3, counts3)
    ax1.plot(energy1, exponential(energy1, *popt1), 'k-')
    ax1.plot(energy2, exponential(energy2, *popt2), 'k-')
    
    ax1.set_title('Energy Distribution for Solid 1 and Solid 2', fontsize = 10)
    ax1.set_ylabel(r'$E_{n}$')
    ax1.set_xlabel('n')
    ax1.plot(energy1,counts1, 'bx', markersize=5)  
    ax1.plot(energy2,counts2, 'ro', markersize=5)
    
    
    
    
    ax2.set_title('Log of Energy Distribution function for Solid 1 and Solid 2 ', fontsize = 10)

    ax2.plot(energy1, exponential(energy1, *popt1), 'k-')
    ax2.plot(energy2, exponential(energy2, *popt2), 'k-')
    ax2.set_ylabel(r'$logE_{n}$')
    ax2.set_xlabel('n')
    ax2.semilogy(energy1,counts1,'bx', markersize=5)
    ax2.semilogy(energy2,counts2,'ro', markersize=5)
    
    
    ax3.set_title('Energy Distribution for Combined solid', fontsize = 10)
    ax3.set_ylabel(r'$E_{n}$')
    ax3.set_xlabel('n')
    ax3.plot(energy3, exponential(energy3, *popt3), 'k-')
    
    ax3.plot(energy3,counts3, 'g+', markersize=7 )
    
    ax4.set_title('Log of Energy Distribution function for combined solid', fontsize = 10)
    ax4.set_ylabel(r'$logE_{n}$')
    ax4.set_xlabel('n')
    ax4.plot(energy3, exponential(energy3, *popt3), 'k-')
    ax4.semilogy(energy3,counts3, 'g+', markersize=7)
    
    ax5.set_title('Log of Energy Distribution for seperate solids and combined solid', fontsize = 10)
    ax5.set_ylabel(r'$logE_{n}$')
    ax5.set_xlabel('n')
    ax5.semilogy(energy3,counts3, 'g+', markersize=7)
    ax5.plot(energy3, exponential(energy3, *popt3), 'k--')
    ax5.semilogy(energy1,counts1,'bx', markersize=5)
    ax5.plot(energy1, exponential(energy1, *popt1), 'k-')
    ax5.semilogy(energy2,counts2,'ro', markersize=5)
    ax5.plot(energy2, exponential(energy2, *popt2), 'k-')
    
    
    ax6.set_title('Energy Distribution for Combined solid', fontsize = 10)
    ax6.set_ylabel(r'$E_{n}$')
    ax6.set_xlabel('n')
    ax6.plot(energy1,counts1, 'bx', markersize=5)  
    ax6.plot(energy1, exponential(energy1, *popt1), 'k-')
    ax6.plot(energy2,counts2, 'ro', markersize=5)
    ax6.plot(energy2, exponential(energy2, *popt2), 'k-')
    ax6.plot(energy3,counts3, 'g+', markersize=7 )
    ax6.plot(energy3, exponential(energy3, *popt3), 'k--')
    
plt.show 


plt.savefig('twoSolids.png', dpi=300)




