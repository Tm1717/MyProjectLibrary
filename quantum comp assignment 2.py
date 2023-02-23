# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 21:36:12 2022

@author: tanzi
"""
'''Importing libraries'''

import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh_tridiagonal

'''Defining Constants'''

N = 500                     #Number of segments   
x = np.zeros(N)             #Empty array for x values    
x = np.linspace(-7,7,N)     #Getting N x values 
dx = np.diff(x)[0]          #Size of segment for x
V0 = 15                     #Finite Well Depth
r = np.zeros(N)             #Empty array for r values(For hydrogen)
r = np.linspace(0.01,25,N)  #Getting N r values 
dr = np.diff(r)[0]          #Size of segment for r



'''FInding Energies and Wavefunctions for Finite Well'''

main_diag = 2*np.ones(N)/dx**2 +V0-V0*((x>=-1)*(x<=1)).astype(float)                    #Defining main diagonal for my matrix as 2/dx^2 + Vx(Vx = potential at x)
off_diag =  -np.ones(N-1)/dx**2                                                         #Defining off diagonal for my matrix as -1/dx^2)

Es1, psis1 = eigh_tridiagonal(main_diag, off_diag, select='i', select_range=(0, 2))     #The built in function returns Eigenvalues (Energies) and Eigenvectors (Wavefunctions) of the matrix defined by the main and off diagonal


'''PLotting Wavefunctions for Finite Well'''

plt.figure(figsize=(24,12))
plt.plot(x, psis1.T[0], label = 'First Energy state')
plt.plot(x, psis1.T[1], label = 'Second Energy state')
plt.plot(x, psis1.T[2], label = 'Third Energy state')
plt.title('Finite Well Potential', fontsize = 30)
plt.ylabel(r'$\Psi$(z)', fontsize = 15)
plt.xlabel(r'z', fontsize = 15)
plt.legend()
plt.savefig('Finite Well Graphs.png')


'''FInding Energies and Wavefunctions for SHO Potential'''
            
main_diag = 2*np.ones(N)/dx**2 +(x**2)                                                  #Defining main diagonal for my matrix as 2/dx^2 + Vx(Vx = potential at x)
off_diag =  -np.ones(N-1)/dx**2                                                         #Defining off diagonal for my matrix as -1/dx^2)


Es2, psis2 = eigh_tridiagonal(main_diag, off_diag, select='i', select_range=(0, 2))     #The built in function returns Eigenvalues (Energies) and Eigenvectors (Wavefunctions) of the matrix defined by the main and off diagonal


'''PLotting Wavefunctions for Finite Well'''

plt.figure(figsize=(24,12))
plt.plot(x, psis2.T[0], label = 'First Energy state', color = 'r')
plt.plot(x, psis2.T[1], label = 'Second Energy state', color = 'm')
plt.plot(x, psis2.T[2], label = 'Third Energy state', color = 'c')     
plt.title('Simple Harmonic Oscillator Potential', fontsize = 30)
plt.ylabel(r'$\Psi$(z)', fontsize = 15)
plt.xlabel(r'z', fontsize = 15)
plt.legend()
plt.savefig('SHO Potential Graphs.png')


'''FInding Energies and Wavefunctions for Hydrogen Potential'''    
    
main_diag = 2*np.ones(N)/dr**2 -(2/r)                                                   #Defining main diagonal for my matrix as 2/dx^2 + Vx(Vx = potential at x)
off_diag =  -np.ones(N-1)/dr**2                                                         #Defining off diagonal for my matrix as -1/dx^2)


Es3, psis3 = eigh_tridiagonal(main_diag, off_diag, select='i', select_range=(0, 2))     #The built in function returns Eigenvalues (Energies) and Eigenvectors (Wavefunctions) of the matrix defined by the main and off diagonal


'''PLotting Wavefunctions for Finite Well'''

plt.figure(figsize=(24, 12))
plt.plot(r, psis3.T[0], color = 'red', label = 'First Energy state')
plt.plot(r, psis3.T[1], color = 'blue', label = 'Second Energy state')
plt.plot(r, psis3.T[2], color = 'Green', label = 'Third Energy state')
plt.axis([0, 24, -0.2, 0.2])      
plt.title('Radial Hydrogen Potential', fontsize = 30)
plt.ylabel(r'$\Psi$(z)', fontsize = 15)
plt.xlabel(r'z', fontsize = 15)
plt.legend()
plt.savefig('Hydrogen Potential Graphs.png')






    
