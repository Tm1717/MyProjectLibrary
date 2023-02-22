# -*- coding: utf-8 -*-
"""
Created on Sun Mar  6 21:33:17 2022

@author: tanzi
"""


""" Libraries Imported here
#---------------------------------------------------------------------------"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# this formats the output of numpy arrays to the console
np.set_printoptions(precision=4)


""" Functions Defined Here 
#---------------------------------------------------------------------------"""


def amplitude(data):
    
    A = (abs(np.min(data)) + np.max(data))/2
    
    

    return A

def cosfit(t, A, P, delta, d):
    return A*np.cos(2*np.pi*t/P - delta) + d

def linOscillator(wo,beta,f,w,Xo,Vo,numPeriods,N):
    
    period = 2*np.pi/w       # define period based on input frequency
    #N = 1000                # number of steps per period
    
    # Create time array using linspace; first point t = 0, then N points 
    # for each period afterwards
    T = np.linspace(0,numPeriods*period,numPeriods*N+1)
    
    dt = T[1]               # time step for Euler-Cromer method

    xEul = np.zeros(T.size)
    vEul = np.zeros(T.size)
    
    xEul[0] = Xo
    vEul[0] = Vo             # initial condition xo = 1, vo = 0
    
    for t in range(T.size-1):
        #accel = -wo**2*xEul[t] -K*xEul[t]**3  -2*beta*vEul[t] + f*wo**2*np.cos(w*T[t])
        accel = -wo**2*np.sin(xEul[t]) -2*beta*vEul[t] + f*wo**2*np.cos(w*T[t])
        vEul[t+1] = vEul[t] + accel*dt
        xEul[t+1] = xEul[t] + vEul[t+1]*dt

    return T,xEul,vEul


w = np.linspace(0, 0.01, 5)
A = np.zeros(w.size)

""" Main body of code
#---------------------------------------------------------------------------"""
# Input parameters for the damped oscillator
m = 1                   # mass
k = 0.6                   # spring constant
b = 0.05                 # damp coefficient

wo = (k/m)**(1/2)       # natural frequency
beta = b/(2*m)          # damping coeff.


# Define parameters for time - can be used as inputs into our function
period = 2*np.pi/wo     # period of oscillation
numPeriods = 100        # calculate motion over # of periods
N = 100                # number of time steps per period

# Define parameters for driving force
F = 0.20                # amplitude of driving force
w = 1.2                 # frequency of driving force
f = F/m                 # normalized force for linOscillator DE

# Initial conditions
initX = 0.0
initV = 0.0

w = np.arange(0.1, 1.8, 0.01)      #create driving frequency values
A = np.zeros(w.size)            
phi = np.zeros(w.size)

for i in range (w.size):                                                        #This loop saves amplitude and phase shift data 
    T, X, V = linOscillator(wo,beta,f,w[i],initX,initV,numPeriods,N)            #It iterates linOscillator (w.size) times 
    A[i] = amplitude(X)
    
    prm1, pcov1 = curve_fit(cosfit, T, X, p0 = [A[i], (2*np.pi/w[i]),1,0] )     #Collect Amplitude and Phase Shift from curve fitting 
    A[i] = abs(prm1[0])
    if prm1[2] < 0:         #This accounts for negative phase shift data 
        prm1[2] += np.pi
    phi[i] =  prm1[2]
    
""" Plotting
#---------------------------------------------------------------------------"""    
    
fig = plt.figure(figsize=(12,9))

ax0 = fig.add_subplot(121)
ax0.set_title('Driving Force vs Amplitude Graph')
ax0.set_ylabel('Amplitude')
ax0.set_xlabel('Driving Frequency')
ax0.plot(w,A, 'r')

ax1 = fig.add_subplot(122) 
ax1.set_title('Driving Force vs Phase Angle Graph')
ax1.set_ylabel('Phase Shift')
ax1.set_xlabel('Driving Frequency')
ax1.plot(w, phi, 'b')

   
plt.savefig('resonance.png', dpi=300)    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    