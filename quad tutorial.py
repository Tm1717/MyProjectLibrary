# Computational Exercise 1 for PHYS*2310 - Winter 2022
# quadDrag.py
# MVM
#------------------------------------------------------------------------------

""" Section 1:  Start by importing relevant libraries
#---------------------------------------------------------------------------"""
import numpy as np
import matplotlib.pyplot as plt


""" Section 2:  Define functions for the script
#---------------------------------------------------------------------------"""
def quadDragCalc(dt, Tf, yo, vo):
    T = np.arange(0,Tf+dt,dt)
    yEul = np.zeros(T.size)
    vEul = np.zeros(T.size)
    
    yEul[0] = yo
    vEul[0] = vo
    
    for t in range(T.size-1):
        accel = g - dragCoeff*vEul[t]**2/m
        vEul[t+1] = vEul[t] + accel*dt
        yEul[t+1] = yEul[t] + vEul[t]*dt
    
    return T, yEul, vEul

""" Section 3a:  Main body of code - define quantities
#---------------------------------------------------------------------------"""
# Input parameters for a falling egg, dropped from an airplane at 10 km alt
# assume all quantities in SI units
kappa=0.25      # dimensionless coeff for sphere = 1/4
rho=1.29        # air density in kg/m^3 at STP
r=0.015          # cotton ball
m=0.0007         # mass of a cotton ball

A=np.pi*r**2    # x-sectional area of sphere (assumption for the ball)
g=9.803         # acceleration due to gravity m/s^2

dragCoeff = kappa*rho*A/2     #In 2310 they use F = -cv^2, c=dragCoeff

# Define the time step (dt) and the total fallTime for the object
dt = 0.1
fallTime = 10.0

# Initial conditions
yo = 0.0
vo = 0.0

# call our function to compute 1D Euler's method
T, Y, V = quadDragCalc(dt,fallTime, yo, vo)

# these equations are the analytical solution ("theoretical") to 1D free fall 
# with quadratic drag, subject to our initial conditions (i.e. from rest)
vTh = np.sqrt(m*g/dragCoeff)*np.tanh(np.sqrt(dragCoeff*g/m)*T)
yTh = (m/dragCoeff)*np.log(np.cosh(np.sqrt(dragCoeff*g/m)*T))



""" Section 4:  Plot your results as a contour plot
#---------------------------------------------------------------------------"""

# fig = plt.figure(figsize=(12,6))
# add axes to this figure using "ax = fig.add_subplot(rci)"
fig = plt.figure(figsize=(12,6))
axY = fig.add_subplot(121)
axV = fig.add_subplot(122)

# plot the computed and theoretical
axY.plot(T,Y,'rx')
axV.plot(T,V,'gx')
axY.plot(T,yTh,'r')
axV.plot(T,vTh,'g')


# NEW PLOTTING FEATURE - NICE!!! 
#plt.tight_layout()

# Save your plot
#plt.savefig('quadDragGraphs.png', dpi=300)
