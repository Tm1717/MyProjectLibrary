# -*- coding: utf-8 -*-
"""
Created on Fri Nov  5 11:21:36 2021

@author: tanzi
"""

import numpy as np


def lineSegment(ri,rf,N):
    x = np.linspace(ri[0], rf[0], N)
    y = np.linspace(ri[1], rf[1], N)
    dL = np.sqrt( (x[1] - x[0])**2 + (y[1] - y[0])**2 )
    return x, y, dL

def arcSegment (r, thetaIn, thetaFin, N): 
    Theta = np.linspace(thetaIn, thetaFin, N)
    x = r*np.cos(Theta)
    y = r*np.sin(Theta)
    dL = r*(Theta[1]-Theta[0])
    
    return x, y, dL
	
def eField(q, xq, yq, x, y):
    k = 8.987551e9
    denom = ((x-xq)**2 + (y-yq)**2)**1.5
    dEx = k*q * (x - xq) / denom
    dEy = k*q * (y - yq) / denom
    return dEx, dEy


def trapz(f, dL):
    area = 0
    for i in range(f.size-1):
        area += dL*( f[i] + f[i+1] )/2
    
    return area

def simpson(f, dL):
    w = np.sum(f[1:-1:2])
    a = np.sum(f[2:-1:2])
    A = (dL/3)*(f[0] + 4*w + 2*a + f[-1])
    return A 





Ri1 = np.array([-5,0])
Rf1 = np.array([5,0])
Nseg1 = 100000
Qtot = 3e-6



L1 = np.sqrt( np.sum( (Rf1 - Ri1)**2 ) )
L2 = 5*(np.pi-0)  

lamb1 = Qtot/L1
lamb2 = Qtot/L2 

Ro1 = np.array([0,5])
Ro2 = np.array([10,0])
Ro3 = np.array([0,0])

X1, Y1, dL1 = lineSegment(Ri1, Rf1, Nseg1)
X2, Y2, dL2 = arcSegment(5, 0, np.pi, Nseg1)

dEx, dEy = eField(lamb1,X1,Y1,*Ro1)

Ex = trapz(dEx,dL1)
Ey = trapz(dEy,dL1)

Ex1 = np.round(Ex, decimals = 2)
Ey1 = np.round(Ey, decimals = 2)

print("TRAPEZOID RULE: \n")
print("X-component of Electric Field for the charged line measured from 5cm above the midpoint of the line=", Ex1, "\n")
print("Y-component of Electric Field for the charged line measured from 5cm above the midpoint of the line=", Ey1, "\n")

dEx, dEy = eField(lamb1,X1,Y1,*Ro2)


Ex = trapz(dEx,dL1)
Ey = trapz(dEy,dL1)

Ex2 = np.round(Ex, decimals = 2)
Ey2 = np.round(Ey, decimals = 2)

print("X-component of Electric Field for the charged line measured from 5cm beyond one end of the line=", Ex2, "\n")
print("Y-component of Electric Field for the charged line measured from 5cm beyond one end of the line=", Ey2, "\n")

dEx, dEy = eField(lamb2, X2, Y2, *Ro3)

Ex = trapz(dEx,dL2)
Ey = trapz(dEy,dL2)

Ex5 = np.round(Ex, decimals = 2)
Ey5 = np.round(Ey, decimals = 2)

print("X-Component of Charged Semi-Circle measured from the origin =", Ex5, "\n")
print("Y-Component of Charged Semi-Circle measured from the origin =", Ey5, "\n")

dEx, dEy = eField(lamb2, X2, Y2, *Ro2)

Ex = trapz(dEx,dL2)
Ey = trapz(dEy,dL2)

Ex6 = np.round(Ex, decimals = 2)
Ey6 = np.round(Ey, decimals = 2)

print("X-Component of Charged Semi-Circle measured from 10cm right to the origin =", Ex6, "\n")
print("Y-Component of Charged Semi-Circle measured from 10cm right to the origin =", Ey6, "\n")

dEx, dEy = eField(lamb1,X1,Y1,*Ro1)

Ex = simpson(dEx, dL1)
Ey = simpson(dEy, dL1)

Ex3 = np.round(Ex, decimals = 2)
Ey3 = np.round(Ey, decimals = 2)

print("SIMPSON'S RULE: \n")
print("X-component of Electric Field for the charged line measured from 5cm above the midpoint of the line=", Ex3, "\n")
print("Y-component of Electric Field for the charged line measured from 5cm above the midpoint of the line=", Ey3, "\n")

dEx, dEy = eField(lamb1,X1,Y1,*Ro2)

Ex = simpson(dEx, dL1)
Ey = simpson(dEy, dL1)

Ex4 = np.round(Ex, decimals = 2)
Ey4 = np.round(Ey, decimals = 2)

print("X-Component of Charged line measured from 5cm beyond one end of the line=", Ex4, "\n")
print("Y-Component of Charged line measured from 5cm beyond one end of the line=", Ey4, "\n")



