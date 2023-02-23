# -*- coding: utf-8 -*-
"""
Created on Sun Oct 16 14:36:12 2022

@author: tanzi
"""

""" Section 1:  Start by importing relevant libraries
#---------------------------------------------------------------------------"""
import matplotlib.pyplot as plt
import numpy as np 

""" Section 2:  Define functions for the script
#---------------------------------------------------------------------------"""

def trapz(f, dL):
    area = 0
    for i in range(f.size-1):
        area += dL*( f[i] + f[i+1] )/2
    
    return area

def Bnvals(x, y, n):
    f = np.zeros(x.size)
    Bn = np.zeros(n+1)
    N = np.zeros(n+1)
    for n in range(0, n+1):    
        for i in range(xdata.size):
            f[i] = (y[i] * np.sin((2 * np.pi * n * x[i]) / L))
        
        Bn[n] =  2/L *  trapz(f, dL)
        
    return Bn    


""" Section 3:  Main body 
#---------------------------------------------------------------------------"""



data = np.loadtxt("wavedata.csv", delimiter = ',')    #Reading data from file.

xdata = data[: , 0]         #seperating x and y data
ydata = data[: , 1]


dL = xdata[1] - xdata[0]
L = 10
f = np.zeros(xdata.size)
B10 = np.zeros(11)
Fx10 = np.zeros(xdata.shape)
Fx4 = np.zeros(xdata.shape)
B10 = np.zeros(11)
N10 = np.zeros(B10.size)
B4 = np.zeros(5)    
N4 = np.zeros(B4.size)

A0 = 1/L * trapz(ydata, dL)             #Finding A0 value 
B10 = Bnvals(xdata, ydata, 10)     #finding Bns for n = 10
B4 = Bnvals(xdata, ydata, 10)      #Finding Bns for n = 4
    

#Used to add up Terms for n = 10 and n = 4
for x in range(xdata.size):
    Fx10[x] = A0 + (B10[0] * np.sin((2 * np.pi* xdata[x] * 0) / L)) + (B10[1] * np.sin((2 * np.pi* xdata[x] * 1) / L)) + (B10[2] * np.sin((2 * np.pi* xdata[x] * 2) / L)) + (B10[3] * np.sin((2 * np.pi* xdata[x] * 3) / L)) + (B10[4] * np.sin((2 * np.pi* xdata[x] * 4) / L)) + (B10[5] * np.sin((2 * np.pi* xdata[x] * 5) / L)) + (B10[6] * np.sin((2 * np.pi* xdata[x] * 6) / L)) + (B10[7] * np.sin((2 * np.pi* xdata[x] * 7) / L)) + (B10[8] * np.sin((2 * np.pi* xdata[x] * 8) / L)) + (B10[9] * np.sin((2 * np.pi* xdata[x] * 9) / L)) + (B10[10] * np.sin((2 * np.pi* xdata[x] * 10) / L))
      
      
for x in range(xdata.size):
    Fx4[x] = A0 + (B10[0] * np.sin((2 * np.pi* xdata[x] * 0) / L)) + (B10[1] * np.sin((2 * np.pi* xdata[x] * 1) / L)) + (B10[2] * np.sin((2 * np.pi* xdata[x] * 2) / L)) + (B10[3] * np.sin((2 * np.pi* xdata[x] * 3) / L)) + (B10[4] * np.sin((2 * np.pi* xdata[x] * 4) / L))

""" Section 4 plotting
#---------------------------------------------------------------------------"""
##Note to future self this is how you add subscript to text output $whateverTakesSubscriptHere_{subscript here}$##
Bnlabels = ['$B_{0}$', '$B_{1}$', '$B_{2}$', '$B_{3}$', '$B_{4}$', '$B_{5}$', '$B_{6}$', '$B_{7}$', '$B_{8}$', '$B_{9}$', '$B_{10}$']

fig = plt.figure(figsize = (24, 8))
ax1 = fig.add_subplot(1,2,1)
ax2 = fig.add_subplot(1,2,2)

ax1.bar(Bnlabels,B10)
ax1.set(xlabel='$B_{n}$', ylabel = 'Contribution to f(x)', title = 'Bargraph showing contribution of each $B_{n}$ coefficient' )
ax2.plot(xdata, ydata, '-k', label = 'Original f(x)')
ax2.plot(xdata, Fx10, '.r', label = 'f(x) for n = 10')
ax2.plot(xdata, Fx4, '--g', label = 'f(x) for n = 4')
ax2.set(xlabel = 'x', ylabel = 'f(x)', title = 'Comparing original f(x) to reconstructed f(x)')
plt.legend()
plt.show()    
    
plt.savefig('Fourier Series Comp Assignment.png')                             

                             
