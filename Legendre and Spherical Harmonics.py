# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 14:04:36 2023

@author: tanzi
"""
import numpy as np 
from scipy.special import factorial2
import math
import matplotlib.pyplot as plt
from scipy.special import lpmv

def legendre(m, n, x):
    # Check if m and n satisfy certain conditions for legendre functions 
    if m > n or n < 0 or abs(m) > n:                
        return 0
    # Check if m is equal to n, in which case we can use Eq. (6) to compute P_n^m(x)
    if m == n:
        return ((-1) ** m) * factorial2(2 * n - 1) * (1 - x ** 2) ** (m / 2)
    # Check if m is equal to n-1, in which case we can use Eq. (7) to compute P_n^m(x)
    elif m == n -   1:
        return x * (2 * n - 1) * legendre(m, n - 1, x)
    #Check if m is negative, in which case we can use Eq. (8) to compute P_n^m(x)
    elif m < 0: 
        return (-1)**m * (math.factorial(n - m)/math.factorial(n + m)) * legendre(m, n, x)
    # If m is not equal to n or n-1, we use Eq. (5) to compute P_n^m(x)
    else:
        return ((2 * n - 1) * x * legendre(m, n - 1, x) - (n + m - 1) * legendre(m, n - 2, x)) / (n - m)
    


x = [-1 + i*(2/99) for i in range(101)]     #used list comprehension to create a list of values from -1 to 1 
del x[-1]                                   #However when I use range(100) i dont get values including 1 so 
                                            #So I used range(101) and deleted the last item in the list in above line

# Compute y values using legendre functions one from my defined function the other for scipi function

y1 = [legendre(2, 4, i) for i in x]  

y2 = [lpmv(2, 4, i) for i in x]

# Plot the results

plt.plot(x, y1)
plt.title("P24(x)")
plt.xlabel("x")
plt.ylabel("P24(x)")
plt.show()    

plt.plot(x, y1, label="Custom function")
plt.plot(x, y2, label="scipy.special.lpmv")
plt.title("P24(x)")
plt.xlabel("x")
plt.ylabel("P24(x)")
plt.legend()
plt.show()




def real_spherical_harmonics(m, n, theta, phi):
    x =  math.sin(theta) * math.cos(phi)
    if m < 0:
        return (-1)**m * math.sqrt(2) * math.sqrt(((2*n + 1) / (4 * math.pi)) * (math.factorial(n - abs(m))/math.factorial(n + abs(m)))) * legendre(abs(m), n, x) * (math.cos(theta)) * math.sin(abs(m) * phi) 
    elif m == 0:
        return math.sqrt((2 * n + 1) / (4 * math.pi)) * legendre(0, n, x) * math.cos(theta)
    else:
        return (-1)**m * math.sqrt(2) * math.sqrt((2 * n + 1) / (4 * math.pi) * (math.factorial(n - m)/math.factorial(n + m))) * math.cos(m * phi) * legendre(m, n, x) * (math.cos(theta))

# Generate a grid of theta and phi values
theta = []
phi = []
n_theta, n_phi = 100, 100
for i in range(n_theta):
    for j in range(n_phi):
        t = (i+1)/n_theta * math.pi
        p = (j)/n_phi * 2 * math.pi
        theta.append(t)
        phi.append(p)

# Compute the function for n=4 and m=2
n, m = 4, 2
r = []
for t, p in zip(theta, phi):
    r.append(abs(real_spherical_harmonics(m, n, t, p)))

# Convert to Cartesian coordinates
x, y, z = [], [], []
for i in range(len(theta)):
    x.append(r[i] * math.sin(theta[i]) * math.cos(phi[i]))
    y.append(r[i] * math.sin(theta[i]) * math.sin(phi[i]))
    z.append(r[i] * math.cos(theta[i]))

# Plot the surface
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(np.array(x).reshape((100,100)), np.array(y).reshape((100,100)), np.array(z).reshape((100,100)), rstride=1, cstride=1, cmap='viridis')
ax.set_title('|Y_42|')
plt.show()