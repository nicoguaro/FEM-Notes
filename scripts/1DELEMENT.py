# -*- coding: utf-8 -*-
"""
--------1D Lagramge interpolation problem---------------
"""
from __future__ import division
import lagrange as la
import matplotlib.pyplot as plt
import numpy as np
from sympy import *
from sympy import init_printing
init_printing()
#
fx = lambda x: x**3 + 4*x**2 - 10.0;
#
# Assign symbols
#
x= symbols('x')
#
npts = 200
xx = np.linspace(-1, 1, npts)
fd = np.array([fx(-1.0), fx(1.0) ,fx(0.0)])
#
# Obtain Lagrange polynomials
#
pol = []
pol.append(simplify(la.LagrangPoly(x, 2, 0, [-1,1,0])))
pol.append(simplify(la.LagrangPoly(x, 2, 1, [-1,1,0])))
pol.append(simplify(la.LagrangPoly(x, 2, 2, [-1,1,0])))
#
#%% Plot P0
#
print "First polynomial", pol[0]
print "Second polynomial", pol[1]
print "Third polynomial", pol[2]

# Create a lambda function for the polynomials
pol_num = lambdify((x), pol, "numpy")

# Plotting the base functions  
plt.subplot(1,2,1)
for k in range(3):
    yy = pol_num(xx)[k]
    plt.plot(xx, yy)

# Plotting the interpolated function
yy = sum(fd[k]*pol_num(xx)[k] for k in range(len(pol)))

plt.subplot(1,2,2)
plt.plot([-1, 1, 0], fd, 'ko')
plt.plot(xx, yy)
plt.show()








































   