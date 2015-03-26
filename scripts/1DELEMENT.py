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
fx = lambda x: 0.5*x + 1.0;
#
# Assign symbols
#
x= symbols('x')
#
npts = 200
xx = np.linspace(-1, 1, npts)
yy = np.zeros((npts))
fd = np.array([0.5, 1.5 ,1.0])
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

# Plotting the base functions
plt.subplot(1,2,1)
for k in range(3):
    for i in range(npts):
        yy[i] = pol[k].subs([(x, xx[i])])
        
    plt.plot(xx, yy)

# Plotting the interpolated function
for i in range(npts):
    yy[i] = fd[0]*pol[0].subs([(x, xx[i])]) + fd[2]*pol[1].subs([(x, xx[i])]) \
            + fd[1]*pol[2].subs([(x, xx[i])])

plt.subplot(1,2,2)
plt.plot([-1, 0, 1], fd, 'ko')
plt.plot(xx, yy)
plt.show()








































   