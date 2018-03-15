# -*- coding: utf-8 -*-
"""
--------2D Lagrange interpolation problem---------------
"""
from __future__ import division
import lagrange as la
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from sympy import *
from sympy import init_printing
init_printing()
#
# Assumed displacements
#
ud = np.array([0.0, -0.1 ,-0.05,0.05])
vd = np.array([0.0,  0.0 , 0.1, 0.1])
#
# Assign symbols
#
x, y= symbols('x y')
#
# Compute 1D Lagrange polynomials 
#
polx = []
poly = []
for i in range(2):
    polx.append(simplify(la.LagrangPoly(x, 1, i, [-1,1])))
    poly.append(simplify(la.LagrangPoly(y, 1, i, [-1,1])))
    
#
# Compute 2D Lagrange polynomials using products of the 1D Pols
#

polxy = [px*py for px in polx for py in poly]
for k in range(len(polxy)):
    print "Polynomial %i: "%k, polxy[k]

#
# Generate sampling points and mesh
#
li=-1.0
ls= 1.1
dl= 0.1
npts=int((ls-li)/dl)
xx, yy = np.mgrid[li:ls:npts*1j, li:ls:npts*1j]
#
# Plot the polynomials and the displaced element
#
fig1 = plt.figure()
fig2 = plt.figure()
fig3 = plt.figure()
U = np.zeros(npts)
V = np.zeros(npts)
for k in range(4):
    f = lambdify((x,y), polxy[k], "numpy")
    z = f(xx, yy)
    U = U + f(xx,yy)*ud[k]
    V = V + f(xx,yy)*vd[k]
    ax =fig1.add_subplot(2, 2, k, projection='3d')
    ax.plot_surface(xx, yy, z, rstride=1, cstride=1, cmap=plt.cm.hot)
ax =fig2.add_subplot(1, 1, 1, projection='3d')
ax.plot_surface(xx, yy, U, rstride=1, cstride=1, cmap=plt.cm.hot)
ax =fig3.add_subplot(1, 1, 1, projection='3d')
ax.plot_surface(xx, yy, V, rstride=1, cstride=1, cmap=plt.cm.hot)
plt.show()









































   