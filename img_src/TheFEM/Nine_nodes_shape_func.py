# -*- coding: utf-8 -*-
"""
Plot the interpolation functions for a 4-node isoparametric
element.

@author: Nicolas Guarin-Zapata
"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import rcParams
import seaborn as sns
rcParams['font.family'] = 'serif'
rcParams['font.size'] = 14

def make_plot(x, y, N):
    coords = [-1, 0, 1]
    pts = np.array([[ptx, pty] for ptx in coords for pty in coords])
    fig = plt.figure(figsize=(8, 5))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot([-1, 1, 1, -1, -1], [-1, -1, 1, 1, -1], "-k", zorder=-10)
    ax.plot(pts[:, 0], pts[:, 1], "ko", zorder=-10)
    ax.plot_surface(x, y, N, cstride=1, rstride=1, cmap="YlGnBu_r",
                alpha=0.6, lw=0.5, zorder=3, vmin=-0.25, vmax=1)
    ax.view_init(azim=-60, elev=30)
    ax.set_xlabel(r"$x$", fontsize=18)
    ax.set_ylabel(r"$y$", fontsize=18)
    ax.set_zlabel(r"$N^%i(x, y)$"%cont, fontsize=18)
    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)
    ax.set_zlim(-0.25, 1)
    plt.savefig("../../img/TheFEM/shape_func-9-nodes-%i.pdf"%cont,
                bbox_inches="tight", pad_inches=0.1, transparent=True)


x, y = np.mgrid[-1:1:21j, -1:1:21j]
L = lambda x: [-0.5*(1 - x)*x, 0.5*(1 + x)*x, 1 - x**2]
shapes = np.array([Lx*Ly for Lx in L(y)for Ly in L(x)])
shapes = shapes[[0, 1, 4, 3, 2, 7, 5, 6, 8], :, :]
cont = 0
for shape in shapes:
    cont = cont + 1
    make_plot(x, y, shape)    

plt.show()