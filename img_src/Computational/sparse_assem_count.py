#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 11:29:17 2017

@author: nguarinz
"""
from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt

anal_count = lambda nx, ny: 9*nx*ny - 9*nx - 3*ny + 4
exponent = 4
nodes_side = np.logspace(0.31, exponent, 100)
nnonzero = anal_count(nodes_side, nodes_side)
total = nodes_side**4


#%% Plots
plt.figure(figsize=(6, 3))

plt.subplot(121)
plt.plot(nodes_side**2, nnonzero)
plt.plot(nodes_side**2, total)
plt.xscale("log")
plt.yscale("log")
plt.xlabel("Number of nodes")
plt.ylabel("Elements to store")
plt.legend(["Sparse matrix", "Dense matrix"],
           frameon=False)

plt.subplot(122)
plt.plot(nodes_side**2, nnonzero/total)
plt.xscale("log")
plt.yscale("log")
plt.xlabel("Number of nodes")
plt.ylabel("Sparsity")

plt.tight_layout()
plt.savefig("../../img/Computational/sparse_storage.pdf",
            bbox_inches="tight")
plt.show()
