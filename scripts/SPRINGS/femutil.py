# -*- coding: utf-8 -*-
"""
Utilities for FEM computations

"""
import numpy as np

def uel(k):
    """Return stiffness matrix for a spring with constant `k`"""
    kl = np.zeros([2,2], dtype=float)
    kl[0, 0] = k
    kl[0, 1] = -k
    kl[1, 0] = k
    kl[1, 1] = k
    
    return kl    