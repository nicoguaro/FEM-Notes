# -*- coding: utf-8 -*-
"""
Created on Tue Feb 17 10:02:49 2015

@author: eafit
"""
import numpy as np
def mvprod(a,b,irow,icol):
    c=np.zeros(irow)
    for i in range(0,irow):
        c[i]=0.0
        for j in range(0,icol):
            c[i]=c[i]+a[i,j]*b[j]
    return c
#
def uel(k):
    kl=np.zeros([2,2],dtype=np.float128)
    kl[0,0]= k
    kl[0,1]=-k
    kl[1,0]=-k
    kl[1,1]= k    
    return kl    