# -*- coding: utf-8 -*-
"""
Computes the displacements and internal forces for a mass-springs
system under static loads.  

Variables
----------
ne : number of elements
nn : number of nodes
nm : number of material profiles
nl : number of point loads
IDN[] : Stores the nodal indentifier
IBC[] : Stores the nodal boundary condition (0:free; 1: restrained)
COR[] : Stores the nodal coordinates
MIE[] : Stores the equation lists for each element.
MAT[] : Stores the material properties.
IM[]  : Stores the material profile for each element.
IP[]  : Stores the nodal identifier for each point load.
PL[]  : Stores points loads vector.
COR[] : Stores nodal coordinates.
KG[]   : Stores the global stiffness matrix.
RHSG[] : Stores the global right hand side vector.
U[]    : Stores the global displacements vector.
    
@author:Juan David Gómez

"""
import numpy as np
import femutil as fem
import matplotlib.pyplot as plt


def read_model(node_file, mater_file, els_file, load_file, verbose=True):
    """Read model data from text files"""
    # rea
    nodes = np.loadtxt(node_file)
    mats = np.loadtxt(mater_file)
    elements = np.loadtxt(els_file, dtype=int)
    loads = np.loadtxt(load_file)
    
    # Generate echo files
    if verbose:
        np.savetxt("KNODES.txt", nodes,    fmt='%5.2f', delimiter=' ')
        np.savetxt("KMATES.txt", mats,     fmt='%5.2f', delimiter=' ')
        np.savetxt("KELEMS.txt", elements, fmt='%5.2f', delimiter=' ')
        np.savetxt("KLOADS.txt", loads,    fmt='%5.2f', delimiter=' ')
    
    return nodes, mats, elements, loads


nodes, mats, elements, loads = read_model('nodes.txt', 'mater.txt',
                                          'eles.txt', 'loads.txt')
#
# Retrieves problem parameters and initializes
# arrays.
#
ne = len(elements[:,0])
nn = len(nodes[:,0])
nm = len(mats[:])
nl = len(loads[:,0])
#
IDN = np.zeros([nn,1],dtype=int)
IBC = np.zeros([nn,1],dtype=int)
COR = np.zeros([nn,1],dtype=float)
MIE = np.zeros([ne,2],dtype=int)
IM  = np.zeros([ne,1],dtype=int)
IP  = np.zeros([nl,1],dtype=int)
PL  = np.zeros([nl,1],dtype=float)
#
# Pre-processing begins
# 
# Counts active equations
#
icount = 0
for i in range(0,nn):
    IDN[i] = int(nodes[i,0])
    COR[i] = nodes[i,1]
    IBC[i] = int(nodes[i,2])
    if IBC[i] == 0:
        IBC[i] = icount
        icount = icount + 1
#
# Assembles MIE in translated form
#        
for i in range(0,ne):
    IM[i] = elements[i,3]
    for j in range(1,3):
        MIE[i,j-1] = IBC[elements[i,j]]
#
# Reads points loads
#
for i in range(0,nl):
    IP[i] = int(loads[i,0])
    PL[i] = loads[i,1]
#
# Starts global system assembly
#
KG = np.zeros([icount,icount], dtype=float)
RHS = np.zeros([icount], dtype=float)
U = np.zeros([icount], dtype=float)
#
# Loops through all the elements
#        
for i in range(0,ne):
    imm = IM[i]
    k = mats[imm]
    lmie = np.zeros([2],dtype=int)
    for j in range(0,2):
        lmie[j] = MIE[i,j]
#
#  Calls UEL
#        
    kloc = fem.uel(k)
    print "Local stiffnes matrix for element %i:\n"%i, kloc
# 
#   Global stiffness
#    
    for ii in range(0,2):
        kk=lmie[ii]
        if kk != -1:
            for jj in range(0,2):
                ll = lmie[jj]
                if ll != -1:
                    KG[kk,ll] = KG[kk,ll] + kloc[ii,jj]
#                    
plt.figure()
plt.spy(KG)
plt.title("Stiff matrix")
plt.ylabel(r"$i$ index", size=14)
plt.xlabel(r"$j$ index", size=14)                   
#                    
# Global RHS
#                    
for i in range(0,icount):
    il = IBC[IP[i]]
    RHS[il] = PL[i]

#%%
#
# Solution begins
#
U = np.linalg.solve(KG, RHS)
print 'Global displacements', U    
equil = np.allclose(np.dot(KG, U), RHS)
if equil:
    print "Equilibrium satisfied."
else:
    print "Equilibrium not satisfied"

#%%
# Post-processing begins
#
for i in range(0,ne):
    imm = IM[i]
    k = mats[imm]
    ul = np.zeros([2],dtype=float)
    rhsl = np.zeros([2],dtype=float)
    for j in range(0,2):
        kk = MIE[i,j]
        if kk == -1:
            ul[j] = 0.0
        else:
            ul[j] = U[IBC[kk]]
    kloc = fem.uel(k)
    rhsl = np.dot(kloc, ul)
    print 'Internal forces for element %i = '  % i, np.round(rhsl, 6)

#
#   Prints element results
#
#
#np.savetxt("UGLOB.txt", U, fmt='%5.2f', delimiter=' ')
plt.figure()
plt.plot(U, 'ro')
plt.ylabel('Global displacements', size=14)
plt.xlabel('Global coordinate', size=14)

print('Program terminated')
plt.show()    

