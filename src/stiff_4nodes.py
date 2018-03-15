"""
Compute the stiffness matrix for a 4-noded square element
"""
from __future__ import division, print_function
from sympy import *


def umat(nu, E):
    """2D Elasticity constitutive matrix"""
    C = zeros(3, 3)
    G = E/(1 - nu**2)
    mnu = (1 - nu)/2.0
    C[0, 0] = G
    C[0, 1] = nu*G
    C[1, 0] = C[0, 1]
    C[1, 1] = G
    C[2, 2] = G*mnu

    return C


def stdm4(x, y):
    """Four noded element strain-displacement matrix"""
    N = zeros(4)
    B = zeros(3, 8)
    N = S(1)/4*Matrix([
         (1 - x)*(1 - y),
         (1 + x)*(1 - y),
         (1 + x)*(1 + y),
         (1 - x)*(1 + y)])
    dhdx=zeros(2, 4)
    for i in range(4):
        dhdx[0,i]=diff(N[i], x)
        dhdx[1,i]=diff(N[i], y)

    for i in range(4):
        B[0, 2*i] = dhdx[0, i]
        B[1, 2*i+1] = dhdx[1, i]
        B[2, 2*i] = dhdx[1, i]
        B[2, 2*i+1] = dhdx[0, i]

    return B


# Assign symbols
x, y = symbols('x y') 
nu, E = symbols('nu E')
h = symbols('h')

K = zeros(8, 8)

# Symbolically compute matrices
C = umat(nu, E)
B = stdm4(x, y)
K_int = B.T * C * B

# Integrate final stiffness
for i in range(8):
    for j in range(8):
        K[i,j] = integrate(K_int[i,j], (x,-h,h), (y,-h,h))

knum = K.subs([(E, S(1)), (nu, S(1)/3.0), (h, S(2))])

print(knum)
