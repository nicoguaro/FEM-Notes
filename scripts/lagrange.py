# -*- coding: utf-8 -*-
"""
Lagrange Interpolation Routines.

@author: Juan Gomez
"""
from __future__ import division
from sympy import *
#
def LagrangPoly(x,order,i,xi=None):
    r"""Compute interpolant Lagrange Polynomials
    
    .. math::
      l_j(x)=\prod_{\substack{0<m\leq k\\ m\neq j}} \frac{x - x_m}{x_j - x_m}
    
    Parameters
    ----------
    x : Sympy symbol
        Variable for the interpolation.
    order : int
        Order of the polynomials.
    i : int
        Number of the polynomial according to `xi`.
    xi : (order + 1) list of float
        Position for the nodes.
        
    Returns
    -------
    poly : Sympy expression
        Interpolant polynomial for the `i`-th node in `xi`.
        
    Examples
    --------
    >>> from lagrange import LagrangPoly as la_poly
    >>> x = symbols('x')
    >>> pol0 = simplify(la_poly(x, 2, 0, [-1,1,0]))
    >>> pol1 = simplify(la_poly(x, 2, 1, [-1,1,0]))
    >>> pol2 = simplify(la_poly(x, 2, 2, [-1,1,0]))
    >>> print(pol0)
    x*(x - 1)/2
    
    >>> print(pol1)
    x*(x + 1)/2
        
    """
    if xi==None:
        xi = symbols('x:%d'%(order + 1))  # No entiendo esta condicion
    index = range(order + 1)
    index.pop(i)
    poly = prod([(x - xi[j])/(xi[i] - xi[j]) for j in index])
    return poly
    
    
# Run examples as tests
if __name__=="__main__":
    import doctest
    doctest.testmod()