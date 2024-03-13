# -*- coding: utf-8 -*-
"""
Created on Tue May  7 21:55:14 2019

@author: Glenn Hammond, Sandia National Laboratories, gehammo@sandia.gov
"""

import math

# this number is platform dependent. it is a large float where exp(-amax)
# avoids underflow and exp(amax) avoids overflow. outside this bound, 
# exp() will throw an error.
amax = 500. # this number is platform dependent. it prevents over/underflow.

def exp_erfc(a,b):
    '''This subroutine is based on Fortran code from:
 
      van Genuchten, M. T. and W.J. Alves, (1982) Analytical solutions of the 
      one-dimensional convective-dispersive solute transport equation, USDA 
      Tech. Bull. 1661.

      where exp(a)*erfc(b) is evaluated at extreme values of a and b.
    '''
    value = 0.
    if abs(a) > amax and b <= 0.:
      return value
    if b == 0.:
      return expd(a)
    c = a-b*b
    if abs(c) > amax and b > 0.:
      return value
    if c < -amax:
      if b < 0.:
        value = 2.*expd(a)
    else:
      x = abs(b)
      if x > 3.:
        y = 0.5641896/(x+0.5/(x+1./(x+1.5/(x+2./(x+2.5/(x+1.))))))
        value = y*expd(c)
        if b < 0.:
          value = 2.*expd(a)-value
      else:
        t = 1./(1.+0.3275911*x)
        y = t*(0.2548296-t*(0.2844967-t*(1.421414-t*(1.453152-1.061405*t))))
        value = y*expd(c)
        if b < 0.:
          value = 2.*expd(a)-value
    return value
      
def expd(x):
    if x < -amax:
      return 0.
    return math.exp(x)
