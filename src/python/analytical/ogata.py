import math
from xerf import exp_erfc
from scipy import special

'''
@author: Glenn Hammond, Pacific Northwest National Laboratory

Anaytical solution for 1D solute transport with a fixed or constant input 
concentration and a first-type boundary condition.

Source: Ogata, A. and R.B. Banks (1961) A Solution of the Differential
        Equation of Longitudinal Dispersion in Porous Media, Fluid
        Movement in Earth Materials, Geological Survey Professional 
        Paper 411-A, US Dept. of Interior.
'''

def ogata(C0,D,v,x,t):
    twosqrtDT = 2.*math.sqrt(D*t)
    return C0*0.5*(special.erfc((x-v*t)/twosqrtDT) + \
                   exp_erfc(v*x/D, \
                            (x+v*t)/twosqrtDT))
