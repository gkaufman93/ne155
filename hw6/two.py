import numpy as np
import matplotlib.pyplot as pp
from math import *

fcof = np.zeros((81,81))
S = np.zeros((81,1))
D = 1
sig = 0.2
h = 0.1
S0 = 8
a = 4
n = 2*a/h + 1
S[:] = S0*pow(h,2)/D
L2 = sig/D #inverse of sqrt(D/sig) more useful for solution
fcof[0][0] = (2+pow(h,2)*L2)
fcof[0][1] = -1
fcof[80][79] = -1
fcof[80][80] = (2+pow(h,2)*L2)
for i in range(1,80):
    fcof[i][i-1] = -1
    fcof[i][i]   = (2+pow(h,2)*L2)
    fcof[i][i+1] = -1
flux = np.linalg.solve(fcof,S)
x = np.linspace(-a,a,81)
L = sqrt(D/sig)
flux_dir = np.zeros((81,1))
for i in range(0,81):
    flux_dir[i] = S0/sig*(1-cosh(x[i]/L)/cosh(a/L))
pp.plot(x,flux,'ro', label = 'Numerical Solution')
pp.plot(x,flux_dir,'g^', label = 'Direct Solution')
pp.xlabel('x (cm)')
pp.ylabel('flux (n/cm^2)')
pp.legend()
pp.show()


