import numpy as np
import matplotlib.pyplot as pp
from math import *

D = 1
sig = 0.2
S0 = 8
a = 4
L2 = sig/D #inverse of sqrt(D/sig) more useful for solution
L  = sqrt(D/sig)

def flux_mesh(h):
    n = int(2*a/h + 1)
    fcof = np.zeros((n,n))
    S = np.zeros((n,1))
    S[:] = S0*pow(h,2)/D
    fcof[0][0] = (2+pow(h,2)*L2)
    fcof[0][1] = -1
    fcof[n-1][n-2] = -1
    fcof[n-1][n-1] = (2+pow(h,2)*L2)
    for i in range(1,n-1):
        fcof[i][i-1] = -1
        fcof[i][i]   = (2+pow(h,2)*L2)
        fcof[i][i+1] = -1
    flux = np.linalg.solve(fcof,S)
    return flux[1:n-1]

def flux_dir(h):
    n = int(2*a/h + 1)
    flux= np.zeros((n,1))
    x = np.linspace(-a,a,n)
    for i in range(0,n):
        flux[i] = S0/sig*(1-cosh(x[i]/L)/cosh(a/L))
    return flux[1:n-1]

err = []
n = []
hs  = [1, 0.5, 0.1, 0.05, 0.01]
for h in hs:
    n.append(2*a/h + 1)
    fm = flux_mesh(h)
    fd = flux_dir(h)
    err.append(1-max(abs((fm-fd)/fd)))
err_arr = np.array(err)
n_arr = np.array(n)
pp.plot(n_arr,err_arr,'ro')
pp.xlabel('Total Number of Meshes')
pp.ylabel('Maximum Relative Error')
pp.show()




