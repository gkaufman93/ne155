import numpy as np
import matplotlib.pyplot as pp
from math import *

D = 1
siga = 0.7
sigf = 0.6
L = sqrt(D/siga)
L2 = siga/D
h = 0.1
a = 4
n = 2*a/h+1

fcof = np.zeros((n,n))
S = np.zeros((n,1))
fcof[0][0] = (2+pow(h,2)*L2)
fcof[0][1] = -1
fcof[n-1][n-2] = -1
fcof[n-1][n-1] = (2+pow(h,2)*L2)
for i in range(1,n-1):
    fcof[i][i-1] = -1
    fcof[i][i]   = (2+pow(h,2)*L2)
    fcof[i][i+1] = -1
flux = np.linalg.solve(fcof,S)
