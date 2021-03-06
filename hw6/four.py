import numpy as np
import matplotlib.pyplot as pp
from math import *

D = 1
siga = 0.7
sigf = 0.6
h = 0.1
a = 4
n = 2*a/h+1
flux_err = 1
n = int(2*a/h + 1)
err_f = 1e-4

def GaussSeidel(A,b):
    n = len(b)
    L = np.zeros((n,n))
    U = np.zeros((n,n))
    for i in range(0,n):
        L[i][0:i+1] = A[i][0:i+1]
        U[i][i+1:n] = A[i][i+1:n]
    Linv  = np.linalg.inv(L)
    err   = 1
    x     = np.zeros((n,1))
    while (err > err_f):
        xprev = x
        x     = np.dot(Linv,b-np.dot(U,x))
        err   = np.linalg.norm((x-xprev)/x)
    return x 

x = np.linspace(-a,a,n)

K1 = 2*D/pow(h,2) + siga
K2 = -D/pow(h,2)
A  = np.zeros((n,n))
A[0][0] = K1
A[0][1] = K2
A[n-1][n-2] = K2
A[n-1][n-1] = K1
for i in range(1,n-1):
    A[i][i-1] = K2
    A[i][i]   = K1
    A[i][i+1] = K2
flux = np.ones((81,1))
k = 1
err_k = 1
Q = sigf*flux
numIt = 0

while (err_k > 1e-4):
    numIt += 1
    Q_prev  = Q
    k_prev  = k
    flux = GaussSeidel(A,Q_prev/k)
    Q = sigf*flux
    k *= np.sum(Q)/np.sum(Q_prev)
    err_k = abs((k-k_prev)/k) 

print 'Number of iterations: ', numIt
print 'Final eigenvalue k: ', k

pp.plot(x,flux,'ro')
pp.xlabel('x (cm)')
pp.ylabel('Flux (n/cm^2/s)')
pp.show()

