import numpy as np
import math

def instAb(n):
    A = np.zeros((n,n))
    b = n*np.ones((n,1))
    A[0][0:2] = 4,-1
    A[n-1][n-2:n] = -1,4
    for i in range(1,n-1):
        A[i][i-1:i+2] = -1,4,-1
    return A,b

def Jacobi(A,b):
    n = len(b)
    D = np.zeros((n,n))
    R = np.zeros((n,n))
    for i in range(0,n):
        D[i][i] = A[i][i]
        if (i-1 >= 0):
            R[i][i-1] = A[i][i-1]
        if (i+1 < n):
            R[i][i+1] = A[i][i+1]
    Dinv  = np.linalg.inv(D)
    numIt = 0
    err   = b
    x     = np.zeros((n,1))
    while ((err > 1e-6).sum() > 0):
        numIt += 1
        x   = np.dot(Dinv,b-np.dot(R,x))
        err = abs(np.dot(A,x)-b)
    return x,numIt

def GaussSeidel(A,b):
    n = len(b)
    L = np.zeros((n,n))
    U = np.zeros((n,n))
    for i in range(0,n):
        L[i][0:i+1] = A[i][0:i+1]
        U[i][i+1:n] = A[i][i+1:n]
    Linv  = np.linalg.inv(L)
    numIt = 0
    err   = b
    x     = np.zeros((n,1))
    while ((err > 1e-6).sum() > 0):
        numIt += 1
        x   = np.dot(Linv,b-np.dot(U,x))
        err = abs(np.dot(A,x)-b)
    return x,numIt

def SOR(A,b,w):
    n = len(b)
    D = np.zeros((n,n))
    L = np.zeros((n,n))
    U = np.zeros((n,n))
    for i in range(0,n):
        D[i][i]     = A[i][i]
        L[i][0:i]   = A[i][0:i]
        U[i][i+1:n] = A[i][i+1:n]
    DLinv  = np.linalg.inv(D+w*L)
    numIt = 0
    err   = b
    x     = np.zeros((n,1))
    while ((err > 1e-6).sum() > 0):
        numIt += 1
        x   = np.dot(DLinv,w*b-np.dot((w*U+(w-1)*D),x))
        err = abs(np.dot(A,x)-b)
    return x,numIt

def main():
    A,b = instAb(5)
    w   = 1.1

    x,numIt = Jacobi(A,b)
    print '\nJacobi method solution:'
    print x
    print 'Number of Iterations: ', numIt

    x,numIt = GaussSeidel(A,b)
    print '\nGauss Seidel method solution:'
    print x
    print 'Number of Iterations: ', numIt

    x,numIt = SOR(A,b,w)
    print '\nSOR method solution:'
    print x
    print 'Number of Iterations: ', numIt
    



