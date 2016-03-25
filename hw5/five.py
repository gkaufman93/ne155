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

def Jacobi(A,b,e):
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
    err   = 1 
    x     = np.zeros((n,1))
    while (err > e):
        numIt += 1
        xprev = x
        x     = np.dot(Dinv,b-np.dot(R,x))
        err   = np.linalg.norm((x-xprev)/x)
    return numIt

def GaussSeidel(A,b,e):
    n = len(b)
    L = np.zeros((n,n))
    U = np.zeros((n,n))
    for i in range(0,n):
        L[i][0:i+1] = A[i][0:i+1]
        U[i][i+1:n] = A[i][i+1:n]
    Linv  = np.linalg.inv(L)
    numIt = 0
    err   = 1
    x     = np.zeros((n,1))
    while (err > e):
        numIt += 1
        xprev = x
        x     = np.dot(Linv,b-np.dot(U,x))
        err   = np.linalg.norm((x-xprev)/x)
    return numIt

def SOR(A,b,w,e):
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
    err   = 1
    x     = np.zeros((n,1))
    while (err > e):
        numIt += 1
        xprev = x
        x     = np.dot(DLinv,w*b-np.dot((w*U+(w-1)*D),x))
        if (np.linalg.norm(x) > 0):
            err = np.linalg.norm(x-xprev)/np.linalg.norm(x)
        else: 
            err = 1
    return numIt

def a():
    A,b = instAb(5)
    w   = 1.1

    e = 1e-6
    print '\nRelative error = ',e
    numIt = Jacobi(A,b,e)
    print 'Jacobi:',numIt
    numIt = GaussSeidel(A,b,e)
    print 'Gauss Seidel:',numIt
    numIt = SOR(A,b,w,e)
    print 'SOR:',numIt

    e = 1e-8
    print '\nRelative error = ',e
    numIt = Jacobi(A,b,e)
    print 'Jacobi:',numIt
    numIt = GaussSeidel(A,b,e)
    print 'Gauss Seidel:',numIt
    numIt = SOR(A,b,w,e)
    print 'SOR:',numIt
    
def b():
    A,b = instAb(5)
    w   = 1.0
    wIt = {}
    print ''
    while (w < 1.1):
        print '\rw = ',w,'        ',
        wIt[w] = SOR(A,b,w,1e-6)
        w += 0.000001
    print '\n'
    print 'w_opt: ',min(wIt, key=wIt.get)



