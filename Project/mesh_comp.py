import numpy as np
from math import *


class comp_flux:
    def __init__(self, d, nx, e, ny, mat_map, mat_dict, S):
        self.d        = d
        self.nx       = nx
        self.e        = e
        self.ny       = ny
        self.mat_dict = mat_dict
        self.mat_map  = mat_map
        self.S        = S


    def calcFlux():
        self.make_coeffs()
        A       = np.array([[]])
        for j in xrange(0,self.ny+1):
            if j > 1:
                row_arr = np.zeros((self.nx+1)*(j-1),self.nx+1)
            else:
                row_arr = np.array([[]])
            row_arr = np.append(row_arr, self.get_sub(j,-1));
            row_arr = np.append(row_arr, self.get_sub(j,0));
            row_arr = np.append(row_arr, self.get_sub(j,1));
            if j < self.nx-1:
                row_arr = np.zeros((self.nx+1)*(self.nx-j-1),self.nx+1)
            A = np.append(A,row_arr,axis=0)
        return GuassSeidel(A,self.S)

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

    def get_sub(j,ind):
        sub_arr = np.zeros(self.nx+1, self.nx+1)
        if ind == 0:
            for i in xrange(0, self.nx+1):
                if (i > 0):
                    sub_arr[i-1,i] = self.get_a(i,j,'l')
                sub_arr[i,i]   = self.get_a(i,j,'c')
                if (i < self.nx):
                    sub_arr[i+1,i] = self.get_a(i,j,'r')
        elif ind == -1:
            if j == 0:
                return np.array([[]])
            else:
                for i in xrange(0,self.nx+1):
                    sub_arr[i,i] = self.get_a(i,j,'b')
        elif ind == 1:
            if j == self.nx:
                return np.array([[]])
            else:
                for i in xrange(0,self.nx+1):
                    sub_arr[i,i] = self.get_a(i,j,'t')
        return sub_arr

    def D(i,j):
        t = mat_dict(mat_map[i,j])
        return t[0]

    def sig(i,j):
        if i != 0 and j != 0:
            t1 = self.mat_dict(self.mat_map[i,j])
        else:
            t1 = (0,0)
        if i != self.nx and j != 0:
            t2 = self.mat_dict(self.mat_map[i+1,j])
        else:
            t2 = (0,0)
        if i != self.nx and j != self.ny:
            t3 = self.mat_dict(self.mat_map[i+1,j+1])
        else:
            t3 = (0,0)
        if i != 0 and j != self.ny
            t4 = self.mat_dict(self.mat_map[i,j+1])
        else:
            t4 = (0,0)
        return 0.25*self.d*self.e*(t1[1]+t2[1]+t3[1]+t4[1])

    def get_a(i,j,flag):
        return a[(i,j,flag)]

    def make_coeffs():
        coeffs        = {} 
        for i in xrange(0, self.nx+1):
            for j in xrange(0, self.ny+1):
                #left
                if i == 0 or j == 0:
                    aL = 0
                elif j == self.ny:
                    aL = -self.D(i,j)*self.e/2/self.d
                else:
                    aL = -(self.D(i,j) + self.D(i,j+1))*self.e/2/self.d 
                coeffs[(i,j,'l')] = aL

                #right
                if i == 0 or i == self.nx or j == 0:
                    aR = 0
                elif j == self.ny:
                    aR = -self.D(i+1,j)*self.e/2/self.d
                else:
                    aR = -(self.D(i+1,j) + self.D(i+1,j+1))*self.e/2/self.d 
                coeffs[(i,j,'r')] = aR

                #bottom
                if i == 0 or j == 0:
                    aB = 0
                elif i == self.nx:
                    aB = -self.D(i,j)*self.d/2/self.e
                else:
                    aB = -(self.D(i,j) + self.D(i+1,j))*self.d/2/self.e 
                coeffs[(i,j,'b')] = aB

                #top
                if i == 0 or j == 0 or j == self.ny:
                    aT = 0
                elif i == self.nx:
                    aT = -self.D(i,j+1)*self.d/2/self.e
                else:
                    aT = -(self.D(i,j+1) + self.D(i+1,j+1))*self.d/2/self.e 
                coeffs[(i,j,'t')] = aT

                #center
                if i == 0 or j == 0:
                    aC = 1
                else:
                    aC = self.sig(i,j)-(aL + aR + aB + aT)
                coeffs[(i,j,'c')] = aC

        self.a = coeffs
        return;
