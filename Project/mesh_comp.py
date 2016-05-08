import numpy as np
from math import *
from sys import stdout


class comp_flux:
    def __init__(self, d, nx, e, ny, mat_map, mat_dict, S):
        self.d        = d
        self.nx       = nx
        self.e        = e
        self.ny       = ny
        self.mat_dict = mat_dict
        self.mat_map  = mat_map
        self.S        = S
        self.a        = {} 


    def calcFlux(self):
        print "Pre-making flux coefficients..."
        self.make_coeffs()
        print "Flux coefficients loaded. Making coefficient matrix..."
        temp_dim = (self.ny+1)*(self.nx+1)
        A = np.array([]).reshape((0,temp_dim))
        for j in xrange(0,self.ny+1):
            stdout.write('\rMaking macro-row %d...' % j)
            if j > 1:
                row_arr = np.zeros((self.nx+1,(self.nx+1)*(j-1)))
            else:
                row_arr = np.array([]).reshape((self.nx+1,0))
            row_arr = np.append(row_arr, self.get_sub(j,-1),axis=1)
            row_arr = np.append(row_arr, self.get_sub(j,0),axis=1)
            row_arr = np.append(row_arr, self.get_sub(j,1),axis=1)
            if j < self.ny-1:
                row_arr = np.append(row_arr,np.zeros\
                                    (((self.nx+1),(self.nx+1)*(self.ny-j-1))), axis=1)
            A = np.append(A,row_arr,axis=0)
            stdout.flush()
        print ''
        print 'Coefficient matrix finished! Running Gauss-Seidel subroutine...'

        L = np.zeros((temp_dim,temp_dim))
        U = np.zeros((temp_dim,temp_dim))
        for i in xrange(0,temp_dim):
            L[i][0:i+1] = A[i][0:i+1]
            U[i][i+1:temp_dim] = A[i][i+1:temp_dim]
        try:
            Linv  = np.linalg.inv(L)
        except np.linalg.linalg.LinAlgError:
            raise Exception('System unsolvable. Check cell definitions and\
                            ensure whole mesh is defined.')
        numIt = 0
        err   = 1 
        x     = np.zeros((temp_dim,1))
        err_cnt = 0
        err_tmp = 1
        while (err > 1e-6 and err_cnt < 500):
            numIt += 1
            err_tmp = err
            stdout.write('\rGauss-Seidel solution iteration: {: >5}, '.format(numIt)\
                         + 'Max error: {:1.4E}'.format(float(err)))
            stdout.flush()
            x   = np.dot(Linv,self.S-np.dot(U,x))
            err = np.linalg.norm(np.dot(A,x)-self.S)/np.linalg.norm(self.S)
            if err_tmp == err:
                err_cnt += 1
            else:
                err_cnt = 0

        print ''
        return x,numIt,err

    def get_sub(self,j,ind):
        if ind == 0:
            sub_arr = np.zeros((self.nx+1, self.nx+1))
            for i in xrange(0, self.nx+1):
                if (i > 0):
                    sub_arr[i-1,i] = self.get_a(i,j,'l')
                sub_arr[i,i]   = self.get_a(i,j,'c')
                if (i < self.nx):
                    sub_arr[i+1,i] = self.get_a(i,j,'r')
        elif ind == -1:
            if j == 0:
                return np.array([]).reshape((self.nx+1,0))
            else:
                sub_arr = np.zeros((self.nx+1, self.nx+1))
                for i in xrange(0,self.nx+1):
                    sub_arr[i,i] = self.get_a(i,j,'b')
        elif ind == 1:
            if j == self.ny:
                return np.array([]).reshape((self.nx+1,0))
            else:
                sub_arr = np.zeros((self.nx+1, self.nx+1))
                for i in xrange(0,self.nx+1):
                    sub_arr[i,i] = self.get_a(i,j,'t')
        return sub_arr

    def D(self,i,j):
        t = self.mat_dict[self.mat_map[i][j]]
        return t[0]

    def sig(self,i,j):
        if i != 0 and j != 0:
            t1 = self.mat_dict[self.mat_map[i][j]]
        else:
            t1 = (0,0)
        if i != self.nx and j != 0:
            t2 = self.mat_dict[self.mat_map[i+1][j]]
        else:
            t2 = (0,0)
        if i != self.nx and j != self.ny:
            t3 = self.mat_dict[self.mat_map[i+1][j+1]]
        else:
            t3 = (0,0)
        if i != 0 and j != self.ny:
            t4 = self.mat_dict[self.mat_map[i][j+1]]
        else:
            t4 = (0,0)
        return 0.25*self.d*self.e*(t1[1]+t2[1]+t3[1]+t4[1])

    def get_a(self,i,j,flag):
        return self.a[(i,j,flag)]

    def make_coeffs(self):
        for i in xrange(0, self.nx+1):
            for j in xrange(0, self.ny+1):
                #left
                if i == 0 or j == 0:
                    aL = 0
                elif j == self.ny:
                    aL = -self.D(i,j)*self.e/2/self.d
                else:
                    aL = -(self.D(i,j) + self.D(i,j+1))*self.e/2/self.d 
                self.a[(i,j,'l')] = aL

                #right
                if i == 0 or i == self.nx or j == 0:
                    aR = 0
                elif j == self.ny:
                    aR = -self.D(i+1,j)*self.e/2/self.d
                else:
                    aR = -(self.D(i+1,j) + self.D(i+1,j+1))*self.e/2/self.d 
                self.a[(i,j,'r')] = aR

                #bottom
                if i == 0 or j == 0:
                    aB = 0
                elif i == self.nx:
                    aB = -self.D(i,j)*self.d/2/self.e
                else:
                    aB = -(self.D(i,j) + self.D(i+1,j))*self.d/2/self.e 
                self.a[(i,j,'b')] = aB

                #top
                if i == 0 or j == 0 or j == self.ny:
                    aT = 0
                elif i == self.nx:
                    aT = -self.D(i,j+1)*self.d/2/self.e
                else:
                    aT = -(self.D(i,j+1) + self.D(i+1,j+1))*self.d/2/self.e 
                self.a[(i,j,'t')] = aT

                #center
                if i == 0 or j == 0:
                    aC = 1
                else:
                    aC = self.sig(i,j)-(aL + aR + aB + aT)
                self.a[(i,j,'c')] = aC
        return;
