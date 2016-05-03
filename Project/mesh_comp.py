import numpy as np
from math import *


class comp_flux:
    def __init__(self, d, nx, e, ny, mat_dict):
        self.d        = d
        self.nx       = nx
        self.e        = e
        self.ny       = ny
        self.mat_dict = mat_dict

    def calcFlux():
        A       = np.array([[]])
        for j in xrange(0,self.ny):
            if j > 0:
                row_arr = np.zeros(self.nx*(j-1),self.ny)
            else:
                row_arr = np.array([[]])
            row_arr = np.append(row_arr, get_sub(j,-1));
            row_arr = np.append(row_arr, get_sub(j,0));
            row_arr = np.append(row_arr, get_sub(j,1));
            if j < (self.nx - 1):
                row_arr = np.zeros(self.nx*(self.nx-j-1),self.ny)
            A = np.append(A,row_arr,axis=0)


    def get_sub(j,ind):
        sub_arr = np.zeros(self.nx, self.ny)
        if ind == 0:

