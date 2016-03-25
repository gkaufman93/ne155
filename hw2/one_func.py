import math as m
import numpy as np

def f(x):
    return m.cos(m.pi*x/2)+x**2/2

def L(x,pts,j):
   prod = 1
   for i in range(0,3):
      if i == j: 
         continue
      else:
         prod *= (x-pts[i])/(pts[j]-pts[i])
   return prod

def P3(pts):
   ln = 100
   x = np.linspace(-0.5,4.5,ln)
   f_temp = np.vectorize(f)
   F = f_temp(x)
   P = np.zeros(ln)
   for j in range(0,ln):
      for i in range(0,4):
         P[j] += f(pts[i])*L(x[j],pts,i)
   return x,F,P 



