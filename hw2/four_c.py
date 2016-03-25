import numpy as np
import math

def main():
   h = [5.00000e-02, 2.50000e-02, 1.25000e-02, 6.25000e-03, 3.12500e-03, 1.56250e-03, 7.81250e-04, 3.90625e-04]
   E = [1.036126e-01, 3.333834e-02, 1.375409e-02, 4.177237e-03, 1.103962e-03, 2.824698e-04, 7.185644e-05, 1.813937e-05]
   lh = []
   lE = []
   for i in h:
      lh.append(math.log(i))
   for i in E:
      lE.append(math.log(i))
   lh = np.array(lh)
   lE = np.array(lE)
   lh_sum = sum(lh)
   lh2_sum = sum(np.multiply(lh,lh))
   lE_sum = sum(lE)
   lElh_sum = sum(np.multiply(lh,lE))
   A = np.array([[8.0, lh_sum],[lh_sum, lh2_sum]])
   b = np.array([[lE_sum],[lElh_sum]])
   x = np.linalg.solve(A,b)
   k = math.exp(x[0])
   p = x[1]
   return A,b,x,k,p
