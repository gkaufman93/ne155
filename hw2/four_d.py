import numpy as np
import math
import matplotlib.pyplot as pp
import four_c

h = np.array([5.00000e-02, 2.50000e-02, 1.25000e-02, 6.25000e-03, 3.12500e-03, 1.56250e-03, 7.81250e-04, 3.90625e-04])
E = np.array([1.036126e-01, 3.333834e-02, 1.375409e-02, 4.177237e-03, 1.103962e-03, 2.824698e-04, 7.185644e-05, 1.813937e-05])
Ef = np.array([])
A,b,x,k,p = four_c.main()
for i in h:
   Ef = np.append(Ef,k*math.pow(i,p))
pp.loglog(h,Ef,'r*',label='Least Squares')
pp.loglog(h,E,'bo',label='Data')
pp.legend(loc='upper right')
pp.show()
