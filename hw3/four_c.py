import math
import matplotlib.pyplot as pp
import numpy as np
from four_b import CompSimp38

h  = np.zeros(100)
E3 = np.zeros(100)
E4 = np.zeros(100)
E5 = np.zeros(100)
n = 3
It = math.sqrt(12)-math.sqrt(5)
for i in range(0,100):
    nt   = 1.0*n
    it   = 1.0*i
    h[i] = 1.0/(it+1.0)/nt
    E3[i] = pow(abs(It - CompSimp38(3.0,4.0,n*(i+1))),1.0/3.0)
    E4[i] = pow(abs(It - CompSimp38(3.0,4.0,n*(i+1))),1.0/4.0)
    E5[i] = pow(abs(It - CompSimp38(3.0,4.0,n*(i+1))),1.0/5.0)

pp.plot(h,E3,'ro', label = 'E^(1/3)')
pp.plot(h,E4,'go', label = 'E^(1/4)')
pp.plot(h,E5,'bo', label = 'E^(1/5)')
pp.legend(loc = 'upper left')
pp.show()

