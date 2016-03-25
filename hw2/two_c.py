import numpy as np
import math

x = np.linspace(-0.5,4.5,1000)
e = np.ones(1000)
for i in range(0,1000):
   e[i] = (np.pi/2.0)**4/math.factorial(4)*x[i]*(x[i]-2.0)*(x[i]-3.0)*(x[i]-4.0)
m = np.amax(e)
print 'Max error: {}'.format(m)
