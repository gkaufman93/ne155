import numpy as np
from scipy.interpolate import interp1d,lagrange,splev,splrep
import matplotlib.pyplot as pp

x = np.array([1.0, 2.0,  3.0,  4.0, 5.0, 6.0])
f = np.array([1.0, 3.0, 15.0, 12.0, 7.0, 3.0])
xp = np.linspace(0.75,6.25,1100)
xp1 = np.linspace(1.0,6.0,1000)
in1 = interp1d(x,f)
in2 = lagrange(x,f)
tck = splrep(x,f,s=0)
in3 = splev(xp,tck,der=0)

temp, axarr = pp.subplots(3, 1)
axarr[0].plot(xp1,in1(xp1),'b--',label='Interpolation')
axarr[0].plot(x,f,'ro',label='Data')
axarr[0].set_title('Piecewise Linear Interpolation')
axarr[0].legend(loc='upper right')
axarr[0].axis([0.75, 6.75,0,20])

axarr[1].plot(xp,in2(xp),'b--',label='Interpolation')
axarr[1].plot(x,f,'ro',label='Data')
axarr[1].set_title('Lagrange Interpolation')
axarr[1].legend(loc='upper right')

axarr[2].plot(xp,in3,'b--',label='Interpolation')
axarr[2].plot(x,f,'ro',label='Data')
axarr[2].set_title('Spline Interpolation')
axarr[2].legend(loc='upper right')

pp.show()
