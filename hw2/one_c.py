import one_func
import numpy as np
import matplotlib.pyplot as pp

v = [0.0,2.0,3.0,4.0]
x,F,P = one_func.P3(v)
f2 = np.vectorize(one_func.f)
x2 = np.array(v)
pts = f2(x2)
pp.plot(x,F,'r-',label='f(x)')
pp.plot(x,P,'g-',label='P(x)')
pp.plot(x2,pts,'bo',label='xj')
pp.legend(loc='upper left')
pp.show()
