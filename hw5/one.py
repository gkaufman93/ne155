import numpy as np
import matplotlib.pyplot as pp
import math

A = np.zeros((100,100))
b = np.zeros((100,1))

A[0][0:2] = 2,-1
b[0]      = 0
for i in range(1,99):
    A[i][i-1:i+2] = -1,2,-1
    b[i]          = i

A[99][98:100] = -1,2
b[99]        = 99


Am = np.asmatrix(A)
bm = np.asmatrix(b)

print 'Condition number: ',np.linalg.norm(Am.getI())*np.linalg.norm(A)

xC = (Am.getI())*bm
xD = np.linalg.solve(A,b)

pp.plot(b,xC,'rs',label='Explicit solution (c)')
pp.plot(b,xD,'g^',label='Built-in solver (d)')
pp.legend(loc = 'top left')
pp.ylabel('x')
pp.xlabel('b')
pp.title('Problem One')
pp.show()

