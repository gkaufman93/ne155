import math

def f(x):
    y = x/math.sqrt(math.pow(x,2.0)-4.0)
    return y

def CompSimp38(a,b,n):
    at = 1.0*a
    h = 1.0*(b-a)/n
    I = 0.0
    for i in range(0,n/3):
        I += (3./8.)*h*(f(at+3.0*i*h)+3.0*f(at+(3.0*i+1.0)*h)+3.0*f(at+(3.0*i+2.0)*h)+f(at+(3.0*i+3.0)*h))
    return I
