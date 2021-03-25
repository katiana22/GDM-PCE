from numpy import *
import pylab as p
from scipy import integrate

def LV(a, b, c, d, t):
    
    X_f0 = array([     0. ,  0.])
    X_f1 = array([ c/(d*b), a/b])

    def dX_dt(X, t=0):
        """ Return the growth rate of fox and rabbit populations. """
        return array([ a*X[0] -   b*X[0]*X[1] ,  
                      -c*X[1] + d*b*X[0]*X[1] ])

    X0 = array([10, 5])                     # initials conditions: 10 rabbits and 5 foxes  

    X, infodict = integrate.odeint(dX_dt, X0, t, full_output=True)
    
    return X, X_f1
