import numpy as np
import math as mt

def GLC_pwts(n):
    """ 
    Gauss-Lobatto-Chebyshev (GLC) points and integration weights over [-1,1]    
    Args: 
      `n`: int, number of nodes in physical space
    Returns 
       `x`: 1D numpy array of size `n`, nodes         
       `w`: 1D numpy array of size `n`, weights
    """
    def c(i,n):
        c_=2.
        if i==0 or i==n-1:
           c_=1.
        return c_
    theta=np.arange(n)*np.pi/(n-1)    
    x=np.cos(theta)
    w=np.zeros(n)    
    for k in range(n):
        tmp_=0.0
        for j in range(1,int((n-1)/2)+1):
            bj=2
            if n%2!=0 and j==int((n-1)/2):
               bj=1 
            tmp_+=bj/(4.*j**2-1)*mt.cos(2*j*theta[k])
        w[k]=(1-tmp_)*c(k,n)/float(n-1)
    return x,w 