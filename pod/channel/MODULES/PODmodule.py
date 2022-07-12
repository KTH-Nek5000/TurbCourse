"""
The module computes POD modes via 
singular value decomposition
"""

import numpy as np
import numpy.linalg as la





def POD(Usnp,mvect,nsnap,ifsym):
    """
    Module to compute POD modes
    Args: 
        - Usnp      = snapshots matrix Usnp
        - mvect     = array with square of mass weights (nGLLe*ncomponents)
        - nsnap      = number of snapshots
        - ifsym      = to mirror data wrt symmetry x axis
    Returns:
        - L    = left singular vector of Usnp(spatial modes in POD expansion)
        - Lam2 = matrix of singular values of Usnp (normalized vector)
        - R    = right singular vector of Usnp
        - A2   = POD coefficients
        
    """    
    
    # Mass matrix weighting
    Usnp = np.multiply(Usnp,np.tile(mvect,(Usnp.shape[1],1)).T)
    
    # SVD: 
    #   - L    = matrix whose columns are the spatial modes in POD expansion
    #   - Lam2 = diagonal matrix of singular values for Usnp
    #   - R    = matrix with temporal coefficients
    L,Lam2,R = la.svd(Usnp,full_matrices=False)
    
    print(L.shape,Lam2.shape,R.shape)

    # Mass matrix weighting
    L = np.divide(L,np.tile(mvect,(L.shape[1],1)).T)
    
    # POD coefficients
    A2= R.T*Lam2
    
    # Normalization of the singular values:
    # to be equal to the eigenvalues of correlation matrix approach
    Lam2=Lam2**2/Usnp.shape[1]
    
    return L,Lam2,R,A2