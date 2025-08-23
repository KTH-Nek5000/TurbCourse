"""
The module computes DMD modes via 
singular value decomposition
"""

import numpy as np
import numpy.linalg as la

def hermitian(arr):
    return arr.conj().T

def DMD(Usnp,mvect,nsnap,r,ifsym,deltaT):
    """
    Module to compute DMD modes
    Args: 
        - Usnp      = snapshots matrix Usnp
        - mvect     = array with square root of mass weights (nGLLe*ncomponents)
        - nsnap      = number of snapshots
        - ifsym      = to mirror data wrt symmetry x axis
        - r          = reduced rank number
        - deltaT     = timestep
    Returns:
        - Phi = DMD Modes
        - Lambdat = DMD Eigenvalues
        - a1 = DMD Coefficients
        
    """    

    # Snapshot pairs
    X1 = Usnp[:,0:-1]
    X2 = Usnp[:,1:]

    # SVD of X1: 
    #   - U    = matrix whose columns are the spatial modes in POD expansion
    #   - Sigma = diagonal matrix of singular values for Usnp
    #   - Vt    = conjugate transpose of matrix with temporal coefficients
    U,Sigma,Vt = la.svd(X1,full_matrices=False)
    V = hermitian(Vt)

    # Rank r Reduction
    U = U[:,:r]
    Sigma = np.diag(Sigma[:r])
    V = V[:,:r]

    # Computation of Atilde
    Atilde = hermitian(U) @ X2 @ V @ la.inv(Sigma)
    
    # Computation of Eigenvalues and Eigenmodes
    Lambdat, W = la.eig(Atilde)
    Phi = X2 @ V @ la.inv(Sigma) @ W
    omega = np.log(Lambdat)/deltaT ### continuous time eigenvalues

    # Normalize the DMD Modes
    Phi = Phi / np.linalg.norm(Phi, axis=0)

    # DMD Coefficients
    a1 = la.lstsq(Phi, Usnp[:,0], rcond=None)[0]

    # Reordering eigenvalues and modes based on mode energies

    mode_energy = np.abs(a1)**2
    idx = np.flip(np.argsort(mode_energy))
    Lambdat = Lambdat[idx]
    omega = omega[idx]
    Phi = Phi[:, idx]

    return Phi, Lambdat, a1, omega