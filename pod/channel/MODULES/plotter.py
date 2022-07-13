"""
The module plots the spectrum of POD/DMD modes
and saves the eigenvalues. 
"""





import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
params = {'legend.fontsize': 15,
          'legend.loc':'best',
          'figure.figsize': (15, 5),
          'lines.markerfacecolor':'none',
          'axes.labelsize': 17,
          'axes.titlesize': 17,
          'xtick.labelsize':15,
          'ytick.labelsize':15,
          'grid.alpha':0.6}
pylab.rcParams.update(params)





def pplot(A2,Lam2,nplt,iff):
    """
    Module to plot eigenvalues spectrum and cumulative sum
    Args: 
        - A2   = POD coefficients
        - Lam2 = matrix of singular values of Usnp (normalized vector)
        - nplt = number of eigenvalues to plot
        - iff  = flag to plot coeff.s or eigen.s
    """        
    
    if (iff=='coeff'):
    
        # Coefficient matrix
        plt.pcolor(A2[:,0:nplt])
        plt.title('Time coefficients')
        plt.colorbar()
        plt.xlabel('mode')
        plt.ylabel('step')
        plt.clim(-0.1,0.1)
        plt.show()

        plt.plot(A2[:,3])
        plt.title('mode 3')
        plt.xlabel('step')
        plt.ylabel('amplitude')
        plt.show()
        
    elif (iff=='eigen'):   
        # Eigenvalues
        print('singular values:',Lam2[0:nplt])
        plt.plot(np.arange(nplt),np.log10((Lam2[0:nplt])/sum(Lam2[0:nplt])),'^--r',label=r'$\lambda_k/\sum_{i=0}^m{\lambda_i}$')
        # plt.plot(np.cumsum(Lam2[0:nplt])/sum(Lam2[0:nplt]),'o-b',label='$\sum_{i=1}^k\lambda_i/\sum_{i=1}^m{\lambda_i}$')
        plt.xlabel(r'$k$')
        plt.ylabel('log10 energy')
        plt.legend(loc='best')
        plt.grid()
        plt.show()

        # here we only show the cummulative sum without the mean
        plt.plot(np.arange(nplt-1)+1,np.cumsum(Lam2[1:nplt])/sum(Lam2[1:nplt]),'o-b',label='$\sum_{i=1}^k\lambda_i/\sum_{i=1}^m{\lambda_i}$')
        plt.xlabel(r'$k$')
        plt.legend(loc='best')
        plt.grid()
        plt.show()
    
    return 



def pplotc(Lambdat):
    """
    Module to plot eigenvalues spectrum 
    Args: 
        - Lambdat = DMD eigenvalues
    """        
    #plot eigenvalues
    plt.figure(figsize=(7,7))
    plt.plot(Lambdat.real,Lambdat.imag,'ob',ms=8)
    thet=np.linspace(0,2*np.pi,100)
    plt.plot(np.cos(thet),np.sin(thet),'-k')
    for i in range(len(Lambdat)):
        if abs(Lambdat[i])<1:
           plt.plot(Lambdat[i].real,Lambdat[i].imag,'ob',mfc='b',label=r'$|\lambda|<1$') 
        elif abs(Lambdat[i])>1:
           plt.plot(Lambdat[i].real,Lambdat[i].imag,'or',mfc='r',label=r'$|\lambda|>1$') 
    
    plt.xlabel(r'$Re(\tilde{\Lambda})$')
    plt.ylabel(r'$Im(\tilde{\Lambda})$')
    plt.grid(alpha=0.4)
    plt.title(r'Eigenvalues $\tilde{\mathbf{\Lambda}}$')
    plt.show()    
  
    return 



def savetxt(Lam2,info):
    """
    Module to save eigenvalues spectrum 
    Args: 
        - Lam2 = POD/DMD eigenvalues
        - info = case info
    """        
    if (info['module']=='POD'):
        
        fileName = info['outputPath']+'PODeigns.txt' #'./OUTPUT/eigns.txt'
        F = open(fileName,'w') 
        F.write("# Snapshots POD result \n") 
        F.write("# Eigenvalues computed via SVD \n") 
        F.write("# According to our convention the mode 0 is the mean value \n") 
        F.write("# ------------------------------------------------------------------\n") 
        F.write("# Sum of all eigenvalues   = %g \n" % sum(Lam2)) 
        F.write("# Number of snapshots      = %g \n" % np.size(Lam2)) 
        F.write("# ------------------------------------------------------------------\n") 
        F.write("# Eigenvalues\t Eigenvalues/Sum\t  Cumulative sum\t  \n") 
        F.write("# ------------------------------------------------------------------\n") 
        for i in range(0,np.size(Lam2)):
           F.write("%g\t%g\t%g\t \n" % \
                   (Lam2[i], Lam2[i]/sum(Lam2), np.cumsum(Lam2)[i]/sum(Lam2) ))
        F.close()
        
    elif (info['module']=='DMD'): 
        
        fileName = info['outputPath']+'DMDeigns.txt' #'./OUTPUT/eigns.txt'
        F = open(fileName,'w') 
        F.write("# Snapshots DMD result \n") 
        F.write("# Eigenvalues of linear operator A: Ax=b \n") 
        F.write("# According to our convention the mode 0 is the mean value \n") 
        F.write("# ------------------------------------------------------------------\n") 
        F.write("# ------------------------------------------------------------------\n") 
        F.write("# Real{Eigenvalues}\t Im{Eigenvalues}\t  \n") 
        F.write("# ------------------------------------------------------------------\n") 
        for i in range(0,np.size(Lam2)):
           F.write("%g\t%g\t \n" % \
                   (Lam2[i].real, Lam2[i].imag ))
        F.close()
            
    return 
