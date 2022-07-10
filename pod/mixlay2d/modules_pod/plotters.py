#################################################
# A set of plotters for POD
#################################################
# Saleh Rezaeiravesh, salehr@mech.kth.se
#------------------------------------------------
#
import numpy as np
import matplotlib.pyplot as plt
#
class contour2d:
    """
    Contour plot of `f` in x-y plane
      `f` is a numpy array of shape (nx*ny,)
      `x`: numpy vector of length nx
      `y`: numpy vector of length ny
    """
    def __init__(self,f,x,y,figsize=(10,5),title='',cmap='seismic'):
        self.f=f
        self.x=x
        self.y=y
        self.nx=self.x.shape[1]
        self.ny=self.x.shape[0]
        self.figsize=figsize
        self.title=title
        self.cmap=cmap
        self._plot()

    def _plot(self):
        f_=np.reshape(self.f,(self.nx,self.ny),order='F')
        plt.figure(figsize=self.figsize)
        plt.contourf(self.x,self.y,f_.T,levels=100,cmap=self.cmap)
        if len(self.title)>0:
           plt.title(self.title) 
        plt.colorbar()
        plt.show()
#
class scatter:
    """
    Scatter plot of `data_x` and `data_y`
    """
    def __init__(self,data_x,data_y,xlab='',ylab='',figsize=(5,5),title=''):
        self.f=data_x
        self.g=data_y
        self.figsize=figsize
        self.title=title
        self.xlab=xlab
        self.ylab=ylab
        self._plot()

    def _plot(self):
        plt.figure(figsize=self.figsize)
        plt.scatter(self.f,self.g,alpha=0.2,s=1)
        min_=min(np.min(self.f),np.min(self.g))*0.9
        max_=max(np.max(self.f),np.max(self.g))*1.1
        plt.plot(np.linspace(min_,max_,100),np.linspace(min_,max_,100),'-r',lw=2)
        plt.xlabel(self.xlab)
        plt.ylabel(self.ylab)
        if len(self.title)>0:
           plt.title(self.title)        
        plt.show()
        acorr2=np.corrcoef(self.f,self.g)[0,1]
        print('Correlation between the observation and POD reconstruction = ',acorr2)
#
