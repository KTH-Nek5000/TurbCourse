##############################################################
# Generate coordinates of the elements for the mixing layer
##############################################################
import numpy as np
import matplotlib.pyplot as plt
from math import pi
#
#--- SETTINGS
# x-dir
Ex=30   #number of elements in x
xSt=0.0
xEnd=10.0
# y-dir
Ey=20   #number of elements in y, !Even number
gam=0.45   #should be <= 0.5 (smaller gamma, higher compression towards the center)
ySt=0.0
yEnd=10
#------------------------------
#
#
# elements in y
tmp=np.linspace(0,1,(Ey+2)/2)
y=1-np.cos(gam*tmp*pi)
yMin=np.min(y)
yMax=np.max(y)
y=(y-yMin)/(yMax-yMin)
#other half
z=-np.flip(y)
z=np.delete(z,-1)
#concatenate both halves
y=np.concatenate((z,y))  #\in [-1,1]
y=0.5*(y+1)*(yEnd-ySt)+ySt

#plot
plt.figure(figsize=(14,2))
plt.plot(y,np.ones(y.shape[0]),'ob',mfc='none')
plt.xlabel('Elements y-coordinates')
plt.grid(alpha=0.2)
plt.show()



#write
F=open('mixlay.box','w')
F.write("-2 \t\t spatial dimension (<0 => .re2 and .rea files)\n")
F.write("1  \t\t number of fields\n")
F.write("box_1 \t\t name of the box\n")
F.write("-%d   %d  1 \t nelx, nely, nelz (<0 => equally divided)\n" %(Ex,y.shape[0]-1))
F.write("%.2F   %.2F  1    xStart, xEnd, gain\n" %(xSt,xEnd))
for y_ in y:
    F.write("%.2F  " %y_)
F.write("\t y1,y2,...\n")    

F.write("v  ,o  ,ON ,ON  bcs (west, east, south, north)")



