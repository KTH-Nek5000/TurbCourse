"""
The code implements relevant data-driven techniques with I/O designed for Nek5000.
It can handle 2D/3D cases and impose symmetry condition, if required.
So far, the following methods are available: 
    - Proper Orthogonal Decomposition (POD)
    - Dynamic Mode Decomposition (DMD)


@authors: Daniele Massaro, Saleh Rezaeiravesh and Philipp Schlatter.
dmassaro@kth.se, salehr@kth.se, pschlatt@mech.kth.se 
SimEx/FLOW, Engineering Mechanics, KTH Royal Institute of Technology, Stockholm, Sweden, 2021
"""

import sys
# Local modules
sys.path.append('$PATH/MODULES/')
from snapMaker   import snpAssembler,snpAssembler_sym
from reader      import read_input
from dbMaker     import dbMan,dbMan_sym
from PODmodule   import POD
from DMDmodule   import DMD
from plotter     import pplot,pplotc,savetxt
from outpWriter  import modes,snaprcn,prdct
from pickManager import pickReader,pickReader_sym


# ------- GENERAL NOTES:
#   - nGLL   = number of GLL points (=n.ofelements * lx1*ly1*lz1)
#   - n      = nGLL * number of variables (e.g. 3 variables if we consider u,v,w vel. comp.s)
#   - nsnap  = number of snapshots 
#   - Usnp   = snapshots matrix with dimensions: {n,snap} (or {n,snap*2} imposing the symmetry)
#   - Mode 0 = according to our convention the mode 0 corresponds to mean value




# ------- READING SETTINGS 
qoiName,nsnap,nplt,r,timeprdc,       \
outMod,outSnp,maxMode,               \
if3D,ifsym,ifPickSave,ifPickRead,    \
info,info_m,info_s  = read_input('$PATH/input.txt')



              
# ------- DATABASE/MATRIX DEFINITION 
if ifsym:
    
    if ifPickRead:    # reading pickle 
        Usnp,db,db_s,mvect = pickReader_sym(info,info_s)
    else:             # building data
        db,db_m,db_s,data_ms = dbMan_sym(info,'field',info_m,'mass',info_s,'field')       # Database generation
        Usnp,mvect = snpAssembler_sym(db,db_s,data_ms,info,nsnap,ifsym,if3D,ifPickSave)   # Snapshots matrix assembly
        
else:
    
    if ifPickRead:    # reading pickle 
        Usnp,db,mvect = pickReader(info)
    else:             # building data
        db,db_m,data_ms = dbMan(info,'field',info_m,'mass')                                # Database generation
        Usnp,mvect = snpAssembler(db,data_ms,info,nsnap,ifsym,if3D,ifPickSave)             # Snapshots matrix assembly     
    



# ------- CALCULATION
if (info['module']=='POD'):
    # POD
    L,Lam2,R,A2 = POD(Usnp,mvect,nsnap,ifsym)
elif (info['module']=='DMD'): 
    # DMD
    Phi,Lambdat,a1 = DMD(Usnp,mvect,nsnap,r,ifsym)




# ------- PLOTS
if (info['module']=='POD'):
    pplot(A2,Lam2,nplt,'coeff')
    pplot(A2,Lam2,nplt,'eigen')
    savetxt(Lam2,info)
elif (info['module']=='DMD'): 
    pplotc(Lambdat)
    savetxt(Lambdat,info)



# ------- OUTPUT
if (info['module']=='POD'):
    modes(outMod,db['data'][0],info,L,if3D)
    snaprcn(outSnp,db['data'][0],info,L,A2,maxMode,if3D)
elif (info['module']=='DMD'): 
    modes(outMod,db['data'][0],info,Phi.real,if3D)
    prdct(timeprdc,db['data'][0],info,Phi,a1,Lambdat,if3D)
    
    
