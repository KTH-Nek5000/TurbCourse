"""
The module assembles the snapshot matrix
"""

from sys import exit
import numpy as np
from pickManager import pickWriter,pickWriter_sym




def snpAssembler_sym(db,db_s,dat_m,info,nsnap,ifsym,if3D,ifPickSave):
    """
    Module to assemble snapshots matrix
    Args: 
        - db         = database which contains elements & GLL points info for each snapshot 
        - db_s       = database which contains elements & GLL points info for each mirrored snapshot 
        - dat_m      = exadata  which contains elements & GLL points info for the mass matrix
        - info       = general info of the case
        - nsnap      = number of snapshots
        - ifsym      = to mirror data wrt symmetry x axis
        - if3D       = to identify 3D or 2D case
        - ifPickSave = to save Pickle file
    Returns:
        - Snapshots matrix Usnp 
        - array with square of mass weights (nGLLe*ncomponents)
        
    """ 

    # Initialisation module
    if (info['variable'] == 'scalar'):    
        Usnp,tmpi,tmpi_s,mvect = init(info,nsnap,db,ifsym,if3D)
    elif (info['variable'] == 'vector'):
        if if3D:
            Usnp,Uu,Vv,Ww,Uu_m,Vv_m,Ww_m,mvect = init(info,nsnap,db,ifsym,if3D)
            # if ifsym=False then Uu_m,Vv_m,Ww_m are empty 
        else:
            Usnp,Uu,Vv,Uu_m,Vv_m,mvect = init(info,nsnap,db,ifsym,if3D)
            # if ifsym=False then Uu_m,Vv_m are empty 
    
    # Number of GLL points in each element
    nGLLe = db['data'][0].lr1[0]*db['data'][0].lr1[1]*db['data'][0].lr1[2]
      
    # Reading data, Mass weighting, Snapshot matrix construction (ELEMENTWISE)
    for j in range(nsnap):                      # Loop snapshots
        # Looking at first snapshot of db
        data = db['data'][j] 
        k = 0            
        #If symmetry is applied, mirrored fields are read   
        data_s = db_s['data'][j]  
        
        for i in range(data.nel):               # Loop elements

            # Mass matrix weights (same for all snapshots)
            if (j==0):
               tmpm = np.sqrt(dat_m.elem[i].vel[0])   
               mvect[k:k+nGLLe] = np.reshape(tmpm,(1,nGLLe),order='F')
            
            # Scalar: temperature
            if (info['qoiName']  == 'temperature'):
                # Temperature
                tmp = data.elem[i].temp    # tmp is a working array
                if if3D:
                    tmpt = tmp[0]
                else:
                    tmpt = tmp[0,0]
                tmpi[k:k+nGLLe]  = np.reshape(tmpt,(1,nGLLe),order='F')
                # Symmetry
                tmp_s = data_s.elem[i].temp    
                if if3D:
                    tmpt_s = tmp_s[0]   
                else:
                    tmpt_s = tmp_s[0,0]
                tmpi_s[k:k+nGLLe]  = np.reshape(tmpt_s,(1,nGLLe),order='F')
                        
            # Scalar: pressure
            elif (info['qoiName']  == 'pressure'):
                
                tmp = data.elem[i].pres    # tmp is a working array
                if if3D:
                    tmpp = tmp[0]
                else:
                    tmpp = tmp[0,0]
                tmpi[k:k+nGLLe] = np.reshape(tmpp,(1,nGLLe),order='F')
                # Symmetry
                tmp_s = data_s.elem[i].pres    
                if if3D:
                    tmpp_s = tmp_s[0]
                else:
                    tmpp_s = tmp_s[0,0]
                tmpi_s[k:k+nGLLe]  = np.reshape(tmpp_s,(1,nGLLe),order='F')
                        
            # Scalar: u velocity component
            elif (info['qoiName']  == 'uvel'):
                
                tmp = data.elem[i].vel 
                tmpu = tmp[0]
                # Symmetry
                tmp_s = data_s.elem[i].vel    
                tmpu_s = tmp_s[0] 
                tmpi_s[k:k+nGLLe] = np.reshape(tmpu_s,(1,nGLLe),order='F')
                    
            # Scalar: v velocity component
            elif (info['qoiName']  == 'vvel'):
                
                tmp = data.elem[i].vel 
                tmpv = tmp[1]
                tmpi[k:k+nGLLe] = np.reshape(tmpv,(1,nGLLe),order='F')         
                # Symmetry                
                tmp_s = data_s.elem[i].vel    
                tmpv_s = tmp_s[1] 
                tmpi_s[k:k+nGLLe] = np.reshape(tmpv_s,(1,nGLLe),order='F')
                    
            # Scalar: w velocity component
            elif (info['qoiName']  == 'wvel'):
                
                if (if3D):               
                    tmp = data.elem[i].vel 
                    tmpw= tmp[2]
                    tmpi[k:k+nGLLe] = np.reshape(tmpw,(1,nGLLe),order='F') 
                    # Symmetry
                    tmp_s = data_s.elem[i].vel    
                    tmpw_s = tmp_s[2] 
                    tmpi_s[k:k+nGLLe] = np.reshape(tmpw_s,(1,nGLLe),order='F')
                else:
                    print('Set if3D=True to read w velocity component')
                    exit(0)
                          

            # Vector: velocity
            elif (info['qoiName']  == 'velocity'):    
                
                # velocity components
                tmp = data.elem[i].vel      
                
                tmpu = tmp[0]     
                tmpv = tmp[1] 

                # Rearrange 
                tmpu = np.reshape(tmpu,(1,nGLLe),order='F')
                tmpv = np.reshape(tmpv,(1,nGLLe),order='F')              
                
                Uu[k:k+nGLLe] = tmpu
                Vv[k:k+nGLLe] = tmpv
                                                   
                if if3D:
                    tmpw = tmp[2] 
                    tmpw = np.reshape(tmpw,(1,nGLLe),order='F')
                    Ww[k:k+nGLLe] = tmpw
                # Symmetry
                tmp_s = data_s.elem[i].vel     
                tmpu_s = tmp_s[0]     
                tmpv_s = tmp_s[1] 

                # Rearrange 
                tmpu_s = np.reshape(tmpu_s,(1,nGLLe),order='F')
                tmpv_s = np.reshape(tmpv_s,(1,nGLLe),order='F')

                tmpu_s=tmpu_s[0]
                tmpv_s=tmpv_s[0]
                Uu_m[k:k+nGLLe] = tmpu_s 
                Vv_m[k:k+nGLLe] = tmpv_s 
                                                   
                if if3D:
                    tmpw_s = tmp_s[2] 
                    tmpw_s = np.reshape(tmpw_s,(1,nGLLe),order='F')
                    Ww_m[k:k+nGLLe] = tmpw_s

            k = k + nGLLe
    
        k1 =  data.nel*nGLLe
        k2 =  k1*2  

        # Mass matrix weights      
        if (j==0):    
            if (info['variable']  == 'vector'):
                mvect[k1:k2]   = mvect[0:k1]
                if if3D:
                    mvect[k2:k1*3] = mvect[0:k1]

            
        if (info['variable']  == 'scalar'):
            Usnp[0:k1,j] = tmpi
            # SYMMETRY
            Usnp[0:k1,j+nsnap] = tmpi_s
                
        elif (info['variable']  == 'vector'):
            # Velocity components concatenation
            Usnp[0:k1,j] = Uu
            Usnp[k1:k2,j] = Vv  
            if if3D:
                Usnp[k2:k1*3,j] = Ww
    
            # SYMMETRY: mirrored fields added as additional fields           
            Usnp[0:k1,j+nsnap] = Uu_m
            Usnp[k1:k2,j+nsnap]  = Vv_m
            if if3D:
                Usnp[k2:k1*3,j+nsnap] = Ww_m
              
    pickWriter_sym(info,Usnp,db,db_s,mvect,ifPickSave)             

    return Usnp,mvect




def snpAssembler(db,dat_m,info,nsnap,ifsym,if3D,ifPickSave):
    """
    Module to assemble snapshots matrix
    Args: 
        - db         = database which contains elements & GLL points info for each snapshot 
        - dat_m      = exadata  which contains elements & GLL points info for the mass matrix
        - info       = general info of the case
        - nsnap      = number of snapshots
        - ifsym      = to mirror data wrt symmetry x axis
        - if3D       = to identify 3D or 2D case
        - ifPickSave = to save Pickle file
    Returns:
        - Snapshots matrix Usnp 
        - array with square of mass weights (nGLLe*ncomponents)
        
    """ 

    # Initialisation module
    if (info['variable'] == 'scalar'):    
        Usnp,tmpi,tmpi_s,mvect = init(info,nsnap,db,ifsym,if3D)
    elif (info['variable'] == 'vector'):
        if if3D:
            Usnp,Uu,Vv,Ww,Uu_m,Vv_m,Ww_m,mvect= init(info,nsnap,db,ifsym,if3D)
            # if ifsym=False then Uu_m,Vv_m,Ww_m are empty 
        else:
            Usnp,Uu,Vv,Uu_m,Vv_m,mvect = init(info,nsnap,db,ifsym,if3D)
            # if ifsym=False then Uu_m,Vv_m are empty 
    
    # Number of GLL points in each element
    nGLLe = db['data'][0].lr1[0]*db['data'][0].lr1[1]*db['data'][0].lr1[2]
    
    
    # Reading data, Mass weighting, Snapshot matrix construction (ELEMENTWISE)
    for j in range(nsnap):                      # Loop snapshots
        # Looking at first snapshot of db
        data = db['data'][j] 
        k = 0            

        for i in range(data.nel):               # Loop elements   

            # Mass matrix weights (same for all snapshots)
            if (j==0):
               tmpm = np.sqrt(dat_m.elem[i].vel[0])   
               mvect[k:k+nGLLe] = np.reshape(tmpm,(1,nGLLe),order='F')        
            
            # Scalar: temperature
            if (info['qoiName']  == 'temperature'):
                # Temperature
                tmp = data.elem[i].temp    # tmp is a working array
                if if3D:
                    tmpt = tmp[0]
                else:
                    tmpt = tmp[0,0]
                tmpi[k:k+nGLLe]  = np.reshape(tmpt,(1,nGLLe),order='F')
                               
            # Scalar: pressure
            elif (info['qoiName']  == 'pressure'):
                
                tmp = data.elem[i].pres    # tmp is a working array
                if if3D:
                    tmpp = tmp[0]
                else:
                    tmpp = tmp[0,0]
                tmpi[k:k+nGLLe] = np.reshape(tmpp,(1,nGLLe),order='F')
                               
            # Scalar: u velocity component
            elif (info['qoiName']  == 'uvel'):
                
                tmp = data.elem[i].vel 
                tmpu = tmp[0]
                tmpi[k:k+nGLLe] = np.reshape(tmpu,(1,nGLLe),order='F')                      
                    
            # Scalar: v velocity component
            elif (info['qoiName']  == 'vvel'):
                
                tmp = data.elem[i].vel 
                tmpv = tmp[1] 
                tmpi[k:k+nGLLe] = np.reshape(tmpv,(1,nGLLe),order='F')                      
                    
            # Scalar: v velocity component
            elif (info['qoiName']  == 'wvel'):
                
                if (if3D):               
                    tmp = data.elem[i].vel 
                    tmpw= tmp[2]
                    tmpi[k:k+nGLLe] = np.reshape(tmpw,(1,nGLLe),order='F') 
                              
                else:
                        print('Set if3D=True to read w velocity component')
                        exit(0)
                          
            # Vector: velocity
            elif (info['qoiName']  == 'velocity'):    
                
                # velocity components
                tmp = data.elem[i].vel      
                
                tmpu = tmp[0]     
                tmpv = tmp[1] 

                # Rearrange 
                tmpu = np.reshape(tmpu,(1,nGLLe),order='F')
                tmpv = np.reshape(tmpv,(1,nGLLe),order='F')
           
                Uu[k:k+nGLLe] = tmpu
                Vv[k:k+nGLLe] = tmpv
                                                   
                if if3D:
                    tmpw = tmp[2] 
                    tmpw = np.reshape(tmpw,(1,nGLLe),order='F')
                    Ww[k:k+nGLLe] = tmpw           

            k = k + nGLLe
    
        k1 =  data.nel*nGLLe
        k2 =  k1*2  
        
        # Mass matrix weights      
        if (j==0):    
            if (info['variable']  == 'vector'):
                mvect[k1:k2]   = mvect[0:k1]
                if if3D:
                    mvect[k2:k1*3] = mvect[0:k1]
        
        if (info['variable']  == 'scalar'):
            Usnp[0:k1,j] = tmpi
                
        elif (info['variable']  == 'vector'):
            # Velocity components concatenation
            Usnp[0:k1,j] = Uu
            Usnp[k1:k2,j] = Vv  
            if if3D:
                Usnp[k2:k1*3,j] = Ww
    
    pickWriter(info,Usnp,db,mvect,ifPickSave)                

    return Usnp,mvect




def init(info,nsnap,db,ifsym,if3D):
    """
    Module to initialise working array 
    Args: 
        - info  = general info of the case
        - nsnap = number of snapshots
        - db    = database which contains elements & GLL points info for each snapshot 
        - mvect = mass weights for each snapshot (array with dims nGLLe*ncomponents)
        - ifsym = to mirror data wrt symmetry x axis
        - if3D  = to identify 3D or 2D case
    Returns:
        - if variable=scalar   --> matrix snapshots (Usnp) + working arrays
        - if variable=vector   --> matrix snapshots (Usnp) + working arrays
        
    """ 

    # Mesh info   
    # Look at first snapshot to extract mesh info
    data = db['data'][0]
    nGLLe = db['data'][0].lr1[0]*db['data'][0].lr1[1]*db['data'][0].lr1[2]
    
    # Snapshots matrix
    if (info['variable'] == 'scalar'):    
        if ifsym:
            Usnp =  np.single(np.zeros((data.nel*nGLLe,nsnap*2))) # np.zeros((data.nel*nGLLe,nsnap))
        else:
            Usnp =  np.single(np.zeros((data.nel*nGLLe,nsnap))) 
        mvect=  np.single(np.zeros((data.nel*nGLLe)))
        # Working vectors to storescalar GLL data
        tmpi   = np.single(np.zeros((data.nel*nGLLe)))
        tmpi_s = np.single(np.zeros((data.nel*nGLLe)))
        
    elif (info['variable'] == 'vector'):
        if ifsym:
            Usnp =  np.single(np.zeros((data.nel*nGLLe*data.ndim,nsnap*2)))
        else:
            Usnp =  np.single(np.zeros((data.nel*nGLLe*data.ndim,nsnap)))
        mvect =  np.single(np.zeros((data.nel*nGLLe*data.ndim)))
        # Working vectors to store velocity GLL data
        Uu = np.single(np.zeros((data.nel*nGLLe)))
        Vv = np.single(np.zeros((data.nel*nGLLe)))
        # SYMMETRY VARIABLES
        # Working vectors to store velocity GLL data for additional mirrored fields
        Uu_m = np.single(np.zeros((data.nel*nGLLe)))
        Vv_m = np.single(np.zeros((data.nel*nGLLe)))
        
        if if3D:
            Ww   = np.single(np.zeros((data.nel*nGLLe)))
            Ww_m = np.single(np.zeros((data.nel*nGLLe)))
        
    # RETURNS
    if (info['variable'] == 'scalar'):    
        return Usnp,tmpi,tmpi_s,mvect
    elif (info['variable'] == 'vector'):
        if if3D:
            return Usnp,Uu,Vv,Ww,Uu_m,Vv_m,Ww_m,mvect
        else:
            return Usnp,Uu,Vv,Uu_m,Vv_m,mvect