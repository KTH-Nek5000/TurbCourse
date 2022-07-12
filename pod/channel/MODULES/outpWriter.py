"""
The module saves a selected number of POD/DMD modes or 
reconstructed fields in Nek5000 data format.
Predicted fields can also be saved.
"""

from pymech import neksuite
import numpy as np





def modes(nmod,db,info,L,if3D):
    """
    Module to save a selected number of POD/DMD modes in Nek5000 data format.
    Use Visit/Paraview to visualise it. 
    Args: 
        - nmod  = largest mode to save. ex: nmod=0-->0th mode, nmod=2-->0th,1st,2nd modes
        - db    = database of a single snapshots which contains elements & GLL points 
        - info  = general info of the case
        - L     = POD modes matrix
        - if3D  = to identify 3D or 2D case
        
    """    
    
    for m in range(nmod+1):
        
        data = db
     
        nGLLe = data.lr1[0]*data.lr1[1]*data.lr1[2]
        # Create data structure to outpost file equal to the original (data_o)
        # Save u,v (and w if 3D) in the u,v,w velocity componets field
        data_o = data           


        k = 0                                                # index to loop on u-block
        k1 = nGLLe*data.nel    # index to loop on v-block
        k2 = k1*2                                            # index to loop on w-block
        for i in range(data.nel): 
            
            if (info['variable'] == 'scalar'):
                tmp1 = L[k:k+nGLLe,m]
                tmp1 = np.reshape(tmp1,(data.lr1[2],data.lr1[1],data.lr1[0]),order='F')
                
            elif (info['variable'] == 'vector'):            
                tmp1 = L[k:k+nGLLe,m]
                tmp1 = np.reshape(tmp1,(data.lr1[2],data.lr1[1],data.lr1[0]),order='F')    
                tmp2 = L[k1:k1+nGLLe,m]
                tmp2 = np.reshape(tmp2,(data.lr1[2],data.lr1[1],data.lr1[0]),order='F')
            
                if if3D:
                    tmp3 = L[k2:k2+nGLLe,m]
                    tmp3 = np.reshape(tmp3,(data.lr1[2],data.lr1[1],data.lr1[0]),order='F')              
    
            if (info['qoiName'] == 'temperature'):   
                data_o.elem[i].temp[0,:,:,:] = tmp1 
            elif (info['qoiName'] == 'pressure'):
                data_o.elem[i].pres[0,:,:,:] = tmp1 
            elif (info['qoiName'] == 'uvel'):
                data_o.elem[i].vel[0,:,:,:] = tmp1  # COMMENT ABOUT OTHER VARIABLES SAVED IN THE FIELDS
            elif (info['qoiName'] == 'vvel'):
                data_o.elem[i].vel[1,:,:,:] = tmp1
            elif (info['qoiName'] == 'wvel'):
                data_o.elem[i].vel[2,:,:,:] = tmp1
            elif (info['qoiName'] == 'velocity'):
                data_o.elem[i].vel[0,:,:,:] = tmp1    
                data_o.elem[i].vel[1,:,:,:] = tmp2 
                if if3D:
                    data_o.elem[i].vel[2,:,:,:] = tmp3    
                
            # Indices update 
            k = k + nGLLe
            k1 = k1 + nGLLe
            k2 = k2 + nGLLe

        if (info['module']=='POD'):
            # Writing data in nek/visit format
            data_o.time = m

            neksuite.writenek(info['outputPath']+'PODmod'+info['caseName']+'0.f'+str(m).zfill(5),data_o)
            print('POD mode number %d has been saved' % (m))
        elif (info['module']=='DMD'): 
            # Writing data in nek/visit format
            neksuite.writenek(info['outputPath']+'DMDmod'+info['caseName']+'0.f'+str(m).zfill(5),data_o)
            print('DMD mode number %d has been saved' % (m))

    print('Writing: '+info['outputPath']+'PODmod'+info['caseName']+'.nek5000')
    with open(info['outputPath']+'PODmod'+info['caseName']+'.nek5000', "w") as f:
       f.write('filetemplate: PODmod' + info['caseName']+'%01d.f%05d\n')
       f.write('firsttimestep: 0\n')
       f.write('numtimesteps: %d\n' % (nmod+1))
    
    return


def snaprcn(nmod,db,info,L,A2,maxMode,if3D):
    """
    Module to save a selected number of POD modes in Nek5000 data format.
    Use Visit/Paraview to visualise it. 
    Args: 
        - nmod    = largest mode to save. ex: nmod=0-->0th mode, nmod=2-->0th,1st,2nd modes
        - db      = database which contains elements & GLL points info for each snapshot 
        - info    = general info of the case
        - L       = POD modes matrix
        - A2      = POD coefficients
        - maxMode = number of modes to use for reconstructing the flow
        - if3D    = to identify 3D or 2D case
    Returns:
        
    """    
    
    uRcn = L[:,:maxMode+1] @ A2[:,:maxMode+1].T
    for m in range(nmod+1):
    
        data = db
        nGLLe = data.lr1[0]*data.lr1[1]*data.lr1[2]
        # Create data structure to outpost file equal to the original (data_o)
        # Save u,v (and w if 3D) in the u,v,w velocity componets field
        data_o = data           
        k = 0                                                # index to loop on u-block
        k1 = nGLLe*data.nel    # index to loop on v-block
        k2 = k1*2                                            # index to loop on w-block
        for i in range(data.nel): 
        
            if (info['qoiName'] == 'temperature'):
                tmp1 = uRcn[k:k+nGLLe,m]
                tmp1 = np.reshape(tmp1,(data.lr1[2],data.lr1[1],data.lr1[0]),order='F')
                
            elif (info['qoiName'] == 'velocity'):            
                tmp1 = uRcn[k:k+nGLLe,m]
                tmp1 = np.reshape(tmp1,(data.lr1[2],data.lr1[1],data.lr1[0]),order='F')    
                tmp2 = uRcn[k1:k1+nGLLe,m]
                tmp2 = np.reshape(tmp2,(data.lr1[2],data.lr1[1],data.lr1[0]),order='F')
            
                if if3D:
                    tmp3 = uRcn[k2:k2+nGLLe,m]
                    tmp3 = np.reshape(tmp3,(data.lr1[2],data.lr1[1],data.lr1[0]),order='F')              
    
            if (info['qoiName'] == 'temperature'):
                data_o.elem[i].temp[0,:,:,:] = tmp1 
            elif (info['qoiName'] == 'velocity'):
                data_o.elem[i].vel[0,:,:,:] = tmp1    
                data_o.elem[i].vel[1,:,:,:] = tmp2 
                if if3D:
                    data_o.elem[i].vel[2,:,:,:] = tmp3    
                
            # Indices update 
            k = k + nGLLe
            k1 = k1 + nGLLe
            k2 = k2 + nGLLe



        # Writing data in nek/visit format
        neksuite.writenek(info['outputPath']+'PODsnaprcn'+info['caseName']+'0.f'+str(m).zfill(5),data_o)
        print('Snapshot reconstructed  number %d has been saved' % (m))
     
    
    return





def prdct(tt,db,info,Phi,a1,Lambdat,if3D):
    """
    Module to save a predicted fields in Nek5000 data format.
    Use Visit/Paraview to visualise it. 
    Args: 
        - tt       = time instant to whom the system state is predicted 
        - db       = database which contains elements & GLL points info for each snapshot 
        - info     = general info of the case
        - Phi      = DMS modes 
        - a1       = DMD coefficients
        - Lambdat  = DMD eigenvalues
        - if3D     = to identify 3D or 2D case
    Returns:
        
    """    
    
    
    #Data driven spectral decomposition
    #k < m
    uPred1=Phi @ np.multiply(Lambdat**tt , a1)   #np.multiply() for elementwise multiplication
    #uPred1=Phi @ np.diag(Lambda**k) @ a1  #alternative implementation
    
    #uPred2=Phi @ np.multiply(Lambda**k , a2)
    #uPred2=Phi @ np.diag(Lambda**k) @ a2   #alternative implementation
    
    Uprd = uPred1.real
    

    m = 0 
    data = db
    nGLLe = data.lr1[0]*data.lr1[1]*data.lr1[2]
    # Create data structure to outpost file equal to the original (data_o)
    # Save u,v (and w if 3D) in the u,v,w velocity componets field
    data_o = data           
    k = 0                                                # index to loop on u-block
    k1 = nGLLe*data.nel    # index to loop on v-block
    k2 = k1*2                                            # index to loop on w-block
    for i in range(data.nel): 
    
        if (info['qoiName'] == 'temperature'):
            tmp1 = Uprd[k:k+nGLLe]
            tmp1 = np.reshape(tmp1,(data.lr1[2],data.lr1[1],data.lr1[0]),order='F')
            
        elif (info['qoiName'] == 'velocity'):            
            tmp1 = Uprd[k:k+nGLLe]
            tmp1 = np.reshape(tmp1,(data.lr1[2],data.lr1[1],data.lr1[0]),order='F')    
            tmp2 = Uprd[k1:k1+nGLLe]
            tmp2 = np.reshape(tmp2,(data.lr1[2],data.lr1[1],data.lr1[0]),order='F')
        
            if if3D:
                tmp3 = Uprd[k2:k2+nGLLe]
                tmp3 = np.reshape(tmp3,(data.lr1[2],data.lr1[1],data.lr1[0]),order='F')              

        if (info['qoiName'] == 'temperature'):
            data_o.elem[i].temp[0,:,:,:] = tmp1 
        elif (info['qoiName'] == 'velocity'):
            data_o.elem[i].vel[0,:,:,:] = tmp1    
            data_o.elem[i].vel[1,:,:,:] = tmp2 
            if if3D:
                data_o.elem[i].vel[2,:,:,:] = tmp3    
            
        # Indices update 
        k = k + nGLLe
        k1 = k1 + nGLLe
        k2 = k2 + nGLLe



    # Writing data in nek/visit format
    neksuite.writenek(info['outputPath']+'DMDsnaprdct'+info['caseName']+'0.f'+str(m).zfill(5),data_o)
    print('System state predicted at time instant %d has been saved' % (tt))
     
    
    return
