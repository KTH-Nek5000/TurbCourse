"""
The module reads the parameters set in input.txt
It assumes that all the read snapshots are used (see start_ and end_ variables)
"""




def read_input(filename):
   params = {}
   with open(filename) as fpar:
       for line in fpar:
          line = line.strip()
          key_val = line.split("= ")
          if len(key_val) == 2:
             params[key_val[0].strip()] = key_val[1]
   # Convert dict
   mod        = params['module']
   path       = params['path_in']
   path_out   = params['path_out']
   caseName   = params['caseName']
   nsnap      = int(params['nsnap'])
   if3D       = eval(params['if3D'])
   qoiName    = params['qoiName']
   variable   = params['variable']
   ifsym      = eval(params['ifsym'])
   path_m     = params['path_m']
   caseName_m = params['caseName_m']
   nplt       = int(params['nplt'])
   outMod     = int(params['outMod'])
   outSnp     = int(params['outSnp'])
   maxMode    = int(params['maxMode'])
   ifPickSave = eval(params['ifPickSave'])
   ifPickRead = eval(params['ifPickRead'])
   r          = int(params['rpar'])
   timeprdc   = int(params['tprdc'])
   
  
   # All snapshots are used.
   start_ = timeprdc         # ID of the starting file, e.g. 5 -> caseName0.f00005
   end_   = timeprdc + nsnap     # ID of the end file

   # Info definition:
   # - snapshots
   info = {'module':mod,
           'dataPath':path,
           'outputPath':path_out,
           'caseName':caseName,
           'startID':start_,
           'endID':end_,
           'variable':variable,
           'qoiName':qoiName}
   # - mass matrix ( X-VELOCITY COMPONENT IS USED, SAVED HERE bm1 MATRIX)
   massName = 'bm1'+caseName
   info_m = {'dataPath':path,
             'caseName':massName}
   
   # Symmetric fields acquisition 
   if ifsym:
       info_s = {'dataPath':path_m,
                 'caseName':caseName_m,
                 'startID':start_,
                 'endID':end_,
                 'variable':variable,
                 'qoiName':qoiName}
   if ifsym:
       return qoiName,nsnap,nplt,r,timeprdc,      \
              outMod,outSnp,maxMode,              \
              if3D,ifsym,ifPickSave,ifPickRead,   \
              info,info_m,info_s     
              #if POD module is selected, r is not used
   else:
       return qoiName,nsnap,nplt,r,timeprdc,      \
              outMod,outSnp,maxMode,              \
              if3D,ifsym,ifPickSave,ifPickRead,   \
              info,info_m,info 
              #if not symm info are back twice
              #if POD module is selected, r is not used

