"""
The module reads and saves database in a pickle format.
The pickle file is read/saved in the OUTPUT path.
"""

import pickle




def pickReader(info):
    pickleFile_=info['outputPath']+info['caseName']+'_'+str(info['startID'])+'_'+str(info['endID'])
    with open(pickleFile_,'rb') as F:
        Usnp=pickle.load(F)
        db=pickle.load(F)
        data_ms=pickle.load(F)
    print('... pickle database read.') 
        
    return Usnp,db,data_ms

def pickReader_sym(info,info_s):
    pickleFile_=info['outputPath']+info['caseName']+'_'+str(info['startID'])+'_'+str(info['endID'])+'_'+'sym'
    with open(pickleFile_,'rb') as F:
        Usnp=pickle.load(F)
        db=pickle.load(F)
        db_s=pickle.load(F)
        data_ms=pickle.load(F)
    print('... pickle database read.')  
    return Usnp,db,db_s,data_ms

def pickWriter(info,Usnp,db,data_ms,iff):
    if (iff):
        pickleFile_=info['outputPath']+info['caseName']+'_'+str(info['startID'])+'_'+str(info['endID'])
        #dump the pickle data    
        with open(pickleFile_,'wb') as F:
            pickle.dump(Usnp,F)
            pickle.dump(db,F)
            pickle.dump(data_ms,F)
        print('... pickle database created and saved.')   
    return 

def pickWriter_sym(info,Usnp,db,db_s,data_ms,iff):
    if (iff):
        pickleFile_=info['outputPath']+info['caseName']+'_'+str(info['startID'])+'_'+str(info['endID'])+'_'+'sym'
        #dump the pickle data    
        with open(pickleFile_,'wb') as F:
            pickle.dump(Usnp,F)
            pickle.dump(db,F)
            pickle.dump(db_s,F)
            pickle.dump(data_ms,F)
        print('... pickle database created and saved.')  
    return 
