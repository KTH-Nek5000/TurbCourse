#!/bin/bash                                                                                         
casename='mixlay'
noPrcs=4    #number of processors

echo $casename > SESSION.NAME
echo $PWD/ >> SESSION.NAME

if [ -f "logfile" ]; then
	   rm logfile
   fi
   mpirun -np $noPrcs ./nek5000 >>logfile&   #use either ./nek5000 or nek5000
