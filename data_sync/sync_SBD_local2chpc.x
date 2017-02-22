#! /bin/tcsh -f

rsync -avzh -e ssh --progress ~/SBD_VAPOR/L[0,1] kingspeak.chpc.utah.edu:/uufs/chpc.utah.edu/common/home/bowen-group1/SBD_VAPOR/L[0,1]
