#! /bin/tcsh -f

rsync -avzh -e ssh --progress --exclude 'to_restore/' kingspeak.chpc.utah.edu:/uufs/chpc.utah.edu/common/home/bowen-group1/SBD_VAPOR/Raw ~/SBD_VAPOR/
