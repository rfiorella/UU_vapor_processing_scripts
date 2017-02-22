#! /bin/tcsh -f

# sync L0 files
#rsync -avzh -e ssh --progress ~/WBB_VAPOR/L0/ kingspeak.chpc.utah.edu:/uufs/chpc.utah.edu/common/home/bowen-group1/WBB_VAPOR/L0/

# sync L1 files
rsync -avzh --update -e ssh --progress ~/WBB_VAPOR/L1 kingspeak.chpc.utah.edu:/uufs/chpc.utah.edu/common/home/bowen-group1/WBB_VAPOR/L1
