# L0_user_specs.R
# rich fiorella, 15mar17.

# this file contains all of the parameters that need to be specified for the L0 data processing script.

# 1. what dates should we process? script works sequentially on dates from startdate to enddate.

startdate <- "2016-12-08" #expects a number in yyyymmdd format - will begin at MIDNIGHT GMT.
enddate   <- "2017-05-01" #expects a number in yyyymmdd format - will end at MIDNIGHT GMT.

# 2. where is the data we're processing? and where should we save output data?
# don't forget trailing slash!
path.to.data <- "~/VaporData/SBD_VAPOR/Raw/"
path.to.output.L0.data <- "~/VaporData/SBD_VAPOR/L0/v1_1_0beta/"

# 3. what do we call the output data? file name will have the format of:
# (path.to.output.L0.data)/(output.file.prefix)_Calib/AmbientData_L0_YYYY-mm-dd_(codeversion).dat
output.file.prefix <- "HDP_Water_Vapor"

# 4. This step takes averages some of the ambient data from its native 1.16 Hz time resolution.
# these parameters specify: (a) how long should each averaging period be (in minutes)?
# (b) how many data points have to be within that averaging period in order for the average
# to be taken/included in the output data? for example, function defaults are set to 1 minute
# and 20 data points, indicating that 20 data points have to be present in each 1 minute period
# to be kept.
# they modify the behavior of the reduce.ambient.data function...

averaging.length.in.minutes <- 1 	# must be a number!
minimum.number.of.datapoints <- 20	# must be a number!

# 5. Is debugging necessary? This parameter will help determine why code is crashing.
debug <- 2 # set to 0 for no diag output, >0 for increased diagnostic levels for debugging.
             
######################################################################################################################################
# SET METADATA
# write out metadata to help data curation - these will be appended to the top of the data files currently. it might be possible in
# later versions to use this to write out an xml file?
######################################################################################################################################

metadata.frame <- read.csv("../internal_metadata/L0_Snowbird_metadata.csv",header=TRUE)