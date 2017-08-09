# L1_user_specs.R
# rich fiorella, 15mar17.

# this file contains all of the parameters that need to be specified for the L1 data processing script.

# 1. what dates should we process? script works sequentially on dates from startdate to enddate.

start.date <- ymd("2016-12-01")
end.date <- ymd("2017-05-01")

# 2. where is the data we're processing? and where should we save output data?
path.to.L0.data <- "~/VaporData/SBD_VAPOR/L0/v1_1_0beta/"
path.to.output.L1.data <-  "~/VaporData/SBD_VAPOR/L1/v110beta/"

# 3. what do we call the output data? file name will have the format of:
# (path.to.output.L0.data)/(output.file.prefix)_Calib/AmbientData_L0_YYYY-mm-dd_(codeversion).dat
output.file.prefix <- "HPD_Water_Vapor"

# 4. Is debugging necessary? This parameter will help determine why code is crashing.
debug <- 0

######################################################################################################################################
# SET METADATA
# write out metadata to help data curation - these will be appended to the top of the data files currently. it might be possible in
# later versions to use this to write out an xml file?
######################################################################################################################################

#metadata.frame <- read.csv("../metadata_templates/L0_WBB_metadata.csv",header=TRUE)