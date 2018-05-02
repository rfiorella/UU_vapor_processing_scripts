# L3_user_specs.R
# rich fiorella, 3apr17.

# this file contains all of the parameters that need to be specified for the L3 data processing script.

start.date <- ymd("2015-12-01")
end.date <- ymd("2016-03-01")

# where is the uncalibrated ambient vapor data?
path.to.L1.data <- "~/VaporData/WBB_VAPOR/L1/v110beta/"

# id calibration data to use:
path.to.L2.calib.data <- "~/VaporData/WBB_VAPOR/L2/v110beta/w1314/"
calib.file <- "WBB_Water_Vapor_CalibrationData_L2_2015-11-01_2016-03-01.rds"

# where should output data go?
path.to.output.L3.data <-  "~/VaporData/WBB_VAPOR/L3/v110beta/"
output.file.prefix <- "WBB_Water_Vapor_AmbientDataVSMOW_L3_"

RUN_PLOTS <- "TRUE"
debug <- 0 # 0 for no debugging, >0 for (increasingly verbose) debugging

# output time frequency (uses xts)
out.interval <- 5 # minutes.

#-------------------------------------------
# How to handle gaps in calibration regressions?
#-------------------------------------------

# a. should we force the first regression period in calibration
# file to correspond to the beginning of the ambient period 
# being processed?
force.to.beginning <- TRUE 

# b. should we force the last regression period in calibration
# file to correspond to the end of the ambient period being
# processed?
force.to.end <- TRUE

# c. how should we deal with gaps in the calibration file?
# three possible options: (a) omit, (b) carry last calibration
# forward, (c) carry next calibration backward.
gap.handling <- "last.carried.forward" 	# omit, 
										# last.carried.forward, 
										# or next.carried.backward

# import fit data from L2 for humidity dependence of delta.
fit.type <- "hyperbolic.offset"
# these parameters have been updated from a humidity cal synthesis
# across all five samplings throughout time. they differ from the values
# in the WBB_July17_humidcal regression - hopefully not by too much!!
Oslope <- -5838.945589
Ointercept <- 0.308163
Hslope <- -9726.7185548
Hintercept <- 0.4081879
