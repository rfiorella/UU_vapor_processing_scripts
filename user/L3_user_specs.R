# L3_user_specs.R
# rich fiorella, 3apr17.

# this file contains all of the parameters that need to be specified for the L3 data processing script.

start.date <- ymd("2013-12-01")
end.date <- ymd("2017-03-01")

# where is the uncalibrated ambient vapor data?
path.to.L1.data <- "~/VaporData/WBB_VAPOR/L1/v110beta/"

# id calibration data to use:
path.to.L2.calib.data <- "~/VaporData/WBB_VAPOR/L2/v110beta/"
calib.file <- "WBB_Water_Vapor_CalibrationRegressionData_L2_2013-05-01_2017-05-01.dat"

# where should output data go?
path.to.output.L3.data <-  "~/VaporData/WBB_VAPOR/L3/v110beta/"
output.file.prefix <- "WBB_Water_Vapor_AmbientDataVSMOW_L3_"

RUN_PLOTS <- "TRUE"
debug <- 0 # 0 for no debugging, >0 for (increasingly verbose) debugging

# output time frequency (uses xts)
out.interval <- 5 # minutes.