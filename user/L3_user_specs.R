# L3_user_specs.R
# rich fiorella, 3apr17.

# this file contains all of the parameters that need to be specified for the L3 data processing script.

start.date <- ymd("2016-10-01")
end.date <- ymd("2017-04-01")

# where is the uncalibrated ambient vapor data?
path.to.L1.data <- "~/VaporData/WBB_VAPOR/L1/v1beta/"

# id calibration data to use:
path.to.L2.calib.data <- "~/SLC_inversion/H2O_iso/raw/wbb/L2/"
calib.file <- "WBB_Water_Vapor_CalibrationRegressionData_L2_2016-10-01_2017-04-01_handclean.dat"

# where should output data go?
path.to.output.L3.data <-  "~/SLC_inversion/H2O_iso/raw/wbb/L3/"
output.file.prefix <- "WBB_Water_Vapor_AmbientDataVSMOW_L3_"

RUN_PLOTS <- "TRUE"
debug <- 0 # 0 for no debugging, >0 for (increasingly verbose) debugging
