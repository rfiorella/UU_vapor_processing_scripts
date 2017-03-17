# L2_user_specs.R
# rich fiorella, 15mar17.

# this file contains all of the parameters that need to be specified for the L2 data processing script.

# 1. what dates should we process? script works sequentially on dates from startdate to enddate.

start.date <- ymd("2016-12-01")
end.date <- ymd("2017-02-01")

# 2. where is the data we're processing? and where should we save output data?
path.to.L1.data <- "~/VaporData/SBD_VAPOR/L1/testing/"
path.to.output.L2.data <-  "~/VaporData/SBD_VAPOR/L2/WBB_SBD_comp/"

# 3. what do we call the output data? file name will have the format of:
# (path.to.output.L0.data)/(output.file.prefix)_Calib/AmbientData_L0_YYYY-mm-dd_(codeversion).dat
output.file.prefix <- "SBD_Water_Vapor"

# 4. should we run diagnostic plots? (logical value)
RUN_PLOTS <- "TRUE"

# 5. Is debugging necessary? This parameter will help determine why code is crashing.
debug <- 1

# 6. What is the dependence of delta on humidity for this analyzer? 
# correction is implemented in point-slope formula in the following way
# delta(@20000ppm) = slope*(1/20000-1/measured.H2O) + measured.delta

# set site...
site <- "WBB"

if (site=="WBB") {
  fit.type <- "hyperbolic" # allowed values: hyperbolic (1/H2O)
  Oslope <- -3.618*10^3
  Hslope <- -1.656*10^4
} else if (site=="Snowbird") {
  fit.type <- "logarithmic" # allowed values: hyberbolic (1/H2O)
  Oslope <- 2.59226
  Hslope <- 15.0711
}

######################################################################################
# User specified information on calibration standards - this example function
# has the information for our two analyzers at WBB and @ Snowbird. this will have
# to be modified for your site!!


## UU Specific Example!
site <- "Snowbird"

assign.standard.names.and.values <- function(data,O18.break=-8.0,D.break=-50.0) {
  if (site=="WBB") {
    # need to assign standard names and d18O, dD values (VSMOW) based on 
    # corrected concentrations.
    # first, break into subsets based on compositions

    std.switch <- as.numeric(as.POSIXct("2013-11-15",tz="GMT",origin="1970-01-01"))

    UD.stds <- which(data$Delta_18_16_bgc < O18.break & data$time.mean < std.switch)
    PZ.stds <- which(data$Delta_18_16_bgc >= O18.break & data$time.mean < std.switch)
    UT.stds <- which(data$Delta_18_16_bgc < O18.break & data$time.mean > std.switch)
    FL.stds <- which(data$Delta_18_16_bgc >= O18.break & data$time.mean > std.switch)

    # assign standard name/values
    PZ.standard <- "PZ"
    PZ.standard.OValue <- 1.65
    PZ.standard.DValue <- 16.9

    UD.standard <- "UD"
    UD.standard.OValue <- -16.52
    UD.standard.DValue <- -123.1

    UT.standard <- "UT"
    UT.standard.OValue <- -16.0
    UT.standard.DValue <- -121.0

    FL.standard <- "FL"
    FL.standard.OValue <- -1.23
    FL.standard.DValue <- -5.51

    # fit standard values into the vectors corresponding to proper indices in standard.data.frame
    add.names <- vector(length=nrow(data))
    add.std.Oval <- vector(length=nrow(data))
    add.std.Dval <- vector(length=nrow(data))

    # if data for the standard exists, loop through and assign it to these rows.
    if (!is.null(UD.stds) & length(UD.stds)>0) {
      add.names[UD.stds] <- rep(UD.standard,length(UD.stds))
      add.std.Oval[UD.stds] <- rep(UD.standard.OValue,length(UD.stds))
      add.std.Dval[UD.stds] <- rep(UD.standard.DValue,length(UD.stds))
    }

    if (!is.null(PZ.stds) & length(PZ.stds)>0) {
      add.names[PZ.stds] <- rep(PZ.standard,length(PZ.stds))
      add.std.Oval[PZ.stds] <- rep(PZ.standard.OValue,length(PZ.stds))
      add.std.Dval[PZ.stds] <- rep(PZ.standard.DValue,length(PZ.stds))
    }

    if (!is.null(UT.stds) & length(UT.stds)>0) {
      add.names[UT.stds] <- rep(UT.standard,length(UT.stds))
      add.std.Oval[UT.stds] <- rep(UT.standard.OValue,length(UT.stds))
      add.std.Dval[UT.stds] <- rep(UT.standard.DValue,length(UT.stds))
    }

    if (!is.null(FL.stds) & length(FL.stds)>0) {
      add.names[FL.stds] <- rep(FL.standard,length(FL.stds))
      add.std.Oval[FL.stds] <- rep(FL.standard.OValue,length(FL.stds))
      add.std.Dval[FL.stds] <- rep(FL.standard.DValue,length(FL.stds))
    }

    # add columns to standard.data.frame
    data <- cbind(data,add.names,add.std.Oval,add.std.Dval)

    # rename column headers to something more useful
    names(data)[names(data) == 'add.names'] <- 'standard.name'
    names(data)[names(data) == 'add.std.Oval'] <- 'standard.O18.VSMOW'
    names(data)[names(data) == 'add.std.Dval'] <- 'standard.H2.VSMOW'

    # return standard.data.frame
    return(data)

  } else if (site=="Snowbird") { 
    # overwrite O break and H break for SBD
    O.break <- -30
    H.break <- -200

    # need to assign standard names and d18O, dD values (VSMOW) based on 
    # corrected concentrations.
    # first, break into subsets based on compositions

    SPS.stds <- which(data$Delta_18_16_bgc > O18.break)
    UD.stds <- which(data$Delta_18_16_bgc <= O18.break)

    # assign standard name/values

    # standard values for snowbird based on SPATIAL samples 21814 (SPS1) 
    # and 21480 (bldg DI) run on 11/2 and 10/12 respectively.

    SPS.standard <- "SPS1"
    SPS.standard.OValue <- -46.81
    SPS.standard.DValue <- -364.31

    UD.standard <- "UD"
    UD.standard.OValue <- -15.32
    UD.standard.DValue <- -116.03

    # fit standard values into the vectors corresponding to proper indices in standard.data.frame
    add.names <- vector(length=nrow(data))
    add.std.Oval <- vector(length=nrow(data))
    add.std.Dval <- vector(length=nrow(data))

    # if data for the standard exists, loop through and assign it to these rows.
    if (!is.null(UD.stds) & length(UD.stds)>0) {
      add.names[UD.stds] <- rep(UD.standard,length(UD.stds))
      add.std.Oval[UD.stds] <- rep(UD.standard.OValue,length(UD.stds))
      add.std.Dval[UD.stds] <- rep(UD.standard.DValue,length(UD.stds))
    }

    if (!is.null(SPS.stds) & length(SPS.stds)>0) {
      add.names[SPS.stds] <- rep(SPS.standard,length(SPS.stds))
      add.std.Oval[SPS.stds] <- rep(SPS.standard.OValue,length(SPS.stds))
      add.std.Dval[SPS.stds] <- rep(SPS.standard.DValue,length(SPS.stds))
    }

    # add columns to standard.data.frame
    data <- cbind(data,add.names,add.std.Oval,add.std.Dval)

    # rename column headers to something more useful
    names(data)[names(data) == 'add.names'] <- 'standard.name'
    names(data)[names(data) == 'add.std.Oval'] <- 'standard.O18.VSMOW'
    names(data)[names(data) == 'add.std.Dval'] <- 'standard.H2.VSMOW'
 	
 	# return standard.data.frame
    return(data)
  }
}

######################################################################################################################################
# SET METADATA
# write out metadata to help data curation - these will be appended to the top of the data files currently. it might be possible in
# later versions to use this to write out an xml file?
######################################################################################################################################

metadata.frame <- read.csv("../metadata_templates/L0_WBB_metadata.csv",header=TRUE)

# Note: several variables left to port here, including: (a) spline parameters, (b) filtering info for standards