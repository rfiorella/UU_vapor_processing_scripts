# L2_user_specs.R
# rich fiorella, 15mar17.

# this file contains all of the parameters that need to be specified for the L2 data processing script.
# NOTE: this is long because there are a lot design choices that needed to be made, and choices that
# I made for my data might not apply for other users!! Please take some time to read through this script 
# and ensure that these choices make sense for your setup as well!

###################################################################
# SECTION 1: Dates to run, data location, include debugging?

# ALL PATHS REQUIRE TRAILING SLASH!! (e.g., all paths below - at least on Unix/os X - 
# require forward slash at end of the path!)

# 1a. what dates should we process? script works sequentially on dates from startdate to enddate.

start.date <- ymd("2016-12-01")
end.date <- ymd("2017-03-01")

# 1b. where is the data we're processing? and where should we save output data?
path.to.L1.data <- "~/Dropbox/SLVhistoricalCDV/data/wateriso/L1/"
path.to.output.L2.data <-  "~/Dropbox/SLVhistoricalCDV/data/wateriso/L2/"

# 1c. what do we call the output data? file name will have the format of:
# (path.to.output.L0.data)/(output.file.prefix)_Calib/AmbientData_L0_YYYY-mm-dd_(codeversion).dat
output.file.prefix <- "WBB_Water_Vapor"

# 1d. should we run diagnostic plots? (logical value) and where should they be written?
RUN_PLOTS <- TRUE
plot.path <- "~/Dropbox/SLVhistoricalCDV/scripts/UU_Vapor_processing_scripts-1.2.0/plots/"

# 1e. Is debugging necessary? This parameter will help determine why code is crashing.
debug <- 2

###################################################################
# SECTION 2: Information specific to the analyzer

# 2a. What is the dependence of delta on humidity for this analyzer? 
# correction is implemented in point-slope formula in the following way
# delta(@20000ppm) = slope*(1/20000-1/measured.H2O) + measured.delta

# set site...
site <- "WBB"

# allowed values for fit.type, and what other variables are expected with
# each type:
# a) hyberbolic, requires: Oslope,Hslope
# b) logarithmic, requires: Oslope,Hslope
# c) hyberbolic.offset, requires: Oslope, Hslope, Ointercept, Hintercept
# note: c treats the dependence on humidity as an offset parameter that is
# to subtracted from measured values - underlying function calculates offset
# based on eqn and then subtracts from analyzer measured value.

if (site=="WBB") {
  fit.type <- "hyperbolic.offset"
  # these parameters have been updated from a humidity cal synthesis
  # across all five samplings throughout time. they differ from the values
  # in the WBB_July17_humidcal regression - hopefully not by too much!!
  Oslope <- -5838.945589
  Ointercept <- 0.308163
  Hslope <- -9726.7185548
  Hintercept <- 0.4081879
} else if (site=="Snowbird") {
  fit.type <- "logarithmic" # allowed values: hyberbolic (1/H2O)
  Oslope <- 2.59226
  Hslope <- 15.0711
}

# 2b. This is a long one - and one that hopefully can be shortened in future.
# User specified information on calibration standards - this example function
# has the information for our two analyzers at WBB and @ Snowbird. this will have
# to be modified for your site!!

assign.standard.names.and.values <- function(data,O18.break=-8.0,D.break=-50.0) {
  if (site=="WBB") {
    # need to assign standard names and d18O, dD values (VSMOW) based on 
    # corrected concentrations.
    # first, break into subsets based on compositions

    # modified 16sep17 to account for accidental change back to PZ in February 2017.
    # switch was done after attempted humidity calibration on 2016-02-21.
    # the light standard - UD - I think remained unchanged.

    #RPF note for ACPD MS - there was an additional standard switch mentioned in Gorski et. al.
    #perhaps this was before Dec13?
    std.switch1 <- as.numeric(as.POSIXct("2017-02-16",tz="GMT",origin="1970-01-01"))
    std.switch2 <- as.numeric(as.POSIXct("2017-02-21",tz="GMT",origin="1970-01-01"))

    UD2.stds <- which(data$Delta_18_16_bgc < O18.break & data$time.mean > std.switch2)
    PZ.stds <- which(data$Delta_18_16_bgc >= O18.break & data$time.mean > std.switch2)
    UD1.stds <- which(data$Delta_18_16_bgc < O18.break & data$time.mean < std.switch1)
    FL.stds <- which(data$Delta_18_16_bgc >= O18.break & data$time.mean < std.switch1)

    # assign standard name/values
    PZ.standard <- "PZ"
    PZ.standard.OValue <- 1.65
    PZ.standard.DValue <- 16.9

    # UD.standard <- "UD"
    # UD.standard.OValue <- -16.52
    # UD.standard.DValue <- -123.1

    UD1.standard <- "UT"
    UD1.standard.OValue <- -16.0
    UD1.standard.DValue <- -121.0

    FL.standard <- "FL"
    FL.standard.OValue <- -1.23
    FL.standard.DValue <- -5.51

    # added second version of UT di water. 
    UD2.standard <- "UTD2"
    UD2.standard.OValue <- -15.88
    UD2.standard.DValue <- -119.66

    # fit standard values into the vectors corresponding to proper indices in standard.data.frame
    add.names <- vector(length=nrow(data))
    add.std.Oval <- vector(length=nrow(data))
    add.std.Dval <- vector(length=nrow(data))

    # if data for the standard exists, loop through and assign it to these rows.
    if (!is.null(UD2.stds) & length(UD2.stds)>0) {
      add.names[UD2.stds] <- rep(UD2.standard,length(UD2.stds))
      add.std.Oval[UD2.stds] <- rep(UD2.standard.OValue,length(UD2.stds))
      add.std.Dval[UD2.stds] <- rep(UD2.standard.DValue,length(UD2.stds))
    }

    if (!is.null(PZ.stds) & length(PZ.stds)>0) {
      add.names[PZ.stds] <- rep(PZ.standard,length(PZ.stds))
      add.std.Oval[PZ.stds] <- rep(PZ.standard.OValue,length(PZ.stds))
      add.std.Dval[PZ.stds] <- rep(PZ.standard.DValue,length(PZ.stds))
    }

    if (!is.null(UD1.stds) & length(UD1.stds)>0) {
      add.names[UD1.stds] <- rep(UD1.standard,length(UD1.stds))
      add.std.Oval[UD1.stds] <- rep(UD1.standard.OValue,length(UD1.stds))
      add.std.Dval[UD1.stds] <- rep(UD1.standard.DValue,length(UD1.stds))
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
    O18.break <- -30
    H2.break <- -200

    # need to assign standard names and d18O, dD values (VSMOW) based on 
    # corrected concentrations.
    # first, break into subsets based on compositions

    SPS.stds <- which(data$Delta_18_16_bgc < O18.break)
    UD.stds <- which(data$Delta_18_16_bgc >= O18.break)

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

###################################################################
# SECTION 3: Script design decisions that determine: what is a good
# standard peak? what exclusion criteria do we use? what amount of
# water vapor do we assume makes it through a drierite column? etc.

# this section is organized by what function in L2_Script.R and L2_functions.R
# is being referenced by individual parameters.

# ID.calib.breakpoints function ------
time.threshold <- 10 # time in seconds.
# meaning: what amount of time passing between two adjacenet
# calibration rows in the data frame should be used to divide
# the data into separate analyses? This number could presumably
# be a bit larger than 10, but there's also no reason it shouldn't
# be this short in most cases...

# fit.calibration.splines function -------
stiff.spline.dfs <- 12 # number of degrees of freedom in "stiff" splines fit to data - 
                       # seeking a value here that will allow for "minor" hiccups in the
                       # data (e.g., a sufficiently small air bubble or the like), but
                       # catch larger excursions in the data.

# extact.stable.calib.indices function ------------
# this function IDs which points in the identified peaks are "stable"
# by selecting only the points where the derivative of the spline created
# in fit.calibration.splines is below each of the following thresholds:
# IMPORTANT NOTE: IF YOU CHANGE STIFF.SPLINE.DFS YOU *MUST* RECHECK THESE
# DERIVATIVES - THEY ARE TUNED BASED ON A VALUE OF 12!

H2O.thres <- 20.0
d18O.thres <- 0.01
d2H.thres <- 0.1

# calculate.standard.averages function -------------

# 1. should we attempt to limit the effect of instrumental memory?
# built in function sequentially throws out 10% of data at the beginning
# of peak until the estimated trend magnitude (permil/min) is less than 
# 0.1 for d18O and 0.5 for d2H, or 80% of the data has been removed. 
memory.filter <- TRUE

# apply.mixingratio.correction ------------
# parameters already included in section 2!

# apply.drygas.correction -----------------

do.correction <- TRUE # set to true if using drierite, can set to false otherwise and keep same data structure
                      # not yet implemented, so this switch currently has no effect.

H2O.bg <- 250 # integer - what ppmv concentration of H2O assumed to get through column?

include.gypsum.fractionation <- FALSE
# additional second-order adjustment to vapor making it through column by 
# including influence of gypsum hydration water fractionation factors onto
# the column. Not yet implemented, so this switch currently has no effect.

# Finally, perform some filtering to the produced average data frame -------
# this section seeks to define additional parameters/checks on what
# constitutes an acceptible standard analysis.

h2o.max.thres <- 35000    # maximum mean H2O concentration allowed
h2o.min.thres <- 2000     # minimum mean H2O concentration allowed
h2o.sdev.thres <- 1000    # maximum standard deviation for H2O allowed
d18O.sdev.thres <- 1      # maximum standard deviation for d18O allowed
d2H.sdev.thres <- 8       # maximum standard deviation for d2H allowed
max.length.thres <- 4200  # maximum number of indices allowed (at 1.16 Hz, 4200 = 1 hour)
min.length.thres <- 35    # minimum number of indices required (at 1.16 Hz, 35 = 30 seconds)