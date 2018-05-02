# L3 functions - need to update the header here...

extract.date.from.L1.files <- function(file.list,dbg.level=0) {
  # print statement announcing the beginning of the function
  if (dbg.level>0) { 
    print("============================================")
    print("Start of extract.date.from.L1.files function") 
  }

  # suppress warnings in this function - as.numeric will create NAs by coercion,
  # and this will be reported as a warning to the screen every time. But since this
  # is precisely what we are hoping to do, warnings are unncessary here.
  oldw <- getOption("warn")
  options(warn = -1)
  # split file name into pieces using *underscores*
  file.name.string.pieces <- strsplit(file.list,split="_")
  # extract dates from resulting string pieces
  file.name.dates <- ymd(unlist(file.name.string.pieces))
  # remove NAs
  file.name.dates.woNAs <- file.name.dates[!is.na(file.name.dates)]
  # set warnings to old status
  options(warn = oldw)

  # print statement announcing the beginning of the function
  if (dbg.level>0) { 
    print("End of extract.date.from.L1.files function") 
    print("==========================================")
  }
  
  # return date list
  return(file.name.dates.woNAs)
}

#---------------------------------------------------------------------------
# attach.L3.Header

attach.L3.Header <- function(output_filename,metadata_dataframe,dbg.level=0) {
  # print statement indicating that this function is starting
  if (dbg.level>0) {
    print("==================================")
    print("starting attach.L3.Header function")
  }

  # initiate writing out to output file...
  sink(output_filename)

  # add L1 information to top of output file...
  cat("# BEGIN L3 HEADER \n")

  # note when L2 script was run
  cat(paste("# date.L3.scripts.run, ",base::date(),"\n"))

  # calibration file used...
  cat(paste("# Calibration file used, ",path.to.L2.calib.data,calib.file,"\n"))

  # close L3 header
  cat("# END L3 HEADER \n")

  # append L0/L1/L3 headers.
  # loop through existing headers and paste below L3 header.
  for (i in 1:length(metadata_dataframe)) {
    cat(paste(metadata_dataframe[i],"\n")) 
  }
  
  # turn off output
  sink()
  
  # print statment denoting the end of this function
  if (dbg.level>0) {
    print("ending attach.L3.Header function")
    print("================================")
  }
}

#---------------------------------------------------------
# reduce.calibrated.ambient.vapor

reduce.calibrated.ambient.vapor <- function(amb.df,interval.in.minutes) {
  use.xts <- require(xts)

  if (use.xts==TRUE) {
      # print statement detailing what this function does...
    print(paste(Sys.time()," Reducing ambient data to ",interval.in.minutes," minute averages"))
    dframe <- as.xts(amb.df,
      order.by=as.POSIXct(amb.df$EPOCH_TIME,tz="UTC",origin="1970-01-01"))

    # get endpoints
    eps <- endpoints(dframe,on="minutes",k=interval.in.minutes)

    # take mean.
    dframe.reduced <- period.apply(dframe,eps,mean)
    
    # return data
    return(as.data.frame(coredata(dframe.reduced)))
  }
}


#-----------------------------------------------------------------------------------
# apply.mixingratio.correction - this function corrects the measured isotopic compositions
# using a linear regression between delta and 1/[H2O] - based on coefficients determined from
# CalibrationStudyRegressionHandClean_gjb. Requires mixing ratio coefficients to be specified 
# in L1_Calibration_Parameters.R. Other regressions are possible for correcting for
# concentration dependence, but have not be implemented in this code base yet.

remove.humidity.dependence <- function(data.frame,fit.type,Oslope,Hslope,dbg.level=0) {
  # print statement indicating that this function is starting
  if (dbg.level>0) {
    print("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&")
    print("starting apply.mixingratio.correction.calibration function")
  }

  #======================================================
  # start work of function

  print(paste(now()," Applying mixing ratio correction to: ",deparse(substitute(data.frame))))

  if (fit.type=="hyperbolic.offset") {
    # save a copy of uncorrected variables.
    data.frame$Delta_18_16_nocorr <- data.frame$Delta_18_16
    data.frame$Delta_D_H_nocorr <- data.frame$Delta_D_H
    # calculate offsets
    O_offset <- Ointercept + Oslope/data.frame$H2O
    H_offset <- Hintercept + Hslope/data.frame$H2O
    # subtract offsets from measurements
    data.frame$Delta_18_16 <- data.frame$Delta_18_16_nocorr - O_offset
    data.frame$Delta_D_H <- data.frame$Delta_D_H_nocorr - H_offset
  }

  # print statment denoting the end of this function
  if (dbg.level>0) {
    print("ending apply.mixingratio.correction.calibration function")
    print("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&")
  }

  # return the whole data frame with the two mixing ratio corrected variables attached.
  return(data.frame)
}

#------------------------------------------------------------