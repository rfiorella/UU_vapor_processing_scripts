# Picarro Water Vapor Isotope Data Reduction 
# L2 Functions

# migrated from WBB-vapor-scripts repo by rich fiorella on 22/2/17.
# code in these functions written by rich fiorella based on work by galen gorski
# (likely some of galen's functions have made it through with no changes!)

# contact info: rich fiorella, rich.fiorella@utah.edu

# this script contains a series of accessory functions for the L1 script.
# as of 22/2/17, this script has been designed to work only with a Picarro L2130.
# other analyzer versions will likely require modification of this script.

extract.date.from.L1.files <- function(file.list) {
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
  # return date list
  return(file.name.dates.woNAs)
}

#------------------------------------------------------------------------------------------
# ID.calib.breakpoints - finds indices in the input cailbration data frame where the
# gap between points exceeds threshold number of seconds (default = 10 seconds)

ID.calib.breakpoints <- function(calibration.data.frame,thres=10) {
  print(paste(Sys.time()," separating calibration periods based on gaps in time..."))
  if (!"EPOCH_TIME" %in% colnames(calibration.data.frame)) {
    stop("Expected time array absent from input data frame, check input data to ID.calib.breakpoints")
  }
  tmp <- diff(calibration.data.frame$EPOCH_TIME)
  # note that diff returns the difference between point i and i+1 at index i
  # ...when pulling out the indices that indicate the start of the subsequent 
  # analysis, these are at index i+1, which has been returned at index i - so
  # need to add 1 to the indices returned by diff 
  out1 <- which(tmp>thres)+1

  # add points to the beginning and end f the array to fully bound first/last plateaus
  out1 <- c(1,out1,nrow(calibration.data.frame)) 

  # in cases of power failure or some other hiccup, sometimes two consecutive points
  # can be identified as break points...this causes havoc in the spline fitting
  # process, and therefore, points that are too close to each other must be filtered out.
  # proposing a threshold here of 30 indices.
  tmp <- diff(out1)

  ini.breaks <- out1[c(1e6,tmp)>30]

  # loop through initial break points, calculate splines, ID points where change is huge.
  fix.break.splines <- vector("list",length(ini.breaks))
  
  # calculate splines

  # NOTE: in this portion of the script, df and the threshold derivative value are TUNED parameters
  # they may not work in all situations...please check your data and ensure this works for your 
  # data...
  for (i in 1:(length(ini.breaks)-1)) {
    tmp.splines <- smooth.spline(calibration.data.frame$EPOCH_TIME[ini.breaks[i]:(ini.breaks[i+1]-1)],
      calibration.data.frame$H2O[ini.breaks[i]:(ini.breaks[i+1]-1)],df=96,tol=1e-4)
    tmp.spline.deriv <- c(NA,diff(tmp.splines$y)) # initial NA required to keep vector same length.
    # find indices where derivative is "exceedingly" far from zero
    tmp.inds <- which(abs(tmp.spline.deriv)>50) # 
    tmp.inds.diff <- c(NA,diff(tmp.inds)) # add initial NA to keep tmp.inds and tmp.inds.diff the same length
    # tmp.inds.diff indicates large index gaps between points with large derivatives.
    # to separate out distinct calibration periods, we seek to find when the difference
    # between points is greater than about 3 minutes. L2130i collects data at ~1.25 Hz,
    # so 3 minutes corresponds to ~225 indices.
    fix.break.splines[[i]] <- ini.breaks[i]+tmp.inds[which(tmp.inds.diff >= 225)]
  }
  
  # unlist indices in fix.break.splines and add to ini.breaks
  addl.breaks <- unlist(fix.break.splines,use.names=FALSE)

  # add additional breaks to initial break column and sort.
  fin.breaks.tmp <- sort(c(ini.breaks,addl.breaks))
  fin.diff.tmp <- c(1e6,diff(fin.breaks.tmp))
  
  # require difference between points to be >3 minutes (e.g., there is no calibration cycle
  # shorter than 3 minutes) - this is likely a conservative threshold (could be higher!) - 
  # but it seems to be sufficient in most cases.
  
  # NOTE: important to be careful about where the additional value goes to make the diff output
  # the same length as the diff input. in many cases it's desirable to know at index i what the
  # difference is between index i and index i-1 (this is the fin.diff.temp construction above).
  # it is most useful here to know if the difference between index i and i+1 is <= a threshold.
  # therefore, the final "spare" value is placed at the end of the array.
  time.diff <- which(c(diff(calibration.data.frame$EPOCH_TIME[fin.breaks.tmp]),1e6)>180)
  array.diff <- which(c(diff(fin.breaks.tmp),1e6)>225)

  fin.breaks <- fin.breaks.tmp[intersect(time.diff,array.diff)]

  return(fin.breaks)
}

#------------------------------------------------------------------------------------------
# fit.calibration.splines - fits a smooth spline to each identified calibration period.

fit.calibration.splines <- function(calibration.data.frame,breakpoints,dfree=12,tolerance=1e-4) {
  print(paste(Sys.time()," fitting splines to calibration periods..."))
  if (!"EPOCH_TIME" %in% colnames(calibration.data.frame)) {
    stop("Expected time array absent from input data frame, check input data to fit.calibration.splines")
  }
  # calculate number of plateaus, set up lists for processing.
  nplat <- length(breakpoints)-1 # number of plateaus - have to add one to account for first plateau...

  # intialize output variable - which will be a list of dataframes
  output.spline.fits <- vector("list",nplat)

  # loop through each ID'ed plateau for the desired variable.
  for (i in 1:nplat) {
    # fit splines.
    #print(i)
    H2O.spline.fits <- smooth.spline(calibration.data.frame$EPOCH_TIME[(breakpoints[i]):(breakpoints[i+1]-1)],
      calibration.data.frame$H2O[(breakpoints[i]):(breakpoints[i+1]-1)],df=dfree,tol=tolerance)
    #print(i)
    d18O.spline.fits <- smooth.spline(calibration.data.frame$EPOCH_TIME[(breakpoints[i]):(breakpoints[i+1]-1)],
      calibration.data.frame$Delta_18_16[(breakpoints[i]):(breakpoints[i+1]-1)],df=dfree,tol=tolerance)
    #print(i)
    d2H.spline.fits <- smooth.spline(calibration.data.frame$EPOCH_TIME[(breakpoints[i]):(breakpoints[i+1]-1)],
      calibration.data.frame$Delta_D_H[(breakpoints[i]):(breakpoints[i+1]-1)],df=dfree,tol=tolerance)
    # extract y value from splines, x value from 1 spline (should all be the same x values!)
    temp1 <- H2O.spline.fits$x
    temp2 <- H2O.spline.fits$y
    temp3 <- d18O.spline.fits$y
    temp4 <- d2H.spline.fits$y
    # make data frame for list element...
    output.spline.fits[[i]] <- data.frame("time"=temp1,"H2O"=temp2,"d18O"=temp3,"d2H"=temp4)
  }

# return spline values.
return(output.spline.fits)
}

#------------------------------------------------------------------------------------------
# calculate.spline.derivatives - simple function taking a spline list generated by fit.calibration.splines
# and calculating the numerical derivative. 

calculate.spline.derivatives <- function(spline.fits,breaks) {
  print(paste(Sys.time()," calculating time derivatives from splines..."))

  # initialize/preallocate output variable
  spline.derivatives <- vector("list",length(spline.fits))

  # loop through plateaus and calculate derivative, output new list including
  # derivatives...
  for (i in 1:length(spline.fits)) {
    temp1 <- diff(spline.fits[[i]]$H2O)
    temp2 <- diff(spline.fits[[i]]$d18O)
    temp3 <- diff(spline.fits[[i]]$d2H)
    temp4 <- spline.fits[[i]]$time[2:(length(spline.fits[[i]]$time))] # keep time associated with the derivatives...
    temp5 <- (breaks[i]+1):(breaks[i+1]-1)
    temp6 <- c(NA,diff(temp1))
    temp7 <- c(NA,diff(temp2))
    temp8 <- c(NA,diff(temp3))
    spline.derivatives[[i]] <- data.frame("d_H2O"=temp1,"d_d18O"=temp2,"d_d2H"=temp3,"time"=temp4,"inds"=temp5,
      "d2_H2O"=temp6,"d2_d18O"=temp7,"d2_d2H"=temp8)
  }

  # return output variable
  return(spline.derivatives)
}

#------------------------------------------------------------------------------------------
# extract.stable.calib.indices -
# using the splines and their derivatives, finds the set of indices
# that should be used to calculate mean/standard deviation values for the analysis. This function
# is somewhat large, but a lot is done under the hood and doesn't need to be stored in memory 
# following calculation. 

extract.stable.calib.indices <- function(spline.derivatives,H2O.thres=5.0,d18O.thres=0.005,d2H.thres=0.01) {
  print(paste(Sys.time()," finding indices corresponding to stable measurements..."))

  # allocate output vector
  inds.to.keep <- vector("list",length(spline.derivatives))

  # loop through and pull out list of inds to keep for each "plateau"
  for (i in 1:length(spline.derivatives)) {
  indvec <- vector("logical",length(spline.derivatives[[i]]$inds))
  #print(indvec)
  indvec <- !(abs(spline.derivatives[[i]]$d_H2O) > H2O.thres |
    abs(spline.derivatives[[i]]$d_d18O) > d18O.thres |
    abs(spline.derivatives[[i]]$d_d2H) > d2H.thres)
  #print(indvec)
  inds.to.keep[[i]] <- spline.derivatives[[i]]$inds[indvec]
  }
  return(inds.to.keep)
}

#------------------------------------------------------------------------------------------
# calculate.standard.averages - 

calculate.standard.averages <- function(calib.data,retained.indices) {

  # calculate mean values and uncertainty from identified plateaus.
  temp1 <- vector("numeric",length(retained.indices))
  temp2 <- vector("numeric",length(retained.indices))
  temp3 <- vector("numeric",length(retained.indices))
  temp4 <- vector("numeric",length(retained.indices))
  temp5 <- vector("numeric",length(retained.indices))
  temp6 <- vector("numeric",length(retained.indices))
  temp7 <- vector("numeric",length(retained.indices))
  temp8 <- vector("numeric",length(retained.indices))
  temp9 <- vector("numeric",length(retained.indices))

  for (i in 1:length(retained.indices)) {
    temp1[i] <- mean(calib.data$H2O[retained.indices[[i]]])
    temp2[i] <- mean(calib.data$Delta_18_16[retained.indices[[i]]])
    temp3[i] <- mean(calib.data$Delta_D_H[retained.indices[[i]]])

    temp4[i] <- sd(calib.data$H2O[retained.indices[[i]]])
    temp5[i] <- sd(calib.data$Delta_18_16[retained.indices[[i]]])
    temp6[i] <- sd(calib.data$Delta_D_H[retained.indices[[i]]])

    temp7[i] <- mean(calib.data$EPOCH_TIME[retained.indices[[i]]])
    temp8[i] <- min(calib.data$EPOCH_TIME[retained.indices[[i]]])
    temp9[i] <- max(calib.data$EPOCH_TIME[retained.indices[[i]]])
  }
  
  std.avgs <- data.frame("H2O.mean"=temp1,"d18O.mean"=temp2,"d2H.mean"=temp3,
    "H2O.sd"=temp4,"d18O.sd"=temp5,"d2H.sd"=temp6,
    "time.mean"=temp7,"time.min"=temp8,"time.max"=temp9)

  print(nrow(std.avgs))

  # filter out some obviously incorrect points.
  h2o.min <- which(std.avgs$H2O.mean > 2000) # must be 2000 ppm
  h2o.max <- which(std.avgs$H2O.mean < 30000)
  h2o.lsd <- which(std.avgs$H2O.sd < 1000)
  d18O.lsd <- which(std.avgs$d18O.sd < 0.5)
  d2H.lsd <- which(std.avgs$d2H.sd < 4)

  retain.inds <- Reduce(intersect,list(h2o.min,h2o.max,h2o.lsd,d18O.lsd,d2H.lsd))

  print(length(retain.inds))
  
  if (length(retain.inds) > 0) {
    std.avgs <- std.avgs[retain.inds,] 
  } else {
    std.avgs <- NULL
  }
  
  return(std.avgs)
}
