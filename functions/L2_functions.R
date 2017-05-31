# Picarro Water Vapor Isotope Data Reduction 
# L2 Functions

# migrated from WBB-vapor-scripts repo by rich fiorella on 22/2/17.
# code in these functions written by rich fiorella based on work by galen gorski
# (likely some of galen's functions have made it through with no changes!)

# contact info: rich fiorella, rich.fiorella@utah.edu

# this script contains a series of accessory functions for the L1 script.
# as of 22/2/17, this script has been designed to work only with a Picarro L2130.
# other analyzer versions will likely require modification of this script.

extract.date.from.L1.files <- function(file.list,dbg.level=0) {
  # print statement indicating that this function is starting
  if (dbg.level>0) {
    print("============================================")
    print("starting extract.date.from.L1.files function")
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
  # print statment denoting the end of this function
  if (dbg.level>0) {
    print("ending extract.date.from.L1.files function")
    print("==========================================")
  }
  # return date list
  return(file.name.dates.woNAs)
}

#------------------------------------------------------------------------------------------
# ID.calib.breakpoints - finds indices in the input calibration data frame where the
# gap between points exceeds threshold number of seconds (default = 10 seconds)

ID.calib.breakpoints <- function(calibration.data.frame,thres=10,dbg.level=0) {
  # print statement indicating that this function is starting
  if (dbg.level>0) {
    print("++++++++++++++++++++++++++++++++++++++")
    print("starting ID.calib.breakpoints function")
  }

  print(paste(Sys.time()," separating calibration periods based on gaps in time..."))
  if (!"EPOCH_TIME" %in% colnames(calibration.data.frame)) {
    stop("Expected time array absent from input data frame, check input data to ID.calib.breakpoints")
  }

  # rebuilding for version 1.1.0.
  # want to find: points with large differences in time or moving average H2O and
  # define as separate peaks.
  # note that diff returns the difference between point i and i+1 at index i
  # ...when pulling out the indices that indicate the start of the subsequent 
  # analysis, these are at index i+1, which has been returned at index i - so
  # need to add 1 to the indices returned by diff 
  #------------------------------------------------------
  # first, find time inds.
  #------------------------------------------------------
  time.diffs <- diff(calibration.data.frame$EPOCH_TIME)

  # find which inds are above time threshold.
  time.inds <- which(time.diffs>thres) + 1

  # caluclate moving average of H2O data frame.
  h2o.smth <- rollapply(calibration.data.frame$H2O,
      width=210, # corresponds to ~30 seconds of 1.16 Hz data
      FUN=median, # running mean
      fill=NA)  # require that output vector have smae length as input vector

  h2o.1st.diff <- c(diff(h2o.smth),NA)

  h2o.inds.tmp <- which(abs(h2o.1st.diff) > 100) + 1

  # need to trim h2o.inds! lots of consecutive points...
  h2o.inds <- h2o.inds.tmp[c(1,which(diff(h2o.inds.tmp)>210)+1,length(h2o.inds.tmp))]

  joined <- unique(sort(c(h2o.inds,time.inds)))

  # remove points that are essentially adjacent - this occurs when both the time and the
  # h2o diff filters detect a peak.
  all.breaks <- c(1,joined[which(diff(joined)>30)],nrow(calibration.data.frame))

  # return output value.
  return(all.breaks)
}

#------------------------------------------------------------------------------------------
# fit.calibration.splines - fits a smooth spline to each identified calibration period.

fit.calibration.splines <- function(calibration.data.frame,breakpoints,dfree=12,dbg.level=0) {
  # print statement indicating that this function is starting
  if (dbg.level>0) {
    print("*****************************************")
    print("starting fit.calibration.splines function")
  }

  #===================================================
  # start work of function
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
    H2O.spline.fits <- smooth.spline(calibration.data.frame$EPOCH_TIME[(breakpoints[i]):(breakpoints[i+1]-1)],
      calibration.data.frame$H2O[(breakpoints[i]):(breakpoints[i+1]-1)],df=dfree)
    
    d18O.spline.fits <- smooth.spline(calibration.data.frame$EPOCH_TIME[(breakpoints[i]):(breakpoints[i+1]-1)],
      calibration.data.frame$Delta_18_16[(breakpoints[i]):(breakpoints[i+1]-1)],df=dfree)
    
    d2H.spline.fits <- smooth.spline(calibration.data.frame$EPOCH_TIME[(breakpoints[i]):(breakpoints[i+1]-1)],
      calibration.data.frame$Delta_D_H[(breakpoints[i]):(breakpoints[i+1]-1)],df=dfree)
    
    # extract y value from splines, x value from 1 spline (should all be the same x values!)
    temp1 <- H2O.spline.fits$x
    temp2 <- H2O.spline.fits$y
    temp3 <- d18O.spline.fits$y
    temp4 <- d2H.spline.fits$y
    
    # make data frame for list element...
    output.spline.fits[[i]] <- data.frame("time"=temp1,"H2O"=temp2,"d18O"=temp3,"d2H"=temp4)
  }

  # print statment denoting the end of this function
  if (dbg.level>0) {
    print("ending fit.calibration.splines function")
    print("***************************************")
  }

  # return spline values.
  return(output.spline.fits)
}

#------------------------------------------------------------------------------------------
# calculate.spline.derivatives - simple function taking a spline list generated by fit.calibration.splines
# and calculating the numerical derivative. 

calculate.spline.derivatives <- function(spline.fits,breaks,dbg.level=0) {
  # print statement indicating that this function is starting
  if (dbg.level>0) {
    print("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&")
    print("starting calculate.spline.derivatives function")
  }

  #======================================================
  # start work of function
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

    if (dbg.level>1) {
      print(i)  
      print(paste(length(temp1),length(temp2),length(temp3),length(temp4),
        length(temp5),length(temp6),length(temp7),length(temp8)))
    }

    spline.derivatives[[i]] <- data.frame("d_H2O"=temp1,"d_d18O"=temp2,"d_d2H"=temp3,"time"=temp4,"inds"=temp5,
      "d2_H2O"=temp6,"d2_d18O"=temp7,"d2_d2H"=temp8)
  }

  # print statment denoting the end of this function
  if (dbg.level>0) {
    print("ending calculate.spline.derivatives function")
    print("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&")
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

# v110 change - a bit of a thorny issue where often at the beginning of an analysis, there's a local
# min/max that gets included in the peak. Most direct way to exclude this is just to restrict the
# values that are returned to be only in the last ~75% (tuned parameter!) of the peak.

extract.stable.calib.indices <- function(spline.derivatives,H2O.thres=5.0,d18O.thres=0.01,d2H.thres=0.1,
  cut.in.span=0.20,cut.out.span=0.05,dbg.level=0) {
  # print statement indicating that this function is starting
  if (dbg.level>0) {
    print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
    print("starting extract.stable.calib.indices function")
  }

  #======================================================
  # start work of function

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

    # set all values in the cutin and cutout periods to false!
    # --------------------------------------------
    # get nearest index to cutin period fraction.
    # floor ensures return value is an integer!
    cutin.max <- floor(cut.in.span*length(spline.derivatives[[i]]$inds))
    # set values in cut in region to false.
    indvec[1:cutin.max] <- FALSE

    # now repeat for cutout span.
    cutout.min <- ceiling((1-cut.out.span)*length(spline.derivatives[[i]]$inds))
    # set values in cut.out region to false.
    indvec[cutout.min:length(indvec)] <- FALSE 

    # store TRUE values in inds.to.keep list for this measurement, 
    # and move on to the next measurement.
    inds.to.keep[[i]] <- spline.derivatives[[i]]$inds[indvec]
  }
  # print statment denoting the end of this function
  if (dbg.level>0) {
    print("ending extract.stable.calib.indices function")
    print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
  }
  # return list of inds to keep for each peak.
  return(inds.to.keep)
}

#------------------------------------------------------------------------------------------
# calculate.standard.averages - 

calculate.standard.averages <- function(calib.data,retained.indices,memory.filter=TRUE,dbg.level=0) {
  # print statement indicating that this function is starting
  if (dbg.level>0) {
    print("=============================================")
    print("starting calculate.standard.averages function")
  }

  #======================================================
  # start work of function

  print(paste(Sys.time()," Calculating average/stdev/ranges of identified peaks..."))

  # calculate mean values and uncertainty from identified plateaus.
  temp1 <- vector("numeric",length(retained.indices)) # mean H2O
  temp2 <- vector("numeric",length(retained.indices)) # mean d18O
  temp3 <- vector("numeric",length(retained.indices)) # mean d2H
  temp4 <- vector("numeric",length(retained.indices)) # std dev H2O
  temp5 <- vector("numeric",length(retained.indices)) # std dev d18O
  temp6 <- vector("numeric",length(retained.indices)) # std dev d2H
  temp7 <- vector("numeric",length(retained.indices)) # mean time
  temp8 <- vector("numeric",length(retained.indices)) # min time
  temp9 <- vector("numeric",length(retained.indices)) # max time
  temp10 <- vector("numeric",length(retained.indices)) # H2O range
  temp11 <- vector("numeric",length(retained.indices)) # d18O range
  temp12 <- vector("numeric",length(retained.indices)) # d2H range
  temp13 <- vector("numeric",length(retained.indices)) # count of points retained.
  temp14 <- vector("numeric",length(retained.indices)) # d18O slope
  temp15 <- vector("numeric",length(retained.indices)) # d18O rsq
  temp16 <- vector("numeric",length(retained.indices)) # d2H slope
  temp17 <- vector("numeric",length(retained.indices)) # d2H rsq
  temp18 <- vector("numeric",length(retained.indices))

  for (i in 1:length(retained.indices)) {
    if (memory.filter==TRUE) {
      print("attempting to filter out for memory effects...")
      # --------------------------
      # NEW MEMORY CORRECTION

      # estimate the slope and strength of the relationship for the peak of interest.
      modtime <- calib.data$EPOCH_TIME[retained.indices[[i]]]-calib.data$EPOCH_TIME[retained.indices[[i]][1]]
      Omod <- lm(calib.data$Delta_18_16[retained.indices[[i]]] ~ modtime)
      Hmod <- lm(calib.data$Delta_D_H[retained.indices[[i]]] ~ modtime)

      Omem.slope <- 60*coef(Omod)[[2]] # convert from permil/sec to permil/min
      Omem.rsq <- summary(Omod)$r.squared

      Hmem.slope <- 60*coef(Hmod)[[2]] # convert from permil/sec to permil/min
      Hmem.rsq <- summary(Hmod)$r.squared

      # how many points are there to begin with in this data frame?
      npts <- length(retained.indices[[i]])
      tnth <- npts %/% 10 # how many points correspond to 1/10 of the data?

      z <- 0  # q will keep track of how many tenths of data have been discarded. algorithm stops
              # when 80% of data has been discarded.

      # set logical conditions that determine whether memory is present in the sample
      conds <- (abs(Omem.slope) > 0.1 | abs(Hmem.slope) > 0.5) & z < 8
      
      print(paste(abs(Omem.slope),abs(Hmem.slope),z))

      while (conds) {
        # cut out first 10% of data in the timeseries.
        new.inds <- retained.indices[[i]][((z+1)*tnth):npts]
        
       # print(new.inds)

        # rerun the models.
        modtime <- calib.data$EPOCH_TIME[new.inds]-calib.data$EPOCH_TIME[new.inds[1]]
        Omod <- lm(calib.data$Delta_18_16[new.inds] ~ modtime)
        Hmod <- lm(calib.data$Delta_D_H[new.inds] ~ modtime)

        Omem.slope <- 60*coef(Omod)[[2]] # convert from permil/sec to permil/min
        Omem.rsq <- summary(Omod)$r.squared

        Hmem.slope <- 60*coef(Hmod)[[2]] # convert from permil/sec to permil/min
        Hmem.rsq <- summary(Hmod)$r.squared
        
        # increment q
        z <- z + 1

        # reëvaluate conds.
        conds <- (abs(Omem.slope) > 0.1 | abs(Hmem.slope) > 0.5) & z < 8

      }
 
      # set mem.retained.inds for what follows.
      if (z>0) {
        mem.retained.inds <- new.inds
      } else {
        mem.retained.inds <- retained.indices[[i]] # no memory correction necessayr.
      }
      
      #---------------------------------------------------------------------
      temp1[i] <- mean(calib.data$H2O[mem.retained.inds],na.rm=TRUE)
      temp2[i] <- mean(calib.data$Delta_18_16[mem.retained.inds],na.rm=TRUE)
      temp3[i] <- mean(calib.data$Delta_D_H[mem.retained.inds],na.rm=TRUE)

      temp4[i] <- sd(calib.data$H2O[mem.retained.inds],na.rm=TRUE)
      temp5[i] <- sd(calib.data$Delta_18_16[mem.retained.inds],na.rm=TRUE)
      temp6[i] <- sd(calib.data$Delta_D_H[mem.retained.inds],na.rm=TRUE)

      temp7[i] <- mean(calib.data$EPOCH_TIME[mem.retained.inds],na.rm=TRUE)
      temp8[i] <- min(calib.data$EPOCH_TIME[mem.retained.inds],na.rm=TRUE)
      temp9[i] <- max(calib.data$EPOCH_TIME[mem.retained.inds],na.rm=TRUE)

      temp10[i] <- max(calib.data$H2O[mem.retained.inds],na.rm=TRUE)-min(calib.data$H2O[mem.retained.inds],na.rm=TRUE)
      temp11[i] <- max(calib.data$Delta_18_16[mem.retained.inds],na.rm=TRUE)-min(calib.data$Delta_18_16[mem.retained.inds],na.rm=TRUE)
      temp12[i] <- max(calib.data$Delta_D_H[mem.retained.inds],na.rm=TRUE)-min(calib.data$Delta_D_H[mem.retained.inds],na.rm=TRUE)

      temp13[i] <- length(mem.retained.inds)

      modtime <- calib.data$EPOCH_TIME[mem.retained.inds]-calib.data$EPOCH_TIME[mem.retained.inds[1]]
      Omod <- lm(calib.data$Delta_18_16[mem.retained.inds] ~ modtime)
      Hmod <- lm(calib.data$Delta_D_H[mem.retained.inds] ~ modtime)
    
      temp14[i] <- 60*coef(Omod)[[2]] # convert from permil/s to permil/min
      temp15[i] <- summary(Omod)$r.squared

      temp16[i] <- 60*coef(Hmod)[[2]]
      temp17[i] <- summary(Hmod)$r.squared

      temp18[i] <- 10*z # percentage of points removed!

    } else { # DO NOT remove any memory issues.
      if (lengths(good.inds)[i] > 0 & !all(is.na(good.inds[[i]]))) {
        #print(i) # debug

        #print("Not attempting to correct for memory effects...")
        temp1[i] <- mean(calib.data$H2O[retained.indices[[i]]],na.rm=TRUE)
        temp2[i] <- mean(calib.data$Delta_18_16[retained.indices[[i]]],na.rm=TRUE)
        temp3[i] <- mean(calib.data$Delta_D_H[retained.indices[[i]]],na.rm=TRUE)

        temp4[i] <- sd(calib.data$H2O[retained.indices[[i]]],na.rm=TRUE)
        temp5[i] <- sd(calib.data$Delta_18_16[retained.indices[[i]]],na.rm=TRUE)
        temp6[i] <- sd(calib.data$Delta_D_H[retained.indices[[i]]],na.rm=TRUE)

        temp7[i] <- mean(calib.data$EPOCH_TIME[retained.indices[[i]]],na.rm=TRUE)
        temp8[i] <- min(calib.data$EPOCH_TIME[retained.indices[[i]]],na.rm=TRUE)
        temp9[i] <- max(calib.data$EPOCH_TIME[retained.indices[[i]]],na.rm=TRUE)

        temp10[i] <- max(calib.data$H2O[retained.indices[[i]]],na.rm=TRUE)-min(calib.data$H2O[retained.indices[[i]]],na.rm=TRUE)
        temp11[i] <- max(calib.data$Delta_18_16[retained.indices[[i]]],na.rm=TRUE)-min(calib.data$Delta_18_16[retained.indices[[i]]],na.rm=TRUE)
        temp12[i] <- max(calib.data$Delta_D_H[retained.indices[[i]]],na.rm=TRUE)-min(calib.data$Delta_D_H[retained.indices[[i]]],na.rm=TRUE)

        temp13[i] <- length(retained.indices[[i]])

        modtime <- calib.data$EPOCH_TIME[retained.indices[[i]]]-calib.data$EPOCH_TIME[retained.indices[[i]][1]]
        Omod <- lm(calib.data$Delta_18_16[retained.indices[[i]]] ~ modtime)
        Hmod <- lm(calib.data$Delta_D_H[retained.indices[[i]]] ~ modtime)
      
        temp14[i] <- 60*coef(Omod)[[2]] # convert from permil/s to permil/min
        temp15[i] <- summary(Omod)$r.squared

        temp16[i] <- 60*coef(Hmod)[[2]]
        temp17[i] <- summary(Hmod)$r.squared

        temp18[i] <- NA # memory not removed!!!
      } else {
        #print("No inds for this period...")
      } # end check for good inds.
    }
  }
  
  std.avgs <- data.frame("H2O.mean"=temp1,"d18O.mean"=temp2,"d2H.mean"=temp3,
    "H2O.sd"=temp4,"d18O.sd"=temp5,"d2H.sd"=temp6,
    "time.mean"=temp7,"time.min"=temp8,"time.max"=temp9,
    "H2O.range"=temp10,"d18O.range"=temp11,"d2H.range"=temp12,
    "ind.count"=temp13,"d18O.trend"=temp14,"d18O.r2.of.trend"=temp15,
    "d2H.trend"=temp16,"d2H.r2.of.trend"=temp17,"pct.removed.memory"=temp18)

  # kludge fix - remove any rows where all the values are zero? not sure why this 
  # is happening.
  pts <- which(std.avgs$time.mean==0)

  if (length(pts)>0) {std.avgs <- std.avgs[-pts,]}

  # print statment denoting the end of this function
  if (dbg.level>0) {
    print("ending calculate.standard.averages function")
    print("===========================================")
  }
  
  # return data frame of average values for the standards
  return(std.avgs)
}

#-----------------------------------------------------------------------------------
# get.ambient.deltas function - 

get.ambient.deltas <- function(calib.averages,ambient.data,dbg.level=0) {
  # print statement indicating that this function is starting
  if (dbg.level>0) {
    print("++++++++++++++++++++++++++++++++++++")
    print("starting get.ambient.deltas function")
  }

  #======================================================
  # start work of function

  print(paste(Sys.time(),"Estimate delta of vapor coming through drierite canister before/after analysis..."))

  # first, check to see if the month has any calibration averages
  # provided to it.
  if (length(calib.averages) == 0 ) {
    warning("No calibration data in data frame passed to get.ambient.deltas...")
    output <- list()
    return(output)
    break
  }

  # get time vars, length of calib averages array
  time.inds <- calib.averages[c("time.min","time.max")]

  # ambient.time
  amb.time <- ambient.data[c("EPOCH_TIME")]

  # loop through rows of time.inds and find the index of
  # nearest point before time.min and the nearest point after time.max
  before.ind <- vector("numeric",nrow(time.inds))
  after.ind <- vector("numeric",nrow(time.inds))

  for (i in 1:nrow(time.inds)) {
      # should be positive for points before calibration point
      tmp1 <- time.inds$time.min[i] - amb.time
      
      # find minimum positive point, assuming any points are positive
      before.ind[i] <- ifelse(any(tmp1>0), # conditional test
        which(tmp1==min(tmp1[tmp1>0],na.rm=TRUE)), # true branch
        NA) # false branch
      # reverse the logic for the point after calibration
      tmp2 <- amb.time - time.inds$time.max[i]
      
      # first point after calibration using eqn above should be 
      # minimum positive number  
      after.ind[i] <- ifelse(any(tmp2>0), # conditional test
        which(tmp2==min(tmp2[tmp2>0],na.rm=TRUE)), # true branch
        NA) # false branch

      # clean up loop
      rm(tmp1)
      rm(tmp2)
  }

  # now need to get characteristics of period before/after calibration
  nmin <- 5 # number of ambient minutes to grab before/after calibration

  before.H2O <- vector("numeric",nrow(time.inds))
  before.18O <- vector("numeric",nrow(time.inds))
  before.2H  <- vector("numeric",nrow(time.inds))

  after.H2O <- vector("numeric",nrow(time.inds))
  after.18O <- vector("numeric",nrow(time.inds))
  after.2H  <- vector("numeric",nrow(time.inds))

  for (i in 1:nrow(time.inds)) {
    # calculate characteristics of ambinet data before calibration
    if (!is.na(before.ind[i]) & before.ind[i]>(nmin-1)) {
      tmp1.inds <- (before.ind[i]-nmin+1):before.ind[i]
      before.H2O[i] <- mean(ambient.data$H2O[tmp1.inds])
      before.18O[i] <- mean(ambient.data$Delta_18_16[tmp1.inds])
      before.2H[i] <- mean(ambient.data$Delta_D_H[tmp1.inds])
    } else {
      before.H2O[i] <- before.18O[i] <- before.2H[i] <- NA
    }
    
    # calculate charateristics of ambient data after calibration
    if (!is.na(after.ind[i]) & (after.ind[i]+nmin-1) <= nrow(ambient.data)) {
      tmp2.inds <- after.ind[i]:(after.ind[i]+nmin-1)
      after.H2O[i] <- mean(ambient.data$H2O[tmp2.inds])
      after.18O[i] <- mean(ambient.data$Delta_18_16[tmp2.inds])
      after.2H[i] <- mean(ambient.data$Delta_D_H[tmp2.inds])
    } else {
      after.H2O[i] <- after.18O[i] <- after.2H[i] <- NA
    }
  }

  # make the data output data frame
  output <- data.frame("before.H2O"=before.H2O,"before.d18O"=before.18O,
    "before.d2H"=before.2H,"after.H2O"=after.H2O,
    "after.d18O"=after.18O,"after.d2H"=after.2H)

  # print statment denoting the end of this function
  if (dbg.level>0) {
    print("ending get.ambient.deltas function")
    print("++++++++++++++++++++++++++++++++++")
  }

  # return output frame
  return(output)
}

#-----------------------------------------------------------------------------------
# apply.mixingratio.correction - this function corrects the measured isotopic compositions
# using a linear regression between delta and 1/[H2O] - based on coefficients determined from
# CalibrationStudyRegressionHandClean_gjb. Requires mixing ratio coefficients to be specified 
# in L1_Calibration_Parameters.R. Other regressions are possible for correcting for
# concentration dependence, but have not be implemented in this code base yet.

apply.mixingratio.correction <- function(avg.data.frame,fit.type,Oslope,Hslope,dbg.level=0) {
  # print statement indicating that this function is starting
  if (dbg.level>0) {
    print("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&")
    print("starting apply.mixingratio.correction.calibration function")
  }

  #======================================================
  # start work of function

  print(paste(now()," Applying mixing ratio correction to: ",deparse(substitute(avg.data.frame))))

  # apply mixing ratio corrections - this takes coeffiecients specified in L1_Calibration_Parameters.R
  if (fit.type=="hyperbolic") {
    Delta_18_16_mrc <- Oslope*(1/20000-1/avg.data.frame$H2O.mean) + avg.data.frame$d18O.mean 
    Delta_D_H_mrc <- Hslope*(1/2000-1/avg.data.frame$H2O.mean) + avg.data.frame$d2H.mean
  } else if (fit.type=="logarithmic") { 
    Delta_18_16_mrc <- Oslope*(log(1/20000,base=10)-log(1/avg.data.frame$H2O.mean,base=10)) +
      avg.data.frame$d18O.mean
    Delta_D_H_mrc <- Hslope*(log(1/20000,base=10)-log(1/avg.data.frame$H2O.mean,base=10)) +
      avg.data.frame$d2H.mean
    }

  # attach new mrc variables to original data frame.
  avg.data.frame <- cbind(avg.data.frame,Delta_18_16_mrc)
  avg.data.frame <- cbind(avg.data.frame,Delta_D_H_mrc)

  # print statment denoting the end of this function
  if (dbg.level>0) {
    print("ending apply.mixingratio.correction.calibration function")
    print("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&")
  }

  # return the whole data frame with the two mixing ratio corrected variables attached.
  return(avg.data.frame)
}

#-----------------------------------------------------------------------------------
# apply.drygas.correction function - corrected now that ambient data is already in data frame.

apply.drygas.correction <- function(data,do.correction=TRUE,H2O.bg,include.gypsum.fractionation=FALSE,dbg.level=0) {
  # print statement indicating that this function is starting
  if (dbg.level>0) {
    print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
    print("starting apply.drygas.correction function")
  }

  #======================================================
  # start work of function
  print(paste(now()," Applying dry gas correction to calibration data..."))

  if (include.gypsum.fractionation==FALSE) {
    # apply corrections using equation S3 of Gorski et al 2014.
    Delta_18_16_bgc <- (data$Delta_18_16_mrc*data$H2O.mean -
      data$before.d18O*H2O.bg)/(data$H2O.mean-H2O.bg)
    Delta_D_H_bgc <- (data$Delta_D_H_mrc*data$H2O.mean -
      data$before.d2H*H2O.bg)/(data$H2O.mean-H2O.bg)

    # temporary kludge fix - some of these data points won't have
    # ambient data available before, but will have data available after.
    # so, for those rows that returned NA previously, go back
    # and fill those using ambient data from after the analysis
    Delta_18_16_bgc[is.na(Delta_18_16_bgc)] <- 
      (data$Delta_18_16_mrc[is.na(Delta_18_16_bgc)]*data$H2O.mean[is.na(Delta_18_16_bgc)] -
      data$after.d18O[is.na(Delta_18_16_bgc)]*H2O.bg)/
      (data$H2O.mean[is.na(Delta_18_16_bgc)]-H2O.bg)
    
    Delta_D_H_bgc[is.na(Delta_D_H_bgc)] <- 
      (data$Delta_D_H_mrc[is.na(Delta_D_H_bgc)]*data$H2O.mean[is.na(Delta_D_H_bgc)] -
      data$after.d2H[is.na(Delta_D_H_bgc)]*H2O.bg)/
      (data$H2O.mean[is.na(Delta_D_H_bgc)]-H2O.bg)

    # add these columns to the original data frame
    data <- cbind(data,Delta_18_16_bgc)
    data <- cbind(data,Delta_D_H_bgc)
    

    # print statment denoting the end of this function
    if (dbg.level>0) {
      print("ending apply.drygas.correction function")
      print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
    }

    # return data frame
    return(data)
  } else if (include.gypsum.fractionation==TRUE) {
    # much more complicated function here than without removing gypsum fractionation

    alpha.CaSO4.hydwater.oxygen <-  1.004 # source: Fontes, Gonfiantini and Fontes     
    alpha.CaSO4.hydwater.hydrogen <- 0.9825 # source: Horita EPSL

    # apply oxygen correction
    #--------------------------
    oxygen.numerator <- data$Delta_18_16_mrc*data$H2O.mean -
      (data$before.d18O*data$before.H2O)-(data$before.H2O-H2O.bg)*
      (alpha.CaSO4.hydwater.oxygen*(data$before.d18O+1000)-1000)

    # kludge fix for missing values.
    oxygen.numerator[is.na(oxygen.numerator)] <- data$Delta_18_16_mrc*data$H2O.mean -
      (data$after.d18O*data$after.H2O)-(data$after.H2O-H2O.bg)*
      (alpha.CaSO4.hydwater.oxygen*(data$after.d18O+1000)-1000)

    # calculate background corrected d18O
    Delta_18_16_bgc <- oxygen.numerator/(data$H2O.mean-H2O.bg)
  
    # apply hydrogen correction
    #----------------------------
    hydrogen.numerator <- data$Delta_D_H_mrc*data$H2O.mean -
      (data$before.d2H*data$before.H2O)-(data$before.H2O-H2O.bg)*
      (alpha.CaSO4.hydwater.hydrogen*(data$before.d2H+1000)-1000)

    # kludge fix for missing values.
    hydrogen.numerator[is.na(hydrogen.numerator)] <- data$Delta_D_H_mrc*data$H2O.mean -
      (data$after.d2H*data$after.H2O)-(data$after.H2O-H2O.bg)*
      (alpha.CaSO4.hydwater.hydrogen*(data$after.d2H+1000)-1000)

    # calculate background corrected d2H
    Delta_D_H_bgc <- oxygen.numerator/(data$H2O.mean-H2O.bg)

    # add these columns to the original data frame
    data <- cbind(data,Delta_18_16_bgc)
    data <- cbind(data,Delta_D_H_bgc)

    # print statment denoting the end of this function
    if (dbg.level>0) {
      print("ending apply.drygas.correction function")
      print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
    }

    # return data frame
    return(data)
  }
}

#-----------------------------------------------------------------------------------
# correct.standards.to.VSMOW function

correct.standards.to.VSMOW <- function(standard.data.frame,method=1,dbg.level=0) {
  # print statement indicating that this function is starting
  if (dbg.level>0) {
    print("============================================")
    print("starting correct.standards.to.VSMOW function")
  }

  #======================================================
  # start work of function

  print(paste(now()," Converting calibration data to VSMOW..."))

  # methods for calibrating:
  # 1 - bracket each ambient period with the standard measurements immediately
  #     before and after ambient period. equally weighted.
  # 2 - bracket each ambient period with the standard measurements immediately
  #     before and after ambient period - weighted by inverse of measurement standard error. (NOT CODED)
  # 3 - larger stencil (NOT CODED)
  # 4 - Bayesian regression? (NOT CODED)

  if (method==1) {
    # ensure that we're working with a data frame...
    stopifnot(is.data.frame(standard.data.frame)) # is this actually necessary?

    # does the data frame actually have data?
    if (nrow(standard.data.frame)==0) {
      return(NULL)
    }

    # find the number of calibration periods in the data frame.
    #----------------------------------------------------------
    # however, cannot simply treat these as distinct periods as
    # the same standard could be analyzed multiple times sequentially
    # either by choice (e.g., as in a concentration calibration)
    # or by error (e.g., one of the SDM needles develops a clog, and is 
    # not immediately noticed and replaced or is otherwise inaccesible)
    # really - want to do something called run length encoding (rle)

    std.names <- as.vector(standard.data.frame$standard.name) # rle requires atomic vector
   
    # check to see if there are any NA values that might interfere with
    # calibrations...
    if (any(is.na(std.names))) { 
      print("There are NAs in assigned standards...check this output very carefully!!!!")
    }

    # run rle - rle gives a list with two outputs: (a) the length of each run and
    # (b) the standard associated with each run (values list).
    std.periods <- rle(std.names) 

    # ok, what we want now is to define individual calibration periods.
    # a calibration period requires there to be measurements of 2 different
    # standards. %/%2 operator requires 2 standards at each period, the -1
    # accounts for the fact we're interested in the periods bracketing 
    # a certain portion of time.
    nperiods <- length(std.periods$lengths)%/% 2 - 1

    # assign each index of std.periods$lengths to a bracket -
    # bracket i will refer to the beginning of period i, bracket i+1
    # will refer to the end of period i
    # (indexing here a little funky, but starts at 2 so that first two
    # elements of brackets both return 1)
    brackets <- (seq(2,length(std.periods$lengths)+1))%/%2 

    # initiate output vectors
    period.begin <- vector("numeric",nperiods)
    period.end   <- vector("numeric",nperiods)
    period.Oslope <- vector("numeric",nperiods)
    period.Ointercept <- vector("numeric",nperiods)
    period.Orsq <- vector("numeric",nperiods)
    period.Hslope <- vector("numeric",nperiods)
    period.Hintercept <- vector("numeric",nperiods)
    period.Hrsq <- vector("numeric",nperiods)
    period.ID <- vector("numeric",nperiods)
    start.points <- vector("numeric",nrow(standard.data.frame))
    end.points <- vector("numeric",nrow(standard.data.frame))

    # loop through the periods and estimate the slope.
    for (i in 1:nperiods) {
      # create an vector the length of std.periods$lengths that matches
      # i to the calibration period.

      bracket.i.lims <- c(min(which(brackets==i)),max(which(brackets==i)))
      bracket.iplus1.lims <- c(min(which(brackets==(i+1))),max(which(brackets==(i+1))))

      # find the indices for the beginning portion of the bracket.
      start.i <- ifelse(i==1,1,sum(std.periods$lengths[1:(bracket.i.lims[1]-1)])+1)
      stop.i <- sum(std.periods$lengths[1:bracket.i.lims[2]])

      # find the indices for the end portion of the brakcet.
      start.iplus1 <- stop.i + 1
      stop.iplus1 <- sum(std.periods$lengths[1:bracket.iplus1.lims[2]])

      # calculate mean times for start and end of period
      period.begin[i] <- mean(standard.data.frame$time.mean[start.i:stop.i])
      period.end[i]   <- mean(standard.data.frame$time.mean[start.iplus1:stop.iplus1])

      # fit model for O slope
      Omod <- lm(standard.O18.VSMOW ~ Delta_18_16_bgc,data=standard.data.frame[start.i:stop.iplus1,])
      # extract O slope, intercept, and rsquared
      period.Oslope[i] <- coef(Omod)[[2]]
      period.Ointercept[i] <- coef(Omod)[[1]]
      period.Orsq[i] <- summary(Omod)$r.squared

      # fit model for H slope
      Hmod <- lm(standard.H2.VSMOW ~ Delta_D_H_bgc,data=standard.data.frame[start.i:stop.iplus1,])
      # extract O slope, intercept, and rsquared
      period.Hslope[i] <- coef(Hmod)[[2]]
      period.Hintercept[i] <- coef(Hmod)[[1]]
      period.Hrsq[i] <- summary(Hmod)$r.squared

      # set period.id as number of period.
      period.ID[i] <- i

      # figure out which points correspond to this period in the standard data frame.
      start.points[start.i:stop.i] <- i
      end.points[start.iplus1:stop.iplus1] <- i

    }
   
    # structure output into a list of data frames.
    #--------------------------------------------------

    # 1. uncalibrated standard data
    calibration.raw.data <- data.frame("standard"=standard.data.frame$standard.name,
      "time.mean"=standard.data.frame$time.mean,
      "time.min"=standard.data.frame$time.mean,
      "time.max"=standard.data.frame$time.min,
      "ind.count"=standard.data.frame$time.max,
      "H2O.mean.uncal"=standard.data.frame$H2O.mean,
      "H2O.stdev.uncal"=standard.data.frame$H2O.sd,
      "H2O.range.uncal"=standard.data.frame$H2O.range,
      "d18O.mean.uncal"=standard.data.frame$d18O.mean,
      "d18O.stdev.uncal"=standard.data.frame$d18O.sd,
      "d18O.range.uncal"=standard.data.frame$d18O.range,
      "d2H.mean.uncal"=standard.data.frame$d2H.mean,
      "d2H.stdev.uncal"=standard.data.frame$d2H.sd,
      "d2H.range.uncal"=standard.data.frame$d2H.range)

    # 2. ambient vapor data just before and after standard analysis
    ambient.correction.data <- data.frame("standard"=standard.data.frame$standard.name,
      "time.mean"=standard.data.frame$time.mean,
      "before.H2O.mean"=standard.data.frame$before.H2O,
      "before.d18O.mean.uncal"=standard.data.frame$before.d18O,
      "before.d2H.mean.uncal"=standard.data.frame$before.d2H,
      "after.H2O.mean"=standard.data.frame$after.H2O,
      "after.d18O.mean.uncal"=standard.data.frame$after.d18O,
      "after.d2H.mean.uncal"=standard.data.frame$after.d2H)

    # 3. diagnostics on the memory correction
    memory.correction.diagnostics <- data.frame("standard"=standard.data.frame$standard.name,
      "time.mean"=standard.data.frame$time.mean,
      "d18O.trend"=standard.data.frame$d18O.trend,
      "d18O.r2.of.trend"=standard.data.frame$d18O.r2.of.trend,
      "d2H.trend"=standard.data.frame$d2H.trend,
      "d2H.r2.of.trend"=standard.data.frame$d2H.r2.of.trend,
      "pct.removed.memory"=standard.data.frame$pct.removed.memory)

    # 4. calibrated data.
    calib.stds <- data.frame("standard"=standard.data.frame$standard.name,
      "time.mean"=standard.data.frame$time.mean,
      "H2O.mean.uncal"=standard.data.frame$H2O.mean,
      "d18O.mrc"=standard.data.frame$Delta_18_16_mrc,
      "d18O.mrcbgc"=standard.data.frame$Delta_18_16_bgc,
      "d2H.mrc"=standard.data.frame$Delta_D_H_mrc,
      "d2H.mrcbgc"=standard.data.frame$Delta_D_H_bgc,
      "starts.period"=start.points,
      "ends.period"=end.points,
      "known.std.d18O"=standard.data.frame$standard.O18.VSMOW,
      "known.std.d2H"=standard.data.frame$standard.H2.VSMOW)

    # 5. regression data
    regression.data <- data.frame("period.id"=period.ID,
      "period.start"=period.begin,
      "period.end"=period.end,
      "O.slope"=period.Oslope,
      "O.intercept"=period.Ointercept,
      "O.r2"=period.Orsq,
      "H.slope"=period.Hslope,
      "H.intercept"=period.Hintercept,
      "H.r2"=period.Hrsq,
      "qflag"=ifelse(period.Orsq < 0.95 | period.Hrsq < 0.95,0,1))

    # create dataframe packaging this info out to return
    output <- list("calibration.raw.data"=calibration.raw.data,
      "ambient.correction.data"=ambient.correction.data,
      "memory.correction.diagnostics"=memory.correction.diagnostics,
      "calib.stds"=calib.stds,
      "regression.data"=regression.data)

    # return the list.
    return(output)

  } # end method
} # end function

#----------------------------------------------------
# create attach.L2.Header function
# this function creates the header for datafiles and attaches the metadata
# loaded into a csv file.

attach.L2.Header <- function(output_filename,metadata_dataframe,dbg.level=0) {
  # print statement indicating that this function is starting
  if (dbg.level>0) {
    print("==================================")
    print("starting attach.L2.Header function")
  }

  # initiate writing out to output file...
  sink(output_filename)

  # add L1 information to top of output file...
  cat("# BEGIN L2 HEADER \n")

  # note when L2 script was run
  cat(paste("# date.L2.scripts.run, ",base::date(),"\n"))

  # what amount of time difference between points qualifies to (roughly) divide analyses?
  cat(paste("# time.threshold.for.separate.standard.runs, ",time.threshold,"\n"))

  # how many degrees of freedom should be used in the stiff spline?
  cat(paste("# stiff.spline.dfs, ",stiff.spline.dfs,"\n"))

  # what threshold value was used to identify an unacceptably large trend in H2O?
  cat(paste("# spline.derivative.H2O.threshold (ppm/sec), ",H2O.thres,"\n"))

  # what threshold value was used to identify an unacceptably large trend in d18O?
  cat(paste("# spline.derivative.d18O.threshold (permil/sec), ",d18O.thres,"\n"))

  # what threshold value was used to identify an unacceptably large trend in d2H?
  cat(paste("# spline.derivative.d2H.threshold (permil/sec), ",d2H.thres,"\n"))

  # did we attempt to remove potential effects of memory?
  cat(paste("# memory.filter.on to calculate averages, ",memory.filter,"\n"))
  # future addition - what thresholds were used????

  # did we correct for incomplete drying of "dry gas" as a mixing equation?
  cat(paste("# dry.gas.correction.applied, ",do.correction,"\n"))
  if (do.correction==TRUE) {
    # if yes, what value of H2O ppmv did we assume made it through the column?
    cat(paste("# background H2O assumed, ",H2O.bg,"\n"))
    # did we account for fractionation during hydration of gypsum? (typically a small effect)
    cat(paste("# gypsum fractionation included in drierite correction, ",include.gypsum.fractionation,"\n"))
  }

  # what type of (linear) function was used to correct for delta dependence on humidity?
  cat(paste("# humidity correction fit type, ",fit.type,"\n"))

  # what was the value of the delta-humidity correction slope for d18O?
  cat(paste("# humidity correction oxygen slope, ",Oslope, "\n"))

  # what was the value of the delta-humidity correction slope for d2H?
  cat(paste("# humidity correction hydrogen slope, ",Hslope,"\n"))

  # summary of filter values used to remove "bad" peaks -------------------------
  # maximum H2O value allowed?
  cat(paste("# maximum mean H2O allowed in peaks, ",h2o.max.thres,"\n"))

  # minimum H2O value allowed?
  cat(paste("# minimum mean H2O allowed in peaks, ",h2o.min.thres,"\n"))

  # maximum sd H2O value allowed?
  cat(paste("# maximum H2O standard deviation allowed in peaks, ",h2o.sdev.thres,"\n"))

  # maximum sd d18O value allowed?
  cat(paste("# maximum d18O standard deviation allowed in peaks, ",d18O.sdev.thres,"\n"))

  # maximum sd d2H value allowed?
  cat(paste("# maximum d2H standard deviation allowed in peaks, ",d2H.sdev.thres,"\n"))

  # maximum number of measurements allowed? (at 1.16 Hz, 4200 indices = 1 hour)
  # this filter is necessary if an analyzer runs out of dry gas source, for example.
  cat(paste("# maximum number of measurements allowed in peaks, ",max.length.thres,"\n"))

  # mainimum number of measurements allowed? (at 1.16 Hz, 35 indices = 30 seconds)
  cat(paste("# minimum number of measurements allowed in peaks, ",min.length.thres,"\n"))

  # close L2 header
  cat("# END L2 HEADER \n")

  # append L0/L1 headers.
  # loop through existing headers and paste below L2 header.
  for (i in 1:length(metadata_dataframe)) {
    cat(paste(metadata_dataframe[i],"\n")) 
  }
  # turn off output
  sink()
  
  # print statment denoting the end of this function
  if (dbg.level>0) {
    print("ending attach.L2.Header function")
    print("================================")
  }
}