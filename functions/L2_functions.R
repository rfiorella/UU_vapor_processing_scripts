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
  # return date list
  return(file.name.dates.woNAs)

  # print statment denoting the end of this function
  if (dbg.level>0) {
    print("ending extract.date.from.L1.files function")
    print("==========================================")
  }

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
  tmp <- diff(calibration.data.frame$EPOCH_TIME)
  # note that diff returns the difference between point i and i+1 at index i
  # ...when pulling out the indices that indicate the start of the subsequent 
  # analysis, these are at index i+1, which has been returned at index i - so
  # need to add 1 to the indices returned by diff 
  out1 <- which(tmp>thres)+1

  # add points to the beginning and end of the array to fully bound first/last plateaus
  out1 <- c(1,out1,nrow(calibration.data.frame)) 

  # in cases of power failure or some other hiccup, sometimes two consecutive points
  # can be identified as break points...this causes havoc in the spline fitting
  # process, and therefore, points that are too close to each other must be filtered out.
  # proposing a threshold here of 30 indices.
  tmp <- diff(out1)

  ini.breaks <- out1[c(1e6,tmp)>10]

  # print previous two vectors if requested by sufficiently high debug level
  if (dbg.level>1) {
    print("Break indices identified solely by time threshold:")
    print(out1)
    print("Time differences between these indices:")
    print(tmp)
  }

  # loop through initial break points, calculate splines, ID points where change is huge.
  fix.break.splines <- vector("list",length(ini.breaks))
  
  # calculate splines

  # NOTE: in this portion of the script, df and the threshold derivative value are TUNED parameters
  # they may not work in all situations...please check your data and ensure this works for your 
  # data...
  for (i in 1:(length(ini.breaks)-1)) {
    tmp.splines <- smooth.spline(calibration.data.frame$EPOCH_TIME[ini.breaks[i]:(ini.breaks[i+1]-1)],
      calibration.data.frame$H2O[ini.breaks[i]:(ini.breaks[i+1]-1)],df=96)
    tmp.spline.deriv <- c(NA,diff(tmp.splines$y)) # initial NA required to keep vector same length.
    # print(paste(i,max(tmp.spline.deriv)))
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

  # print statment denoting the end of this function
  if (dbg.level>0) {
    print("ending ID.calib.breakpoints function")
    print("++++++++++++++++++++++++++++++++++++")
  }

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
    #print(i)
    H2O.spline.fits <- smooth.spline(calibration.data.frame$EPOCH_TIME[(breakpoints[i]):(breakpoints[i+1]-1)],
      calibration.data.frame$H2O[(breakpoints[i]):(breakpoints[i+1]-1)],df=dfree)
    #print(i)
    d18O.spline.fits <- smooth.spline(calibration.data.frame$EPOCH_TIME[(breakpoints[i]):(breakpoints[i+1]-1)],
      calibration.data.frame$Delta_18_16[(breakpoints[i]):(breakpoints[i+1]-1)],df=dfree)
    #print(i)
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

  # return spline values.
  return(output.spline.fits)

  # print statment denoting the end of this function
  if (dbg.level>0) {
    print("ending fit.calibration.splines function")
    print("***************************************")
  }

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
    #temp5 <- temp4 # momentary kludge fix for debugging.
    temp6 <- c(NA,diff(temp1))
    temp7 <- c(NA,diff(temp2))
    temp8 <- c(NA,diff(temp3))

    #if (dbg.level>1) {
      print(i)  
      print(paste(length(temp1),length(temp2),length(temp3),length(temp4),
        length(temp5),length(temp6),length(temp7),length(temp8)))
    #}

    spline.derivatives[[i]] <- data.frame("d_H2O"=temp1,"d_d18O"=temp2,"d_d2H"=temp3,"time"=temp4,"inds"=temp5,
      "d2_H2O"=temp6,"d2_d18O"=temp7,"d2_d2H"=temp8)
  }

  # return output variable
  return(spline.derivatives)

  # print statment denoting the end of this function
  if (dbg.level>0) {
    print("ending calculate.spline.derivatives function")
    print("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&")
  }

}

#------------------------------------------------------------------------------------------
# extract.stable.calib.indices -
# using the splines and their derivatives, finds the set of indices
# that should be used to calculate mean/standard deviation values for the analysis. This function
# is somewhat large, but a lot is done under the hood and doesn't need to be stored in memory 
# following calculation. 

extract.stable.calib.indices <- function(spline.derivatives,H2O.thres=5.0,d18O.thres=0.01,d2H.thres=0.1,dbg.level=0) {
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
  inds.to.keep[[i]] <- spline.derivatives[[i]]$inds[indvec]
  }
  return(inds.to.keep)

  # print statment denoting the end of this function
  if (dbg.level>0) {
    print("ending extract.stable.calib.indices function")
    print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
  }

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

        #print(paste(abs(Omem.slope),abs(Hmem.slope),z))
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
      } else {
        #print("No inds for this period...")
      } # end check for good inds.
    }
  }
  
  std.avgs <- data.frame("H2O.mean"=temp1,"d18O.mean"=temp2,"d2H.mean"=temp3,
    "H2O.sd"=temp4,"d18O.sd"=temp5,"d2H.sd"=temp6,
    "time.mean"=temp7,"time.min"=temp8,"time.max"=temp9,
    "H2O.range"=temp10,"d18O.range"=temp11,"d2H.range"=temp12,
    "ind.count"=temp13,"d18O.trend"=temp14,"d18O.r2trend"=temp15,
    "d2H.trend"=temp16,"d2H.r2trend"=temp17)

  # filter out some obviously incorrect points.
  h2o.min <- which(std.avgs$H2O.mean > 2000) # must be at least 2000 ppm
  h2o.max <- which(std.avgs$H2O.mean < 40000) # must be no greater than 30000 ppm
  h2o.lsd <- which(std.avgs$H2O.sd < 1000)  # enforce a measure of stability
  d18O.lsd <- which(std.avgs$d18O.sd < 1) # enforce a measure of isotopic stability
  d2H.lsd <- which(std.avgs$d2H.sd < 8) # enforce a measure of isotopic stability
  # enforce a max length on the calibration data set - calibration periods identified
  # are unlikely to be longer than a half hour. at 1.16 Hz: 60 sec*1.16Hz*30 min = 2088 inds.
  max.length <- which(std.avgs$ind.count <= 4200)
  # also, require at least 2 minutes of data for data points 1.16 Hz*60sec*3 min = 140 inds
  min.length <- which(std.avgs$ind.count >= 35)

  retain.inds <- Reduce(intersect,list(h2o.min,h2o.max,h2o.lsd,d18O.lsd,d2H.lsd,max.length,min.length))
  
  #retain.inds <- seq(1,nrow(std.avgs),1) # for debugging - keeps all data.

  #print(retain.inds)
  #if (memory.filter) {print(mem.count)}
  # print some diagnostics. 
  print(paste(nrow(std.avgs)-length(retain.inds),
     " of ",nrow(std.avgs)," (",round(100*(nrow(std.avgs)-length(retain.inds))/nrow(std.avgs),2),
     "%) calibration periods failed filter...removing..."))
  
  if (length(retain.inds) > 0) {
    std.avgs <- std.avgs[retain.inds,] 
  } else {
    std.avgs <- list() # try returning empty list to avoid this from not existing...
  }
  
  return(std.avgs)

  # print statment denoting the end of this function
  if (dbg.level>0) {
    print("ending calculate.standard.averages function")
    print("===========================================")
  }

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

  # print(head(time.inds))
  # print(nrow(time.inds))

  # loop through rows of time.inds and find the index of
  # nearest point before time.min and the nearest point after time.max
  before.ind <- vector("numeric",nrow(time.inds))
  after.ind <- vector("numeric",nrow(time.inds))

  for (i in 1:nrow(time.inds)) {
      # should be positive for points before calibration point
      tmp1 <- time.inds$time.min[i] - amb.time
      #print(tmp1)  
      # find minimum positive point, assuming any points are positive
      before.ind[i] <- ifelse(any(tmp1>0), # conditional test
        which(tmp1==min(tmp1[tmp1>0],na.rm=TRUE)), # true branch
        NA) # false branch
      # reverse the logic for the point after calibration
      tmp2 <- amb.time - time.inds$time.max[i]
      #print(tmp2)
      # first point after calibration using eqn above should be 
      # minimum positive number  
      after.ind[i] <- ifelse(any(tmp2>0), # conditional test
        which(tmp2==min(tmp2[tmp2>0],na.rm=TRUE)), # true branch
        NA) # false branch

      # clean up loop
      rm(tmp1)
      rm(tmp2)
  }

  # print(before.ind)
  # print(after.ind)

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

  # print some diagnostics
  # print(before.H2O)
  # print(after.H2O)
  # print(before.18O)
  # print(after.18O)

  # make the data output data frame
  output <- data.frame("before.H2O"=before.H2O,"before.d18O"=before.18O,
    "before.d2H"=before.2H,"after.H2O"=after.H2O,
    "after.d18O"=after.18O,"after.d2H"=after.2H)

  # return output frame
  return(output)

  # print statment denoting the end of this function
  if (dbg.level>0) {
    print("ending get.ambient.deltas function")
    print("++++++++++++++++++++++++++++++++++")
  }

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

  # return the whole data frame with the two mixing ratio corrected variables attached.
  return(avg.data.frame)

  # print statment denoting the end of this function
  if (dbg.level>0) {
    print("ending apply.mixingratio.correction.calibration function")
    print("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&")
  }

}

#-----------------------------------------------------------------------------------
# apply.drygas.correction function - corrected now that ambient data is already in data frame.

apply.drygas.correction <- function(data,do.correction=TRUE,H2O.bg,include.gypsum.fractionation,dbg.level=0) {
  # print statement indicating that this function is starting
  if (dbg.level>0) {
    print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
    print("starting apply.drygas.correction function")
  }

  #======================================================
  # start work of function

  print(paste(now()," Applying dry gas correction to calibration data..."))
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
  
  # return data frame
  return(data)

  # print statment denoting the end of this function
  if (dbg.level>0) {
    print("ending apply.drygas.correction function")
    print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
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

  if (method==1) {
    # ensure that we're working with a data frame...
    stopifnot(is.data.frame(standard.data.frame))

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
    }
   
    # create dataframe packaging this info out to return
    output <- data.frame("start.time"=period.begin,"stop.time"=period.end,
      "O.slope"=period.Oslope,"O.intercept"=period.Ointercept,"O.r2"=period.Orsq,
      "H.slope"=period.Hslope,"H.intercept"=period.Hintercept,"H.r2"=period.Hrsq)

    # we should cut out any obviously bad calibrations - with 2 points, this can be any period where
    # r2 is less than 0.95.
    bad.Hinds <- which(output$H.r2 < 0.9)
    bad.Oinds <- which(output$O.r2 < 0.9)

    all.bad.inds <- intersect(bad.Hinds,bad.Oinds)
    print(paste(bad.Hinds,bad.Oinds,all.bad.inds))

    # fix time array to avoid losing coverage (just assume previous calibration period extends to current one for now...)
    output$stop.time[(all.bad.inds-1)] <- output$stop.time[all.bad.inds]

    # now remove the bad inds...
    if (length(all.bad.inds)>0) {
      output <- output[-all.bad.inds,]
    }
    # return the variables.
    return(output)

  }

  # print statment denoting the end of this function
  if (dbg.level>0) {
    print("ending correct.standards.to.VSMOW function")
    print("==========================================")
  }
}

#----------------------------------------------------
# create attach.L2.Header function
# this function creates the header for datafiles and attaches the metadata
# loaded into a csv file.

attach.L2.Header <- function(output_filename,metadata_dataframe,dbg.level=0) {
  # print statement indicating that this function is starting
  if (dbg.level>0) {
    print("==================================")
    print("starting attach.L1.Header function")
  }

  # initiate writing out to output file...
  sink(output_filename)

  # add L1 information to top of output file...
  cat("# BEGIN L2 HEADER \n")

  # note when L1 script was run
  cat(paste("# date.L2.scripts.run, ",base::date(),"\n"))
  cat(paste("# time.threshold.for.separate.standard.runs, ",time.threshold,"\n"))
  cat(paste("# stiff.spline.dfs, ",stiff.spline.dfs,"\n"))
  cat(paste("# spline.derivative.H2O.threshold (ppm/sec), ",H2O.thres,"\n"))
  cat(paste("# spline.derivative.d18O.threshold (permil/sec), ",d18O.thres,"\n"))
  cat(paste("# spline.derivative.d2H.threshold (permil/sec), ",d2H.thres,"\n"))
  cat(paste("# memory.filter.on to calculate averages, ",memory.filter,"\n"))
  cat(paste("# dry.gas.correction.applied, ",do.correction,"\n"))
  if (do.correction==TRUE) {
    cat(paste("# background H2O assumed, ",H2O.bg,"\n"))
    cat(paste("# gypsum fractionation included in drierite correction, ",include.gypsum.fractionation,"\n"))
  }
  cat(paste("# humidity correction fit type, ",fit.type,"\n"))
  cat(paste("# humidity correction oxygen slope, ",Oslope, "\n"))
  cat(paste("# humidity correction hydrogen slope, ",Hslope,"\n"))

  # future dev: include information about standard? probably not necessary since this is included
  # in output variables.

  # close L1 header
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

