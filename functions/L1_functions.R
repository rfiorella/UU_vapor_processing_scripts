# Picarro Water Vapor Isotope Data Reduction 
# L1 Functions

# migrated from WBB-vapor-scripts repo by rich fiorella on 21/2/17.
# code in these functions written by rich fiorella based on work by galen gorski
# (likely some of galen's functions have made it through with no changes!)

# contact info: rich fiorella, rich.fiorella@utah.edu

# this script contains a series of accessory functions for the L1 script.
# as of 21/2/17, this script has been designed to work only with a Picarro L2130.
# other analyzer versions will likely require modification of this script.

# load any required packages for these functions
library(parallel) # allow for parallel computing
library(compiler) # allow functions to be precompiled.
library(lubridate) # include useful time functions
##################################################################################
# determine the number of cores available for parallel computing

no_cores <- detectCores() - 1 # but leave on thread available so as not to crash computer
							  # this is important if runnning on a local computer, but -1 
							  # can be removed if running on a dedicated cluster
print(paste(Sys.time(),"Number of available cores for parallel processing:",no_cores))
#----------------------------------------------------------------------------------
# L1 functions follow....
#------------------------------------------------------------------------------------
# extract.date.from.L0.files - break down file names to provide just the date
# back as the return value. input is a list of file names - can have a full directory
# path, output is a list of just dates. from this, the indices matching each month can 
# easily be extracted.

extract.date.from.L0.files <- function(file.list,dbg.level=0) {
  # print statement indicating that this function is starting
  if (dbg.level>0) {
    print("============================================")
    print("starting extract.date.from.L0.files function")
  }

  # suppress warnings in this function - as.numeric will create NAs by coercion,
  # and this will be reported as a warning to the screen every time. But since this
  # is precisely what we are hoping to do, warnings are unncessary here.
  oldw <- getOption("warn")
  options(warn = -1)
  # split file name into pieces using *underscores*
  file.name.string.pieces <- strsplit(file.list,split="_")
  # extract dates from resulting string pieces
  file.name.dates <- ymd(strptime(unlist(file.name.string.pieces),format="%Y-%m-%d"))
  # remove NAs
  file.name.dates.woNAs <- file.name.dates[!is.na(file.name.dates)]
  # set warnings to old status
  options(warn = oldw)
  # print statment denoting the end of this function
  if (dbg.level>0) {
    print("ending extract.date.from.L0.files function")
    print("==========================================")
  }
  # return date list
  return(file.name.dates.woNAs)
}

#------------------------------------------------------------------------------------
# concatenate.to.monthly - load all files that are within a single month, and
# concatenate to a single monthly data frame. requires as an input a list of files to 
# concatenate, and returns the monthly data frame.

concatenate.to.monthly.uncompiled <- function(month.list,useParallel=TRUE,dbg.level=0) {
  # print statement indicating that this function is starting
  if (dbg.level>0) {
    print("++++++++++++++++++++++++++++++++++++++++")
    print("starting concatenate.to.monthly function")
  }

  print(paste(Sys.time()," Concatenating daily L0 files to a single monthly file..."))
  # require that indices be of type list...
  stopifnot(is.list(month.list))

  # make a vector of files
  files <- unlist(month.list)

  if (useParallel==TRUE) { # run function in parallel
  	cl <- makeCluster(no_cores)
  	d <- parLapply(cl,files,function(x){read.table(x,header=TRUE,stringsAsFactors=FALSE,sep=",")})
  	stopCluster(cl)

	output <- do.call(rbind,d)

	# clean up memory 
	rm(d)
	gc()
  } else { # run function in serial...
 	# concatenate all files together
  	output <- read.table(files[1],header=TRUE,stringsAsFactors=FALSE,sep=",")
  	# if month has more than one day of data, 
  	if (length(files)>1){
    	for (i in (2:length(files))) {
    	temp <- read.table(files[i],header=TRUE,stringsAsFactors=FALSE,sep=",")
    	output <- rbind(output,temp)
    	} # end for
    } # end length(files) if statement
  } # end (isParallel) if statement

  # print statment denoting the end of this function
  if (dbg.level>0) {
    print("ending concatenate.to.monthly function")
    print("++++++++++++++++++++++++++++++++++++++")
  }

  # return concatenated data
  return(output)
}

concatenate.to.monthly <- cmpfun(concatenate.to.monthly.uncompiled)

#-----------------------------------------------------------------------------------
# qflag function - read in data table, run a quick check. this varies from 
# older versions of qflag in the tolerance - reduced digits parameter from 
# 2 to 0 - changes tolerance from ±0.005°C or ±0.005 torr to ±0.5°C and ±0.5 torr
# in quick tests, this reduced the amount of data rejected from ~70% to <1%.
# old tolerances seem to be unnecessarily strict.

# input is the monthly data frame, output is the monthly data frame with flagged
# points removed.

qflag.uncompiled <- function(x,dbg.level=0){
  # print statement indicating that this function is starting
  if (dbg.level>0) {
    print("&&&&&&&&&&&&&&&&&&&&&&&")
    print("starting qflag function")
  }

  print(paste(now()," Removing points that fail instrument quality checks..."))
  # save number of rows initially to calculate how many data points are removed.
  begin <- nrow(x)
  # use the round function to identify points that will round to 80°C or 50 torr
  x = x[(round(x$CavityTemp, digits = 0) == 80),]
  x = x[(round(x$CavityPressure, digits = 0) == 50),]
  # find how many rows remain that pass the rounding test
  end <- nrow(x)
  # calculate the percentage of points that bave been removed, and print to the screen
  print(paste(100*(1-end/begin),"% of data points removed"))
  # return the filtered data frame and the percentage of points removed
  output <- list("data"=x,"pct.discarded"=round(100*(1-end/begin),2))
  # print statment denoting the end of this function
  if (dbg.level>0) {
    print("ending qflag function")
    print("&&&&&&&&&&&&&&&&&&&&&")
  }
  # return list
  return(output)
}

qflag <- cmpfun(qflag.uncompiled)
#-----------------------------------------------------------------------
# data.sanity.check function
# common sense data flag - pull out values where the output makes no sense
# (e.g., impossible delta or H2o values). This function is explicitly conservative
# to only remove values that are extremely improbably to be correct...

data.sanity.check.uncompiled <- function(x,dbg.level=0) {
  # print statement indicating that this function is starting
  if (dbg.level>0) {
    print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
    print("starting data.sanity.check function")
  }

  print(paste(now()," Removing points that fail data sanity checks..."))
  # save number of initial rows for future calculations
  begin <- nrow(x)
  # apply data sanity standards - at a minimum, values should be:
  # (a) for H2O, values should be positive.
  # (b) for d18O, values MUST be above -1000, and ALMOST CERTAINLY will be below 50
  # (c) for d2H, values MUST be above -1000, and ALMOST CERTAINLY will be below 400
  x = x[(x$H2O>=0),] 
  x = x[(x$Delta_18_16 > -1000),]
  x = x[(x$Delta_D_H > -1000),]
  x = x[(x$Delta_18_16 < 50),]
  x = x[(x$Delta_D_H < 400),]
  # how many points pass these tests?
  end <- nrow(x)
  # calculate the percentage of points that have been removed, and print to stdout
  print(paste(100*(1-end/begin),"% of data points removed"))
  # return the filtered data frame and the percentage of points removed
  output <- list("data"=x,"pct.discarded"=round(100*(1-end/begin),2))
  # print statment denoting the end of this function
  if (dbg.level>0) {
    print("ending data.sanity.check function")
    print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
  }
  # return output list
  return(output)
}

data.sanity.check <- cmpfun(data.sanity.check.uncompiled)

#----------------------------------------------------
# create attach.L0.Header function
# this function creates the header for datafiles and attaches the metadata
# loaded into a csv file.

attach.L1.Header.uncompiled <- function(output_filename,metadata_dataframe,dbg.level=0) {
  # print statement indicating that this function is starting
  if (dbg.level>0) {
    print("==================================")
    print("starting attach.L1.Header function")
  }

  # initiate writing out to output file...
  sink(output_filename)

  # add L1 information to top of output file...
  cat("# BEGIN L1 HEADER \n")

  # note when L1 script was run
  cat(paste("# date.L1.scripts.run, ",base::date(),"\n"))
  
  # close L1 header
  cat("# END L1 HEADER \n")

  # append L0 header.
  # loop through existing L0 header and paste below L1 header.
  for (i in 1:length(metadata_dataframe)) {
    cat(paste(metadata_dataframe[i],"\n")) 
  }

  # turn off output
  sink()
  
  # print statment denoting the end of this function
  if (dbg.level>0) {
    print("ending attach.L1.Header function")
    print("================================")
  }
}

# compile the function
attach.L1.Header <- cmpfun(attach.L1.Header.uncompiled)

#----------------------------------------------------
# create post.calibration.filter function
# often after a calibratino period, there are apparent discontinuities in the 
# data that make it appear is if there are very rapid h2o/isotope changes
# in the ambient data that are wholly artifacts of the SDM/Picarro valve sequences
# this function is designed to filter out these discontinuities.

post.calibration.filter <- function(qcchecked.data.frame,dbg.level=0) {
  # print statement indicating that this function is starting
  if (dbg.level>0) {
    print("+++++++++++++++++++++++++++++++++++++++++")
    print("starting post.calibration.filter function")
  }

  # first, find where the calibration periods are...
  ind.after.calib <- which(diff(qcchecked.data.frame$EPOCH_TIME)>120) + 1

  # identify point 15 minutes later...
  end.calib.memory.period <- ind.after.calib + 30 # note: THIS ASSUMES USING 1 MINUTE AVERAGES!!!!

  # loop through create a list of indices corresponding to this and remove...
  inds.to.remove.list <- vector("list",length(ind.after.calib))

  for (i in 1:length(ind.after.calib)) { 
    inds.to.remove.list[[i]] <- seq(ind.after.calib[[i]],end.calib.memory.period[[i]],1)
  }

  inds.to.remove <- unlist(inds.to.remove.list)

  print(length(inds.to.remove)/nrow(qcchecked.data.frame))
  # omit these rows and return the data frame...
  output <- qcchecked.data.frame[-inds.to.remove,]

  # print statment denoting the end of this function
  if (dbg.level>0) {
    print("ending post.calibration.filter function")
    print("+++++++++++++++++++++++++++++++++++++++")
  }

  # return the output
  return(output)
}