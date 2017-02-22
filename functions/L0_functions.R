# Picarro Water Vapor Isotope Data Reduction 
# L0 Functions

# migrated from WBB-vapor-scripts repo by rich fiorella on 21/2/17.
# code in these functions written by rich fiorella based on work by galen gorski
# (likely some of galen's functions have made it through with no changes!)

# contact info: rich fiorella, rich.fiorella@utah.edu

# this script contains a series of accessory functions for the L0 script.
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

##################################################################################
# Here are the functions used in the L0 script!
##################################################################################
# create function to extract dates from log file output
# NOTE: this so far likely only works for Picarro instruments and will
# have to be modified for different instrument systems!

extract.date.from.files.uncompiled <- function(file.list) {
	# suppress warnings in this function - as.numeric will create NAs by coercion,
	# and this will be reported as a warning to the screen every time. But since this
	# is precisely what we are hoping to do, warnings are unncessary here.
	oldw <- getOption("warn")
	options(warn = -1)
	# split the list of files into pieces using folder structure
	file.name.string.pieces <- strsplit(file.list,split="/") 
	# force string pieces to be numeric - should only leave the yyy, mm, and dd pieces
	file.name.string.vector <- as.numeric(unlist(file.name.string.pieces)) 
	# remove NAs from prior vector, leaving only the components associated with the date
	date.vector <- file.name.string.vector[!is.na(file.name.string.vector)]
	# restructure to a matrix with dimensions daysX3 (for yyyy,mm,dd) and fill by row, not by column
	date.matrix <- matrix(date.vector,ncol=3,byrow=TRUE)
	# finally, concatenate these as strings again and format as a date object so R can use later
	# to specify whether certain days are in or out.
	date.list <- as.Date(paste(date.matrix[,1],"-",date.matrix[,2],"-",date.matrix[,3],sep=""))
	# return warnings to old status
	options(warn = oldw)
	# return date list as variable of interest
	return(date.list)
}

# compile function
extract.date.from.files <- cmpfun(extract.date.from.files.uncompiled)
#------------------------------------------------------------------------
# pad.file.indices function
# adds "padding" to the identified file indices to ensure that each day 
# has all of the data corresponding to that particular day, which may be 
# in "neighboring day folders"

pad.file.indices.uncompiled <- function(indices,file.list) {
	minval <- min(indices)-2
  	maxval <- max(indices)+2
  	# set limits on minval and max val to the boundaries of file.list
  	# to ensure that they don't go below limits of file.list
  	if (minval<1) {
    	minval <- 1
  	}
  	if (maxval>length(file.list)) {
    	maxval <- length(file.list)
  	}
  	# reset array to proper values
  	inds <- seq(minval,maxval,1)
  	return(inds)
}

# compile function
pad.file.indices <- cmpfun(pad.file.indices.uncompiled)

#----------------------------------------------------
# concatenate.to.daily function
# combines all of the data from the files corresponding to 
# a day, cuts to the limits of the day, and returns the proper array

# NOTE TO USERS: This code is definitely faster than the non-parallel 
# version, but it is terribly inefficient - speeds up overall script by
# ~40%, but triples the computing resources required.

concatenate.to.daily.uncompiled <- function(indices,date,file.list,useParallel=TRUE) {
	
	if (useParallel==TRUE) { # if running function in parallel...
		cl <- makeCluster(no_cores)
		d <- parLapply(cl,file.list[indices],function(x){read.table(x,header=TRUE,stringsAsFactors=FALSE,fill=TRUE)})
		stopCluster(cl)

		# concatenate into same data frame
		output <- do.call(rbind,d)

		# clean up memory 
		rm(d)
	} else { # if require serial...
		output <- read.table(file.list[indices[1]],header=TRUE,stringsAsFactors=FALSE,fill=TRUE)
		for (i in (2:length(indices))) {
			temp <- read.table(file.list[indices[i]],header=TRUE,stringsAsFactors=FALSE,fill=TRUE)
			# add to data frame
			output <- rbind(output,temp)
		} # end for loop
	} # end parallel if statement

	# set time limits using seconds
	min.seconds <- as.numeric(as.POSIXct(date,tz="gmt"))
	max.seconds <- as.numeric(as.POSIXct(date+1,tz="gmt"))

	# trim output to time limits
	output.trimmed <- output[output$EPOCH_TIME >= min.seconds & output$EPOCH_TIME < max.seconds,]

	# return the data frame.
	return(output.trimmed)
}

# compile the function
concatenate.to.daily <- cmpfun(concatenate.to.daily.uncompiled)


#----------------------------------------------------
# restructure.time.variables function
# takes DATE and TIME columns and restructures to yyyy, month, dd, hh, minute, ss variables

restructure.time.variables.uncompiled <- function(data.frame) {
	# requires that dataframe has 2 variables called DATE and TIME
	# probably only works for Picarro analyzers!
	
	# restructure date variable
	date.fix <- strsplit(data.frame$DATE,split="-")
	date.vector <- as.numeric(unlist(date.fix))
	date.matrix <- as.data.frame(matrix(date.vector,ncol=3,byrow=TRUE))
	names(date.matrix) <- c("YYYY","MM","DD")

	# restructure time variable
	time.fix <- strsplit(data.frame$TIME,split=":")
	time.vector <- as.numeric(unlist(time.fix))
	time.matrix <- as.data.frame(matrix(time.vector,ncol=3,byrow=TRUE))
	names(time.matrix) <- c("hh","mm","ss")

	# bind new time variable structures together
	newtime <- cbind(date.matrix,time.matrix)

	# remove old variables from input data frame.
	data.frame$DATE <- data.frame$TIME <- NULL

	# combine new time variables with input data frame
	new.data.frame <- cbind(newtime,data.frame)

	# return new data frame
	return(new.data.frame)
}

# compile the function
restructure.time.variables <- cmpfun(restructure.time.variables.uncompiled)

#----------------------------------------------------
# create attach.L0.Header function
# this function creates the header for datafiles and attaches the metadata
# loaded into a csv file.

attach.L0.Header.uncompiled <- function(output_filename,metadata_dataframe) {
  sink(output_filename)
  cat(paste("# ",nrow(metadata_dataframe)+4,"\n",sep=""))
  cat("# BEGIN HEADER \n")
  for (i in 1:nrow(metadata_dataframe)) {
    	cat(paste("# ",metadata_dataframe[i,1],", ",metadata_dataframe[i,2],"\n",sep=""))   	
  }
  cat(paste("# date.scripts.run, ",base::date(),"\n"))
  # date.scripts.run flag modified on 2feb17 from date() to base:date() to avoid conflicts
  # if the lubridate package is active.
  cat("# END HEADER \n")
  sink()
}

# compile the function
attach.L0.Header <- cmpfun(attach.L0.Header.uncompiled)

#-----------------------------------------------------------------------------------
# extract.calib.periods function - identifies the potential periods that could be calibration
# data. As of 24aug16, this function will work *ONLY* with Picarro analyzers - should
# work at least for the l2120 and 2130 analyzers.

# valve code 6 indicates valves 2 and 3 in the vaporizer are open and is delivering
# standard to the analyzer cavity.

extract.calib.periods <- function(data) {
  # new strategy - identify all points where vm is 6, then find point just before or just after 
  CalibPeriod <- data[data$ValveMask == 6 & is.na(data$ValveMask) == FALSE,]
  # print maximum number of calibration observations
  print(paste(Sys.time()," Maximum number of calibration points: ",nrow(CalibPeriod)))
  # return calibration only data frame
  return(CalibPeriod)
}

#-----------------------------------------------------------------------------------
# extract.ambient.periods function - identifies the potential periods that could be ambient
# data. As of 24aug16, this function will work *ONLY* with Picarro analyzers - should
# work at least for the l2120 and 2130 analyzers.

# valve code 0 indicates all valves in the vaporizer are closed, defaulting to ambient intake.

extract.ambient.periods <- function(data) {
  # identify the data points where: valvemask = 0, and valvemask exists.  
  AmbPeriod <- data[data$ValveMask == 0 & is.na(data$ValveMask) == FALSE,]
  # print number of ambient observations  
  print(paste(Sys.time()," Number of ambient points: ",nrow(AmbPeriod)))
  # return ambient only data frame
  return(AmbPeriod)
}

#-----------------------------------------------------------------------------------
# reduce.ambient.data function - this function reduces the 1.16 Hz analyzer data to 
# minute or longer averages. requires at least 20 data points to be in a minute to
# to be recorded.

# 2 accessory functions embedded in reduce.ambient.data and printed below. 
# get.time.indices determines the indices within each period that will be averaged
# get.time.mean takes averages over the period requested. 

#----------------------------------------
# start reduce.ambient.data function

reduce.ambient.data.uncompiled <- function(ambient.data.frame,time.length.average=1,
	minimum.points.to.average=20,useParallel=TRUE) {
  print(paste(Sys.time()," Reducing ambient data to ",time.length.average," minute averages"))
  # find month we're working in and how many days are in it.
  mon <- month(ambient.data.frame$MM[1]) # probably can do better...
  year <- ambient.data.frame$YYYY[1]
  dinm <- days_in_month(mon)

  # determine number of potential data points in output data.frame
  nperiods <- dinm*(24)*(60/time.length.average)

  # print(nperiods)
  # create output array
  startdate <- ymd(paste(year,mon,"01"))
  if (!is.na(nperiods) & nperiods > 0) {
  	  time.averages.array <- as.numeric(as.POSIXct(startdate + c(0:nperiods) * minutes(time.length.average),origin="1970-01-01"))
	  
	  # data.time.array
	  data.time.array <- as.numeric(as.POSIXct(ymd_hms(paste(ambient.data.frame$YYYY,ambient.data.frame$MM,
	    ambient.data.frame$DD,ambient.data.frame$hh,ambient.data.frame$mm,ambient.data.frame$ss)),origin="1970-01-01"))

	  # define accessory functions 

	  # this function - get.time.indices - is obnoxiously slow.
	  get.time.indices <- function(i) {
		# attempts to speed up this function - it's the worst!
		# tsubset <- data.time.array[seq(1,length(data.time.array),500)]
		# print(tsubset)
		# closest.ind <- which(tsubset-time.averages.array[i] == min(tsubset-time.averages.array[i]))
		# print(closest.ind)
		# # narrows the next step down from ~2 million lines to ~20,000
		# ind.subset <- max(c(1,500*(closest.ind-20))):min(c(length(data.time.array),500*(closest.ind-20)))

	  	# print(ind.subset)
	  	# 	output <- which(data.time.array[ind.subset] >= time.averages.array[i-1] & 
	  	# 		data.time.array[ind.subset] < time.averages.array[i])
	  	# 	return(output)
	   
		## old, slow version
	  	output <- which(data.time.array >= time.averages.array[i-1] & 
	  		data.time.array < time.averages.array[i])

	  	# return list of indices that correspond to desired time interval
	  	return(output)
	  }

	  get.time.mean <- function(i) {
	  	if (length(time.inds[[i]]) >= minimum.points.to.average) {
	  		output <- colMeans(ambient.data.frame[time.inds[[i]],],na.rm=TRUE)
	  	}
	  }
	  
	  if (useParallel==TRUE) {
		# loop through times, find indices within each time window, and average variables over them
		cl <- makeCluster(no_cores,"FORK")
		print(paste(Sys.time()," finding indices associated with each averaging period..."))
		time.inds <- parLapply(cl,2:(length(time.averages.array)+1),get.time.indices)

		# loop through index lists and make column averages
		print(paste(Sys.time()," averaging columns in data frame..."))
		tmp.output <- parLapply(cl,1:(length(time.averages.array)),get.time.mean)

		stopCluster(cl)

	  } else {
	  	stop("No serial version of reduce.ambient.data has been written...")
	  }

	  # bind rows together
	  reduced.ambient.data <- do.call(rbind,tmp.output)

	  # return the vector
	  return(reduced.ambient.data)

  } else {
  	print("this date has no ambient data - appears to be all calibration data or no data")
  	return(NULL)
  } # end if statement re: nperiods
}

# compile function
reduce.ambient.data <- cmpfun(reduce.ambient.data.uncompiled)
