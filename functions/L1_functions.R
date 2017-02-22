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

extract.date.from.L0.files <- function(file.list) {
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
  # return date list
  return(file.name.dates.woNAs)
}

#------------------------------------------------------------------------------------
# concatenate.to.monthly - load all files that are within a single month, and
# concatenate to a single monthly data frame. requires as an input a list of files to 
# concatenate, and returns the monthly data frame.

concatenate.to.monthly.uncompiled <- function(month.list,useParallel=TRUE) {
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

qflag.uncompiled <- function(x){
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
  # return list
  return(output)
}

qflag <- cmpfun(qflag.uncompiled)
#-----------------------------------------------------------------------
# data.sanity.check function
# common sense data flag - pull out values where the output makes no sense
# (e.g., impossible delta or H2o values). This function is explicitly conservative
# to only remove values that are extremely improbably to be correct...

data.sanity.check.uncompiled <- function(x) {
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
  # return output list
  return(output)
}

data.sanity.check <- cmpfun(data.sanity.check.uncompiled)

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
