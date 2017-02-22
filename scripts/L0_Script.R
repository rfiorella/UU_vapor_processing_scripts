#Picarro Water Vapor Isotope Data Reduction 
# L0 Script, v1.0.2
#This is the first of 3 Script/Function pairs that take raw 1Hz Picarro water isotope data and work it up to five minute calibrated data
#
# v0 written 04/15 by GAG
# algorithms optimized in 05/16 by RPF. 
#
#General Explanation:
#The idea of this first level is to take the raw 1 hz data, stored in hourly files and place the rows in daily data frames. 
#The data frames have one row per second and each data frame is the exact same size 86400 rows and 26 columns, some entries 
#will contain only zeros as no measurement was taken that second, some rows will be overwritten as multiple measurements were 
#taken in that second. Some large blocks of time will be zeros as well, due gaps when the instrument was off or not taking data. 
#All the rows will be placed in the data frames based on the count of the seconds since midnight.

#All times are in GMT, although in the raw .dat files produced by the intstrument, the DATE, TIME and EPOCH TIME are GMT while the FRAC_DAYS_SINCE_JAN1, 
#FRAC_HRS_SINCE_JAN_1, and JULIAN_DAYS are all based on the internal computer time (set to PST). This code and the following layers of script only use
#DATE, TIME and EPOCH_TIME columns so all input and outputs are in GMT. More detailed comments can be found in the word document 
#"Comments and Explanation for L0 Script and L0 Functions", each [#] in the script refers to an explanation in that document.

# dev notes:
# 30jan17 - in an attempt to reduce the processing time of this script, added
#   defined sizes to variables that are called (e.g., vectors are "preallocated" now)
#   this change sped up processing time by ~14%.


rm(list=ls())

set.seed(1) # reproducability.

#Rprof("L0profile_30jan17")
######################################################################################################################################
# load associated function file(s) 

source("../functions/L0_Functions.R")
print(paste(Sys.time(),"Loading and compiling L0 functions..."))

######################################################################################################################################
# SET USER OPTIONS
######################################################################################################################################

startdate <- "2015-10-12" #expects a number in yyyymmdd format - will begin at MIDNIGHT GMT.
enddate   <- "2016-01-07" #expects a number in yyyymmdd format - will end at MIDNIGHT GMT.

path.to.data <- "~/WBB_VAPOR/Raw/"
path.to.output.L0.data <- "~/WBB_VAPOR/L0/testing/"
             
######################################################################################################################################
# SET METADATA
# write out metadata to help data curation - these will be appended to the top of the data files currently. it might be possible in
# later versions to use this to write out an xml file?
######################################################################################################################################

metadata.frame <- read.csv("../metadata_templates/L0_WBB_metadata.csv",header=TRUE)

######################################################################################################################################
# SET DATE LIMITS

DateLimits <- as.Date(c(startdate,enddate))
ndays <- as.numeric(DateLimits[2]-DateLimits[1],units="days") # requires end date to be after start date...

# sanity check - stop script if enddate is before startdate
if (ndays<0) {stop("Start date is after end date - fix input dates")}
#####################################################################################################################################
# create lists of files required to process the data

# make an empty list to hold all raw file names
raw.file.list <- list()

# list of all raw files
raw.file.list <- list.files(path = path.to.data, full.names=TRUE, recursive=TRUE)

# extract a vector of dates from the file names corresponding to each file.
dates <- extract.date.from.files(raw.file.list)

################################################################################################
# loop through days and create daily files.
#-----------------------------------------------------------
print(c("Maximum number of days to process:",ndays))

# pre-allocate vectors for logging dataframe.
log.date <- vector(mode="double",length=ndays)
log.nfil <- vector(mode="integer",length=ndays)
log.npts <- vector(mode="integer",length=ndays)
log.time <- vector(mode="double",length=ndays)

for (i in 1:ndays) {
  ptm <- proc.time()

  # looking for the index of the first file corresponding to the date of interest...
  date.to.process <- as.Date(startdate)+(i-1)

  # print out the day that is about to be processed:
  print(paste(Sys.time(),"  Processing: ",date.to.process, sep="")) 
 
  # find files that correspond to this day
  files.within.day <- which(dates == date.to.process)

  #--------------------------------------------------------------------
  # LOG FILE - INFORMATION ABOUT HOW MUCH DATA CORRESPONDS TO EACH DAY
  # log.date[i] <- as.Date(startdate)+(i-1)
  # log.nfil[i] <- length(files.within.day)
  #--------------------------------------------------------------------

  # some days will be missing, so need a function that will cause script to move to the 
  # next day if there are no files corresponding to that day
  if (length(files.within.day)==0) {
  	print(paste(Sys.time(),"  No data for : ",date.to.process, sep=""))
  	next
  } # skip to next index if the day is missing...

  # add padding to these numbers - look to the ~two files before and after to ensure that each
  # file has a full 24 hours of data where available. This is required since if a log file begins
  # after 23:00 the previous day, some data for day x will be held in the folder of day x-1.
  file.inds.to.include <- pad.file.indices(files.within.day,raw.file.list)

  # loop through log files for that particular day and concatenate
  daily.data.frame <- concatenate.to.daily(file.inds.to.include,date.to.process,raw.file.list,useParallel=FALSE)

  # fix time variables for output - 2 fixes need to be made:
  # (1) there are several redundant time variables. remove extra ones.
  vars.to.remove <- names(daily.data.frame) %in% c("FRAC_DAYS_SINCE_JAN1","FRAC_HRS_SINCE_JAN1","JULIAN_DAYS")
  daily.data.frame <- daily.data.frame[!vars.to.remove]

  # (2) restructure DATE and TIME columns to give yyyy, month, dd, hh, minute, ss vars
  daily.data.wtimefix <- restructure.time.variables(daily.data.frame)

  #--------------------------------------------------------------------
  # LOG FILE - INFORMATION ABOUT HOW MUCH DATA CORRESPONDS TO EACH DAY
  # log.npts[i] <- nrow(daily.data.wtimefix)
  #--------------------------------------------------------------------

  # divide into ambient and calibration data files
  calib.data <- extract.calib.periods(daily.data.wtimefix)
  amb.data <- extract.ambient.periods(daily.data.wtimefix)

  # clean up original data
  rm(daily.data.wtimefix)
  rm(daily.data.frame)
  gc()

  # reduce ambient data.
  amb.data.avgd <- reduce.ambient.data(amb.data)
  
  # WRITE OUT DAILY FILES FOR BOTH AMBIENT DATA, CALIBRATION DATA.
  #-----------------------------------
  # write out ambient data
  #---------------------------
  if (!is.null(amb.data.avgd)) {
    if (nrow(amb.data.avgd) > 0) {
      print(paste(Sys.time()," Writing out ambient data file..."))
      aoutput.fname <- paste(path.to.output.L0.data,"WBB_Water_Vapor_AmbientData_L0_",
        date.to.process,"_",metadata.frame$Value[metadata.frame$Variable=="code.version"],".dat",sep="")

      # attach metadata
      attach.L0.Header(aoutput.fname,metadata.frame)

      # write out data portion of the data file.
      write.table(amb.data.avgd,file=aoutput.fname,sep=",",append = TRUE,row.names=FALSE)
    } else {
      print(paste(Sys.time(), "No data passed filters for this day...."))
    }
  }

  # write out calibration data
  #---------------------------
  if (nrow(calib.data) > 0) {
    print(paste(Sys.time()," Writing out calibration data frame..."))
    # generate output filename
    coutput.fname <- paste(path.to.output.L0.data,"WBB_Water_Vapor_CalibData_L0_",
        date.to.process,"_",metadata.frame$Value[metadata.frame$Variable=="code.version"],".dat",sep="")

    # attach metadata
    attach.L0.Header(coutput.fname,metadata.frame)

    # write out data portion of the data file.
    write.table(calib.data,file=coutput.fname,sep=",",append = TRUE,row.names=FALSE)
  } else {
    print(paste(Sys.time(), "No data passed filters for this day...."))
  }
 }

# WRITE OUT LOG FILE
#log.dataframe <- data.frame(log.date,log.nfil,log.npts,log.time)
#log.fname <- paste("log_",Sys.time(),"_",startdate,"_to_",enddate,sep="")
#write.table(log.dataframe,file=log.fname,sep=",")

#Rprof(NULL)
