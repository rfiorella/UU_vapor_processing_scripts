# L1_Script.R
# L1 script, richf 
# heavily borrowing script and algorithms from L1 scripts by Galen Gorski.

# starting from a fresh code here since: (a) forces me to go through and understand all old code,
# (b) I think there are some major coding-philosophy differences to be sorted out
# (c) structure of L0 files changed in ways that might not be directly compatible.

# there are a few *GIGANTIC* functions in the old L1 functions file. This code will likely be more
# readable if broken into smaller functions.

# clear workspace
rm(list=ls())

# load libraries
library(lubridate)

# load associated functions
source("../functions/L1_functions.R")

# Set user variables
start.date <- ymd("2013-04-01")
end.date <- ymd("2017-02-01")

path.to.L0.data <- "~/WBB_VAPOR/L0/testing/"
path.to.output.L1.data <-  "~/WBB_VAPOR/L1/testing/"
######################################################################################################################################
# SET METADATA
# write out metadata to help data curation - these will be appended to the top of the data files currently. it might be possible in
# later versions to use this to write out an xml file?
######################################################################################################################################

metadata.frame <- read.csv("../metadata_templates/L0_WBB_metadata.csv",header=TRUE)

##########################################################################
# LOAD DAILY DATA, CONCATENATE TO MONTHLY, THEN SEPARATE INTO 
# CALIBRATION AND AMBIENT MEASUREMENT PERIODS
##########################################################################
period <- interval(start.date,end.date)

nmonths <- period %/% months(1)

print("==================================================================")
print(paste("Starting L1 Processing : ",now()))
print("==================================================================")

# Make lists of relevant *CALIBRATION* data files
#---------------------------------------------------
# make an empty list to hold all raw file names
cal.file.list <- list()

# list of all raw files
cal.file.list <- list.files(path=path.to.L0.data,pattern="CalibData",full.names=TRUE,recursive=TRUE)

# find month for each file
cal.file.dates <- extract.date.from.L0.files(cal.file.list)	

# subset list based on period of interest
cal.subset <- cal.file.list[cal.file.dates %within% period]

# Make lists of relevant *CALIBRATION* data files
#---------------------------------------------------
# make an empty list to hold all raw file names
amb.file.list <- list()

# list of all raw files
amb.file.list <- list.files(path=path.to.L0.data,pattern="AmbientData",full.names=TRUE,recursive=TRUE)

# find month for each file
amb.file.dates <- extract.date.from.L0.files(amb.file.list) 

# subset list based on period of interest
amb.subset <- amb.file.list[amb.file.dates %within% period]

# initiate log - need to keep track of how much data has been removed.
log.yyyy <- vector()
log.mm 	 <- vector()
log.pct  <- vector()

# initiate vector to hold data
file.list.by.month <- vector("list",nmonths)
#-----------------------------------------------------------------------------
# loop through months for calibration data.
#-----------------------------------------------------------------------------

for (i in 1:nmonths) {
	# print status header.
	print("==================================================================")
	print(paste("Processing L1 calibration data for month:",month(start.date %m+% months(i-1)),"/",year(start.date %m+% months(i-1))))
	print("==================================================================")

	# CONCATENATE CALIBRATION DATA FIRST.
  	file.list.by.month[[i]] <- as.list(cal.file.list[month(cal.file.dates)==month(start.date %m+% months(i-1)) &
  		year(cal.file.dates)==year(start.date %m+% months(i-1))])

	# check to see if data exists for month in question
	if (length(file.list.by.month[[i]])==0) {
		print(paste("No data for :",month(start.date %m+% months(i-1)),"/",year(start.date %m+% months(i-1))))
		next
	}

	# concatenate to monthly array
	mondata <- concatenate.to.monthly(file.list.by.month[[i]])

	# rough instrument status check and data sanity checks.
	mondata.qflag <- qflag(mondata)
	# clean up memory and gc
	rm(mondata)
	gc()

	# NOTE: $data required here since output of qflag is a list!
	mondata.filtered <- data.sanity.check(mondata.qflag$data) 

	# get data row for log file - 
	log.yyyy[i] <- year(start.date %m+% months(i-1))
	log.mm[i] <- month(start.date %m+% months(i-1))
	# small bug here to fix at some point - the reference points for these are different,
	# so ultimately, this value can exceed 100% in rare cases when all of the data fails
	# one check or the other.
	log.pct[i] <- mondata.filtered$pct.discarded + mondata.qflag$pct.discarded

	# clean up memory and gc
	rm(mondata.qflag)
	gc()

 	# Write calibration data file.
 	##########################################################################

 	if (nrow(mondata.filtered$data) > 0) {
 		print(paste(now()," Writing out calibration data frame..."))
 		# generate output filename
 		coutput.fname <- paste(path.to.output.L1.data,"WBB_Water_Vapor_CalibData_L1_",
 	  	  (start.date %m+% months(i-1)),"_",metadata.frame$Value[metadata.frame$Variable=="code.version"],".dat",sep="")

 		# attach metadata
 		attach.L0.Header(coutput.fname,metadata.frame)

 		# write out data portion of the data file.
 		write.table(mondata.filtered$data,file=coutput.fname,sep=",",append = TRUE,row.names=FALSE)
 	} else {
 		print(paste(now(), "No data passed filters for this month...."))
 	}

 	# clean up calibration data and gc
 	rm(mondata.filtered)
 	gc()
}

# clean up remaining variables corresponding to calibration data. 
rm(file.list.by.month)
gc()

# re-initiate vector to hold data
file.list.by.month <- vector("list",nmonths)

#-----------------------------------------------------------------------------
# loop through months for ambient data.
#-----------------------------------------------------------------------------

for (i in 1:nmonths) {

	# print status header.
	print("==================================================================")
	print(paste("Processing L1 ambient data for month:",month(start.date %m+% months(i-1)),"/",year(start.date %m+% months(i-1))))
	print("==================================================================")

	# CONCATENATE AMIENT DATA SECOND.
  	file.list.by.month[[i]] <- as.list(amb.file.list[month(amb.file.dates)==month(start.date %m+% months(i-1)) &
  		year(amb.file.dates)==year(start.date %m+% months(i-1))])

	# check to see if data exists for month in question
	if (length(file.list.by.month[[i]])==0) {
		print(paste("No data for :",month(start.date %m+% months(i-1)),"/",year(start.date %m+% months(i-1))))
		next
	}

	# concatenate to monthly array
	mondata <- concatenate.to.monthly(file.list.by.month[[i]])

	# rough instrument status check and data sanity checks.
	mondata.qflag <- qflag(mondata)
	# clean up memory and gc
	rm(mondata)
	gc()

	# NOTE: $data required here since output of qflag is a list!
	mondata.filtered <- data.sanity.check(mondata.qflag$data) 

	# get data row for log file - 
	log.yyyy[i] <- year(start.date %m+% months(i-1))
	log.mm[i] <- month(start.date %m+% months(i-1))
	# small bug here to fix at some point - the reference points for these are different,
	# so ultimately, this value can exceed 100% in rare cases when all of the data fails
	# one check or the other.
	log.pct[i] <- mondata.filtered$pct.discarded + mondata.qflag$pct.discarded

	# clean up memory and gc
	rm(mondata.qflag)
	gc()

 	# Write ambient data file.
 	##########################################################################

 	if (nrow(mondata.filtered$data) > 0) {
 		print(paste(now()," Writing out ambient data frame..."))
 		# generate output filename
 		coutput.fname <- paste(path.to.output.L1.data,"WBB_Water_Vapor_AmbientData_L1_",
 	  	  (start.date %m+% months(i-1)),"_",metadata.frame$Value[metadata.frame$Variable=="code.version"],".dat",sep="")

 		# attach metadata
 		attach.L0.Header(coutput.fname,metadata.frame)

 		# write out data portion of the data file.
 		write.table(mondata.filtered$data,file=coutput.fname,sep=",",append = TRUE,row.names=FALSE)
 	} else {
 		print(paste(now(), "No data passed filters for this month...."))
 	}

 	# clean up calibration data and gc
 	rm(mondata.filtered)
 	gc()
}

Rprof(NULL)
# log.fname <- paste("L1_log",".dat",sep="")
# log.output <- data.frame(log.yyyy,log.mm,log.pct)
# write.table(log.output,file=log.fname,sep=",")
