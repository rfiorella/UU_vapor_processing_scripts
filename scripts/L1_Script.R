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

# load user parameters
source("../user/L1_user_specs.R")

# profile code.
Rprof("L1_v1_0_0test")
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
cal.file.list <- list.files(path=path.to.L0.data,pattern="CalibData",full.names=TRUE,recursive=FALSE)

# find month for each file
cal.file.dates <- extract.date.from.L0.files(cal.file.list,dbg.level=debug)	

# subset list based on period of interest
cal.subset <- cal.file.list[cal.file.dates %within% period]

# Make lists of relevant *CALIBRATION* data files
#---------------------------------------------------
# make an empty list to hold all raw file names
amb.file.list <- list()

# list of all raw files
amb.file.list <- list.files(path=path.to.L0.data,pattern="AmbientData",full.names=TRUE,recursive=FALSE)

# find month for each file
amb.file.dates <- extract.date.from.L0.files(amb.file.list,dbg.level=debug) 

# subset list based on period of interest
amb.subset <- amb.file.list[amb.file.dates %within% period]

# initiate log - need to keep track of how much data has been removed.
log.yyyy <- vector()
log.mm 	 <- vector()
log.pctc <- vector()
log.pcta <- vector()
log.nptsc <- vector()
log.nptsa <- vector()
log.timec <- vector()
log.timea <- vector()

# initiate vector to hold data
file.list.by.month <- vector("list",nmonths)

# Get metadata from an L0 file in here to pass along to L1 header...
tmp <- readLines(amb.file.list[[1]])

md.frame <- tmp[grep("#",tmp)] # pull out lines starting with comment character...

#-----------------------------------------------------------------------------
# loop through months for calibration data.
#-----------------------------------------------------------------------------

for (i in 1:nmonths) {

	# set up timer.
	ptm <- proc.time()

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
	mondata <- concatenate.to.monthly(file.list.by.month[[i]],dbg.level=debug)

	# rough instrument status check and data sanity checks.
	mondata.qflag <- qflag(mondata,dbg.level=debug)
	# clean up memory and gc
	rm(mondata)
	gc()

	# NOTE: $data required here since output of qflag is a list!
	mondata.filtered <- data.sanity.check(mondata.qflag$data,dbg.level=debug) 

	# get data row for log file - 
	log.yyyy[i] <- year(start.date %m+% months(i-1))
	log.mm[i] <- month(start.date %m+% months(i-1))
	# small bug here to fix at some point - the reference points for these are different,
	# so ultimately, this value can exceed 100% in rare cases when all of the data fails
	# one check or the other.
	log.pctc[i] <- mondata.filtered$pct.discarded + mondata.qflag$pct.discarded
	log.nptsc[i] <- nrow(mondata.filtered$data)

	# clean up memory and gc
	rm(mondata.qflag)
	gc()

 	# Write calibration data file.
 	##########################################################################

 	if (nrow(mondata.filtered$data) > 0) {
 		print(paste(now()," Writing out calibration data frame..."))
 		# generate output filename
 		coutput.fname <- paste(path.to.output.L1.data,output.file.prefix,"_CalibData_L1_",
 	  	  (start.date %m+% months(i-1)),".dat",sep="")

 		# attach metadata
 		attach.L1.Header(coutput.fname,md.frame,dbg.level=debug)

 		# write out data portion of the data file.
 		write.table(mondata.filtered$data,file=coutput.fname,sep=",",append = TRUE,row.names=FALSE)
 	} else {
 		print(paste(now(), "No data passed filters for this month...."))
 	}

 	# clean up calibration data and gc
 	rm(mondata.filtered)
 	gc()

 	# write out time for calibration data. 
 	tmp <- proc.time() - ptm
 	log.timec[i] <- tmp[3] # take out elapsed time.
 	rm(ptm) # remove prior time marker
 	rm(tmp) # remove tmp time holder
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
	# set up timer.
	ptm <- proc.time()

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
	mondata <- concatenate.to.monthly(file.list.by.month[[i]],dbg.level=debug)

	# rough instrument status check and data sanity checks.
	mondata.qflag <- qflag(mondata,dbg.level=debug)
	# clean up memory and gc
	rm(mondata)
	gc()

	# NOTE: $data required here since output of qflag is a list!
	mondata.qccheck <- data.sanity.check(mondata.qflag$data,dbg.level=debug) 

	# small bug here to fix at some point - the reference points for these are different,
	# so ultimately, this value can exceed 100% in rare cases when all of the data fails
	# one check or the other.
	log.pcta[i] <- mondata.qccheck$pct.discarded + mondata.qccheck$pct.discarded

	# ambient data needs a bit more filtering here - namely, after a calibration period
	# there are often apparent discontinuities in the data that are problematic...

	mondata.pcfilter <- post.calibration.filter(mondata.qccheck$data,dbg.level=debug)

	# get number of rows.
	log.nptsa[i] <- nrow(mondata.pcfilter)

	# clean up memory and gc
	rm(mondata.qflag)
	rm(mondata.qccheck)
	gc()

 	# Write ambient data file.
 	##########################################################################

 	if (nrow(mondata.pcfilter) > 0) {
 		print(paste(now()," Writing out ambient data frame..."))
 		# generate output filename
 		coutput.fname <- paste(path.to.output.L1.data,output.file.prefix,"_AmbientData_L1_",
 	  	  (start.date %m+% months(i-1)),".dat",sep="")

 		# attach metadata
 		attach.L1.Header(coutput.fname,md.frame,dbg.level=debug)

 		# write out data portion of the data file.
 		write.table(mondata.pcfilter,file=coutput.fname,sep=",",append = TRUE,row.names=FALSE)
 	} else {
 		print(paste(now(), "No data passed filters for this month..."))
 	}

 	# clean up calibration data and gc
 	rm(mondata.pcfilter)
 	gc()

 	 # write out time for calibration data. 
 	tmp <- proc.time() - ptm
 	log.timea[i] <- tmp[3] # take out elapsed time.
 	rm(ptm) # remove prior time marker
 	rm(tmp) # remove prior temporary placeholder
 }

Rprof(NULL)
log.fname <- paste("L1_log_v110beta",".dat",sep="")
log.output <- data.frame(log.yyyy,log.mm,log.pctc,log.pcta,log.nptsc,log.nptsa,log.timec,log.timea)
write.table(log.output,file=log.fname,sep=",")
