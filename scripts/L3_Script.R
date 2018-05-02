# L3_ambient_processing.R
# richf, 27jan17

# combined script meant to calibrate Picarro vapor data - relies on the water vapor calibration functions
rm(list=ls())

# run all calibration data at once.
# pros: can enforce a consistent fit.
# cons: - hard to pull out nearest ambient value. 

library(lubridate)
library(xts)

# load calibration parameters 
source("../functions/L3_functions.R")

# load user variables
source("../user/L3_user_specs.R")

# LOAD CALIBRATION DATA
period <- interval(start.date,end.date)

nmonths <- period %/% months(1)

print("==================================================================")
print(paste("Starting L3 Processing : ",now()))
print("==================================================================")

# make an empty list to hold all raw file names
raw.file.list <- list()

# list of all raw files
raw.file.list <- list.files(path=path.to.L1.data,pattern="AmbientData",full.names=TRUE,recursive=TRUE)

# find month for each file
raw.file.dates <- extract.date.from.L1.files(raw.file.list,dbg.level=debug)	

# subset list based on 
subset.list <- raw.file.list[raw.file.dates %within% period]

#------------------------------------------------------------------------------
# GET CALIBRATION FILE - note: this is done differently than in previous steps!
calib.data <- readRDS(paste0(path.to.L2.calib.data,calib.file))

# get choice from user file - how should gaps be filled?

# cut out bad indices.
calib.data$regression.data <- 
    calib.data$regression.data[calib.data$regression.data$qflag==1,]

# find 
if (gap.handling == "last.carried.forward") {
    for (i in 2:nrow(calib.data$regression.data)) {
        if (calib.data$regression.data$period.end[i-1] !=
            calib.data$regression.data$period.start[i]) {
            # if current period start doesn't equal last period end,
            # fill gap by changing start of next period to end of last
            # period.
            calib.data$regression.data$period.start[i] <- 
                calib.data$regression.data$period.end[i-1]
            # print report:
            print(paste("Changed index:",i))
        }
    }
} else if (gap.handling == "next.carried.backward") {
    for (i in 2:nrow(calib.data$regression.data)) {
        if (calib.data$regression.data$period.end[i-1] !=
            calib.data$regression.data$period.start[i]) {
            # if current period start doesn't equal last period end,
            # fill gap by changing start of next period to end of last
            # period.
            calib.data$regression.data$period.end[i-1] <- 
                calib.data$regression.data$period.start[i]
            # print report:
            print(paste("Changed index:",i))
        }
    }
} else if (gap.handling == "omit") {
    # I actually think nothing is needed here...
}

# get choice from user file - force first regression to beginning of ambient file?
if (force.to.beginning == TRUE) {
    print(paste("Changing beginning of first period in regression file to match",
        "beginning of ambient period. Here's what these values are in UNIX time:"))
    print(paste(as.POSIXct(start.date,origin="1970-01-01",tz="UTC"),
        as.POSIXct(calib.data$regression.data$period.start[1],
        origin="1970-01-01",tz="UTC")))
    # set first period.start to value of process start time.
    calib.data$regression.data$period.start[1] <- as.numeric(
        as.POSIXct(start.date,origin="1970-01-01",tz="UTC"))
}

# get choice from user file - force last regression to end of ambient file?
if (force.to.end == TRUE) {
    print(paste("Changing end of last period in regression file to match",
        "end of ambient period. Here's what these values are in UNIX time:"))
    print(paste(as.POSIXct(end.date,origin="1970-01-01",tz="UTC"),
        as.POSIXct(calib.data$regression.data$period.end[nrow(calib.data$regression.data)],
        origin="1970-01-01",tz="UTC")))
    # set first period.start to value of process start time.
    calib.data$regression.data$period.end[nrow(calib.data$regression.data)] <- 
        as.numeric(as.POSIXct(end.date,origin="1970-01-01",tz="UTC"))
}

# initiate log - need to keep track of how much data has been removed.
# log.yyyy <- vector()
# log.mm 	 <- vector()
# log.pct  <- vector()

# SET UP METADATA HEADERS
# Get metadata from an L1 file in here to pass along to L2 header...
tmp <- readLines(raw.file.list[[1]]) # get first calibration file.
md.frame <- tmp[grep("#",tmp)] # pull out lines starting with comment character...

# BEGIN LOOPING THROUGH MONTHS.

# start file counter
    fcount <- 1

for (j in 1:nmonths) {
	# print status header.
	print("==================================================================")
	print(paste("Processing L3 ambient data for month:",
        month(start.date %m+% months(j-1)),"/",
        year(start.date %m+% months(j-1))))
	print("==================================================================")

	# ensure that the file we're loading matches the month we think we're processing...
    desired.month <- month(start.date %m+% months(j-1))
    desired.year <- year(start.date %m+% months(j-1))

    tmp <- strsplit(subset.list[[fcount]],split="_")
    file.month <- month(tmp[[1]][7])
    file.year <- year(tmp[[1]][7])

    if (desired.month != file.month | desired.year != file.year) {
        print("File month/year doesn't match expected month/year")
        print("Skip to next file...")
        next
    } else {
        print("File month/year matches expected month/year")
        print(paste("Loading: ",subset.list[[fcount]]))
        # load associated data file.
        ambient.data <- read.table(subset.list[[fcount]],header=TRUE,stringsAsFactors=FALSE,sep=",")
        # increment fcount
        fcount <- fcount+1
    }

    # apply humidity calibration.
    ambient.data <- remove.humidity.dependence(ambient.data)

    # create list of indices.
    amb.inds.in.calib.row <- vector("list",nrow(calib.data$regression.data))

    # find where in the calibration file we're hoping to analyze
    for (i in 1:nrow(calib.data$regression.data)) {
        # for this row, does any data in this day fit within this calib period?
        if (any(ambient.data$EPOCH_TIME > calib.data$regression.data$period.start[i] & 
            ambient.data$EPOCH_TIME <= calib.data$regression.data$period.end[i])) {
            amb.inds.in.calib.row[[i]] <- which(ambient.data$EPOCH_TIME > calib.data$regression.data$period.start[i] & 
            ambient.data$EPOCH_TIME <= calib.data$regression.data$period.end[i])
            print(paste(length(amb.inds.in.calib.row[[i]])," inds in calib row:", i))
        } else {
            print(paste("No ambient inds in calib row:",i))
        }
    }

    Delta_18_16_vsmow <- vector("numeric",nrow(ambient.data))
    Delta_D_H_vsmow <- vector("numeric",nrow(ambient.data))

    # loop through amb.inds.in.calib.row list and use to correct to VSMOW
    for (i in 1:length(amb.inds.in.calib.row)) {
        if (!is.null(amb.inds.in.calib.row[[i]])) { # if there is actually ambient data...
            # calibrate d18O
            Delta_18_16_vsmow[amb.inds.in.calib.row[[i]]] <-
                calib.data$regression.data$O.slope[i]*ambient.data$Delta_18_16[amb.inds.in.calib.row[[i]]] + 
                calib.data$regression.data$O.intercept[i]

            # calibrate d2H
            Delta_D_H_vsmow[amb.inds.in.calib.row[[i]]] <-
                calib.data$regression.data$H.slope[i]*ambient.data$Delta_D_H[amb.inds.in.calib.row[[i]]] + 
                calib.data$regression.data$H.intercept[i]
        }
    }

    # bind calibrated variables to ambient df
    ambient.data <- cbind(ambient.data,Delta_18_16_vsmow)
    ambient.data <- cbind(ambient.data,Delta_D_H_vsmow)

    # average down to desired time interval.
    ambient.data <- reduce.calibrated.ambient.vapor(ambient.data,out.interval)

    # print a quick diagnostic plot
    pdf(paste("diag_plots/ambientData_",file.year,"_",file.month,".pdf",sep=""),width=11,height=8)
    plot(ambient.data$EPOCH_TIME,ambient.data$H2O,type="l")

    plot(ambient.data$EPOCH_TIME,ambient.data$Delta_18_16,type="l")
    points(ambient.data$EPOCH_TIME,ambient.data$Delta_18_16_vsmow,pch=".",col="red")

    plot(ambient.data$EPOCH_TIME,ambient.data$Delta_D_H,type="l")
    points(ambient.data$EPOCH_TIME,ambient.data$Delta_D_H_vsmow,pch=".",col="red")

    dev.off()

    # pdf(paste("diag_plots/ambientData_diff_",file.year,"_",file.month,".pdf",sep=""),width=11,height=8)
    # plot(Delta_18_16_vsmow-ambient.data$Delta_18_16,type="l")
    # dev.off()
    # Save L3 file

    # add uncertainty to the data based on the humidity calibration
    ambient.data <- attach.humidcal.sigmas(ambient.data,humidcal.filename)


    print(paste("Saving calibrated ambient data"))
   
    if (desired.month < 10) {
        calib.dt.name <- paste(path.to.output.L3.data,output.file.prefix,
            desired.year,"_0",desired.month,".dat",sep="")
    } else {
        calib.dt.name <- paste(path.to.output.L3.data,output.file.prefix,
            desired.year,"_",desired.month,".dat",sep="")
    }

    # Attach metadata header.
    attach.L3.Header(calib.dt.name,md.frame,dbg.level=debug)
    
    # write out data table of calib.averages.wamb.mrc.bgc.wstds
    write.table(ambient.data,file=calib.dt.name,
        sep=",",row.names=FALSE)

}