# L3_ambient_processing.R
# richf, 27jan17

# combined script meant to calibrate Picarro vapor data - relies on the water vapor calibration functions
rm(list=ls())

# run all calibration data at once.
# pros: can enforce a consistent fit.
# cons: - hard to pull out nearest ambient value. 

library(lubridate)
library(zoo)

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

# GET CALIBRATION FILE - note: this is done differently than in previous steps!!!!!
calib.data <- read.table(paste(path.to.L2.calib.data,calib.file,sep=""),
    sep=",",stringsAsFactors=FALSE,header=TRUE)

# initiate log - need to keep track of how much data has been removed.
log.yyyy <- vector()
log.mm 	 <- vector()
log.pct  <- vector()

# BEGIN LOOPING THROUGH MONTHS.

# start file counter
    fcount <- 1

for (j in 1:nmonths) {
	# print status header.
	print("==================================================================")
	print(paste("Processing L3 ambient data for month:",month(start.date %m+% months(j-1)),"/",year(start.date %m+% months(j-1))))
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

    # print some diagnostics
    # print(head(ambient.data))
    # print(colnames(calib.data))
    # print(nrow(calib.data))

    amb.inds.in.calib.row <- vector("list",nrow(calib.data))

    # find where in the calibration file we're hoping to analyze
    for (i in 1:nrow(calib.data)) {
        # for this row, does any data in this day fit within this calib period?
        if (any(ambient.data$EPOCH_TIME > calib.data$start.time[i] & 
            ambient.data$EPOCH_TIME <= calib.data$stop.time[i])) {
            amb.inds.in.calib.row[[i]] <- which(ambient.data$EPOCH_TIME > calib.data$start.time[i] & 
            ambient.data$EPOCH_TIME <= calib.data$stop.time[i])
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
                calib.data$O.slope[i]*ambient.data$Delta_18_16[amb.inds.in.calib.row[[i]]] + 
                calib.data$O.intercept[i]

            # calibrate d2H
            Delta_D_H_vsmow[amb.inds.in.calib.row[[i]]] <-
                calib.data$H.slope[i]*ambient.data$Delta_D_H[amb.inds.in.calib.row[[i]]] + 
                calib.data$H.intercept[i]
        }
    }

    # print a quick diagnostic plot
    pdf(paste("diag_plots/ambientData_",file.year,"_",file.month,".pdf",sep=""),width=11,height=8)
    plot(ambient.data$EPOCH_TIME,ambient.data$H2O,pch=".")

    plot(ambient.data$EPOCH_TIME,ambient.data$Delta_18_16,pch=".")
    points(ambient.data$EPOCH_TIME,Delta_18_16_vsmow,pch=".",col="red")

    plot(ambient.data$EPOCH_TIME,ambient.data$Delta_D_H,type="l")
    points(ambient.data$EPOCH_TIME,Delta_D_H_vsmow,pch=".",col="red")

    dev.off()


    pdf(paste("diag_plots/ambientData_diff_",file.year,"_",file.month,".pdf",sep=""),width=11,height=8)
    plot(Delta_18_16_vsmow-ambient.data$Delta_18_16,type="l")
    dev.off()
    # Save L3 file

    # bind calibrated variables to ambient df
    ambient.data <- cbind(ambient.data,Delta_18_16_vsmow)
    ambient.data <- cbind(ambient.data,Delta_D_H_vsmow)

    print(paste("Saving calibrated ambient data"))
   
    if (desired.month < 10) {
        calib.dt.name <- paste(path.to.output.L3.data,output.file.prefix,
            desired.year,"_0",desired.month,".dat",sep="")
    } else {
        calib.dt.name <- paste(path.to.output.L3.data,output.file.prefix,
            desired.year,"_",desired.month,".dat",sep="")
    }

    # write out data table of calib.averages.wamb.mrc.bgc.wstds
    write.table(ambient.data,file=calib.dt.name,
        sep=",",row.names=FALSE)


}