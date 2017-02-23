# L2_combined_script.R
# richf, 27jan17

# combined script meant to calibrate Picarro vapor data - relies on the water vapor calibration functions
rm(list=ls())

# run all calibration data at once.
# pros: can enforce a consistent fit.
# cons: - hard to pull out nearest ambient value. 

# load associated functions
#source("../src/L1_Functions_osx.R")
source("../src/data_structure.R")
#source("../src/quality_checks.R")
source("../src/wv_standard_calibration_functions.R")

library(data.table)
library(lubridate)
library(zoo)
library(RColorBrewer)

# load calibration parameters 
# source("L1_Calibration_Parameters_osx.R")

# Set user variables

start.date <- ymd("2013-05-01")
end.date <- ymd("2017-02-01")

path.to.L1.data <- "~/WBB_VAPOR/L1/testing/"
path.to.output.L2.data <-  "~/WBB_VAPOR/L2/testing/"

RUN_PLOTS <- "TRUE"

# LOAD CALIBRATION DATA
period <- interval(start.date,end.date)

nmonths <- period %/% months(1)

print("==================================================================")
print(paste("Starting L2 Processing : ",now()))
print("==================================================================")

# make an empty list to hold all raw file names
raw.file.list <- list()

# list of all raw files
raw.file.list <- list.files(path=path.to.L1.data,pattern="CalibData",full.names=TRUE,recursive=TRUE)

# find month for each file
raw.file.dates <- extract.date.from.L1.files(raw.file.list)	

# subset list based on 
subset.list <- raw.file.list[raw.file.dates %within% period]

# initiate log - need to keep track of how much data has been removed.
log.yyyy <- vector()
log.mm 	 <- vector()
log.pct  <- vector()

# BEGIN LOOPING THROUGH MONTHS.

# preallocate variables
calib.avgs <- vector("list",nmonths)

# start file counter
    fcount <- 1

for (i in 1:nmonths) {
	# print status header.
	print("==================================================================")
	print(paste("Processing L2 data for month:",month(start.date %m+% months(i-1)),"/",year(start.date %m+% months(i-1))))
	print("==================================================================")
    
    # ensure that the file we're loading matches the month we think we're processing...
    desired.month <- month(start.date %m+% months(i-1))
    desired.year <- year(start.date %m+% months(i-1))

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
        calib.data <- read.table(subset.list[[fcount]],header=TRUE,stringsAsFactors=FALSE,sep=",")
        # increment fcount
        fcount <- fcount+1
    }

	# ##########################################################################
	# # PROCESS CALIBRATION DATA
	# ##########################################################################
    # set month for diagnostic plots
    time.suffix <- paste(year(start.date %m+% months(i-1)),"_",month(start.date %m+% months(i-1)),collapse="",sep="")

	# identify breakpoints
	breaks <- ID.calib.breakpoints(calib.data)

	# calculate initial set of splines as a check on breaks
	spline.fits <- fit.calibration.splines(calib.data,breaks)
	

    if (RUN_PLOTS) {
        # make plots of the identified peaks with associated splines.
        print("Making plot showing calibration data and calculated splines...")
        pdf(paste("diag_plots/calib_raw_peakssplines_",time.suffix,".pdf",collapse="",sep=""),width=11,height=8)
        par(mfrow=c(3,1),mar=c(4,4,1,1),oma=c(0,0,3,0))
        for (k in 1:length(spline.fits)) {
            plot(calib.data$EPOCH_TIME[breaks[k]:(breaks[k+1]-1)],calib.data$H2O[breaks[k]:(breaks[k+1]-1)])
            lines(spline.fits[[k]]$time,spline.fits[[k]]$H2O,col="red",lwd=2)

            plot(calib.data$EPOCH_TIME[breaks[k]:(breaks[k+1]-1)],calib.data$Delta_18_16[breaks[k]:(breaks[k+1]-1)])
            lines(spline.fits[[k]]$time,spline.fits[[k]]$d18O,col="red",lwd=2)

            plot(calib.data$EPOCH_TIME[breaks[k]:(breaks[k+1]-1)],calib.data$Delta_D_H[breaks[k]:(breaks[k+1]-1)])
            lines(spline.fits[[k]]$time,spline.fits[[k]]$d2H,col="red",lwd=2)
        }
        dev.off() 
    }

	# calculate derivatives of splines...
	spline.derivatives <- calculate.spline.derivatives(spline.fits,breaks)

    if (RUN_PLOTS) {
        print("Making plot of spline derivatives...")
        # make plots of the spline derivatives.        
        pdf(paste("diag_plots/raw_data_splinederivs_",time.suffix,".pdf",collapse="",sep=""),width=11,height=8)
        par(mfrow=c(3,1),mar=c(4,4,1,1),oma=c(0,0,3,0))
        for (k in 1:length(spline.derivatives)) { 
                plot(spline.derivatives[[k]]$d_H2O,type="l",ylim=c(-50,50))
                #lines(spline.derivatives[[k]]$d2_H2O,col="red")
                plot(spline.derivatives[[k]]$d_d18O,type="l",ylim=c(-0.01,0.01))
                #lines(spline.derivatives[[k]]$d2_d18O,col="red")
                plot(spline.derivatives[[k]]$d_d2H,type="l",ylim=c(-0.1,0.1))
                #lines(spline.derivatives[[k]]$d2_d2H,col="red")
        }
        dev.off() 
    }
    
	# calculate list of indices to keep...
	good.inds <- extract.stable.calib.indices(spline.derivatives)

    # quick filter here - should there be a minimum number of "good" points
    # required to define an analysis as a "standard" and retain the averages?
    # this can be done in ~3-4 lines, so no need for a function here.
    for (j in 1:length(good.inds)) {
        if (length(good.inds[[j]])<60) {
            good.inds[[j]] <- rep(NA,length(good.inds[[j]]))
        }
    }

	# calculate averages for each plateau
	calib.avgs[[i]] <- calculate.standard.averages(calib.data,good.inds)

    print(length(good.inds))
	if (RUN_PLOTS) {
        # plot the data after identifying plateaus using splines
        pdf(paste("diag_plots/spline_selected_points_",time.suffix,".pdf",collapse="",sep=""),width=11,height=8)
        par(mfrow=c(3,1),mar=c(4,4,1,1),oma=c(0,0,3,0))
        mycols <- rainbow(8)
        if (!is.null(calib.avgs[[i]])) {
            for (k in 1:nrow(calib.avgs[[i]])) { 
                plot(spline.fits[[k]]$time,spline.fits[[k]]$H2O)
                points(calib.data$EPOCH_TIME[good.inds[[k]]],
                    calib.data$H2O[good.inds[[k]]],col=ifelse(k%%8 != 0, mycols[k%%8], mycols[8]))
                mtext(paste("Mean: ",round(calib.avgs[[i]]$H2O.mean[k],2)," stdev: ",round(calib.avgs[[i]]$H2O.sd[k],2)),side=1,line=-2)

                plot(spline.fits[[k]]$time,spline.fits[[k]]$d18O)
                points(calib.data$EPOCH_TIME[good.inds[[k]]],
                    calib.data$Delta_18_16[good.inds[[k]]],col=ifelse(k%%8 != 0, mycols[k%%8], mycols[8]))
                mtext(paste("Mean: ",round(calib.avgs[[i]]$d18O.mean[k],2)," stdev: ",round(calib.avgs[[i]]$d18O.sd[k],2)),side=1,line=-2)

                plot(spline.fits[[k]]$time,spline.fits[[k]]$d2H)
                points(calib.data$EPOCH_TIME[good.inds[[k]]],
                    calib.data$Delta_D_H[good.inds[[k]]],col=ifelse(k%%8 != 0, mycols[k%%8], mycols[8]))
                mtext(paste("Mean: ",round(calib.avgs[[i]]$d2H.mean[k],2)," stdev: ",round(calib.avgs[[i]]$d2H.sd[k],2)),side=1,line=-2)

                title(paste("Spline selected points (red) for : ",time.suffix),outer=TRUE)
            }
            dev.off()
        }
	}

    if (RUN_PLOTS) {
        if (!is.null(calib.avgs[[i]])) {
            # plot monthly average data
            pdf(paste("diag_plots/monavgs_",time.suffix,".pdf",collapse="",sep=""),width=11,height=8)
            par(mfrow=c(3,2))
            mycols <- brewer.pal(6,"Paired")
            inds <- vector("list",6)

            inds[[1]] <- which(calib.avgs[[i]]$H2O.mean <= 5000)
            inds[[2]] <- which(calib.avgs[[i]]$H2O.mean > 5000 & calib.avgs[[i]]$H2O.mean <= 10000)
            inds[[3]] <- which(calib.avgs[[i]]$H2O.mean > 10000 & calib.avgs[[i]]$H2O.mean <= 15000)
            inds[[4]] <- which(calib.avgs[[i]]$H2O.mean > 15000 & calib.avgs[[i]]$H2O.mean <= 20000)
            inds[[5]] <- which(calib.avgs[[i]]$H2O.mean > 20000 & calib.avgs[[i]]$H2O.mean <= 25000)
            inds[[6]] <- which(calib.avgs[[i]]$H2O.mean > 25000)


            plot(calib.avgs[[i]]$time.mean[inds[[1]]],calib.avgs[[i]]$H2O.mean[inds[[1]]],col=mycols[1],
                ylim=c(min(calib.avgs[[i]]$H2O.mean,na.rm=TRUE),max(calib.avgs[[i]]$H2O.mean,na.rm=TRUE)),
                xlim=c(min(calib.avgs[[i]]$time.mean,na.rm=TRUE),max(calib.avgs[[i]]$time.mean,na.rm=TRUE)))
            points(calib.avgs[[i]]$time.mean[inds[[2]]],calib.avgs[[i]]$H2O.mean[inds[[2]]],col=mycols[2])
            points(calib.avgs[[i]]$time.mean[inds[[3]]],calib.avgs[[i]]$H2O.mean[inds[[3]]],col=mycols[3])
            points(calib.avgs[[i]]$time.mean[inds[[4]]],calib.avgs[[i]]$H2O.mean[inds[[4]]],col=mycols[4])
            points(calib.avgs[[i]]$time.mean[inds[[5]]],calib.avgs[[i]]$H2O.mean[inds[[5]]],col=mycols[5])
            points(calib.avgs[[i]]$time.mean[inds[[6]]],calib.avgs[[i]]$H2O.mean[inds[[6]]],col=mycols[6])

            plot(calib.avgs[[i]]$time.mean[inds[[1]]],calib.avgs[[i]]$H2O.sd[inds[[1]]],col=mycols[1],
                ylim=c(min(calib.avgs[[i]]$H2O.sd,na.rm=TRUE),max(calib.avgs[[i]]$H2O.sd,na.rm=TRUE)),
                xlim=c(min(calib.avgs[[i]]$time.mean,na.rm=TRUE),max(calib.avgs[[i]]$time.mean,na.rm=TRUE)))
            points(calib.avgs[[i]]$time.mean[inds[[2]]],calib.avgs[[i]]$H2O.sd[inds[[2]]],col=mycols[2])
            points(calib.avgs[[i]]$time.mean[inds[[3]]],calib.avgs[[i]]$H2O.sd[inds[[3]]],col=mycols[3])
            points(calib.avgs[[i]]$time.mean[inds[[4]]],calib.avgs[[i]]$H2O.sd[inds[[4]]],col=mycols[4])
            points(calib.avgs[[i]]$time.mean[inds[[5]]],calib.avgs[[i]]$H2O.sd[inds[[5]]],col=mycols[5])
            points(calib.avgs[[i]]$time.mean[inds[[6]]],calib.avgs[[i]]$H2O.sd[inds[[6]]],col=mycols[6])

            plot(calib.avgs[[i]]$time.mean[inds[[1]]],calib.avgs[[i]]$d18O.mean[inds[[1]]],col=mycols[1],
                ylim=c(min(calib.avgs[[i]]$d18O.mean,na.rm=TRUE),max(calib.avgs[[i]]$d18O.mean,na.rm=TRUE)),
                xlim=c(min(calib.avgs[[i]]$time.mean,na.rm=TRUE),max(calib.avgs[[i]]$time.mean,na.rm=TRUE)))
            points(calib.avgs[[i]]$time.mean[inds[[2]]],calib.avgs[[i]]$d18O.mean[inds[[2]]],col=mycols[2])
            points(calib.avgs[[i]]$time.mean[inds[[3]]],calib.avgs[[i]]$d18O.mean[inds[[3]]],col=mycols[3])
            points(calib.avgs[[i]]$time.mean[inds[[4]]],calib.avgs[[i]]$d18O.mean[inds[[4]]],col=mycols[4])
            points(calib.avgs[[i]]$time.mean[inds[[5]]],calib.avgs[[i]]$d18O.mean[inds[[5]]],col=mycols[5])
            points(calib.avgs[[i]]$time.mean[inds[[6]]],calib.avgs[[i]]$d18O.mean[inds[[6]]],col=mycols[6])

            plot(calib.avgs[[i]]$time.mean[inds[[1]]],calib.avgs[[i]]$d18O.sd[inds[[1]]],col=mycols[1],
                ylim=c(min(calib.avgs[[i]]$d18O.sd,na.rm=TRUE),max(calib.avgs[[i]]$d18O.sd,na.rm=TRUE)),
                xlim=c(min(calib.avgs[[i]]$time.mean,na.rm=TRUE),max(calib.avgs[[i]]$time.mean,na.rm=TRUE)))
            points(calib.avgs[[i]]$time.mean[inds[[2]]],calib.avgs[[i]]$d18O.sd[inds[[2]]],col=mycols[2])
            points(calib.avgs[[i]]$time.mean[inds[[3]]],calib.avgs[[i]]$d18O.sd[inds[[3]]],col=mycols[3])
            points(calib.avgs[[i]]$time.mean[inds[[4]]],calib.avgs[[i]]$d18O.sd[inds[[4]]],col=mycols[4])
            points(calib.avgs[[i]]$time.mean[inds[[5]]],calib.avgs[[i]]$d18O.sd[inds[[5]]],col=mycols[5])
            points(calib.avgs[[i]]$time.mean[inds[[6]]],calib.avgs[[i]]$d18O.sd[inds[[6]]],col=mycols[6])

            plot(calib.avgs[[i]]$time.mean[inds[[1]]],calib.avgs[[i]]$d2H.mean[inds[[1]]],col=mycols[1],
                ylim=c(min(calib.avgs[[i]]$d2H.mean,na.rm=TRUE),max(calib.avgs[[i]]$d2H.mean,na.rm=TRUE)),
                xlim=c(min(calib.avgs[[i]]$time.mean,na.rm=TRUE),max(calib.avgs[[i]]$time.mean,na.rm=TRUE)))
            points(calib.avgs[[i]]$time.mean[inds[[2]]],calib.avgs[[i]]$d2H.mean[inds[[2]]],col=mycols[2])
            points(calib.avgs[[i]]$time.mean[inds[[3]]],calib.avgs[[i]]$d2H.mean[inds[[3]]],col=mycols[3])
            points(calib.avgs[[i]]$time.mean[inds[[4]]],calib.avgs[[i]]$d2H.mean[inds[[4]]],col=mycols[4])
            points(calib.avgs[[i]]$time.mean[inds[[5]]],calib.avgs[[i]]$d2H.mean[inds[[5]]],col=mycols[5])
            points(calib.avgs[[i]]$time.mean[inds[[6]]],calib.avgs[[i]]$d2H.mean[inds[[6]]],col=mycols[6])

            plot(calib.avgs[[i]]$time.mean[inds[[1]]],calib.avgs[[i]]$d2H.sd[inds[[1]]],col=mycols[1],
                ylim=c(min(calib.avgs[[i]]$d2H.sd,na.rm=TRUE),max(calib.avgs[[i]]$d2H.sd,na.rm=TRUE)),
                xlim=c(min(calib.avgs[[i]]$time.mean,na.rm=TRUE),max(calib.avgs[[i]]$time.mean,na.rm=TRUE)))
            points(calib.avgs[[i]]$time.mean[inds[[2]]],calib.avgs[[i]]$d2H.sd[inds[[2]]],col=mycols[2])
            points(calib.avgs[[i]]$time.mean[inds[[3]]],calib.avgs[[i]]$d2H.sd[inds[[3]]],col=mycols[3])
            points(calib.avgs[[i]]$time.mean[inds[[4]]],calib.avgs[[i]]$d2H.sd[inds[[4]]],col=mycols[4])
            points(calib.avgs[[i]]$time.mean[inds[[5]]],calib.avgs[[i]]$d2H.sd[inds[[5]]],col=mycols[5])
            points(calib.avgs[[i]]$time.mean[inds[[6]]],calib.avgs[[i]]$d2H.sd[inds[[6]]],col=mycols[6])

            title(paste("Monthly data averages: ",time.suffix),outer=TRUE)
            dev.off()
        }
    }
}

# make a line plot of how many calibrations were performed versus retained
# from month to month.

ncal <- vector("numeric",nmonths)
ngcal <- vector("numeric",nmonths)

for (i in 1:nmonths) {
    #ncal[i] <-
    if (!is.null(calib.avgs[[i]])) {
        ngcal[i] <- nrow(calib.avgs[[i]])
    } else {
        ngcal[i] <- NA
    }
}

pdf("calibration_number.pdf")
plot(ngcal,type="l")
dev.off()


# # filter averages, and assign a standard name and concentration to the table
# # filter = flag values that have an unusual H2O concentration
# # bind all elements of the list together
calib.averages.all <- do.call(rbind,calib.avgs)

# clean up some of the averages - remove values with:
# low mean H2O, high sd H2O, d18O, and d2H

# h2o.min <- which(calib.averages.all$H2O.mean > 2000) # must be 2000 ppm
# h2o.max <- which(calib.averages.all$H2O.mean < 25000)
# h2o.lsd <- which(calib.averages.all$H2O.sd < 1000)
# d18O.lsd <- which(calib.averages.all$d18O.sd < 0.5)
# d2H.lsd <- which(calib.averages.all$d2H.sd < 4)

# common.inds <- Reduce(intersect,list(h2o.min,h2o.max,h2o.lsd,d18O.lsd,d2H.lsd))

# #print(common.inds)

# quartz()
# par(mfrow=c(3,2))
# plot(calib.averages.all$time.mean,calib.averages.all$H2O.mean)
# plot(calib.averages.all$time.mean,calib.averages.all$H2O.sd)
# plot(calib.averages.all$time.mean,calib.averages.all$d18O.mean)
# plot(calib.averages.all$time.mean,calib.averages.all$d18O.sd)
# plot(calib.averages.all$time.mean,calib.averages.all$d2H.mean)
# plot(calib.averages.all$time.mean,calib.averages.all$d2H.sd)

# quartz()
# par(mfrow=c(3,1))
# plot(calib.averages.all$time.mean[common.inds],
#     calib.averages.all$H2O.mean[common.inds])
# plot(calib.averages.all$time.mean[common.inds],
#     calib.averages.all$d18O.mean[common.inds])
# plot(calib.averages.all$time.mean[common.inds],
#     calib.averages.all$d2H.mean[common.inds])
