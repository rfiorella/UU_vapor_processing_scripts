# L2_combined_script.R
# richf, 27jan17

# combined script meant to calibrate Picarro vapor data - relies on the water vapor calibration functions
rm(list=ls())

# run all calibration data at once.
# pros: can enforce a consistent fit.
# cons: - hard to pull out nearest ambient value. 

# load associated functions
source("../functions/L2_functions.R")

# LOAD REQUIRED LIBRARIES
#==========================================================================
#library(data.table)
library(lubridate)
library(zoo)
library(RColorBrewer)


# Set user variables
#--------------------------------------------------------------------------
# set initial and final dates for processing.
start.date <- ymd("2013-05-01")
end.date <- ymd("2014-06-01")

# where is the input data and where should we write the output data?
path.to.L1.data <- "~/WBB_VAPOR/L1/testing/"
path.to.output.L2.data <-  "~/WBB_VAPOR/L2/testing/"

# run several pdf diagnostic plots? (logical value)
RUN_PLOTS <- "TRUE"

#--------------------------------------------------------------------------
# Begin processing calibration data...

print("==================================================================")
print(paste("Starting L2 Processing : ",now()))
print("==================================================================")

# define the period given by start, end dates
period <- interval(start.date,end.date)

# calculate the number of months represented by start,end dates
# (script is designed to work on complete months) 
nmonths <- period %/% months(1)

# Make lists of relevant *CALIBRATION* data files
#---------------------------------------------------
# make an empty list to hold all raw file names
raw.file.list <- list()

# list of all raw files
raw.file.list <- list.files(path=path.to.L1.data,pattern="CalibData",full.names=TRUE,recursive=TRUE)

# find month for each file
raw.file.dates <- extract.date.from.L1.files(raw.file.list) 

# subset list based on period of interest
subset.list <- raw.file.list[raw.file.dates %within% period]

# Make lists of relevant *AMBIENT* data files
#---------------------------------------------------
# make an empty list to hold all raw file names
amb.file.list <- list()

# list of all raw files
amb.file.list <- list.files(path=path.to.L1.data,pattern="AmbientData",full.names=TRUE,recursive=TRUE)

# find month for each file
amb.file.dates <- extract.date.from.L1.files(amb.file.list) 

# subset list based on period of interest
amb.subset <- amb.file.list[amb.file.dates %within% period]

# # initiate log - need to keep track of how much data has been removed.
# log.yyyy <- vector()
# log.mm 	 <- vector()
# log.pct  <- vector()

# preallocate variables
calib.avgs.nomem <- vector("list",nmonths)
calib.avgs <- vector("list",nmonths)
ambient.buffers <- vector("list",nmonths)

# start file counter
fcount <- 1
gcount <- 1 # perhaps not necessary, but coding defensively here in a rush...

# OK, after setting up everything, start looping through the months requested...
for (i in 1:nmonths) {
	# print status header.
	print("==================================================================")
	print(paste("Processing L2 data for month:",month(start.date %m+% months(i-1)),"/",year(start.date %m+% months(i-1))))
	print("==================================================================")
    
    # set month for diagnostic plots
    time.suffix <- paste(year(start.date %m+% months(i-1)),"_",month(start.date %m+% months(i-1)),collapse="",sep="")

    #-------------------------------------------------------------------------
    # FILE LOAD CHECK
    # ensure that the file we're loading matches the month we think we're processing...
    #-------------------------------------------------------------------------

    desired.month <- month(start.date %m+% months(i-1))
    desired.year <- year(start.date %m+% months(i-1))

    # split the file name to extract the date
    tmp <- strsplit(subset.list[[fcount]],split="_")
    file.month <- month(tmp[[1]][7]) # index 7 of resulting list has the date info
    file.year <- year(tmp[[1]][7])   # index 7 of resulting list has the date info

    # logical test - does the mm/yy of the file match our expectations?
    if (desired.month != file.month | desired.year != file.year) {
        print("File month/year doesn't match expected month/year")
        print("Skip to next file...")
        next # break loop, skip to next index of i
    } else { # e.g., mm/yy matches expectations
        print("File month/year matches expected month/year")
        print(paste("Loading: ",subset.list[[fcount]]))
        # load associated data file.
        calib.data <- read.table(subset.list[[fcount]],header=TRUE,stringsAsFactors=FALSE,sep=",")
        # increment fcount so that next file is attempted the next time through the loop
        fcount <- fcount+1
    }

	##########################################################################
	# PROCESS CALIBRATION DATA
	##########################################################################
    
	# identify breakpoints
	breaks <- ID.calib.breakpoints(calib.data)

	# calculate initial set of splines as a check on breaks
	spline.fits <- fit.calibration.splines(calib.data,breaks)
	
    #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    # DIAGNOSTIC PLOT 1: Plot each peak period identified by ID.calib.breakpoints
    # and the splines that are fit to each peak period.

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
    #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

	# now that we have nice smooth splines, estimate derivatives...
	spline.derivatives <- calculate.spline.derivatives(spline.fits,breaks)

    #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    # DIAGNOSTIC PLOT 2: plot the first derivative of each identified peak period
    # to help identify suitable thresholds for identifying stable portions of peaks.

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
    #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    
	# based on threshold values calculate list of indices to keep from each peak
	good.inds <- extract.stable.calib.indices(spline.derivatives)

    # quick filter here - should there be a minimum number of "good" points
    # required to define an analysis as a "standard" and retain the averages?
    # this can be done in ~3-4 lines, so no need for a function here.
    for (j in 1:length(good.inds)) {
        if (length(good.inds[[j]])<60) {
            good.inds[[j]] <- rep(NA,length(good.inds[[j]]))
        }
    }


	# based on the good points identified, calculate averages for each plateau
	calib.avgs[[i]] <- calculate.standard.averages(calib.data,good.inds)

    # attach nearby ambient data to calibration data frame.
    #-------------------------------------------------------

    # load relevant ambient data file
    #---------------------------------
    tmp <- strsplit(amb.subset[[gcount]],split="_")
    file.month <- month(tmp[[1]][7])
    file.year <- year(tmp[[1]][7])

    #print(paste(desired.month,desired.year,file.month,file.year))

    while (desired.month != file.month | desired.year != file.year) {
        print(paste(Sys.time(),"File month/year doesn't match expected month/year"))
        print("Likely indicates that ambient data was collected for this month, but no calibration data...")
        print("Skip to next file...")
        # increment to next file
        gcount <- gcount + 1
        # reset file.month, file.year based on new gcount.
        tmp <- strsplit(amb.subset[[gcount]],split="_")
        file.month <- month(tmp[[1]][7])
        file.year <- year(tmp[[1]][7])
    } 

        print("Found the file where month/year matches expected month/year")
        print(paste("Loading: ",amb.subset[[gcount]]))
        # load associated data file.
        amb.data <- read.table(amb.subset[[gcount]],header=TRUE,stringsAsFactors=FALSE,sep=",")
        # increment gcount
        gcount <- gcount+1    

    ambient.buffers[[i]] <-  get.ambient.deltas(calib.avgs[[i]],amb.data)

    #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    # DIAGNOSTIC PLOT 3: plot the points selected by the spline thresholds 
    # for each peak. raw data are shown as points, splines are shown as lines.

    if (RUN_PLOTS) {
        # plot the data after identifying plateaus using splines
        pdf(paste("diag_plots/spline_selected_points_",time.suffix,".pdf",
            collapse="",sep=""),width=11,height=8)
        par(mfrow=c(3,1),mar=c(4,4,1,1),oma=c(0,0,3,0))
        mycols <- rainbow(8)
        if (length(calib.avgs) > 0) { # i : months
            if (!is.null(calib.avgs[[i]]) == TRUE & length(calib.avgs[[i]]) > 0) {
                # start separate counter
                j <- 1
                z <- 1
                for (k in 1:length(spline.fits)) { # k loops through spline.fits in month i
                    # problem with straightaway plotting this: length of spline.fits[[k]] and good.inds[[k]]
                    # is NOT the same if points are removed through the filter.
                    # keep a separate counter for good.inds rather than spline fits.
                   
                    # want to line up all of the data with the proper times - there are three different
                    # lists that have different lengths and therefore different time variables
                    # this is a bit of an obnoxious task. there are three possible outcomes here: 
                    # 1) no points in the spline fit passed the derivatives based test 
                    # identifying a stable point, and no data is returned (e.g., spline.has.data = FALSE)
                    # 2) points returned in spline fit, but does not pass filter in averaging
                    # (e.g., spline.has.data = TRUE but spline.has.good.data = FALSE)
                    # 3) points returned in spline fit, and does pass filter in averaging
                    # (e.g., spline.has.data = TRUE and spline.has.good.data = TRUE)

                    # four different indices sprinkled through here: i refers to the month,
                    # k refers to the peak in spline fits, z refers to the peak with good
                    # inds returned, and j refers to the peaks with good inds that pass the filter.
                    # in general, the maximums of these indices should follow this rule: k >= z >= j.

                    tmpmin <- min(spline.fits[[k]]$time)
                    tmpmax <- max(spline.fits[[k]]$time)
                    gimean <- mean(calib.data$EPOCH_TIME[good.inds[[z]]],na.rm=TRUE)
                    camean <- calib.avgs[[i]]$time.mean[j]

                    # print(paste(j,z,k))
                    # print(paste(tmpmin,tmpmax,calib.avgs[[i]]$time.mean[j]))
                    # print(paste(tmpmin,tmpmax,gimean))

                    # first test: was a peak found?
                    if (gimean < tmpmax & gimean >= tmpmin) {
                        # yes, a peak was found.
                        spline.has.data <- "TRUE"
                        # second test: did any data make it through the filter
                        if (camean < tmpmax & camean >= tmpmin) { # yes, data made it through filter
                            spline.has.good.data <- "TRUE"
                        } else { # no, data did not make it through filter.
                            spline.has.good.data <- "FALSE"
                        }
                    } else { # no peak was found.
                        spline.has.data <- "FALSE"
                        spline.has.good.data <- "FALSE"
                    }
                        
                    # print(paste(j,z,k,length(calib.avgs[[i]]$time.mean),length(good.inds),
                    #       length(spline.fits),spline.has.data,spline.has.good.data))

                    # finally make plots.
                    if (spline.has.data==TRUE & spline.has.good.data==TRUE) {
                        plot(spline.fits[[k]]$time,spline.fits[[k]]$H2O,
                            ylim=c(min(c(min(spline.fits[[k]]$H2O,na.rm=TRUE),min(calib.data$H2O[good.inds[[z]]],na.rm=TRUE))),
                            max(c(max(spline.fits[[k]]$H2O,na.rm=TRUE),max(calib.data$H2O[good.inds[[z]]],na.rm=TRUE)))))
                        points(calib.data$EPOCH_TIME[good.inds[[z]]],
                            calib.data$H2O[good.inds[[z]]],col=ifelse(k%%8 != 0, mycols[k%%8], mycols[8]))
                         mtext(paste("Mean: ",round(calib.avgs[[i]]$H2O.mean[j],2),
                             " stdev: ",round(calib.avgs[[i]]$H2O.sd[j],2)),side=1,line=-2)

                        plot(spline.fits[[k]]$time,spline.fits[[k]]$d18O,
                            ylim=c(min(c(min(spline.fits[[k]]$d18O,na.rm=TRUE),min(calib.data$Delta_18_16[good.inds[[z]]],na.rm=TRUE))),
                            max(c(max(spline.fits[[k]]$d18O,na.rm=TRUE),max(calib.data$Delta_18_16[good.inds[[z]]],na.rm=TRUE)))))
                        points(calib.data$EPOCH_TIME[good.inds[[z]]],
                            calib.data$Delta_18_16[good.inds[[z]]],col=ifelse(k%%8 != 0, mycols[k%%8], mycols[8]))
                        mtext(paste("Mean: ",round(calib.avgs[[i]]$d18O.mean[j],2)," stdev: ",round(calib.avgs[[i]]$d18O.sd[j],2)),side=1,line=-2)

                        plot(spline.fits[[k]]$time,spline.fits[[k]]$d2H,
                            ylim=c(min(c(min(spline.fits[[k]]$d2H,na.rm=TRUE),min(calib.data$Delta_D_H[good.inds[[z]]],na.rm=TRUE))),
                            max(c(max(spline.fits[[k]]$d2H,na.rm=TRUE),max(calib.data$Delta_D_H[good.inds[[z]]],na.rm=TRUE)))))
                        points(calib.data$EPOCH_TIME[good.inds[[z]]],
                            calib.data$Delta_D_H[good.inds[[z]]],col=ifelse(k%%8 != 0, mycols[k%%8], mycols[8]))
                        mtext(paste("Mean: ",round(calib.avgs[[i]]$d2H.mean[j],2)," stdev: ",round(calib.avgs[[i]]$d2H.sd[j],2)),side=1,line=-2)
                    } else if (spline.has.data==TRUE & spline.has.good.data==FALSE) { # spline has data, but fails filter tests.
                        plot(spline.fits[[k]]$time,spline.fits[[k]]$H2O,
                            ylim=c(min(c(min(spline.fits[[k]]$H2O,na.rm=TRUE),min(calib.data$H2O[good.inds[[z]]],na.rm=TRUE))),
                            max(c(max(spline.fits[[k]]$H2O,na.rm=TRUE),max(calib.data$H2O[good.inds[[z]]],na.rm=TRUE)))))
                        points(calib.data$EPOCH_TIME[good.inds[[z]]],
                            calib.data$H2O[good.inds[[z]]],col=ifelse(k%%8 != 0, mycols[k%%8], mycols[8]))
                         mtext("Peak fails test in calculate.standard.averages fn",side=1,line=-2)

                        plot(spline.fits[[k]]$time,spline.fits[[k]]$d18O,
                            ylim=c(min(c(min(spline.fits[[k]]$d18O,na.rm=TRUE),min(calib.data$Delta_18_16[good.inds[[z]]],na.rm=TRUE))),
                            max(c(max(spline.fits[[k]]$d18O,na.rm=TRUE),max(calib.data$Delta_18_16[good.inds[[z]]],na.rm=TRUE)))))
                        points(calib.data$EPOCH_TIME[good.inds[[z]]],
                            calib.data$Delta_18_16[good.inds[[z]]],col=ifelse(k%%8 != 0, mycols[k%%8], mycols[8]))
                        mtext("Peak fails test in calculate.standard.averages fn",side=1,line=-2)

                        plot(spline.fits[[k]]$time,spline.fits[[k]]$d2H,
                            ylim=c(min(c(min(spline.fits[[k]]$d2H,na.rm=TRUE),min(calib.data$Delta_D_H[good.inds[[z]]],na.rm=TRUE))),
                            max(c(max(spline.fits[[k]]$d2H,na.rm=TRUE),max(calib.data$Delta_D_H[good.inds[[z]]],na.rm=TRUE)))))
                        points(calib.data$EPOCH_TIME[good.inds[[z]]],
                            calib.data$Delta_D_H[good.inds[[z]]],col=ifelse(k%%8 != 0, mycols[k%%8], mycols[8]))
                        mtext("Peak fails test in calculate.standard.averages fn",side=1,line=-2)
                    } else { # if spline returns no good data!
                        plot(spline.fits[[k]]$time,spline.fits[[k]]$H2O)
                        mtext("No acceptable peak found!",side=1,line=-2)

                        plot(spline.fits[[k]]$time,spline.fits[[k]]$d18O)
                        mtext("No acceptable peak found!",side=1,line=-2)

                        plot(spline.fits[[k]]$time,spline.fits[[k]]$d2H)
                        mtext("No acceptable peak found!",side=1,line=-2)
                    }
                        title(paste("Spline selected points for : ",time.suffix),outer=TRUE)
                # increment z and j as necessary
                if (spline.has.data==TRUE) {z <- min(c(z+1,length(good.inds)))}
                if (spline.has.good.data==TRUE) {j <- min(c(j+1,nrow(calib.avgs[[i]])))}
                }
                dev.off()
            } else { print("skipping plot because data frame is empty...")}
        }
    }

    #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    # DIAGNOSTIC PLOT 4: plots mean/sd values for H2O, delta18O, and delta2H 
    # by month. 
    if (RUN_PLOTS) {
        if (length(calib.avgs) > 0) {
            if (!is.null(calib.avgs[[i]]) == TRUE & length(calib.avgs[[i]]) > 0) {
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


                plot(as.POSIXct(calib.avgs[[i]]$time.mean[inds[[1]]],origin="1970-01-01"),calib.avgs[[i]]$H2O.mean[inds[[1]]],col=mycols[1],
                    ylim=c(min(calib.avgs[[i]]$H2O.mean,na.rm=TRUE),max(calib.avgs[[i]]$H2O.mean,na.rm=TRUE)),
                    xlim=c(min(calib.avgs[[i]]$time.mean,na.rm=TRUE),max(calib.avgs[[i]]$time.mean,na.rm=TRUE)))
                points(as.POSIXct(calib.avgs[[i]]$time.mean[inds[[2]]],origin="1970-01-01"),calib.avgs[[i]]$H2O.mean[inds[[2]]],col=mycols[2])
                points(as.POSIXct(calib.avgs[[i]]$time.mean[inds[[3]]],origin="1970-01-01"),calib.avgs[[i]]$H2O.mean[inds[[3]]],col=mycols[3])
                points(as.POSIXct(calib.avgs[[i]]$time.mean[inds[[4]]],origin="1970-01-01"),calib.avgs[[i]]$H2O.mean[inds[[4]]],col=mycols[4])
                points(as.POSIXct(calib.avgs[[i]]$time.mean[inds[[5]]],origin="1970-01-01"),calib.avgs[[i]]$H2O.mean[inds[[5]]],col=mycols[5])
                points(as.POSIXct(calib.avgs[[i]]$time.mean[inds[[6]]],origin="1970-01-01"),calib.avgs[[i]]$H2O.mean[inds[[6]]],col=mycols[6])


                plot(as.POSIXct(calib.avgs[[i]]$time.mean[inds[[1]]],origin="1970-01-01"),calib.avgs[[i]]$H2O.sd[inds[[1]]],col=mycols[1],
                    ylim=c(min(calib.avgs[[i]]$H2O.sd,na.rm=TRUE),max(calib.avgs[[i]]$H2O.sd,na.rm=TRUE)),
                    xlim=c(min(calib.avgs[[i]]$time.mean,na.rm=TRUE),max(calib.avgs[[i]]$time.mean,na.rm=TRUE)))
                points(as.POSIXct(calib.avgs[[i]]$time.mean[inds[[2]]],origin="1970-01-01"),calib.avgs[[i]]$H2O.sd[inds[[2]]],col=mycols[2])
                points(as.POSIXct(calib.avgs[[i]]$time.mean[inds[[3]]],origin="1970-01-01"),calib.avgs[[i]]$H2O.sd[inds[[3]]],col=mycols[3])
                points(as.POSIXct(calib.avgs[[i]]$time.mean[inds[[4]]],origin="1970-01-01"),calib.avgs[[i]]$H2O.sd[inds[[4]]],col=mycols[4])
                points(as.POSIXct(calib.avgs[[i]]$time.mean[inds[[5]]],origin="1970-01-01"),calib.avgs[[i]]$H2O.sd[inds[[5]]],col=mycols[5])
                points(as.POSIXct(calib.avgs[[i]]$time.mean[inds[[6]]],origin="1970-01-01"),calib.avgs[[i]]$H2O.sd[inds[[6]]],col=mycols[6])

                plot(as.POSIXct(calib.avgs[[i]]$time.mean[inds[[1]]],origin="1970-01-01"),calib.avgs[[i]]$d18O.mean[inds[[1]]],col=mycols[1],
                    ylim=c(min(calib.avgs[[i]]$d18O.mean,na.rm=TRUE),max(calib.avgs[[i]]$d18O.mean,na.rm=TRUE)),
                    xlim=c(min(calib.avgs[[i]]$time.mean,na.rm=TRUE),max(calib.avgs[[i]]$time.mean,na.rm=TRUE)))
                points(as.POSIXct(calib.avgs[[i]]$time.mean[inds[[2]]],origin="1970-01-01"),calib.avgs[[i]]$d18O.mean[inds[[2]]],col=mycols[2])
                points(as.POSIXct(calib.avgs[[i]]$time.mean[inds[[3]]],origin="1970-01-01"),calib.avgs[[i]]$d18O.mean[inds[[3]]],col=mycols[3])
                points(as.POSIXct(calib.avgs[[i]]$time.mean[inds[[4]]],origin="1970-01-01"),calib.avgs[[i]]$d18O.mean[inds[[4]]],col=mycols[4])
                points(as.POSIXct(calib.avgs[[i]]$time.mean[inds[[5]]],origin="1970-01-01"),calib.avgs[[i]]$d18O.mean[inds[[5]]],col=mycols[5])
                points(as.POSIXct(calib.avgs[[i]]$time.mean[inds[[6]]],origin="1970-01-01"),calib.avgs[[i]]$d18O.mean[inds[[6]]],col=mycols[6])

                plot(as.POSIXct(calib.avgs[[i]]$time.mean[inds[[1]]],origin="1970-01-01"),calib.avgs[[i]]$d18O.sd[inds[[1]]],col=mycols[1],
                    ylim=c(min(calib.avgs[[i]]$d18O.sd,na.rm=TRUE),max(calib.avgs[[i]]$d18O.sd,na.rm=TRUE)),
                    xlim=c(min(calib.avgs[[i]]$time.mean,na.rm=TRUE),max(calib.avgs[[i]]$time.mean,na.rm=TRUE)))
                points(as.POSIXct(calib.avgs[[i]]$time.mean[inds[[2]]],origin="1970-01-01"),calib.avgs[[i]]$d18O.sd[inds[[2]]],col=mycols[2])
                points(as.POSIXct(calib.avgs[[i]]$time.mean[inds[[3]]],origin="1970-01-01"),calib.avgs[[i]]$d18O.sd[inds[[3]]],col=mycols[3])
                points(as.POSIXct(calib.avgs[[i]]$time.mean[inds[[4]]],origin="1970-01-01"),calib.avgs[[i]]$d18O.sd[inds[[4]]],col=mycols[4])
                points(as.POSIXct(calib.avgs[[i]]$time.mean[inds[[5]]],origin="1970-01-01"),calib.avgs[[i]]$d18O.sd[inds[[5]]],col=mycols[5])
                points(as.POSIXct(calib.avgs[[i]]$time.mean[inds[[6]]],origin="1970-01-01"),calib.avgs[[i]]$d18O.sd[inds[[6]]],col=mycols[6])

                plot(as.POSIXct(calib.avgs[[i]]$time.mean[inds[[1]]],origin="1970-01-01"),calib.avgs[[i]]$d2H.mean[inds[[1]]],col=mycols[1],
                    ylim=c(min(calib.avgs[[i]]$d2H.mean,na.rm=TRUE),max(calib.avgs[[i]]$d2H.mean,na.rm=TRUE)),
                    xlim=c(min(calib.avgs[[i]]$time.mean,na.rm=TRUE),max(calib.avgs[[i]]$time.mean,na.rm=TRUE)))
                points(as.POSIXct(calib.avgs[[i]]$time.mean[inds[[2]]],origin="1970-01-01"),calib.avgs[[i]]$d2H.mean[inds[[2]]],col=mycols[2])
                points(as.POSIXct(calib.avgs[[i]]$time.mean[inds[[3]]],origin="1970-01-01"),calib.avgs[[i]]$d2H.mean[inds[[3]]],col=mycols[3])
                points(as.POSIXct(calib.avgs[[i]]$time.mean[inds[[4]]],origin="1970-01-01"),calib.avgs[[i]]$d2H.mean[inds[[4]]],col=mycols[4])
                points(as.POSIXct(calib.avgs[[i]]$time.mean[inds[[5]]],origin="1970-01-01"),calib.avgs[[i]]$d2H.mean[inds[[5]]],col=mycols[5])
                points(as.POSIXct(calib.avgs[[i]]$time.mean[inds[[6]]],origin="1970-01-01"),calib.avgs[[i]]$d2H.mean[inds[[6]]],col=mycols[6])

                plot(as.POSIXct(calib.avgs[[i]]$time.mean[inds[[1]]],origin="1970-01-01"),calib.avgs[[i]]$d2H.sd[inds[[1]]],col=mycols[1],
                    ylim=c(min(calib.avgs[[i]]$d2H.sd,na.rm=TRUE),max(calib.avgs[[i]]$d2H.sd,na.rm=TRUE)),
                    xlim=c(min(calib.avgs[[i]]$time.mean,na.rm=TRUE),max(calib.avgs[[i]]$time.mean,na.rm=TRUE)))
                points(as.POSIXct(calib.avgs[[i]]$time.mean[inds[[2]]],origin="1970-01-01"),calib.avgs[[i]]$d2H.sd[inds[[2]]],col=mycols[2])
                points(as.POSIXct(calib.avgs[[i]]$time.mean[inds[[3]]],origin="1970-01-01"),calib.avgs[[i]]$d2H.sd[inds[[3]]],col=mycols[3])
                points(as.POSIXct(calib.avgs[[i]]$time.mean[inds[[4]]],origin="1970-01-01"),calib.avgs[[i]]$d2H.sd[inds[[4]]],col=mycols[4])
                points(as.POSIXct(calib.avgs[[i]]$time.mean[inds[[5]]],origin="1970-01-01"),calib.avgs[[i]]$d2H.sd[inds[[5]]],col=mycols[5])
                points(as.POSIXct(calib.avgs[[i]]$time.mean[inds[[6]]],origin="1970-01-01"),calib.avgs[[i]]$d2H.sd[inds[[6]]],col=mycols[6])

                title(paste("Monthly data averages: ",time.suffix),outer=TRUE)
                dev.off()
            } else { print("skipping plot because data frame is empty...")}
        }
    }
    #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
} # end looping through months.

######################################################################################################
# COMBINE ALL CALIBRATION DATA POINTS, 
# # filter averages, and assign a standard name and concentration to the table
# # filter = flag values that have an unusual H2O concentration
# # bind all elements of the list together
calib.averages.all <- do.call(rbind,calib.avgs)

ambient.bracket.all <- do.call(rbind,ambient.buffers)

# attach ambient data columns to the calibration period averages data frame
calib.averages.all2 <- cbind(calib.averages.all,ambient.bracket.all)

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
