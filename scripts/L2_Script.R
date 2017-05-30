# L2_combined_script.R
# richf, 27jan17

# combined script meant to calibrate Picarro vapor data - relies on the water vapor calibration functions
rm(list=ls())

# LOAD REQUIRED LIBRARIES
#==========================================================================
library(lubridate)
library(zoo)
library(RColorBrewer)

# LOAD ASSOCIATED R CODE
#==========================================================================
# load associated functions
source("../functions/L2_functions.R")

# load user specifications
source("../user/L2_user_specs.R")

#--------------------------------------------------------------------------
# Begin L2 Data processing
#--------------------------------------------------------------------------

print("==================================================================")
print(paste("Starting L2 Processing : ",now()))
print("==================================================================")

# define the period given by start, end dates
period <- interval(start.date,end.date)

# calculate the number of months represented by start,end dates
# (script is designed to work on complete months) 
nmonths <- period %/% months(1)

#---------------------------------------------------
# Make lists of relevant *CALIBRATION* data files
#---------------------------------------------------
# make an empty list to hold all raw file names
raw.file.list <- list()

# list of all raw files
raw.file.list <- list.files(path=path.to.L1.data,pattern="CalibData",full.names=TRUE,recursive=TRUE)

# find month for each file
raw.file.dates <- extract.date.from.L1.files(raw.file.list,dbg.level=debug) 

# subset list based on period of interest
subset.list <- raw.file.list[raw.file.dates %within% period]

# Make lists of relevant *AMBIENT* data files
#---------------------------------------------------
# make an empty list to hold all raw file names
amb.file.list <- list()

# list of all raw files
amb.file.list <- list.files(path=path.to.L1.data,pattern="AmbientData",full.names=TRUE,recursive=TRUE)

# find month for each file
amb.file.dates <- extract.date.from.L1.files(amb.file.list,dbg.level=debug) 

# subset list based on period of interest
amb.subset <- amb.file.list[amb.file.dates %within% period]

# preallocate variables
calib.avgs.nomem <- vector("list",nmonths)
calib.avgs <- vector("list",nmonths)
ambient.buffers <- vector("list",nmonths)

# initiate log - need to keep track of how much data has been removed.
log.yyyy <- vector("numeric",nmonths)
log.mm   <- vector("numeric",nmonths)
log.nosplines <- vector("numeric",nmonths)
log.noderivsp <- vector("numeric",nmonths)
log.nogoodind <- vector("numeric",nmonths)
log.calibavgs <- vector("numeric",nmonths)

# start file counter
fcount <- 1
gcount <- 1 # perhaps not necessary, but coding defensively here in a rush...

#---------------------------------------------------------------------
# SET UP METADATA HEADERS
# Get metadata from an L1 file in here to pass along to L2 header...
tmp <- readLines(raw.file.list[[1]]) # get first calibration file.

md.frame <- tmp[grep("#",tmp)] # pull out lines starting with comment character...

#`````````````````````````````````````````````````````````````````````
# SET UP PLOTTING RESOURCES
#`````````````````````````````````````````````````````````````````````

if (RUN_PLOTS==TRUE) { # if we are running the diagnostic plots,
    # create a directory in plotting diag plot outpud directory that captures
    # script runtime.

    plot.path <- paste(plot.path,toString(base::date()),"/",sep="") # append current time to plot path
    # note: call to base::date is required as lubridate is loaded!!!

    # create output directory file path.
    dir.create(plot.path,recursive=TRUE)
}

#-------------------------------------------------------------------------------
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
    file.month <- month(tmp[[1]][7]) # index 7 of resulting list has the date info with current directory structure
    file.year <- year(tmp[[1]][7])   # index 7 of resulting list has the date info with current directory structure

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
	breaks <- ID.calib.breakpoints(calib.data,dbg.level=debug)

	# # calculate initial set of splines as a check on breaks
	spline.fits <- fit.calibration.splines(calib.data,breaks,dbg.level=debug)
	
    #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    # DIAGNOSTIC PLOT 1: Plot each peak period identified by ID.calib.breakpoints
    # and the splines that are fit to each peak period.

    if (RUN_PLOTS) {
        # make plots of the identified peaks with associated splines.
        print("Making plot showing calibration data and calculated splines...")
        pdf(paste(plot.path,"calib_raw_peakssplines_",time.suffix,".pdf",collapse="",sep=""),width=11,height=8)
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
 	spline.derivatives <- calculate.spline.derivatives(spline.fits,breaks,dbg.level=debug)

    #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    # DIAGNOSTIC PLOT 2: plot the first derivative of each identified peak period
    # to help identify suitable thresholds for identifying stable portions of peaks.

    if (RUN_PLOTS) {
        print("Making plot of spline derivatives...")
        # make plots of the spline derivatives.        
        pdf(paste(plot.path,"raw_data_splinederivs_",time.suffix,".pdf",collapse="",sep=""),width=11,height=8)
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
	good.inds <- extract.stable.calib.indices(spline.derivatives,dbg.level=debug)

    # quick filter here - should there be a minimum number of "good" points
    # required to define an analysis as a "standard" and retain the averages?
    # this can be done in ~3-4 lines, so no need for a function here.
    for (j in 1:length(good.inds)) {
        if (length(good.inds[[j]])<35) {
            good.inds[[j]] <- rep(NA,length(good.inds[[j]]))
        }
    }


	# based on the good points identified, calculate averages for each plateau
	calib.avgs[[i]] <- calculate.standard.averages(calib.data,good.inds,memory.filter=FALSE,dbg.level=debug)

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

    ambient.buffers[[i]] <-  get.ambient.deltas(calib.avgs[[i]],amb.data,dbg.level=debug)

    #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    # DIAGNOSTIC PLOT 3: plot the points selected by the spline thresholds 
    # for each peak. raw data are shown as points, splines are shown as lines.

    if (RUN_PLOTS) {
        print("Running diagnostic plot 3...")
        # plot the data after identifying plateaus using splines
        pdf(paste(plot.path,"spline_selected_points_",time.suffix,".pdf",
            collapse="",sep=""),width=11,height=8)
        par(mfrow=c(3,1),mar=c(4,4,1,1),oma=c(0,0,3,0))
        mycols <- rainbow(8)
        if (length(calib.avgs) > 0) { # i : months
            if (!is.null(calib.avgs[[i]]) == TRUE & length(calib.avgs[[i]]) > 0) {
                # start separate counter
                yy <- 1
                zz <- 1
                for (xx in 1:length(spline.fits)) { # k loops through spline.fits in month i
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

                    tmpmin <- min(spline.fits[[xx]]$time)
                    tmpmax <- max(spline.fits[[xx]]$time)
                    # screen for NA in gimean...
                    gimean <- ifelse(all(is.na(good.inds[[zz]])),NA,
                        mean(calib.data$EPOCH_TIME[good.inds[[zz]]],na.rm=TRUE))
                    camean <- calib.avgs[[i]]$time.mean[yy]

                    # print(paste(xx,yy,zz))
                    # print(paste(tmpmin,tmpmax,calib.avgs[[i]]$time.mean[yy]))
                    # print(paste(tmpmin,tmpmax,gimean))

                    # first test: was a peak found?
                    #print(paste(j,z,k,tmpmin,tmpmax,gimean,camean))

                    if (!is.na(gimean) & gimean < tmpmax & gimean >= tmpmin) {
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
                        plot(spline.fits[[xx]]$time,spline.fits[[xx]]$H2O,
                            ylim=c(min(c(min(spline.fits[[xx]]$H2O,na.rm=TRUE),min(calib.data$H2O[good.inds[[zz]]],na.rm=TRUE))),
                            max(c(max(spline.fits[[xx]]$H2O,na.rm=TRUE),max(calib.data$H2O[good.inds[[zz]]],na.rm=TRUE)))))
                        points(calib.data$EPOCH_TIME[good.inds[[zz]]],
                            calib.data$H2O[good.inds[[zz]]],col=ifelse(xx%%8 != 0, mycols[xx%%8], mycols[8]))
                         mtext(paste("Mean: ",round(calib.avgs[[i]]$H2O.mean[yy],2),
                             " stdev: ",round(calib.avgs[[i]]$H2O.sd[yy],2)),side=1,line=-2)

                        plot(spline.fits[[xx]]$time,spline.fits[[xx]]$d18O,
                            ylim=c(min(c(min(spline.fits[[xx]]$d18O,na.rm=TRUE),min(calib.data$Delta_18_16[good.inds[[zz]]],na.rm=TRUE))),
                            max(c(max(spline.fits[[xx]]$d18O,na.rm=TRUE),max(calib.data$Delta_18_16[good.inds[[zz]]],na.rm=TRUE)))))
                        points(calib.data$EPOCH_TIME[good.inds[[zz]]],
                            calib.data$Delta_18_16[good.inds[[zz]]],col=ifelse(xx%%8 != 0, mycols[xx%%8], mycols[8]))
                        mtext(paste("Mean: ",round(calib.avgs[[i]]$d18O.mean[yy],2)," stdev: ",round(calib.avgs[[i]]$d18O.sd[yy],2)),side=1,line=-2)

                        plot(spline.fits[[xx]]$time,spline.fits[[xx]]$d2H,
                            ylim=c(min(c(min(spline.fits[[xx]]$d2H,na.rm=TRUE),min(calib.data$Delta_D_H[good.inds[[zz]]],na.rm=TRUE))),
                            max(c(max(spline.fits[[xx]]$d2H,na.rm=TRUE),max(calib.data$Delta_D_H[good.inds[[zz]]],na.rm=TRUE)))))
                        points(calib.data$EPOCH_TIME[good.inds[[zz]]],
                            calib.data$Delta_D_H[good.inds[[zz]]],col=ifelse(xx%%8 != 0, mycols[xx%%8], mycols[8]))
                        mtext(paste("Mean: ",round(calib.avgs[[i]]$d2H.mean[yy],2)," stdev: ",round(calib.avgs[[i]]$d2H.sd[yy],2)),side=1,line=-2)
                    
                    } else if (spline.has.data==TRUE & spline.has.good.data==FALSE) { # spline has data, but fails filter tests.
                        plot(spline.fits[[xx]]$time,spline.fits[[xx]]$H2O,
                            ylim=c(min(c(min(spline.fits[[xx]]$H2O,na.rm=TRUE),min(calib.data$H2O[good.inds[[zz]]],na.rm=TRUE))),
                            max(c(max(spline.fits[[xx]]$H2O,na.rm=TRUE),max(calib.data$H2O[good.inds[[zz]]],na.rm=TRUE)))))
                        points(calib.data$EPOCH_TIME[good.inds[[zz]]],
                            calib.data$H2O[good.inds[[zz]]],col=ifelse(xx%%8 != 0, mycols[xx%%8], mycols[8]))
                         mtext("Peak fails test in calculate.standard.averages fn",side=1,line=-2)

                        plot(spline.fits[[xx]]$time,spline.fits[[xx]]$d18O,
                            ylim=c(min(c(min(spline.fits[[xx]]$d18O,na.rm=TRUE),min(calib.data$Delta_18_16[good.inds[[zz]]],na.rm=TRUE))),
                            max(c(max(spline.fits[[xx]]$d18O,na.rm=TRUE),max(calib.data$Delta_18_16[good.inds[[zz]]],na.rm=TRUE)))))
                        points(calib.data$EPOCH_TIME[good.inds[[zz]]],
                            calib.data$Delta_18_16[good.inds[[zz]]],col=ifelse(xx%%8 != 0, mycols[xx%%8], mycols[8]))
                        mtext("Peak fails test in calculate.standard.averages fn",side=1,line=-2)

                        plot(spline.fits[[xx]]$time,spline.fits[[xx]]$d2H,
                            ylim=c(min(c(min(spline.fits[[xx]]$d2H,na.rm=TRUE),min(calib.data$Delta_D_H[good.inds[[zz]]],na.rm=TRUE))),
                            max(c(max(spline.fits[[xx]]$d2H,na.rm=TRUE),max(calib.data$Delta_D_H[good.inds[[zz]]],na.rm=TRUE)))))
                        points(calib.data$EPOCH_TIME[good.inds[[zz]]],
                            calib.data$Delta_D_H[good.inds[[zz]]],col=ifelse(xx%%8 != 0, mycols[xx%%8], mycols[8]))
                        mtext("Peak fails test in calculate.standard.averages fn",side=1,line=-2)
                    } else { # if spline returns no good data!
                        plot(spline.fits[[xx]]$time,spline.fits[[xx]]$H2O)
                        mtext("No acceptable peak found!",side=1,line=-2)

                        plot(spline.fits[[xx]]$time,spline.fits[[xx]]$d18O)
                        mtext("No acceptable peak found!",side=1,line=-2)

                        plot(spline.fits[[xx]]$time,spline.fits[[xx]]$d2H)
                        mtext("No acceptable peak found!",side=1,line=-2)
                    }
                        title(paste("Spline selected points for : ",time.suffix),outer=TRUE)
                # increment z and j as necessary
                if (spline.has.data==TRUE | (is.na(gimean))) {zz <- min(c(zz+1,length(good.inds)))}
                if (spline.has.good.data==TRUE) {yy <- min(c(yy+1,nrow(calib.avgs[[i]])))}

                }
                dev.off()
            } else { print("skipping plot because data frame is empty...")}
        }
    }

    #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    # DIAGNOSTIC PLOT 4: plots mean/sd values for H2O, delta18O, and delta2H 
    # by month. 
    if (RUN_PLOTS) {
        print("Running diagnostic plot 4...")
        if (length(calib.avgs) > 0) {
            if (!is.null(calib.avgs[[i]]) == TRUE & length(calib.avgs[[i]]) > 0) {
                # plot monthly average data
                pdf(paste(plot.path,"monavgs_",time.suffix,".pdf",collapse="",sep=""),width=11,height=8)
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

    # Add some logging data to help track what's going on in this code...
    # print(paste(i,nmonths,length(spline.fits)))
    # log.yyyy[i] <- file.year
    # log.mm[i] <- file.month
    # log.nosplines[i] <- nrow(spline.fits[[i]])
    # log.noderivsp[i] <- nrow(spline.derivatives[[i]])
    # log.nogoodind[i] <- length(good.inds)
    # log.calibavgs[i] <- ifelse(is.list(calib.avgs[[i]]) & length(calib.avgs[[i]])!=0,nrow(calib.avgs[[i]]),0)

    # clean up variables that don't need to be carried to next i
    # rm(spline.fits)
    # rm(spline.derivatives)
    # rm(good.inds)

} # end looping through months.

# create log data frame
# log.df <- data.frame("yyyy"=log.yyyy,"mm"=log.mm,"spline.count"=log.nosplines,
#     "spline.derivative.count"=log.noderivsp,"std.average.count"=log.calibavgs)

######################################################################################################
# COMBINE ALL CALIBRATION DATA POINTS 

# bind all elements of the calibration and ambient lists together
calib.averages.all <- do.call(rbind,calib.avgs)
ambient.bracket.all <- do.call(rbind,ambient.buffers)

# attach ambient data columns to the calibration period averages data frame
calib.averages.wamb <- cbind(calib.averages.all,ambient.bracket.all)

# correct for delta dependence on concentration
calib.averages.wamb.mrc <- apply.mixingratio.correction(calib.averages.wamb,fit.type,Oslope,Hslope,
    dbg.level=debug)

# correct standard values for vapor bleeding through drierite canister...
calib.averages.wamb.mrc.bgc <- apply.drygas.correction(calib.averages.wamb.mrc,H2O.bg,
    include.gypsum.fractionation,dbg.level=debug)

# assign standard names to each calibration period
calib.averages.wamb.mrc.bgc.wstds <- assign.standard.names.and.values(calib.averages.wamb.mrc.bgc)

# do some filtering! --------------------------------------
h2o.min <- which(calib.averages.wamb.mrc.bgc.wstds$H2O.mean > h2o.min.thres) 
h2o.max <- which(calib.averages.wamb.mrc.bgc.wstds$H2O.mean < h2o.max.thres) 
h2o.sd <- which(calib.averages.wamb.mrc.bgc.wstds$H2O.sd < h2o.sdev.thres) 
d18O.sd <- which(calib.averages.wamb.mrc.bgc.wstds$d18O.sd < d18O.sdev.thres) 
d2H.sd <- which(calib.averages.wamb.mrc.bgc.wstds$d2H.sd < d2H.sdev.thres) 
max.length <- which(calib.averages.wamb.mrc.bgc.wstds$ind.count <= max.length.thres)
min.length <- which(calib.averages.wamb.mrc.bgc.wstds$ind.count >= min.length.thres)

retain.inds <- Reduce(intersect,list(h2o.min,h2o.max,h2o.sd,d18O.sd,d2H.sd,max.length,min.length))

calib.averages.wamb.mrc.bgc.wstds.filtered <- calib.averages.wamb.mrc.bgc.wstds[retain.inds,]

# calculate the correction slopes/intercepts.
calibration.regressions <- correct.standards.to.VSMOW(calib.averages.wamb.mrc.bgc.wstds.filtered,dbg.level=debug)

# make diagnostic plots of the calibration periods.
if (TRUE==FALSE) { # replace me when you actually are ready to run this code...
    print("Running diagnostic plot 5...")
    # plot the data after identifying plateaus using splines
    pdf(paste(plot.path,"regression_data_",time.suffix,".pdf",
        collapse="",sep=""),width=11,height=8)
    dev.off()
}


#------------------------------------------------------------------------
# Write out calibration parameters into a data file.

print(paste("Saving calibration data in two separate files"))
print(paste("The first - save all the information about each identified calibration data point..."))

calib.dt.name <- paste(path.to.output.L2.data,output.file.prefix,"_CalibrationAverages_L2_",
    start.date,"_",end.date,".dat",sep="")

# attach metadata
attach.L2.Header(calib.dt.name,md.frame,dbg.level=debug)

# write out data table of calib.averages.wamb.mrc.bgc.wstds
write.table(calib.averages.wamb.mrc.bgc.wstds,file=calib.dt.name,
    sep=",",row.names=FALSE,append=TRUE)

# ok, now write out the data file that contains the regression parameters...
print(paste("The second - save the regression parameters for each period..."))

regress.dt.name <- paste(path.to.output.L2.data,output.file.prefix,"_CalibrationRegressionData_L2_",
    start.date,"_",end.date,".dat",sep="")

# attach metadata
attach.L2.Header(regress.dt.name,md.frame,dbg.level=debug)

# write out datatable.
write.table(calibration.regressions,file=regress.dt.name,sep=",",row.names=FALSE,append=TRUE)

# make a couple quick diagnostic plots of regression parameters
quartz()
par(mfrow=c(3,1),mar=c(4,4,0.5,0.5))
plot(calibration.regressions$start.time,calibration.regressions$O.slope)
plot(calibration.regressions$start.time,calibration.regressions$O.intercept,col="red")
plot(calibration.regressions$start.time,calibration.regressions$O.r2,col="blue")

quartz()
par(mfrow=c(3,1),mar=c(4,4,0.5,0.5))
plot(calibration.regressions$start.time,calibration.regressions$H.slope)
plot(calibration.regressions$start.time,calibration.regressions$H.intercept,col="green")
plot(calibration.regressions$start.time,calibration.regressions$H.r2,col="blue")
