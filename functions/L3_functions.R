# L3 functions - need to update the header here...

extract.date.from.L1.files <- function(file.list,dbg.level=0) {
  # print statement announcing the beginning of the function
  if (dbg.level>0) { 
    print("============================================")
    print("Start of extract.date.from.L1.files function") 
  }

  # suppress warnings in this function - as.numeric will create NAs by coercion,
  # and this will be reported as a warning to the screen every time. But since this
  # is precisely what we are hoping to do, warnings are unncessary here.
  oldw <- getOption("warn")
  options(warn = -1)
  # split file name into pieces using *underscores*
  file.name.string.pieces <- strsplit(file.list,split="_")
  # extract dates from resulting string pieces
  file.name.dates <- ymd(unlist(file.name.string.pieces))
  # remove NAs
  file.name.dates.woNAs <- file.name.dates[!is.na(file.name.dates)]
  # set warnings to old status
  options(warn = oldw)

  # print statement announcing the beginning of the function
  if (dbg.level>0) { 
    print("End of extract.date.from.L1.files function") 
    print("==========================================")
  }
  
  # return date list
  return(file.name.dates.woNAs)
}

#---------------------------------------------------------------------------
# attach.L3.Header

attach.L3.Header <- function(output_filename,metadata_dataframe,dbg.level=0) {
  # print statement indicating that this function is starting
  if (dbg.level>0) {
    print("==================================")
    print("starting attach.L3.Header function")
  }

  # initiate writing out to output file...
  sink(output_filename)

  # add L1 information to top of output file...
  cat("# BEGIN L3 HEADER \n")

  # note when L2 script was run
  cat(paste("# date.L3.scripts.run, ",base::date(),"\n"))

  # calibration file used...
  cat(paste("# Calibration file used, ",path.to.L2.calib.data,calib.file,"\n"))

  # close L3 header
  cat("# END L3 HEADER \n")

  # append L0/L1/L3 headers.
  # loop through existing headers and paste below L3 header.
  for (i in 1:length(metadata_dataframe)) {
    cat(paste(metadata_dataframe[i],"\n")) 
  }
  
  # turn off output
  sink()
  
  # print statment denoting the end of this function
  if (dbg.level>0) {
    print("ending attach.L3.Header function")
    print("================================")
  }
}






