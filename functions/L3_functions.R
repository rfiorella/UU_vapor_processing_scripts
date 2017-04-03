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

