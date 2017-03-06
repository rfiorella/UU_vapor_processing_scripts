# UU_vapor_processing_scripts

This repository contains scripts to process raw log file data from a Picarro L2130 analyzer and process them to a common calibrated format. Written by Rich Fiorella (@rfiorella), Gabriel Bowen (@bumbanian), and Galen Gorski. As of version 1.0, these scripts are designed to work on the .dat log files in the Picarro User Directory. We hope to alter these to work on HDF files (either as a transition from .dat files or in addition to), in accordance with the upcoming/anticipated NEON data streams.

**The location of this repository is likely to change once development settles down! Its final home will likely be in the SPATIAL Lab Group GitHub page**

While this currently works on L2130-i analyzers only, we are interested in broadening this to more analyzer platforms. Please contact Rich if you are interested in discussing this. This document will be updated as the code is updated (and also provides a few hints at current development priorities)

### Usage:
This repository contains three directories: (1) scripts, which contains the executable scripts to process the data, (2) functions, which contains functions called by the associated script files - these can be quite long and are removed from the script file to help with file readability, and (3) metadata_templates, which contains examples of the data that is expected to go into file headers at each stage of data processing. 

There is also a powerpoint presentation/pdf of this presentation in the repository that was initially presented at the ICWEBS workshop in February 2017. Ideally, this will be updated to reflect changes in the code, but failing this, there will be a timestamp on the first slide to indicate the last time these files were updated.

A brief description of what each file in the scripts directory does, and what it produces, is provided below. 

##### Brief description of files in the scripts directory:
1. The _L0_ script concatenate these data files into daily files (starting/ending 0000 GMT), adds metadata,
restructures the time variables to make them a little more useful, and excludes a few redundant variables. Separate _L0_ files are written for calibration and ambient data for each day - calibration data files retain the same data resolution as the log files, ambient data files are averaged to 1 minute resolution (reducing space taken by these files by ~70x without averaging out "resolvable" high frequency variability).


2. The _L1_ script concatenates the daily L0 files to monthly files, and adds a few additional data quality checks. First, data points where the temperature and pressure of the analyzer cavity deviate too far out of spec are removed. Second, a basic data sanity check enforces that only data points where: (a) humidity is positive, (b) delta values are above -1000, (c) delta values are not unreasonably high, defined conservatively as +50 for d18O and +400 for dD. The reason for these permissive bounds is that this sieve is designed to remove only periods where the analyzer is not operating properly and recording values that are not physically possible. These filters were developed based on output from the 1175-HIDS2045 analyzer in the Browning Building, which had some issues with the laser, which has had to be repaired. It is highly likely that expanding these scripts to other Picarro analyzers and analyzer systems will require these filters to be updated.

3. The _L2_ script processes the calibration data, and returns files detailing the calibration analyses and providing the linear regression coefficients necessary to process the ambient data to a humidity-corrected VSMOW scale. By far, this is the most involved data processing step, and the one that could vary the most between research groups. Current functionality is outlined below, but other methods could be proposed and tested by those who wish to contribute code or ideas (ideally through issue reporting on GitHub or through collaborative push requests) (NOTE: everything in quotation marks below represents a tuned parameter or design choice - therefore, your results may be sensitive to the choices of these values. _Therefore, the value used at each choice will be saved as metadata in the file header!_). First, looping through each month of data:
* rough pass at defining time breakpoints in the calibration dataset (e.g., the rows in the matrix/data frame that identify where one peak ends and the next begins) by using the EPOCH_TIME variable in the log files. (contained in the ID.calib.breakpoints function!)
* refine breakpoint identification by fitting a "flexible" smoothing spline to each candidate peak identified in the prior step. Smoothing splines are ideal because they have smoothly varying derivatives, unlike the raw (noisy!) data. This step helps identify peaks that are closely spaced in time, but have very different humidities (e.g., you have been using the Picarro SDM and measure standard waters at different injection rates, and these different injection rates haven't been picked up by the prior step) (contained in the ID.calib.breakpoints function!)
* fit a "stiff" smoothing spline to each peak identified in the ID.calib.breakpoints function (can be modified, but in general, a stiffer spline requires larger deviations for points to be excluded!)
* calculate the time derivative of the smoothing splines (the calculate.spline.derivatives function)
* find the indices in each peak corresponding to a "stable, acceptable" measurement period based on **spline** derivatives (the extract.stable.calib.indices function, default assumes derivative thresholds of (absolute value) 5.0 for H2O, 0.01 for d18O, and 0.1 for d2H - these can be changed by specifying function arguments H2O.thres, d18O.thres, and d2H.thres respectively! _Values used will be written into file header._)
* test to ensure that peaks have "enough" data (default is 35 points, which corresponds to ~30 seconds of data on the L2130i.)
* calculate means and statistics on data calibration peaks (this data frame forms the "CalibrationAverages" data file described for L2 below) - please be aware this function (calculate.standard.averages) is more complex than it sounds! There are two major QA/QC issues embedded in this function:
	* there is a filter included in this function that excludes peaks based on a variety of quality control measures! Defaults are hard-coded in the function currently, but will be moved to be function arguments so that they can easily be changed. Current thresholds are: (a) 2000 ppm < peak H2O.mean < 30,000 ppm; (b) H2O standard deviation < 1000 ppm; (c) d18O standard deviation < 0.5; (d) d2H standard deviation < 4; (e) peak cannot be longer than 4200 indices (1 hour of data at 1.16 Hz, or 70/min - Picarro L2130-i log frequency) - this particular filter is required for the SLC data because it appears there is a brief period in the 2013 data where the drying column may have been exhausted and discrete peaks shorter than an hour could not be found!; (f) peak must be longer than 35 indices.
	* there is an embedded algorithm that seeks to minimize the memory effect of each peak while maximizing the amount of data kept in the peak (this is an alternative to simply taking the last e.g., 2 minutes of a peak and only using that.) this algorithm estimates the trend within the peak, and sequentially throws out 10% of the data from the beginning of the peak until both (a) the slope of the 18O peak is < 0.1 permil/min; (b) the slope of the 2H peak is < 0.5 permil/min, or (c) 80% of the peak is removed.
* estimate the isotopic composition of ambient vapor measured immediately before and after each calibration period (get.ambient.deltas function)

After this process has been completed for each month of data:
* calibration averages have been stored as a list in R - at this point, the list is flattened into a single data frame
* likewise, ambient values before and after calibration have been stored as a list. this is also flattened into a single data frame
* concatenate ambient values data frame to the standard averages data frame.
* _apply correction accounting for apparent delta dependence on cavity humidity (NB outline more as wiki entry?)_ (apply.mixingratio.calibration.correction function)
* correct standard values for incomplete drying of background gas (definitely - if using a drierite canister; NB - calculate isotopic difference required to have measureable effect on dry gas canister if present at 3 ppm.) (apply.drygas.correction function)
* assign names and known isotopic values of standards (assign.standard.names.and.values function) - this requires the user to modify the assign.standard.names.and.values function to fit their standards!!!)
* perform linear regression on calibration data to determine how to correct data on the analyzer scale to VSMOW (correct.standards.to.VSMOW function) - currently, this function identifies how many distinct periods of ambient data measurement are bracketed by standard analyses and performs separate linear regressions for 18O and 2H on the data preceeding and following each ambient measurement period between known 18O and 2H values and analyzer measured mean values (does not currently include uncertainty in these values - future update will include this as a weighted linear regression using 1/variance as regression weights). 


4. The _L3_ script takes uncalibrated monthly _L1_ ambient data files and calibrates them to a humidity-corrected VSMOW scale using the _L2_ calibration outputs. Each row in the _L1_ data is corrected using analyzer-specific cavity humidity-delta corrections (_see same wiki article mentioned above?_) and corrected to VSMOW by identifying the row in the _L2_ CalibrationRegressionData file that corresponds to that ambient measurement period.

5. _Future development: L4-Lx? scripts that attach important environmental data to the data file - this may include: meteorological data, other trace gas concentrations (e.g., CO2), other isotopic variables (13C of CO2?). I likely will develop one that attaches data from MesoWest or the AmeriFlux network; integrating other data sources will likely be responsibility of end users._

##### Description of data hierarchies used - what is contained in the output at each data level?
_Raw_ picarro logfiles are at a data resolution of ~1 Hz and are written hourly.
_L0_ files are daily split into calibration and ambient data - calibration data remains at the same data frequency as analyzer log files, ambient data has been averaged to 1 minute averages to optimize space usage.
_L1_ files are concatenated L0 files to separate monthly files for calibration and ambient data.

_L2_ output files include two different files containing the calibration information: 
1. File names containing "CalibrationAverages" have information about each peak value: 
* average, standard deviation, and range of H2O, d18O, d2H values through the peak; 
* min, max, and mean time corresponding to the peak; 
* count of analyzer measurements contributing to inferred peak;
* estimates of the trend (slope of linear regression) and strength of trend (r2 of linear regression) for peak d2H and d18O (required to help remove significant memory effects from peaks - there is a threshold limit in the calibration routines!)
* estimates of the H2O, d18O, and d2H of ambient vapor immediately before and after the peak (uncalibrated, as measured by analyzer); 
* mixing ratio corrected d18O and d2H (requires estimates of analyzer dependence of delta measurements on cavity humidity); 
* mixing ratio and background corrected d18O and d2H (accounts for incomplete drying of ambient vapor by drierite column); 
* known standard name, d18O, and d2H (assigned based on user supplied information on standards used in analysis, and when they were used if they changed through the time series). 

2. File names containing "CalibrationRegressionData" include:
* start and end times of intervals of ambient data corresponding to each linear regression.
* slope, intercept, and r2 of oxygen regression based on calibration data immediately preceeding and following each ambient data measurement period.
* slope, intercept, and r2 of hydrogen regression based on calibration data immediately preceeding and following each ambient data measurement period.

_L3_ output are essentially the calibrated version of the _L1_ ambient files. They take the _L1_ ambient file and use the _L2_ CalibrationRegressionData file to adjust the _L1_ ambient data to a humidity-corrected VSMOW scale.
