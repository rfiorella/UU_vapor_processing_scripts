# WBB-vapor-scripts
Scripts for postprocessing water vapor isotope data from the William B Browning Building, U Utah.
The same scripts will be used to process the data from Hidden Peak @ Snowbird.

This data package is currently designed for Picarro L2130 analyzers, though there are plans to expand this to other analyzers in the future.

### Usage:
Code for processing and plotting the data are divided between three folders: (1) scripts, which contains the executable scripts to process the data, (2) src, which contains the accessory functions to do the data processing as "source code", (c) plots, which contains a few basic scripts to visualize the data:

##### Brief description of files in the scripts directory:
_Raw_ picarro logfiles are at a data resolution of ~1 Hz and are written hourly.
_L0_ scripts concatenate these data files into daily files (starting/ending 0000 GMT), adds metadata,
alters the time variables, and excludes a few redundant variables without any changes to underlying
data resolution.
_L1_ concatenates the L0 data to monthly files, and adds a few additional data quality checks. First, data points where the temperature and pressure of the analyzer cavity deviate too far out of spec are removed. Second, a basic data sanity check enforces that only data points where: (a) humidity is positive, (b) delta values are above -1000, (c) delta values are not unreasonably high, defined conservatively as +50 for d18O and +400 for dD. The reason for these extremely permissive bounds is that this sieve is designed to remove only periods where the analyzer is not operating properly and recording values that are not physically possible. The 1175-HIDS2045 analyzer in the Browning Building has had some issues with the laser, which has had to be repaired - this sieve filters out these periods.
_L2_ will calibrate the data, and is in the process of being written.

##### Brief description of files in the plots directory:
_L1 plotting_ makes plots of the monthly-scale data.