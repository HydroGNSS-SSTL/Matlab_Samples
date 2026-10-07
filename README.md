# HydroGNSS L1 Data Extraction and Example Plotting
For further information on data structure for received satellite data, please see the documents at https://www.hydrognss.org/:
HydroGNSS PDGS Product Manual
HydroGNSS L1 Algorithm Theoretical Baseline Document (ATBD)

### Scripts and Functions:
 *ExtractDataExample.mlx*
 *readL1DDM.m*
 *readL1MetaData.m*
 *readL1Global.m*

## ExtractDataExample.mlx
Scripted example of data selection, filtering and plotting. Using the provided functions.

## readL1DDM.m
Reads the scaled L1 DDM file and generates a MATLAB Struct.

## readL1MetaData.m
Read data from a HydroGNSS L1 merged metadata with H5library (low level), generates a MATLAB Struct of L1 parameters.

## readL1Global.m
Reads the global parameters from an L1 NetCDF file, works with blackbody and direct file types. Generates a MATLAB Struct of parameters.