MRSWAT Forecast System
User Manual

Version 26Feb2026

Abstract:
The objective of this manual is to provide user instructions to use MRSWAT in a simple and replicable manner. The MRSWAT Forecast System (MFS) is a Windows batch script that should automatically run the MRSWAT executable along with preprocessing and postprocessing scripts and log a standardized saltwater intrusion forecast.

Summary Methodology:
STEP 1.  Set up a working environment.
Set up a folder to run MFS. The folder should include the following elements:
[1] MRSWAT_Launcher.bat – This is the Winsdows batch script that will run one iteration of MFS. 
[2] The Template folder. This folder contains the MRSWAT executable and the required static (i.e., generally non-varying between simulations) input files. MFS will copy these files when it generates a new simulation folder. This folder also contains a subfolder called ‘.venv’. This folder contains the python environment required by MFS python scripts. Outdated, the .venv folder can likely be ignored, see [4]
[3] The Code Storage Main folder. This folder contains the master set of python scripts used by MFS; they will be copied into new simulation folders as needed.
[4] python_env. Complete transportable python environment for all scripts in the MFS.

STEP 2. Run the MRSWAT_Launcher.bat script.
This script should run one iteration of the MRSWAT FS. This includes running the MRSWAT so that it provides predictions of the 9 ppt and 0.25 ppt isohalines within the Lower Mississippi River at 10 and 28 days into the future. The predicted values will be logged into an excel spreadsheet generated in the immediate working environment called “MRSWAT_Data.xlsx”. As the MFS uses a batch script (.bat), Windows Task Scheduler can be used to repeatedly run MFS at a predefined schedule; all results should be logged within the single MRSWAT_Data spreadsheet.





Detailed Description of Files:

Please see the in-file annotation for each code file, it should be detailed.

MRSWAT_Launcher.bat – This is the master batch code that runs one iteration of the MFS. Most of the MFS actions are programmed in individual python code that are called by this script; this script also runs the MRSWAT executable. This code controls organization of the simulation folders creating new folder for each new simulation. A simulation folder is organized by the folder name, which is named in the format MRSWAT_YYYYMMDD_##, where YYYYMMDD is the date in which the folder was created and, presumably, when the simulation was run; ## is an integer representing the order in which the simulation was run for a given date if more than one simulation was run on thatdate.

Create_BC.py – This code creates the river discharge and stage boundary conditions for each simulation, which are saved in the new simulation folder as Stage.txt and Discharge.txt. The boundary condition time-series should span 180 back from the current day (USGS observations) to 30 days in the future (interpreted NOAA predictions). 

Create_RunFile*.py – These [three] files create the main.inp input file for MRSWAT based on criteria. If there are no recent simulations in the working environment from which a hotstart file may be based, the input file is modified so that full 180-day historical time series is simulated after a cold start. If a recent simulation folder is available, a hotstart file is generated from that simulation’s prediction of salinity dynamics for the current day. If it is determined that the simulation predictions of saltwater intrusion have significantly drifted  from observations, a simple synthetic hotstart file is created for the current day that assumes the saltwater wedge toe is located where last measured.

Postprocess_Results.py – This file provides various alternative outputs from the MRSWAT output file once the simulation is complete. Possible outputs include figures, animations, and logging the forecasted salinity intrusion metrics to the MRSWAT_Data.xlsx file.

Observations_check.py – This file checks to see if there are any recent vessel-based CTD cast measurements of saltwater wedge toe position at  \\mvd\mvn\Data_EDR\02 Stream Gaging Section\Gaging Stations\Gages_RandomProjects\Saltwater Wedge\Field Data. If so, the observed 9ppt and 0.25 ppt upstream-most locations are added to the ‘Observations’ sheet of the MRSWAT_Data.xlsx spreadsheet.

 HS_QC.py – At the conclusion of a MRSWAT simulation that utilized a hotstart file based on the results of a previous simulation results file

River_miles.csv – Text file with information defining Latitude and Longitude for Mississippi mile markers. Originally downloaded from USACE online sources. This version has 11 columns with headings: OBJECTID, name, LONGITUDE1, LATITUDE1, MILE, RIVER_CODE, RIVER_NAME, RIVER_NUMB, SOURCE, x, y.

Detailed Description of MRSWAT Forecast System workflow:

When the MRSWAT_Launcher.bat script is executed, the following processes should occur (Figure 1).
[Step 1] A new simulation folder is created named after the current date. Required MRSWAT model files are copied in from the Template folder.
[Step 2] A Python script is run to generate the MRSWAT boundary conditions. The observational and forecasted values are downloaded from online sources.
[Step 3] A Python script prepares the model run/master definition file that controls simulation duration and hotstart/coldstart dependent if recent simulation results are available.
[Step 4] The MRSWAT executable is run with the new boundary conditions run file.
[Step 5] The MRSWAT results are post-processed. If the user requested figure and animation output they are saved in a output folder within the simulation folder. The calculated saltwater wedge toe and 0.25 near-surface isohaline are recorded in the MRSWAT_Data spreadsheet.
[Step 6] If the simulation was initiated with a hotstart file, the MRSWAT FS checks to see if recent observational toe measurements are available. If so, the toe position in the hotstart file is checked against the observed values. If significant drift (offset) is observed, a synthetic hotstart is generated with the toe located at the observed position. An alternative MRSWAT simulation is then executed with the new hotstart.
