# Weak-radio-signal-data-processing
This repository contains the Bachelors thesis and algorithms of Weak radioastronomical object observation data processing - calibration, filtration and result analysis. 


In /docs/ the Latex source code for document is found along side of higher resolution photos.
While the Bachelors degree realization resulted in a lot of different scripts, only the main ones are chosen for display in /scripts/ as they perform functions described in thesis. Primary scripts MPI_calibr_v1.py and MPI_calibr_v1.py are fully commented.

**Scripts do actions, as follows:**
ACU_parser.py - Takes provided observation log file and directory of ACU logs to calculate the difference between log files.

doppler_reshift.py - Takes a directory of log files and result files from calibration and performs doppler compensation with the reshifting function described at the end of chapter 1.4.

doppler_v2.py - Takes a directory of log files and result files from calibration and performs doppler compensation using the old method of using vlsr.py script in helper directory.

MPI_calibr_v1.py - First iteration of MPI algorithms for data processing, takes log file and data directory.

MPI_calibr_v2.py - Final iteration of MPI algorithms for data processing, takes log file and data directory.

plot_visibility.py - Takes start and end dates, along side of object name to plot object visibility in that time period. Uses get_ra_dec_horizons.py helper.


**scripts/helper/ codes** are helper functions for main codes or not significant enough to be put in "main code" directory:

calc_measure_info.py - Reads through all of log files in provided directory to calculate the total observation time.

find_element.py - Used when parsing ACU log files to find location of certain headers.

get_brightness.py - Given start and end time, along side of object name connects to HORIZONS and fetches brightness values for object in that time period

get_ra_dec_horizons.py - Used by other scripts to fetch data from HORIZONS.

vlsr.py - Helper script for velocity calculation for objects outside solar system. Not developed as a part of this Bachelors degree, but rather is used by doppler_v2.py to calculate velocity.



As an addition, the log file for atlas-r-f1666-ir-1 observation is added for display purposes, since it is one of the observations analyzed more in the document.

