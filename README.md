# Zircon-LA-SP-ICP-ToF-MS-Data-Reduction-and-Calculations
Scripts used for "A mechanistic study of matrix effects.." Including Python scripts used to filter particles and calculate particle diameters from counts.


Python Scripts to use on raw data after being exported from NuQuant. The data used in this study are included in the supplementary files. 
File structure needed:
- 1 folder per ablation pit that has a CSV for each shot as exported by NuQuants' Single Particle method. 

Code Application Order:
1. "Concatenate 1.py"
      - Before running, make sure that only the raw CSVs from NuQuant are included within the file. If additional CSVs are included from any other step, the code will try to concatenate them as well. 
      - Will concatenate each of the CSVs from NuQuant
      - Blank the cells that are below the LOD
      - Record the Shot number, reference zircon name, and fluence
      - Calculates Peak Width into a new column
- Input: Path to folder that has all of the CSVs of a single hole (shots 1-100)
- Output 1: ___Folder Name___ Combined.csv
       - Saves to the input path
       - Concatenated file without any calculations or data filtering
- Output 2: ___Folder Name___ Combined LOD.csv
       - Saves to the input path
       - Concatenated file with calculations and cell filtering
       - File to be used for the next step
2. "FalsePositive Determination.py"
       - Uses the blank data file, after the first concatenate step, to determine the minimum peak width and sum that is usable to filter out potential false positives
       - Calculates the rate of potential false positives allowed after filtering
       - Creates a figure showing the loose and strict thresholds to filter out 99% of the false positive transient signals
       - This script is only intended for the blank data, and information from this step (strict and loose threshold peak widths and sum) is used in step 3
- Folder: Path to the blank runs folder that has the folder with all of the CSVs for the blank 20 shots
- File: the folder path plus the file name + Combined LOD.csv (From step 1)
- Output Flagged: ___Folder Name___ Flagged.csv
- Output Threshold: ___Folder Name___ Thresholds.csv
3. "Density and Diameter1.py"
       - Uses the __ Combined LOD.csvs from each parameter, the strict and loose thresholds, and ablation pit diameter and depth per shot to calculate particle diameter using the composition-derived model of zircon versus oxide particles
       - Will 
- Input:
- Output Threshold:
- Output Density:
       - File to be used for the next step
5. "Downhole Evolution averages.py"
- Input:
- Output:
       - File to be used for the next step
7. "Concatenate2.py"
- Input:
- Output:
       - Each summary file can be used and plotted against each other to compare ablation parameters and PSCD with depth.
