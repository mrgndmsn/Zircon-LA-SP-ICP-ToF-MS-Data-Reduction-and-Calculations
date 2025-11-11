# Zircon-LA-SP-ICP-ToF-MS-Data-Reduction-and-Calculations
Scripts used for "A mechanistic study of matrix effects in zircon U-Th/Pb geochronology using single-particle laser ablation time-of-flight mass spectrometry," Including Python scripts used to filter particles and calculate particle diameters from counts.


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
      - Will use isotope sensitivity and long-term U counts over the depth of the hole to root the counts per shot to the number of atoms
      - The number of atoms will be used to reconstruct the Molar mass of each particle using our two-path model into zircon or oxide particles using subsitution from the amount of atoms of each element within each identified particle
      - The Molar mass of each element within each particle, based on the particle composition, will be used to calculate the total mass from each particle
      - Using the calculated mass of each particle, zircon particles will use a zircon effective density, oxide particles will use their oxide's effective density, and the particle diameter is calculated
      - Calculates and records a new column for each of these calculated parameters per particle per shot
- Input: ___Folder Name___ Combined LOD.csv
- Output Threshold: ____Folder Name___ Combined LOD Threshold.csv
      - Filters suspected false positive particles that have too small peak widths combined with too low peak sum counts, as determined by FalsePositive Determination.py
- Output Density: ____Folder Name___ Combined LOD Threshold Density.csv
      - Main file calculated from sensitivity and density models to determine particle volume
      - File to be used for the next step
4. "Downhole Evolution averages.py"
      - The script will filter through all folders within a major day folder to average particles per shot per hole across multiple holes
      - Concatenates all CSVs in a parameter folder within a larger COMBINED folder within each day (it is looking for the ____Folder Name___ Combined LOD Threshold Density.csv that were created in the previous step, and for each run that you want to be included, it needs to be copied inside a smaller parameter folder within the COMBINED folder)
      - Will calculate the number of particles, sum element intensities etc, downhole, and average element intensities per particle downhole, each including  1 SD uncertainties
      - Will calculate Pb/Hf, U/Th, Pb/U, and U/Hf for each downhole and per particle
      - Will do this for an average particle and total volume of the hole across all runs, and the same for each of the 15 different identified particle types to better understand the size and intensities of each particle type per shot and per hole across all runs
- Input: Day folder, with a COMBINED folder within it that hosts a folder for each ablation parameter, which inside has the ____Folder Name___ Combined LOD Threshold Density.csv for each run of that one ablation parameter
- Output: Within each Day/COMBINED/Parameter folder,
  1. ___Source___ Per_shot_summary.csv
      - Concatenates all of the runs within that parameter with all of their particle values from the Combined LOF Threshold Densities
  2. ___Run_Label___ Combined LOD Threshold Density.csv
      - Averages for each shot per run create one summary file like this for each run
  3. ___Run_Label___ parameter_averaged_per_shot.csv
      - Averages across all runs into one 1-100 shot average file
      - File to be used for the next step
  4. parameter_runs_included.csv
      - tells which runs were averaged together to record what happened
5. "Concatenate2.py"
      - Will average across all shots into a representative value of particle size and composition for each ablation parameter and matrix
- Input: Day folder, with a COMBINED folder within it that hosts a folder for each ablation parameter, which inside has the ___Run_Label___ parameter_averaged_per_shot.csv for each ablation parameter 
- Output: ___Run_Label___ hole_averaged.csv
      - Each summary file for the average and total average counts per particle and per hole for a given ablation parameter and matrix
      - Can be used and plotted against eachother to compare ablation parameters and PSCD with depth

IN TOTAL: 
- particle number, intensities, and diameters are weighted averages per particle and per shot, and finally per hole for each ablation parameter and each matrix
- 1SD is propagated using pooled averaging in quadrature and includes internal uncertainty from the variance between particles generated per shot from a single run through to averaging and adding uncertainties to one another from multiple runs into a single mean value for total and single particle types
- values from ___Run_Label___ parameter_averaged_per_shot.csv's were used to plot and calculate downhole fractionation indexes
- values from ___Run_Label___ hole_averaged.csv were used to plot and calculate representative particle and total particle means, etc, and laser-induced fractionation indexes 
