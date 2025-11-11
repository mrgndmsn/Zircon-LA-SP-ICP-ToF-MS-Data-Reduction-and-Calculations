#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  2 00:27:46 2025

@author: morganadamson
"""

##This script 
    #concatenates each of the csvs exported from NuQuant ~10000 particles per shot initially
    #blanks the below LOD values for to reduce error in further calculations
    #Records the csv it came from in a column called CSV
    #Records the shot number it came from in a column called Shot Number, calculated from the last numerical digit in the CSV column
    #Records what standard is being ablatied in a column called Standard
    #records the fluence being used in a column named Fluence
    #Calculates Peak Width as an additional column

import os, re
import pandas as pd 
import numpy as np 


# ---- paths ----
folder  = '/Users/morganadamson/Desktop/RESEARCH /ICPMS/TOF ICPMS/TOF DATA/Particles/082725/ABLATION/91500 6J Spot 2'
output  = '91500 6J Spot 2 Combined'
output2 = '91500 6J Spot 2 Combined LOD' #
path    = os.path.join(folder, output  + '.csv')
path2   = os.path.join(folder, output2 + '.csv')

files = sorted([f for f in os.listdir(folder) if f.lower().endswith(".csv")])

folder_name = os.path.basename(folder)

standard_folder      = re.split(r"[ _\-]+", folder_name.strip())[0]
ablation_type_folder = "dusting" if "dust" in folder_name.lower() else "ablation"
m_folder             = re.search(r"(\d+(?:\.\d+)?)\s*[Jj]\b", folder_name)
folder_fluence       = float(m_folder.group(1)) if m_folder else (0.6 if ablation_type_folder == "dusting" else np.nan)

# normalize shots 1–100 based on first name ending with ..01
shot_offset = None
for f0 in files:
    name0 = os.path.splitext(f0)[0]
    ints0 = re.findall(r"\d+", name0)
    if not ints0:
        continue
    n0 = int(ints0[-1])
    if n0 % 100 == 1:
        shot_offset = n0 - 1
        break

dfs = []
for f in files:
    filepath = os.path.join(folder, f)
    header = pd.read_csv(filepath, skiprows=8, nrows=0).columns
    df = pd.read_csv(filepath, skiprows=14, names=header)
    df.columns = [c.strip() for c in df.columns]

    name = os.path.splitext(f)[0]
    
    df["input_file"] = f

    # shot number (normalize to 1–100)
    ints = re.findall(r"\d+", name)
    if ints:
        n = int(ints[-1])
        if shot_offset is not None:
            k = n - shot_offset
        else:
            k = n % 100
            if k == 0:
                k = 100
    else:
        k = np.nan
    df["shot_number"]   = k
    df["standard"]      = standard_folder
    df["ablation_type"] = ablation_type_folder
    df["fluence"]       = folder_fluence

    # peak width = Peak End - Peak Start (in samples)
    if "Peak End" in df.columns and "Peak Start" in df.columns:
        df["Peak End"]   = pd.to_numeric(df["Peak End"], errors="coerce")
        df["Peak Start"] = pd.to_numeric(df["Peak Start"], errors="coerce")
        df["peak_width"] = df["Peak End"] - df["Peak Start"]
    else:
        df["peak_width"] = np.nan

    dfs.append(df)
    
# first output (raw combined)
combined = pd.concat(dfs, ignore_index=True)

combined.drop(columns=["197U", "237U"], inplace=True, errors="ignore")

combined.to_csv(path, index=False)
print("Wrote:", path, "shape:", combined.shape)

# second output: blank any cell that contains '<='
combined_lod = combined.copy()

# Only scan object (string) columns; set matching cells to ""
obj_cols = combined_lod.select_dtypes(include="object").columns
for c in obj_cols:
    s = combined_lod[c].astype(str)
    mask = s.str.contains(r"(<=)", na=False)
    combined_lod.loc[mask, c] = ""


combined_lod.to_csv(path2, index=False)
print("Wrote:", path2, "shape:", combined_lod.shape)




