#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  4 00:31:34 2025

@author: morganadamson
"""

##This script 
    #takes the 1000s of blank particles from the NuQuant CSVs that were concatenated into one CSV from Concatenate1
    #calcualtes the number of particles identified per shot without further filtering and saves it as variable FalseP1
    #determines the upper boundary of peak width that filters out 95% of the identified particles and reports that value as variable UpperPW
    #determines the lower boundary of peak width that filters out 50% of the identified particles and reports that value as variable LowerPW
    #determines the minimum Peak Sum count that when combined with LowerPW filters out 95% of the identified particles and reports that value as MinCountPS
    #Plots a graph Peak Width vs Peak Sum with each identified particle as a data point
        #Plots the UpperPW, LowerPW, and MinCountPS as lines overlain
        #Calculates the number of particles identified after filtering based on these thresholds as variable FalseP2


import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# ======= USER PATHS =======
folder = '/Users/morganadamson/Desktop/RESEARCH /ICPMS/TOF ICPMS/TOF DATA/Particles/082725/Plesovice 20 shots Blank'
file   = folder + '/Plesovice 20 shots Blank Combined LOD.csv'
out_flagged    = folder + '/Plesovice 20 shots Blank Flagged.csv'
out_thresholds = folder + '/Blank Calculated Thresholds.csv'

# ======= SETTINGS =======
# percentiles (update comments to match values!)
P_STRICT = 99.9   # strict width percentile (e.g., 95.0 or 99.0 or 99.9)
P_LOOSE  = 50.0    # loose width floor (median)
P_COUNTS = 99.9    # peak-sum percentile within [LowerPW, UpperPW)

# ======= LOAD =======
df = pd.read_csv(file)
df.columns = [c.strip() for c in df.columns]

df["peak_width"] = pd.to_numeric(df["peak_width"], errors='coerce')
ps_col = "Peak sum" 
df[ps_col] = pd.to_numeric(df[ps_col], errors='coerce')

w = df["peak_width"].to_numpy(dtype=float)
s = df[ps_col].to_numpy(dtype=float)

FalseP1 = int(np.isfinite(w).sum())
if FalseP1 == 0:
    raise ValueError("No valid peak_width values in the blank dataset.")

# --- thresholds from blank ---
UpperPW = float(np.nanpercentile(w, P_STRICT))
LowerPW = float(np.nanpercentile(w, P_LOOSE))

if not (LowerPW < UpperPW):
    wuniq = np.unique(w[np.isfinite(w)])
    idx = np.searchsorted(wuniq, UpperPW, side="left")
    LowerPW = float(wuniq[idx-1]) if idx > 0 else max(0.0, np.nextafter(UpperPW, -np.inf))

band_mask = (w >= LowerPW) & (w < UpperPW)
band_s = s[band_mask]
if np.isfinite(band_s).sum() >= 5:
    MinCountPS = float(np.nanpercentile(band_s, P_COUNTS))
else:
    MinCountPS = float(np.nanpercentile(s, P_COUNTS))

# --- two-lane rule ---
keep_strict = (w >= UpperPW)
keep_loose  = (w >= LowerPW) & (s >= MinCountPS)
keep_any    = keep_strict | keep_loose

FalseP2  = int(np.nansum(keep_any))
pass_rate = (FalseP2 / FalseP1) if FalseP1 > 0 else np.nan

# ======= SAVE CSVs =======
thr = pd.DataFrame([{
    "P_STRICT": P_STRICT,
    "P_LOOSE":  P_LOOSE,
    "P_COUNTS": P_COUNTS,
    "UpperPW": UpperPW,
    "LowerPW": LowerPW,
    "MinCountPS": MinCountPS,
    "FalseP1": FalseP1,
    "FalseP2": FalseP2,
    "PassRate_after_filters": pass_rate
}])
thr.to_csv(out_thresholds, index=False)

df2 = df.copy()
df2["FP_keep_strict"] = keep_strict
df2["FP_keep_loose"]  = keep_loose
df2["FP_keep_any"]    = keep_any
df2.to_csv(out_flagged, index=False)

print("Thresholds written ->", out_thresholds)
print("Flagged blank table ->", out_flagged)
print(f"FalseP1={FalseP1}  FalseP2={FalseP2}  pass_rate={pass_rate:.3f}")

# ======= PLOT =======
plt.figure(figsize=(6.5,4.5))
plt.scatter(w, s, s=25, alpha=1, label="False Positives", color = 'black')
# vertical strict line
plt.axvline(UpperPW, linewidth=1.5, linestyle="--",
            label=f"Strict Width (p{P_STRICT})", color='r')
plt.text(UpperPW, plt.ylim()[1]*0.9, f"x = {UpperPW:.3f}",
         rotation=90, color='black', ha='right', va='top', fontsize=12)

# vertical loose line
plt.axvline(LowerPW, linewidth=1.5, linestyle=":",
            label=f"Loose Width (p{P_LOOSE})", color='b')
plt.text(LowerPW, plt.ylim()[1]*0.9, f"x = {LowerPW:.3f}",
         rotation=90, color='black', ha='right', va='top', fontsize=12)

# horizontal line
plt.axhline(MinCountPS, linewidth=1.5, linestyle="--",
            label=f"Min Peak Sum (p{P_COUNTS})", color='r')
plt.text(plt.xlim()[1]*0.95, MinCountPS, f"y = {MinCountPS:.2f}",
         color='black', ha='right', va='bottom', fontsize=12)

plt.xlabel("Peak Width (ms)", fontsize=12)  
plt.ylabel("Peak Sum (counts)", fontsize=12)
plt.title("False Positives: Peak Width vs Peak Sum", fontsize=14)

handles, labels = plt.gca().get_legend_handles_labels()

handles.append(Line2D([], [], linestyle='None'))
labels.append(f"Pass rate after filters: {pass_rate:.2%}  ({FalseP2}/{FalseP1})")

plt.legend(handles, labels, loc='best', frameon=True)

plt.rcParams['font.size'] = 12

plt.xlim(0.14, 0.32)
plt.ylim(0, 18)

plt.legend()
plt.tight_layout()


plot_png = os.path.join(folder, "FalsePositives.png")
plot_pdf = os.path.join(folder, "FalsePositives.pdf")
plot_svg = os.path.join(folder, "FalsePositives.svg")

plt.savefig(plot_png, dpi=1000)
plt.savefig(plot_pdf)   
plt.savefig(plot_svg)   

print("saved:", plot_png, plot_pdf, plot_svg)

plt.show()


