#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  16 10:28:02 2025

@author: morganadamson
"""

##This script 
    # Creates a pooled mean and 1SD file for each of the headers calculating the average composition of each particle both total, for zircon and oxide particles 
    # uses the file from Downhole Evolution averages1.py ___Run_Label___ parameter_averaged_per_shot.csv to average across all 1-100 shots from that averaged file across all ablation parameters
    

import os
import numpy as np
import pandas as pd

# ---------------- USER SETTINGS ----------------
ROOT_DAY_DIR  = r'/Users/morganadamson/Desktop/RESEARCH /ICPMS/TOF ICPMS/TOF DATA/Particles/082725'
COMBINED_ROOT = os.path.join(ROOT_DAY_DIR, 'COMBINED')

# Filenames used by the previous script
PARAM_AVG_NAME = "__parameter_averaged_per_shot.csv"

# Columns we expect from the parameter_averaged_per_shot stage
SUM_COLS = ["Mean Sum I Hf","Mean Sum I Th","Mean Sum I Pb","Mean Sum I U",
            "Mean Sum PeakSum","Mean Number of Particles"]
SUM_SD_COLS = ["SD Sum I Hf","SD Sum I Th","SD Sum I Pb","SD Sum I U",
               "SD Sum PeakSum","SD Number of Particles"]  # last may be missing/NaN

AVG_COLS = ["Pooled Avg I Hf","Pooled Avg I Th","Pooled Avg I Pb","Pooled Avg I U","Pooled Avg Diameter"]
AVG_SD_COLS = ["Pooled Stdev I Hf","Pooled Stdev I Th","Pooled Stdev I Pb","Pooled Stdev I U","Pooled Stdev Diameter"]

# ---------------- Helpers ----------------
def _num(s): return pd.to_numeric(s, errors='coerce')

def _pooled_mean_sd(means, sds, ns):
    means = np.asarray(means, dtype=float)
    sds   = np.asarray(sds,   dtype=float)
    ns    = np.asarray(ns,    dtype=float)
    mask  = (ns > 0) & np.isfinite(means) & np.isfinite(sds)
    if mask.sum() == 0:
        return (np.nan, np.nan)
    m, s, n = means[mask], sds[mask], ns[mask]
    N = n.sum()

    if mask.sum() == 1:
        # exactly one contributing group → pass-through its mean and SD
        return (float(m[0]), float(s[0]))

    if N <= 1:
        # degenerate weights; fall back to n-weighted mean, SD undefined
        return (float(np.average(m, weights=n)), np.nan)

    mean_p = float(np.sum(n * m) / N)
    within  = np.sum((n - 1.0) * (s ** 2.0))
    between = np.sum(n * ((m - mean_p) ** 2.0))
    var_p   = (within + between) / (N - 1.0)
    return (mean_p, float(np.sqrt(max(var_p, 0.0))))

def _safe_mean_sd(vec):
    x = _num(vec)
    m = float(np.nanmean(x)) if np.isfinite(x).any() else np.nan
    s = float(np.nanstd(x, ddof=1)) if np.isfinite(x).sum() > 1 else np.nan
    return m, s

def _pooled_total(vals_mean, vals_sd):
    """Sum of means and quadrature of SDs across shots."""
    vm = _num(vals_mean)
    vs = _num(vals_sd)
    S  = float(np.nansum(vm)) if np.isfinite(vm).any() else np.nan
    Sd = float(np.sqrt(np.nansum(np.square(vs)))) if np.isfinite(vs).any() else np.nan
    return S, Sd

# ---------------- Core ----------------
def process_parameter_folder(folder):
    run_label = os.path.basename(folder)
    in_path   = os.path.join(folder, f"{run_label}{PARAM_AVG_NAME}")
    if not os.path.isfile(in_path):
        print(f"[skip] {folder}: no parameter_averaged_per_shot file.")
        return None

    df = pd.read_csv(in_path)
    if df.empty:
        print(f"[skip] {folder}: parameter file is empty.")
        return None

    # Coerce key numeric columns once
    for c in SUM_COLS + SUM_SD_COLS + AVG_COLS + AVG_SD_COLS:
        if c in df.columns:
            df[c] = _num(df[c])

    rows_out = []
    for comp, g in df.groupby("Composition", dropna=False):
        comp_name = "TOTAL" if str(comp).upper() == "TOTAL" else ("none" if pd.isna(comp) else str(comp))
        shots_in_comp = g["Shot"].nunique() if "Shot" in g.columns else int(g.shape[0])
        reps_sum = int(_num(g.get("Replicates", np.nan)).dropna().sum()) if "Replicates" in g.columns else np.nan

        row = {"Composition": comp_name, "Shots Included": int(shots_in_comp), "Sum of Replicates": reps_sum}

        # ---- 1) Sums & counts: mean across shots and SD across shots ----
        for c in SUM_COLS:
            if c in g.columns:
                m, s = _safe_mean_sd(g[c])
                row[f"Shot-Mean {c}"] = m
                row[f"Shot-SD {c}"]   = s

        # ---- 2) Pooled totals across shots for sums/counts ----
        pairs = [
            ("Mean Sum I Hf", "SD Sum I Hf"),
            ("Mean Sum I Th", "SD Sum I Th"),
            ("Mean Sum I Pb", "SD Sum I Pb"),
            ("Mean Sum I U",  "SD Sum I U"),
            ("Mean Sum PeakSum", "SD Sum PeakSum"),
            ("Mean Number of Particles", "SD Number of Particles"),
        ]
        for mean_col, sd_col in pairs:
            if (mean_col in g.columns) and (sd_col in g.columns):
                S, Sd = _pooled_total(g[mean_col], g[sd_col])
                row[f"Pooled Total across shots: {mean_col}"] = S
                row[f"Pooled Total SD across shots: {mean_col}"] = Sd

        # ---- 3) Averages: pooled again across shots using each shot’s (mean, sd, n) ----
        n_for_pooled = g["Mean Number of Particles"] if "Mean Number of Particles" in g.columns else np.full(len(g), np.nan)
        for mean_col, sd_col in zip(AVG_COLS, AVG_SD_COLS):
            if (mean_col in g.columns) and (sd_col in g.columns):
                mean_p, sd_p = _pooled_mean_sd(g[mean_col], g[sd_col], n_for_pooled)
                row[f"Shot-Pooled {mean_col}"] = mean_p
                row[f"Shot-Pooled {sd_col}"]   = sd_p

        rows_out.append(row)

    out = pd.DataFrame(rows_out)

    # Column ordering (no ratio columns)
    lead_cols = ["Composition","Shots Included","Sum of Replicates"]
    sum_mean_cols = [c for c in out.columns if c.startswith("Shot-Mean ")]
    sum_sd_cols   = [c for c in out.columns if c.startswith("Shot-SD ")]
    pooled_total_cols    = [c for c in out.columns if c.startswith("Pooled Total across shots: ")]
    pooled_total_sd_cols = [c for c in out.columns if c.startswith("Pooled Total SD across shots: ")]
    shot_pooled_mean_cols = [c for c in out.columns if c.startswith("Shot-Pooled Pooled Avg")]
    shot_pooled_sd_cols   = [c for c in out.columns if c.startswith("Shot-Pooled Pooled Stdev")]

    ordered = (lead_cols + sum_mean_cols + sum_sd_cols +
               pooled_total_cols + pooled_total_sd_cols +
               shot_pooled_mean_cols + shot_pooled_sd_cols)
    tail = [c for c in out.columns if c not in ordered]
    out = out[ordered + tail]

    out_path = os.path.join(folder, f"{run_label}__hole_averaged.csv")
    out.to_csv(out_path, index=False)
    print(f"[hole-averaged] {out_path}  rows={len(out)}")
    return out_path

# ---------------- MAIN ----------------
if __name__ == "__main__":
    if not os.path.isdir(COMBINED_ROOT):
        raise RuntimeError(f"COMBINED folder not found: {COMBINED_ROOT}")

    any_ok = False
    for entry in sorted(os.scandir(COMBINED_ROOT), key=lambda e: e.name):
        if not entry.is_dir():
            continue
        folder = entry.path
        print(f"\n=== Hole averaging: {folder} ===")
        try:
            outp = process_parameter_folder(folder)
            any_ok = any_ok or (outp is not None)
        except Exception as e:
            print(f"[error] {folder}: {e}")

    if not any_ok:
        print("No subfolders produced hole-averaged outputs. Check that each contains __parameter_averaged_per_shot.csv")
