#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  15 02:51:06 2025

@author: morganadamson
"""

##This script 
    # uses the density csv to combine all particles identified per hole into a shot by shot average
    # Calculates mean and 1 sd deviations of particle counts per particle and then the sum values per shot
    # outputs a summary for each hole 
    # creates files used to calculate LIFI, DHF, and particle size boxplots to compare against other ablation parameters and matrices


import os, glob
import numpy as np
import pandas as pd

# ---------------- USER SETTINGS ----------------
ROOT_DAY_DIR  = r'/Users/morganadamson/Desktop/RESEARCH /ICPMS/TOF ICPMS/TOF DATA/Particles/082725'
COMBINED_ROOT = os.path.join(ROOT_DAY_DIR, 'COMBINED')

# Column preferences / fallbacks
SHOT_CANDIDATES   = ['shot_number','Shot Number','Shot','Shot #','ShotID','Shot Id']
DIAM_CANDIDATES   = ['d_particle_um', 'd_particle_um_atoms']
PEAKSUM_CANDIDATE = 'Peak sum'
COMP_CANDIDATES   = ['Composition', 'Elements']

# Element → isotopes (intensity column names are "I_<iso>")
ELEM_ISOS = {
    'Hf': [176,177,178,179,180],
    'Pb': [206,207,208],
    'Th': [232],
    'U' : [234,235,238], 
}

# ---------------- UTILITIES ----------------
def _drop_unnamed(df):
    return df.loc[:, ~df.columns.str.startswith('Unnamed')]

def _find_first_col(df, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    return None

def _to_num(s):
    return pd.to_numeric(s, errors='coerce')

def _summ_stats(series):
    """Return (n, sum, mean, sd) for numeric series (ignore NaN)."""
    x = _to_num(series).dropna()
    n = int(x.size)
    if n == 0:
        return 0, 0.0, np.nan, np.nan
    s  = float(x.sum())
    mu = float(x.mean())
    sd = float(x.std(ddof=1)) if n > 1 else 0.0
    return n, s, mu, sd

def _sd_of_sum_from_sample(sd_sample, n):
    """If x_i iid with sample sd s, then sd(sum x_i) = s * sqrt(n)."""
    if n <= 1 or not np.isfinite(sd_sample):
        return 0.0
    return float(sd_sample) * np.sqrt(float(n))

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


def _element_intensity_series(df, element):
    """
    Return a per-row Series of element INTENSITY by summing available I_<iso> columns.
    Works whether there's 0, 1, or many isotope columns present.
    """
    cols = [f'I_{iso}' for iso in ELEM_ISOS[element] if f'I_{iso}' in df.columns]
    if not cols:
        return pd.Series(0.0, index=df.index)

    if len(cols) == 1:
        # single Series → safe for to_numeric
        s = pd.to_numeric(df[cols[0]], errors='coerce').fillna(0.0)
        return s

    # multiple columns → convert each col to numeric, then sum row-wise
    block = df[cols].apply(pd.to_numeric, errors='coerce').fillna(0.0)
    return block.sum(axis=1)

def _zircon_mask(df):
    """
    Classify zircon-present rows for TOTAL zircon/oxide splits.
    Uses (N_sites>0) OR (zircon_present True) OR (Composition contains 'zircon').
    """
    m = pd.Series(False, index=df.index)
    if 'N_sites' in df.columns:
        m = m | (_to_num(df['N_sites']) > 0)
    if 'zircon_present' in df.columns:
        m = m | df['zircon_present'].astype(str).str.lower().isin(['true','1','yes'])
    if 'Composition' in df.columns:
        m = m | df['Composition'].astype(str).str.contains('zircon', case=False, na=False)
    return m.fillna(False)

def _ratio_with_sigma(A, sA, B, sB):
    """
    R = A/B ;  σ_R = |R| * sqrt( (σA/A)^2 + (σB/B)^2 ), independent errors
    Returns (R, σ_R) with NaN if invalid.
    """
    if not np.isfinite(A) or not np.isfinite(B) or B <= 0:
        return (np.nan, np.nan)
    R = A / B
    termA = 0.0 if (A <= 0 or not np.isfinite(sA)) else (sA / A)**2
    termB = 0.0 if (B <= 0 or not np.isfinite(sB)) else (sB / B)**2
    sig = abs(R) * np.sqrt(termA + termB)
    return (float(R), float(sig))

def _find_shot_column(df):
    shot_col = _find_first_col(df, SHOT_CANDIDATES)
    if shot_col is None:
        return None, None
    raw = df[shot_col].astype(str).str.strip()
    digits = raw.str.extract(r'(\d+)')[0]
    shot_int = pd.to_numeric(digits, errors='coerce').astype('Int64')
    return shot_col, shot_int

# ---------------- CORE STEPS ----------------
def concatenate_folder_nonrecursive(folder):
    """
    Concatenate all CSVs directly in a parameter folder.
    Adds 'source_csv' column to track which input file each row came from.
    """
    run_label    = os.path.basename(folder)
    combined_name = f"{run_label} Combined LOD Threshold Density.csv"

    all_csvs = sorted(glob.glob(os.path.join(folder, "*.csv")))
    # avoid recursively reusing our combined output on re-runs
    all_csvs = [p for p in all_csvs if os.path.basename(p) != combined_name]
    if not all_csvs:
        raise FileNotFoundError(f"No CSVs found in {folder}")

    frames = []
    for p in all_csvs:
        try:
            d = pd.read_csv(p)
            d = _drop_unnamed(d)
            d["source_csv"] = os.path.basename(p)
            frames.append(d)
        except Exception as e:
            print(f"[skip] {os.path.basename(p)}  error={e}")

    if not frames:
        raise RuntimeError(f"No readable CSVs in {folder}")

    combined = pd.concat(frames, ignore_index=True, sort=False)
    combined = _drop_unnamed(combined)
    out_path = os.path.join(folder, combined_name)
    combined.to_csv(out_path, index=False)
    print(f"[combined] {out_path}  rows={len(combined)}  cols={len(combined.columns)}")
    return combined, out_path

def summarize_one_run(df_run, run_name, out_folder):
    """
    Build per-shot summary for a single run file (one spot).
    - Element intensities are computed from available I_<iso> columns.
    - TOTAL rows get zircon-vs-oxide splits + ratios with 1σ propagation.
    """
    df_run = _drop_unnamed(df_run).copy()

    shot_col, shot_int = _find_shot_column(df_run)
    if shot_col is None or shot_int.isna().all():
        print(f"[warn] {run_name}: no shot column found; writing empty summary.")
        return _write_empty_summary(run_name, out_folder)

    df_run['_shot_int'] = shot_int

    diam_col = _find_first_col(df_run, DIAM_CANDIDATES)
    peak_col = PEAKSUM_CANDIDATE if PEAKSUM_CANDIDATE in df_run.columns else None
    comp_col = _find_first_col(df_run, COMP_CANDIDATES)

    if diam_col: df_run[diam_col] = _to_num(df_run[diam_col])
    if peak_col: df_run[peak_col] = _to_num(df_run[peak_col])

    rows = []

    # Pre-compute per-row element intensity series (so we don't recompute inside loops)
    elem_series = {el: _element_intensity_series(df_run, el) for el in ELEM_ISOS}
    # Attach them (temporary) to sub-dataframes via columns
    for el, s in elem_series.items():
        df_run[f'__I_{el}__'] = s

    for shot_key, g_shot in df_run.groupby('_shot_int', dropna=False):
        if pd.isna(shot_key) or g_shot.empty:
            continue
        shot_out = int(shot_key)

        def compute_row(subdf, comp_label, enable_zr_ox=False):
            # ---- per-element INTENSITY stats on this subset ----
            per_el = {}
            for el in ['Hf','Th','Pb','U']:
                col = f'__I_{el}__'
                if col in subdf.columns:
                    n, s, mu, sd = _summ_stats(subdf[col])
                    sd_sum = _sd_of_sum_from_sample(sd, n)
                else:
                    n, s, mu, sd, sd_sum = (0, 0.0, np.nan, np.nan, np.nan)
                per_el[el] = dict(n=n, s=s, mu=mu, sd=sd, sd_sum=sd_sum)

            # PeakSum total and diameter stats
            if peak_col:
                _, s_peak, _, _ = _summ_stats(subdf[peak_col])
            else:
                s_peak = 0.0
            if diam_col:
                _, _, mu_d, sd_d = _summ_stats(subdf[diam_col])
            else:
                mu_d, sd_d = (np.nan, np.nan)

            row = {
                "Shot": shot_out,
                "Composition": comp_label,
                "Number of Particles": int(len(subdf)),

                "Sum I Hf": per_el['Hf']['s'], "SD Sum I Hf": per_el['Hf']['sd_sum'],
                "Sum I Th": per_el['Th']['s'], "SD Sum I Th": per_el['Th']['sd_sum'],
                "Sum I Pb": per_el['Pb']['s'], "SD Sum I Pb": per_el['Pb']['sd_sum'],
                "Sum I U":  per_el['U']['s'],  "SD Sum I U":  per_el['U']['sd_sum'],

                "Avg I Hf": per_el['Hf']['mu'], "Stdev I Hf": per_el['Hf']['sd'],
                "Avg I Th": per_el['Th']['mu'], "Stdev I Th": per_el['Th']['sd'],
                "Avg I Pb": per_el['Pb']['mu'], "Stdev I Pb": per_el['Pb']['sd'],
                "Avg I U":  per_el['U']['mu'],  "Stdev I U":  per_el['U']['sd'],

                "Sum PeakSum": s_peak,
                "Avg Diameter": mu_d,
                "Stdev Diameter": sd_d,
                "Run": run_name,
            }

            # Zircon vs Oxide split (TOTAL only)
            if enable_zr_ox:
                zm = _zircon_mask(subdf)
                om = ~zm

                def sums_for(mask, tag):
                    res = {}
                    for el in ['Hf','Th','Pb','U']:
                        col = f'__I_{el}__'
                        n, s, mu, sd = _summ_stats(subdf.loc[mask, col]) if col in subdf.columns else (0,0.0,np.nan,np.nan)
                        sd_sum = _sd_of_sum_from_sample(sd, n)
                        res[f"{tag} N {el}"] = n
                        res[f"{tag} Sum I {el}"] = s
                        res[f"{tag} SD Sum I {el}"] = sd_sum
                    return res

                zd = sums_for(zm, "Zirc")
                od = sums_for(om, "Oxide")
                row.update(zd); row.update(od)

                # Ratios (Pb/U, Pb/Hf, U/Th, Pb/Th, U/Hf) with error propagation
                def add_ratios(tag):
                    pairs = [("Pb","U"), ("Pb","Hf"), ("U","Th"), ("Pb","Th"), ("U","Hf")]
                    for num, den in pairs:
                        A  = row.get(f"{tag} Sum I {num}", np.nan)
                        sA = row.get(f"{tag} SD Sum I {num}", np.nan)
                        B  = row.get(f"{tag} Sum I {den}", np.nan)
                        sB = row.get(f"{tag} SD Sum I {den}", np.nan)
                        R, Rsd = _ratio_with_sigma(A, sA, B, sB)
                        row[f"{tag} R {num}/{den}"]    = R
                        row[f"{tag} R {num}/{den} SD"] = Rsd

                add_ratios("Zirc")
                add_ratios("Oxide")

            return row

        # TOTAL row (with Zirc/Oxide splits + ratios)
        rows.append(compute_row(g_shot, "TOTAL", enable_zr_ox=True))

        # Per-composition rows (no Zirc/Oxide columns to keep schema clear)
        if comp_col:
            for comp, g_comp in g_shot.groupby(comp_col, dropna=False):
                comp_name = "none" if pd.isna(comp) else str(comp)
                rows.append(compute_row(g_comp, comp_name, enable_zr_ox=False))

    # Build & save
    out = pd.DataFrame(rows)

    # Ensure stable column order (missing columns will be created as NaN if needed)
    base_cols = [
        "Shot","Composition","Number of Particles",
        "Sum I Hf","SD Sum I Hf","Sum I Th","SD Sum I Th","Sum I Pb","SD Sum I Pb","Sum I U","SD Sum I U","Sum PeakSum",
        "Avg I Hf","Stdev I Hf","Avg I Th","Stdev I Th","Avg I Pb","Stdev I Pb","Avg I U","Stdev I U",
        "Avg Diameter","Stdev Diameter","Run",
        # Zirc/Oxide totals (TOTAL rows only; NaN otherwise)
        "Zirc N Hf","Zirc N Th","Zirc N Pb","Zirc N U",
        "Zirc Sum I Hf","Zirc SD Sum I Hf","Zirc Sum I Th","Zirc SD Sum I Th",
        "Zirc Sum I Pb","Zirc SD Sum I Pb","Zirc Sum I U","Zirc SD Sum I U",
        "Oxide N Hf","Oxide N Th","Oxide N Pb","Oxide N U",
        "Oxide Sum I Hf","Oxide SD Sum I Hf","Oxide Sum I Th","Oxide SD Sum I Th",
        "Oxide Sum I Pb","Oxide SD Sum I Pb","Oxide Sum I U","Oxide SD Sum I U",
        "Zirc R Pb/U","Zirc R Pb/U SD","Zirc R Pb/Hf","Zirc R Pb/Hf SD",
        "Zirc R U/Th","Zirc R U/Th SD","Zirc R Pb/Th","Zirc R Pb/Th SD","Zirc R U/Hf","Zirc R U/Hf SD",
        "Oxide R Pb/U","Oxide R Pb/U SD","Oxide R Pb/Hf","Oxide R Pb/Hf SD",
        "Oxide R U/Th","Oxide R U/Th SD","Oxide R Pb/Th","Oxide R Pb/Th SD","Oxide R U/Hf","Oxide R U/Hf SD",
    ]
    for c in base_cols:
        if c not in out.columns:
            out[c] = np.nan
    out = out[base_cols].sort_values(["Shot","Composition"], kind="mergesort")

    base = os.path.splitext(os.path.basename(run_name))[0]
    out_name = f"{base}__per_shot_summary.csv"
    out_path = os.path.join(out_folder, out_name)
    out.to_csv(out_path, index=False)
    print(f"[per-run summary] {out_path}  rows={len(out)}")
    return out_path

def _write_empty_summary(run_name, out_folder):
    cols = [
        "Shot","Composition","Number of Particles",
        "Sum I Hf","SD Sum I Hf","Sum I Th","SD Sum I Th","Sum I Pb","SD Sum I Pb","Sum I U","SD Sum I U","Sum PeakSum",
        "Avg I Hf","Stdev I Hf","Avg I Th","Stdev I Th","Avg I Pb","Stdev I Pb","Avg I U","Stdev I U",
        "Avg Diameter","Stdev Diameter","Run"
    ]
    out = pd.DataFrame(columns=cols)
    base = os.path.splitext(os.path.basename(run_name))[0]
    out_name = f"{base}__per_shot_summary.csv"
    out_path = os.path.join(out_folder, out_name)
    out.to_csv(out_path, index=False)
    print(f"[per-run summary] {out_path}  rows=0")
    return out_path

def average_parameter(folder, per_run_paths):
    """
    Parameter-averaged per-shot summary (across runs) keyed by (Shot, Composition).

    - For "Sum ..." and "Number of Particles": mean across runs + SD across runs
    - For "Avg ..." (intensities, diameter): pooled mean & pooled SD using within-run SDs and per-row counts
    - For TOTAL rows: compute zircon/oxide ratios from pooled sums with uncertainty
    """
    if not per_run_paths:
        return None

    dfs = []
    for p in per_run_paths:
        try:
            d = pd.read_csv(p)
            d = _drop_unnamed(d)
            dfs.append(d)
        except Exception as e:
            print(f"[skip avg] {os.path.basename(p)}  error={e}")
    if not dfs:
        return None

    big = pd.concat(dfs, ignore_index=True, sort=False)

    # Prepare numeric conversions
    def numcol(c):
        if c in big.columns:
            big[c] = _to_num(big[c])

    sum_cols = ["Sum I Hf","Sum I Th","Sum I Pb","Sum I U","Sum PeakSum","Number of Particles"]
    avg_cols = ["Avg I Hf","Avg I Th","Avg I Pb","Avg I U","Avg Diameter"]
    sd_cols  = ["Stdev I Hf","Stdev I Th","Stdev I Pb","Stdev I U","Stdev Diameter"]

    for c in sum_cols + avg_cols + sd_cols:
        numcol(c)

    # Pooled-sum (TOTAL-only) columns for zircon/oxide
    for tag in ["Zirc","Oxide"]:
        for el in ["Hf","Th","Pb","U"]:
            numcol(f"{tag} Sum I {el}")
            numcol(f"{tag} SD Sum I {el}")

    def ratio_sigma(A, sA, B, sB):
        return _ratio_with_sigma(A, sA, B, sB)

    rows = []
    for (shot, comp), g in big.groupby(["Shot","Composition"], dropna=False):
        row = {"Shot": int(shot), "Composition": comp, "Replicates": int(g.shape[0])}

        # Sums/counts: mean & SD across runs
        for c in sum_cols:
            if c in g.columns:
                vals = g[c].to_numpy(dtype=float)
                row[f"Mean {c}"] = float(np.nanmean(vals)) if np.isfinite(vals).any() else np.nan
                row[f"SD {c}"]   = float(np.nanstd(vals, ddof=1)) if np.isfinite(vals).sum() > 1 else np.nan
            else:
                row[f"Mean {c}"] = np.nan
                row[f"SD {c}"]   = np.nan

        # Averages: pooled mean & pooled SD using within-run SDs and per-row counts
        n_i = g["Number of Particles"].to_numpy(dtype=float) if "Number of Particles" in g.columns else np.full(len(g), np.nan)
        for i, a_col in enumerate(avg_cols):
            sd_col = sd_cols[i]
            if (a_col in g.columns) and (sd_col in g.columns):
                means = g[a_col].to_numpy(dtype=float)
                sds   = g[sd_col].to_numpy(dtype=float)
                mean_p, sd_p = _pooled_mean_sd(means, sds, n_i)
                row[f"Pooled {a_col}"] = mean_p
                row[f"Pooled {sd_col}"] = sd_p
            else:
                row[f"Pooled {a_col}"] = np.nan
                row[f"Pooled {sd_col}"] = np.nan

        # TOTAL-only: ratios from pooled sums (sum across runs; SD via quadrature)
        if str(comp).upper() == "TOTAL":
            for tag in ["Zirc","Oxide"]:
                S = {}; Sd = {}
                for el in ["Hf","Th","Pb","U"]:
                    cS  = f"{tag} Sum I {el}"
                    cSd = f"{tag} SD Sum I {el}"
                    if (cS in g.columns) and (cSd in g.columns):
                        vals  = g[cS].to_numpy(dtype=float)
                        svals = g[cSd].to_numpy(dtype=float)
                        S[el]  = float(np.nansum(vals)) if np.isfinite(vals).any() else np.nan
                        Sd[el] = float(np.sqrt(np.nansum(np.square(svals)))) if np.isfinite(svals).any() else np.nan
                    else:
                        S[el], Sd[el] = (np.nan, np.nan)

                pairs = [("Pb","U"), ("Pb","Hf"), ("U","Th"), ("Pb","Th"), ("U","Hf")]
                for num, den in pairs:
                    R, Rsd = ratio_sigma(S.get(num, np.nan), Sd.get(num, np.nan),
                                         S.get(den, np.nan), Sd.get(den, np.nan))
                    row[f"{tag} R {num}/{den}"]     = R
                    row[f"{tag} R {num}/{den} SD"] = Rsd

        rows.append(row)

    out_avg = pd.DataFrame(rows).sort_values(["Shot","Composition"], kind='mergesort')
    run_label = os.path.basename(folder)
    out_path  = os.path.join(folder, f"{run_label}__parameter_averaged_per_shot.csv")
    out_avg.to_csv(out_path, index=False)
    print(f"[parameter average] {out_path}  rows={len(out_avg)}")
    return out_path

# ---------------- MAIN ----------------
if __name__ == '__main__':
    if not os.path.isdir(COMBINED_ROOT):
        raise RuntimeError(f"COMBINED folder not found: {COMBINED_ROOT}")

    any_ok = False
    for entry in sorted(os.scandir(COMBINED_ROOT), key=lambda e: e.name):
        if not entry.is_dir():
            continue
        folder = entry.path
        print(f"\n=== Processing parameter folder: {folder} ===")
        try:
            # 1) Concatenate (non-recursive) → keeps 'source_csv' (file name) so we can summarize per-spot
            df_all, combined_path = concatenate_folder_nonrecursive(folder)

            # 2) Per-run summaries — one per distinct source_csv (Spot 1, Spot 2, …)
            source_list = df_all['source_csv'].dropna().unique().tolist()
            per_run_paths = []
            for run_name in source_list:
                sub = df_all.loc[df_all['source_csv'] == run_name].copy()
                per_run_paths.append(summarize_one_run(sub, run_name, folder))

            # index of included runs
            pd.DataFrame({"Included runs": source_list}).to_csv(
                os.path.join(folder, "parameter_runs_included.csv"), index=False
            )

            # 3) Averaged per-shot summary across runs (same parameter)
            average_parameter(folder, per_run_paths)

            any_ok = True
        except Exception as e:
            print(f"[error] {folder}: {e}")

    if not any_ok:
        print("No subfolders processed. Make sure COMBINED contains run folders with the two spot CSVs.")
