#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  5 00:44:23 2025

@author: morganadamson
"""

##This script 
    #applies a user defined filter based on UpperPW, LowerPW, and MinCountPS
    #uses the long term Hf and U to create a sensitivity calibration curve so we can measure N from I
    #Adds a column based on the calibrated N for each isotope named N_** for each isotope
    #Adds a column for the summed counts of each element SI_** for each element
    #Adds a column for the summed atoms of each element SN_** for each element\
    #Adds a column for the calibrated Peak Sum for each particle named S_Peak Sum
    #adds a column based on the elements present in each particle in a variable named Composition
    ##For density calculation
        #Adds a column for the approximation of NZr based on the NHf and fHf and Nsites
        #Adds a column for the approximation of Munit based on stoichiometric proportions
        #Adds a column for the approximation of mzircon based on Munit
        #Adds a column for the approximation of mPbO2, mThO2, and mUO2 based on the molar mass M of each oxide
        #Adds a column for mp by adding together all of the calculated m
        #Adds a column for each element calculating the mass fractions w for each element
        #adds a column calculating peff based on the mass fractionatons and masses calculated and respective partiaal densities
        #adds a column approximating the particle dimeter based on peff and spherical shape
    #Exports new filtered and built upon CSV as its own thing named Density and diameter calcuated


import os
import numpy as np
import pandas as pd
from math import pi

# ---------------- USER SETTINGS ----------------
INPUT_CSV  = '/Users/morganadamson/Desktop/RESEARCH /ICPMS/TOF ICPMS/TOF DATA/Particles/082725/DUSTING/91500 6J Dusting Spot 2/91500 6J Dusting Spot 2 Combined LOD.csv'

# Outputs mirroring original behavior
OUTPUT_THRESHOLD = os.path.splitext(INPUT_CSV)[0] + ' Threshold.csv'
OUTPUT_DENSITY   = os.path.splitext(INPUT_CSV)[0] + ' Threshold Density.csv'

# Drop rows with Composition=="none" from the density file
DROP_NONE_FROM_DENSITY = True

# Filtering thresholds (tune from FalsePositiveDetermination script)
STRICT_PW = 0.2054400000015450
LOOSE_PW  = 0.1540800000002490
MIN_PS    = 12.8990250536495

# Geometry & densities
CRATER_DIAM_UM    = 50.0           # spot size (µm)
DEPTH_PER_SHOT_UM = 0.1        # set from your profilometry (µm/shot)
RHO = {
    'zircon': 4.65,   # g/cm3
    'PbO2'  : 9.38,
    'ThO2'  : 10.00,
    'UO2'   : 10.97
}

# Standards & constants
STD_PPM_U = {'Plesovice': 755.0, 'Plešovice': 755.0, '91500': 80.0}
NA = 6.02214076e23

# Natural isotope fractions (unitless)
f_Hf = {176:0.05206, 177:0.18606, 178:0.27297, 179:0.13629, 180:0.35100}
f_U  = {235:0.00720, 238:0.99280}
f_Pb = {206:0.241, 207:0.221, 208:0.524}  

# Accurate masses (amu ≈ g/mol)
M = {
    176:175.941, 177:176.943, 178:177.943, 179:178.945, 180:179.947,
    232:232.038, 235:235.044, 238:238.051,
    206:205.974, 207:206.976, 208:207.977
}
# Element molar masses
M_Zr=91.224; M_Hf=178.49; M_Si=28.0855; M_O=15.999
M_Pb=207.2;  M_Th=232.038; M_U=238.0289
M_PbO2 = M_Pb + 2*16.0
M_ThO2 = M_Th + 2*16.0
M_UO2  = M_U  + 2*16.0

# Hf fraction on the Zr site (Hf/(Zr+Hf)) to infer site count from Hf
F_HF_REF = 0.012   

# CSV column mapping (edit if your headers differ)
COL = {
    176:'176Hf', 177:'177Hf', 178:'(T) 178Hf', 179:'179Hf', 180:'(T) 180Hf',
    232:'(T) 232Th',
    235:'(T) 235U', 238:'(T) 238U',
    206:'(T) 206Pb', 207:'(T) 207Pb', 208:'(T) 208Pb'
}
SHOT_COL_CANDIDATES = ['Shot Number','shot_number','Shot','Shot #','ShotID','Shot Id']

# ---------------- HELPERS ----------------
def cm_from_um(x_um): return x_um*1e-4

def crater_mass_removed_g(n_shots, rho_g_cm3, D_um, depth_per_shot_um):
    r_cm = cm_from_um(D_um)/2
    depth_cm = cm_from_um(n_shots*depth_per_shot_um)
    vol_cm3 = pi*r_cm**2*depth_cm
    return rho_g_cm3*vol_cm3

def atoms_from_ppm(total_mass_g, ppm, molar_mass):
    m_el_g = total_mass_g*(ppm/1e6)
    moles  = m_el_g/molar_mass
    return moles*NA

def safe_sum(df, col):
    return pd.to_numeric(df.get(col, 0.0), errors='coerce').fillna(0.0).sum()

# ---------------- LOAD ----------------
df = pd.read_csv(INPUT_CSV)
df.columns = [c.strip() for c in df.columns]

# ---------------- FILTER ----------------
pw = pd.to_numeric(df['peak_width'], errors='coerce') if 'peak_width' in df else pd.Series(np.inf, index=df.index)
ps = pd.to_numeric(df.get('Peak sum', 0.0), errors='coerce').fillna(0.0)
mask = (pw >= STRICT_PW) | ((pw >= LOOSE_PW) & (ps >= MIN_PS))
df = df[mask].copy()
df.reset_index(drop=True, inplace=True)

df.to_csv(OUTPUT_THRESHOLD, index=False)
print(f"[OK] Saved filtered-only: {OUTPUT_THRESHOLD}")


assert 'standard' in df.columns, "Need a 'standard' column."
std_candidates = df['standard'].dropna().unique().tolist()
std_key = next((k for k in STD_PPM_U if any(k.lower() in str(s).lower() for s in std_candidates)), None)
assert std_key, f"No supported standard in {std_candidates}"
U_ppm = STD_PPM_U[std_key]
df_std = df[df['standard'].astype(str).str.lower().str.contains(std_key.lower())].copy()
assert len(df_std) > 0, "No rows matched the standard subset."

# Shot count for geometric anchor
SHOT_COL = next((c for c in SHOT_COL_CANDIDATES if c in df.columns), None)
if SHOT_COL is None:
    SHOT_COL='__oneshot__'; df[SHOT_COL]=1
N_SHOTS_STD = df_std[SHOT_COL].nunique()

# ---------------- BUILD RELATIVE SENSITIVITY ----------------
pts_m, pts_Srel = [], []

def add_isotope_family(df_std, nat_fracs, anchor_iso):
    if len(nat_fracs)==1:
        iso = list(nat_fracs.keys())[0]
        c = COL.get(iso)
        if c in df_std.columns:
            I = safe_sum(df_std, c)
            if I>0:
                pts_m.append(M[iso]); pts_Srel.append(1.0)
        return
    anchor_f = nat_fracs[anchor_iso]
    anchor_c = COL.get(anchor_iso)
    if anchor_c not in df_std.columns:
        return
    I_anchor = safe_sum(df_std, anchor_c)
    if I_anchor<=0:
        return
    for iso, f in nat_fracs.items():
        c = COL.get(iso)
        if c in df_std.columns:
            I = safe_sum(df_std, c)
            if I<=0: 
                continue
            Srel = (I/I_anchor) / (f/anchor_f)
            pts_m.append(M[iso]); pts_Srel.append(Srel)

add_isotope_family(df_std, f_Hf, 180)
add_isotope_family(df_std, f_U,  238)
add_isotope_family(df_std, f_Pb, 208)
add_isotope_family(df_std, {232:1.0}, 232)  # Th single isotope

assert len(pts_m)>=3, "Not enough isotopic points to fit a mass law."
x = np.log(np.array(pts_m))
y = np.log(np.array(pts_Srel))
beta, logk_rel = np.polyfit(x, y, 1)
k_rel = np.exp(logk_rel)

def S_shape(m):  # relative counts/atom vs mass
    return k_rel*(m**beta)

# ---------------- ABSOLUTE ANCHOR FROM 238U PPM ----------------
mass_std_g = crater_mass_removed_g(N_SHOTS_STD, RHO['zircon'], CRATER_DIAM_UM, DEPTH_PER_SHOT_UM)
N_U_total  = atoms_from_ppm(mass_std_g, U_ppm, M_U)
N_238      = N_U_total * f_U[238]
I_238      = safe_sum(df_std, COL[238])
assert I_238>0, "238U on standard is zero."
S238_abs   = I_238 / N_238
C_abs      = S238_abs / S_shape(M[238])

def S_abs(m):
    return C_abs * S_shape(m)

def counts_to_atoms(I_counts, iso):
    if I_counts is None or np.isnan(I_counts): return 0.0
    s = S_abs(M[iso])
    return 0.0 if s<=0 else float(I_counts)/s

# ---------------- PER-PARTICLE COUNTS TO ATOMS ----------------
for iso, col in COL.items():
    if col in df.columns:
        df[f'I_{iso}'] = pd.to_numeric(df[col], errors='coerce').fillna(0.0)
        df[f'N_{iso}'] = df[f'I_{iso}'].apply(lambda v: counts_to_atoms(v, iso))

df['N_Hf'] = df[[c for c in df.columns if c.startswith('N_') and any(k in c for k in ['176','177','178','179','180'])]].sum(axis=1)
df['N_U']  = df[[c for c in df.columns if c.startswith('N_') and any(k in c for k in ['235','238'])]].sum(axis=1)
df['N_Pb'] = df[[c for c in df.columns if c.startswith('N_') and any(k in c for k in ['206','207','208'])]].sum(axis=1)
df['N_Th'] = df[[c for c in df.columns if c.startswith('N_') and '232' in c]].sum(axis=1)

# ---------------- ZIRCON SITE COUNT FROM Hf ----------------
EPS = 1e-12
N_sites = np.where(df['N_Hf']>EPS, df['N_Hf']/max(F_HF_REF,EPS), 0.0)
zircon_present = N_sites > 0

# ---------------- PARTICLE MASS MODEL ----------------
N_Uv  = df['N_U'].values
N_Thv = df['N_Th'].values
N_Pbv = df['N_Pb'].values
N_Hfv = df['N_Hf'].values

N_Zr_sub = np.maximum(N_sites - (N_Hfv + N_Uv + N_Thv + N_Pbv), 0.0)

M_unit_sub = ( (N_Zr_sub*M_Zr) + (N_Hfv*M_Hf) + (N_Uv*M_U) + (N_Thv*M_Th) + (N_Pbv*M_Pb) ) / np.maximum(N_sites,1) \
             + M_Si + 4*M_O
m_sub_g = (N_sites/NA) * M_unit_sub

m_PbO2_g = (N_Pbv/NA) * M_PbO2
m_ThO2_g = (N_Thv/NA) * M_ThO2
m_UO2_g  = (N_Uv /NA) * M_UO2
m_ox_g   = m_PbO2_g + m_ThO2_g + m_UO2_g

m_particle_g = np.where(zircon_present, m_sub_g, m_ox_g)

# ---------------- EFFECTIVE DENSITY AND DIAMETER ----------------
rho_eff = np.full(len(df), np.nan)
rho_eff[zircon_present] = RHO['zircon']

w_PbO2 = np.zeros(len(df)); w_ThO2=np.zeros(len(df)); w_UO2=np.zeros(len(df))
mask_ox = ~zircon_present & (m_particle_g>0)
w_PbO2[mask_ox] = m_PbO2_g[mask_ox]/m_particle_g[mask_ox]
w_ThO2[mask_ox] = m_ThO2_g[mask_ox]/m_particle_g[mask_ox]
w_UO2 [mask_ox] = m_UO2_g [mask_ox]/m_particle_g[mask_ox]
with np.errstate(divide='ignore', invalid='ignore'):
    inv_rho = (w_PbO2/RHO['PbO2']) + (w_ThO2/RHO['ThO2']) + (w_UO2/RHO['UO2'])
rho_eff[mask_ox] = np.where(inv_rho[mask_ox]>0, 1.0/inv_rho[mask_ox], np.nan)

valid = (m_particle_g>0) & np.isfinite(rho_eff) & (rho_eff>0)
d_cm = np.full(len(df), np.nan)
d_cm[valid] = (6.0*m_particle_g[valid]/(pi*rho_eff[valid]))**(1.0/3.0)
d_um = d_cm*1e4

# ---------------- OUTPUT ----------------
df_out = df.copy()
df_out['N_sites']          = N_sites
df_out['m_particle_g']     = m_particle_g
df_out['rho_eff_g_cm3']    = rho_eff
df_out['d_particle_um']    = d_um
df_out['zircon_present']   = zircon_present

def label(i):
    if zircon_present[i]:
        parts = ['zircon(sub cations)']
        if N_Pbv[i]>0: parts.append('Pb(sub)')
        if N_Thv[i]>0: parts.append('Th(sub)')
        if N_Uv[i] >0: parts.append('U(sub)')
        return ' + '.join(parts)
    else:
        parts=[]
        if N_Pbv[i]>0: parts.append('PbO2')
        if N_Thv[i]>0: parts.append('ThO2')
        if N_Uv[i] >0: parts.append('UO2')
        return ' + '.join(parts) if parts else 'none'

df_out['Composition'] = [label(i) for i in range(len(df_out))]

if DROP_NONE_FROM_DENSITY:
    n_before = len(df_out)
    df_out = df_out[df_out['Composition'].str.lower() != 'none'].copy()
    df_out.reset_index(drop=True, inplace=True)
    print(f"[info] Dropped {n_before - len(df_out)} rows with Composition=='none' for density output.")

df_out.to_csv(OUTPUT_DENSITY, index=False)
print(f"[OK] Saved filtered + physics: {OUTPUT_DENSITY}")
print("Fit beta (mass law) =", round(float(beta),4))

