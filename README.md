# NCP and NPP from BGC-Argo

Estimates Net Community Production (NCP) and Net Primary Production (NPP)
in the open ocean from BGC-Argo float data. Nitrate is predicted from
hydrographic observations using the CANYON-B neural network, and NCP is
derived from the temporal change in depth-integrated nitrate. NPP is
independently estimated using the Carbon-based Productivity Model (CbPM,
Westberry et al. 2008) applied to phytoplankton carbon (Cphyto) profiles
derived from BGC-Argo backscatter measurements, with surface PAR from
daily MODIS-Aqua L3 composites.

---

## Table of contents

1. [Scientific background](#1-scientific-background)
2. [Requirements](#2-requirements)
3. [Repository layout](#3-repository-layout)
4. [Quick start](#4-quick-start)
5. [Configuration reference](#5-configuration-reference)
6. [Pipeline steps in detail](#6-pipeline-steps-in-detail)
7. [Uncertainty estimation](#7-uncertainty-estimation)
8. [Outputs](#8-outputs)
9. [Running individual scripts](#9-running-individual-scripts)

---

## 1. Scientific background

### What is NCP?

Net Community Production is the difference between gross photosynthesis and
total community respiration in the surface ocean. A positive NCP indicates
that the ecosystem is fixing more carbon than it is consuming, and represents
the potential for organic carbon export to depth. It is expressed here in
units of mmol C m-2 d-1.

### Nitrate as a proxy for NCP

In the absence of nitrogen fixation and denitrification, the drawdown of
nitrate in the surface ocean tracks biological carbon fixation via the
Redfield ratio (C:N = 106:16 = 6.625). This pipeline uses nitrate as its
primary tracer. Because BGC-Argo floats measure temperature, salinity, and
dissolved oxygen but not always nitrate directly, nitrate is predicted at
every depth level using CANYON-B, a Bayesian neural network trained on
global ocean data (Bittig et al. 2018).

### NCP calculation

For each time window of width dt (e.g. 10 days), NCP is computed as:

    NCP = dN/dt

where dN/dt is the change in depth-integrated nitrate between consecutive
time windows, converted to carbon units via the Redfield ratio (C:N = 6.625)
and normalised by dt.

The integration depth for each window is the maximum of the current MLD,
the previous MLD, the euphotic-zone depth (Zeu), and the previous Zeu,
ensuring that the full productive layer is always included.

MLD is computed from the density profiles of each float using a threshold
criterion of 0.03 kg m-3 relative to the 0-5 m average. A LOESS smoother
is applied to the MLD time series before integration to reduce noise from
individual profiles.

### What is NPP and how is it estimated?

Net Primary Production is the rate of organic carbon fixation by
phytoplankton, before accounting for community respiration. It is estimated
here using the Carbon-based Productivity Model (CbPM, Westberry et al.
2008), which is driven by three inputs:

- **Cphyto** -- phytoplankton carbon concentration (mg C m-3), derived
  from BGC-Argo particulate backscatter at 700 nm (bbp700). bbp700 is
  despiked, spectrally converted to bbp470, and converted to Cphyto using
  the empirical relationship of Graff et al. (2015):
  Cphyto = 12128 × bbp470 + 0.59.
- **PAR** -- photosynthetically active radiation at the surface (mol
  photons m-2 d-1), from daily MODIS-Aqua L3m 9 km composites downloaded
  from NASA Earthdata.
- **mu** -- phytoplankton growth rate, computed internally by CbPM from
  Cphyto and chlorophyll-a (CHLA_ADJUSTED from BGC-Argo).

CbPM is applied profile by profile on a 0-200 m depth grid. Volumetric
NPP (mg C m-3 d-1) is integrated to the euphotic zone depth to produce
column NPP (mg C m-2 d-1). These are then averaged within each basin and
binned in time to produce a NPP time series comparable to NCP.

---

## 2. Requirements

### Software

- Python >= 3.9
- R >= 4.1
- Snakemake >= 7.0 (the workflow manager)
- A NASA Earthdata account (free) with a `~/.netrc` entry for
  `urs.earthdata.nasa.gov` (required for PAR downloads; run
  `earthaccess.login(persist=True)` once to create it)

Install Snakemake into a conda or pip environment:

    pip install snakemake
    # or
    conda install -c bioconda snakemake

The recommended approach is to use the provided `environment.yml`:

    conda env create -f environment.yml
    conda activate argo-ncp

### Python packages

    pip install pandas pyyaml xarray earthaccess matplotlib shapely

The float downloader uses `argo_gdac` from the `argopy` ecosystem. Ensure
`scripts/utils/float_download.py` is present (it is included in this repo).

### R packages

    install.packages(c(
      "tidyverse", "ncdf4", "zoo", "lubridate", "readr",
      "gsw", "mgcv", "data.table", "seacarb", "here"
    ))

    # castr for MLD calculation (GitHub)
    remotes::install_github("cassidymwagner/castr")

The CANYON-B weights are included in `scripts/utils/canyon_weights/` and
are loaded automatically.

---

## 3. Repository layout

    .
    +-- Snakefile                      Workflow definition
    +-- config.yaml                    All user-facing settings
    +-- environment.yml                Conda environment specification
    +-- scripts/
    |   +-- download_sprof.py          Step 1: download BGC-Argo Sprof files
    |   +-- process_sprof.R            Step 2: parse and QC doxy profiles
    |   +-- process_cphyto.py          Step 2b: extract Cphyto from CHLA/BBP
    |   +-- download_par.py            Step 3 (NPP): download MODIS PAR subsets
    |   +-- match_par.py               Step 4 (NPP): match PAR pixel to each profile
    |   +-- predict_nitrate.R          Step 3: predict nitrate with CANYON-B
    |   +-- merge_nitrate_files.R      Step 4: merge floats, compute MLD
    |   +-- compute_ncp.R              Step 5: compute NCP per basin
    |   +-- compute_ncp_uncertainty.R  Step 5b: Monte Carlo uncertainty (optional)
    |   +-- compare_basins.R           Step 6: cross-basin NCP comparison plot
    |   +-- compute_npp.py             Step 5 (NPP): run CbPM per profile
    |   +-- npp_timeseries.py          Step 6 (NPP): bin NPP, overlay with NCP
    |   +-- compare_basins_npp.py      Step 7 (NPP): cross-basin NPP comparison
    |   +-- utils/
    |       +-- fastr_canyon.R         CANYON-B forward pass
    |       +-- float_download.py      GDAC download helper
    |       +-- cbpm_argo.py           CbPM NPP model (Westberry et al. 2008)
    |       +-- cphyto_graff2015.py    Cphyto from bbp (Graff et al. 2015)
    |       +-- bbp_spectral.py        bbp spectral conversion (700 → 470 nm)
    |       +-- bbp_from_beta.py       VSF → bbp conversion
    |       +-- despike.py             Briggs (2011) despiking
    |       +-- daylength.py           Photoperiod for CbPM
    |       +-- AustinPetzold_1986.py  PAR attenuation coefficients
    |       +-- canyon_weights/        Pre-trained CANYON-B weight files
    +-- data/
    |   +-- raw/
    |   |   +-- modis_par/             Shared bbox-subsetted MODIS PAR .nc files
    |   +-- {run_name}/                Per-run intermediate files
    +-- output/
        +-- {run_name}/                Per-run results and figures

The separation between `data/raw/` (shared) and `data/{run_name}/`
(per-run) means that changing spatial or temporal bounds triggers
re-processing but never re-downloads Sprof or PAR files already on disk.

---

## 4. Quick start

### Step 1 -- edit the configuration

Open `config.yaml` and define your basins and period:

    run_name: "MyRegion_2022"
    basins:
      IcelandBasin:
        polygon:
          - [-32, 57.4]
          - [-22.4, 63]
          - [-12, 63]
          - [-23, 55.7]
    date_start: "2022-01-01"
    date_end:   "2023-01-01"

To enable NPP estimation, add:

    compute_npp: true
    par_resolution_km: 9   # 4 or 9 km MODIS grid

### Step 2 -- authenticate with NASA Earthdata (NPP only)

If `compute_npp: true`, you need a free NASA Earthdata account. Run once
in a Python session:

    import earthaccess
    earthaccess.login(persist=True)

This writes your credentials to `~/.netrc` and is read automatically on
subsequent runs.

### Step 3 -- run the pipeline

From the repository root (inside the conda environment):

    snakemake --cores 4

Snakemake will print a list of the steps it intends to run and then execute
them in order. You do not need to call the individual scripts yourself.

### Step 4 -- inspect the outputs

Results are written to `output/{run_name}/`:

    output/MyRegion_2022/
    +-- float_map.png
    +-- ncp/
    |   +-- {basin}/ncp_results.csv
    |   +-- {basin}/ncp_sensitivity.png
    |   +-- {basin}/ncp_uncertainty.csv      (if uncertainty: true)
    |   +-- {basin}/ncp_uncertainty.png
    |   +-- basins_comparison.png
    +-- npp/                                  (if compute_npp: true)
        +-- {basin}/npp_timeseries.csv
        +-- {basin}/npp_vs_ncp.png
        +-- {basin}/npp_map.png
        +-- basins_npp_comparison.png
        +-- basins_npp_map.png

### What Snakemake does

Snakemake is a workflow manager. You declare rules (each rule is one
processing step with defined inputs and outputs), and Snakemake figures out
the order of execution and only reruns steps whose inputs have changed or
whose outputs are missing. You do not need prior experience with Snakemake
to use this pipeline -- the single command above is all that is required.

If you change `config.yaml` and rerun, Snakemake will rerun only the
affected downstream steps. To force a single rule to rerun regardless:

    snakemake --cores 4 --forcerun compute_ncp

---

## 5. Configuration reference

| Parameter | Type | Default | Description |
|---|---|---|---|
| `run_name` | string | — | Name of this run; all outputs go under `data/{run_name}/` and `output/{run_name}/`. |
| `data_dir` | string | `data` | Root directory for intermediate data. |
| `raw_dir` | string | `data/raw` | Shared directory for downloaded Sprof and PAR files. |
| `output_dir` | string | `output` | Root directory for final outputs. |
| `basins` | mapping | — | Named basins, each with a `polygon` key containing a list of `[lon, lat]` vertex pairs. The pipeline runs once per basin. |
| `date_start` | string | — | Start of the analysis period (YYYY-MM-DD). |
| `date_end` | string | — | End of the analysis period (YYYY-MM-DD). |
| `ncp_time_steps` | list | — | Time step widths for NCP (e.g. `["10 days", "15 days"]`). Multiple values produce a sensitivity plot. |
| `mld_spar` | float | `0.3` | Smoothing span for the LOESS fit to the MLD time series (0–1). |
| `zeu_default` | float | `40` | Euphotic zone depth in metres when not computed dynamically. |
| `uncertainty` | bool | `false` | Set to `true` to run the Monte Carlo NCP uncertainty step. |
| `n_mc` | integer | `200` | Number of Monte Carlo iterations for uncertainty estimation. |
| `canyon_rmse` | float | `1.2` | CANYON-B nitrate prediction uncertainty (mmol m-3), applied as correlated per-profile noise. |
| `compute_npp` | bool | `false` | Set to `true` to activate the full NPP branch (Cphyto + PAR download + CbPM). |
| `par_resolution_km` | integer | `9` | Spatial resolution of MODIS-Aqua PAR product to download (4 or 9 km). |

---

## 6. Pipeline steps in detail

### Step 1 -- download_sprof.py

Queries the Argo GDAC index for BGC-Argo floats carrying a dissolved oxygen
sensor (DOXY) within the union bounding box of all configured basins.
Downloads `*_Sprof.nc` files to `data/raw/`. Floats already on disk are
skipped, so rerunning or changing temporal bounds does not trigger
unnecessary downloads.

Outputs written to `data/{run_name}/raw/`:
- `download_manifest.csv` -- list of WMO numbers and file paths
- `wmo_list.txt` -- plain list of WMO numbers

### Step 2 -- process_sprof.R

Reads each `*_Sprof.nc` file. Extracts temperature, salinity, and
DOXY_ADJUSTED (QC flags 1 and 2 only). Interpolates each profile to a
regular 1 m depth grid from 0 to 2000 m using linear interpolation, with
a minimum of 10 valid observations required. Saves one CSV per float.

Output directory: `data/{run_name}/intermediate/doxy_profiles/`

### Step 2b -- process_cphyto.py *(NPP branch)*

Reads the same `*_Sprof.nc` files and extracts CHLA_ADJUSTED and
BBP700_ADJUSTED profiles (QC flags 1, 2, 5, 8 accepted). For each profile:

1. Interpolates each variable to a regular 1 m grid (0–200 m). NaN gaps
   ≤ 5 m are filled by linear interpolation; profiles with any larger gap
   are dropped entirely.
2. Despiked BBP700 using the Briggs (2011) minmax method (window = 11 m).
3. Converts bbp700 to bbp470 using a power-law spectral slope (γ = 1.0).
4. Converts bbp470 to Cphyto (mg C m-3) using Graff et al. (2015):
   Cphyto = 12128 × bbp470 + 0.59.

Outputs one CSV per float to `data/{run_name}/intermediate/cphyto_profiles/`
and a manifest `cphyto_manifest.csv`.

### Step 3 (NPP) -- download_par.py *(NPP branch)*

Downloads daily MODIS-Aqua L3m PAR granules (`MODISA_L3m_PAR` v2022.0)
from NASA Earthdata Cloud via `earthaccess`. Rather than downloading full
global files (~8–25 MB each), the script opens each granule lazily over
HTTPS and slices to the union bounding box of all configured basins,
producing small (~100 kB) local NetCDF subsets. Files already on disk are
cached and reused across runs.

Requires a valid `~/.netrc` entry for `urs.earthdata.nasa.gov`. Run
`earthaccess.login(persist=True)` once to create it.

PAR files are saved to `data/raw/modis_par/` (shared across runs).
A manifest `data/{run_name}/raw/par_download_manifest.csv` lists
the date, filename, and path of every downloaded granule.

### Step 4 (NPP) -- match_par.py *(NPP branch)*

For each Cphyto profile, finds the nearest PAR pixel in the corresponding
daily granule using the profile's date, longitude, and latitude. Handles
the descending latitude axis of MODIS L3m grids automatically.

Output: `data/{run_name}/intermediate/par_matched/par_matched.csv`
with columns `float_wmo, prof_number, date, lon, lat, par`.

### Step 3 -- predict_nitrate.R

Applies CANYON-B to each interpolated doxy profile to predict nitrate
concentration at every depth level. CANYON-B takes date, latitude,
longitude, pressure, temperature, salinity, and dissolved oxygen as inputs
and returns nitrate in mmol m-3. The pre-trained weights in
`scripts/utils/canyon_weights/` are used directly without any online
component.

Output directory: `data/{run_name}/intermediate/nitrate_profiles/`

### Step 4 -- merge_nitrate_files.R

Combines all per-float CSV files into a single data frame. Computes
Absolute Salinity, Conservative Temperature, and potential density
(sigma0) using the TEOS-10 equations via the `gsw` package. Computes MLD
for each profile using a density threshold criterion of 0.03 kg m-3
relative to the 0-5 m reference layer. Saves the merged dataset and a
float track map.

Output: `data/{run_name}/intermediate/merged/merged_ncp.csv`

### Step 5 -- compute_ncp.R

Applies the NCP calculation described in Section 1 to the merged dataset.
Runs for each time step listed in `ncp_time_steps`. For each time step:

1. Filters the data to the configured spatial polygon and temporal bounds.
2. Smooths the MLD time series with LOESS (span = `mld_spar`).
3. Builds a regular time grid and assigns each profile to a bin.
4. Averages nitrate profiles within each bin and interpolates across
   empty bins with linear interpolation.
5. Integrates nitrate from the surface to the integration depth for each
   bin using trapezoidal integration.
6. Computes dN/dt and converts to carbon units via the Redfield ratio.
7. Applies a final LOESS smooth to the NCP time series.

The smoothing span for the final NCP LOESS is set adaptively as
`max(0.4, 6 / n)` to avoid numerical failure when a long time step
produces few bins.

Outputs per basin: `output/{run_name}/ncp/{basin}/ncp_results.csv` and
`output/{run_name}/ncp/{basin}/ncp_sensitivity.png`

### Step 5 (NPP) -- compute_npp.py *(NPP branch)*

Runs CbPM (Westberry et al. 2008) on every Cphyto profile that has a
matched surface PAR value. For each profile:

1. Loads the 0-200 m Cphyto and CHLA profiles (1 m resolution).
2. Computes the euphotic zone depth (Zeu) from PAR attenuation using
   Austin & Petzold (1986) Kd coefficients based on chlorophyll.
3. Runs CbPM to produce volumetric NPP (mg C m-3 d-1) at each depth.
4. Integrates volumetric NPP from 0 to Zeu (and 0 to 200 m) using a 1 m
   step to yield column NPP (mg C m-2 d-1).

Outputs:
- `data/{run_name}/intermediate/npp/npp_profiles.csv` -- volumetric NPP
  at each depth for every profile
- `data/{run_name}/intermediate/npp/npp_integrated.csv` -- one row per
  profile with `par_surface, zeu_m, npp_int_0_200, npp_int_0_zeu`

### Step 6 (NPP) -- npp_timeseries.py *(NPP branch)*

Filters integrated NPP to profiles inside the basin polygon (using
Shapely), bins them in time using the first entry in `ncp_time_steps`,
and computes mean ± sd per bin. Produces two outputs per basin:

- **npp_vs_ncp.png** -- NPP mean ± 1σ on the left axis, NCP (converted
  from mmol to mg C via × 12.011) on the right axis.
- **npp_map.png** -- scatter plot of profile locations coloured by NPP
  with the basin polygon overlaid.

Outputs per basin under `output/{run_name}/npp/{basin}/`.

### Step 7 (NPP) -- compare_basins_npp.py *(NPP branch)*

Produces two cross-basin summary figures:

- **basins_npp_comparison.png** -- NPP mean time series for all basins
  overlaid on a single axes, with ± 1σ shading.
- **basins_npp_map.png** -- all profile locations across all basins
  coloured by NPP with every basin polygon drawn.

Outputs: `output/{run_name}/npp/basins_npp_comparison.png` and
`output/{run_name}/npp/basins_npp_map.png`

---

## 7. Uncertainty estimation

Set `uncertainty: true` in `config.yaml` to activate
`compute_ncp_uncertainty.R`. This runs after the main NCP step and writes
its results alongside the standard outputs.

### Method: profile-level bootstrap with CANYON noise

The uncertainty has two sources.

**Spatial sampling uncertainty.** Each time bin typically contains
observations from 5 to 10 floats. The bin average is sensitive to which
profiles happened to be sampled. For each Monte Carlo iteration, profiles
are resampled with replacement within each bin (bootstrap). This preserves
the full vertical structure of each profile and correctly propagates the
covariance between depth levels through the integration step. Sampling
each depth level independently would cause errors to cancel during
integration and severely underestimate the true spread.

**CANYON-B prediction uncertainty.** CANYON-B has a root-mean-square
error of approximately 1.2 mmol m-3 for nitrate. This error is correlated
within a profile because it reflects a systematic offset for a given water
mass rather than random noise at each depth. It is modelled as a single
offset drawn from N(0, canyon_rmse) and added uniformly to all depth
levels of each resampled profile before integration.

For each of the `n_mc` iterations and each time step, the full integration
pipeline is run on the perturbed nitrate fields. The output summarises the
5th percentile, mean, and 95th percentile of NCP across iterations at each
time bin.

Outputs per basin: `output/{run_name}/ncp/{basin}/ncp_uncertainty.csv`
and `output/{run_name}/ncp/{basin}/ncp_uncertainty.png`

---

## 8. Outputs

### NCP outputs

| File | Description |
|---|---|
| `output/{run_name}/float_map.png` | Map of float tracks within the study region |
| `output/{run_name}/ncp/{basin}/ncp_results.csv` | NCP time series (one row per time bin per time step) |
| `output/{run_name}/ncp/{basin}/ncp_sensitivity.png` | NCP time series overlaid for all configured time steps |
| `output/{run_name}/ncp/{basin}/ncp_uncertainty.csv` | MC uncertainty summary: mean, sd, q05, q95 per bin |
| `output/{run_name}/ncp/{basin}/ncp_uncertainty.png` | NCP mean with 5–95% confidence ribbon per time step |
| `output/{run_name}/ncp/basins_comparison.png` | Cross-basin NCP uncertainty comparison |

The main results CSV (`ncp_results.csv`) contains the following columns:

| Column | Description |
|---|---|
| `date_grid` | Centre date of the time bin |
| `mld` | Mean smoothed MLD for the bin (m) |
| `NCP` | Net Community Production (mmol C m-2 d-1) |
| `time_step_label` | Which time step this row belongs to |

### NPP outputs *(if `compute_npp: true`)*

| File | Description |
|---|---|
| `output/{run_name}/npp/{basin}/npp_timeseries.csv` | Binned NPP mean, sd, and profile count per time bin |
| `output/{run_name}/npp/{basin}/npp_vs_ncp.png` | NPP time series overlaid with NCP |
| `output/{run_name}/npp/{basin}/npp_map.png` | Profile locations coloured by NPP within the basin |
| `output/{run_name}/npp/basins_npp_comparison.png` | Cross-basin NPP mean time series |
| `output/{run_name}/npp/basins_npp_map.png` | All-basin profile map coloured by NPP |

The NPP time series CSV (`npp_timeseries.csv`) contains:

| Column | Description |
|---|---|
| `date_bin` | Start date of the time bin |
| `npp_mean` | Mean column NPP across profiles in the bin (mg C m-2 d-1) |
| `npp_sd` | Standard deviation of column NPP |
| `n_profiles` | Number of profiles in the bin |

---

## 9. Running individual scripts

Each Python script can be run standalone without Snakemake. The scripts
read positional arguments from the command line and fall back to default
paths (defined near the top of each file) when no arguments are provided.

    # Cphyto extraction
    python scripts/process_cphyto.py

    # PAR download (reads config.yaml for bbox)
    python scripts/download_par.py config.yaml \
        data/{run_name}/intermediate/cphyto_profiles \
        data/raw/modis_par \
        data/{run_name}/raw/par_download_manifest.csv

    # NPP computation
    python scripts/compute_npp.py

    # NPP time series for a specific basin
    python scripts/npp_timeseries.py IcelandBasin

The R scripts similarly fall back to defaults:

    Rscript scripts/compute_ncp.R \
        data/MyRegion_2022/intermediate/merged/merged_ncp.csv \
        output/MyRegion_2022/ncp

    Rscript scripts/compute_ncp_uncertainty.R \
        data/MyRegion_2022/intermediate/merged/merged_ncp.csv \
        output/MyRegion_2022/ncp/ncp_uncertainty.csv \
        output/MyRegion_2022/ncp/ncp_uncertainty.png

---

## References

Bittig, H. C. et al. (2018). An alternative to static climatologies:
robust estimation of open ocean CO2 variables and nutrient concentrations
from T, S, and O2 data using Bayesian neural networks. Frontiers in Marine
Science, 5, 328.

Graff, J. R. et al. (2015). Analytical phytoplankton carbon measurements
spanning diverse ecosystems. Deep-Sea Research I, 102, 16-25.

Redfield, A. C. (1958). The biological control of chemical factors in the
environment. American Scientist, 46(3), 205-221.

Westberry, T. et al. (2008). Carbon-based primary productivity modeling
with vertically resolved photoacclimation. Global Biogeochemical Cycles,
22, GB2024.
