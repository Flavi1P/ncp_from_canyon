# NCP from CANYON-B

Estimates Net Community Production (NCP) in the open ocean from BGC-Argo
float data. Nitrate is predicted from hydrographic observations using the
CANYON-B neural network, and NCP is derived from the temporal change in
depth-integrated nitrate within a defined region and period.

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

### Nitrate as a proxy

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

---

## 2. Requirements

### Software

- Python >= 3.9
- R >= 4.1
- Snakemake >= 7.0 (the workflow manager)

Install Snakemake into a conda or pip environment:

    pip install snakemake
    # or
    conda install -c bioconda snakemake

### Python packages

    pip install pandas pyyaml

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
    +-- scripts/
    |   +-- download_sprof.py          Step 1: download BGC-Argo Sprof files
    |   +-- process_sprof.R            Step 2: parse and interpolate profiles
    |   +-- predict_nitrate.R          Step 3: predict nitrate with CANYON-B
    |   +-- merge_nitrate_files.R      Step 4: merge floats, compute MLD
    |   +-- compute_ncp.R              Step 5: compute NCP
    |   +-- compute_ncp_uncertainty.R  Step 5b: Monte Carlo uncertainty (optional)
    |   +-- utils/
    |       +-- fastr_canyon.R         CANYON-B forward pass
    |       +-- float_download.py      GDAC download helper
    |       +-- canyon_weights/        Pre-trained CANYON-B weight files
    +-- data/
    |   +-- raw/                       Shared Sprof .nc files (reused across runs)
    |   +-- {run_name}/                Per-run intermediate files
    +-- output/
        +-- {run_name}/                Per-run results and figures

The separation between `data/raw/` (shared) and `data/{run_name}/`
(per-run) means that changing spatial or temporal bounds triggers
re-processing but never re-downloads float files that are already on disk.

---

## 4. Quick start

### Step 1 -- edit the configuration

Open `config.yaml` and set your region and period:

    run_name:   "MyRegion_2022"
    lon_min: -45
    lon_max: -10
    lat_min:  54
    lat_max:  63
    date_start: "2022-01-01"
    date_end:   "2023-01-01"

### Step 2 -- run the pipeline

From the repository root:

    snakemake --cores 1

Snakemake will print a list of the steps it intends to run and then execute
them in order. You do not need to call the individual scripts yourself.

To use multiple cores (useful for the CANYON-B prediction step):

    snakemake --cores 4

### Step 3 -- inspect the outputs

Results are written to `output/{run_name}/`:

    output/MyRegion_2022/
    +-- float_map.png           Tracks of all floats used
    +-- ncp/
        +-- ncp_results.csv     NCP time series for each time step
        +-- ncp_sensitivity.png Sensitivity of NCP to the choice of time step

### What Snakemake does

Snakemake is a workflow manager. You declare rules (each rule is one
processing step with defined inputs and outputs), and Snakemake figures out
the order of execution and only reruns steps whose inputs have changed or
whose outputs are missing. You do not need prior experience with Snakemake
to use this pipeline -- the single command above is all that is required.

If you change `config.yaml` and rerun, Snakemake will rerun only the
affected downstream steps. To force a single rule to rerun regardless:

    snakemake --cores 1 --forcerun compute_ncp

---

## 5. Configuration reference

| Parameter | Type | Description |
|---|---|---|
| `run_name` | string | Name of this run. All outputs go under `data/{run_name}/` and `output/{run_name}/`. |
| `data_dir` | string | Root directory for intermediate data (default: `data`). |
| `raw_dir` | string | Shared directory for downloaded Sprof .nc files (default: `data/raw`). |
| `output_dir` | string | Root directory for final outputs (default: `output`). |
| `basins` | mapping | Named basins, each with a `polygon` key containing a list of `[lon, lat]` vertex pairs. Add as many basins as needed; the pipeline runs once per basin in parallel. |
| `date_start` | string | Start of the analysis period (YYYY-MM-DD, exclusive). |
| `date_end` | string | End of the analysis period (YYYY-MM-DD, exclusive). |
| `ncp_time_steps` | list | Time step widths for the NCP calculation (e.g. `["10 days", "15 days"]`). Multiple values produce a sensitivity plot. |
| `mld_spar` | float | Smoothing span for the LOESS fit to the MLD time series (0 to 1, default 0.3). |
| `zeu_default` | float | Euphotic zone depth in metres when not computed dynamically (default 40). |
| `uncertainty` | bool | Set to `true` to run the Monte Carlo uncertainty step. |
| `n_mc` | integer | Number of Monte Carlo iterations (default 200). |
| `canyon_rmse` | float | CANYON-B prediction uncertainty in mmol/m3, applied as correlated per-profile noise (default 1.2). |

---

## 6. Pipeline steps in detail

### Step 1 -- download_sprof.py

Queries the Argo GDAC index for BGC-Argo floats carrying a dissolved oxygen
sensor (DOXY) within the bounding box. Downloads `*_Sprof.nc` files to
`data/raw/`. Floats already on disk are skipped, so rerunning or changing
the temporal bounds does not trigger unnecessary downloads.

Outputs written to `data/{run_name}/raw/`:
- `download_manifest.csv` -- list of WMO numbers and file paths
- `wmo_list.txt` -- plain list of WMO numbers

### Step 2 -- process_sprof.R

Reads each `*_Sprof.nc` file. Extracts temperature, salinity, and
DOXY_ADJUSTED (QC flags 1 and 2 only). Interpolates each profile to a
regular 1 m depth grid from 0 to 2000 m using linear interpolation, with
a minimum of 10 valid observations required. Saves one CSV per float.

Output directory: `data/{run_name}/intermediate/doxy_profiles/`

### Step 3 -- predict_nitrate.R

Applies CANYON-B to each interpolated profile to predict nitrate
concentration at every depth level. CANYON-B takes date, latitude,
longitude, pressure, temperature, salinity, and dissolved oxygen as inputs
and returns nitrate in mmol/m3. The pre-trained weights in
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

1. Filters the data to the configured spatial and temporal bounds.
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

Outputs: `output/{run_name}/ncp/ncp_results.csv` and
`output/{run_name}/ncp/ncp_sensitivity.png`

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
error of approximately 1.2 mmol/m3 for nitrate. This error is correlated
within a profile because it reflects a systematic offset for a given water
mass rather than random noise at each depth. It is modelled as a single
offset drawn from N(0, canyon_rmse) and added uniformly to all depth
levels of each resampled profile before integration.

For each of the `n_mc` iterations and each time step, the full integration
pipeline is run on the perturbed nitrate fields. The output summarises the
5th percentile, mean, and 95th percentile of NCP across iterations at each
time bin.

Outputs: `output/{run_name}/ncp/ncp_uncertainty.csv` and
`output/{run_name}/ncp/ncp_uncertainty.png`

---

## 8. Outputs

| File | Description |
|---|---|
| `output/{run_name}/float_map.png` | Map of float tracks within the study region |
| `output/{run_name}/ncp/ncp_results.csv` | NCP time series (one row per time bin per time step) |
| `output/{run_name}/ncp/ncp_sensitivity.png` | NCP time series overlaid for all configured time steps |
| `output/{run_name}/ncp/ncp_uncertainty.csv` | MC uncertainty summary: mean, sd, q05, q95 per bin |
| `output/{run_name}/ncp/ncp_uncertainty.png` | NCP mean with 5-95% confidence ribbon per time step |
| `data/{run_name}/intermediate/merged/mld_timeseries.png` | Mixed layer depth across all profiles |

The main results CSV (`ncp_results.csv`) contains the following columns:

| Column | Description |
|---|---|
| `date_grid` | Centre date of the time bin |
| `mld` | Mean smoothed MLD for the bin (m) |
| `NCP` | Net Community Production (mmol C m-2 d-1) |
| `time_step_label` | Which time step this row belongs to |

---

## 9. Running individual scripts

Each script can be run from the command line without Snakemake, which is
useful for testing or modifying a single step. Pass input and output paths
as positional arguments:

    Rscript scripts/compute_ncp.R \
        data/MyRegion_2022/intermediate/merged/merged_ncp.csv \
        output/MyRegion_2022/ncp

    Rscript scripts/compute_ncp_uncertainty.R \
        data/MyRegion_2022/intermediate/merged/merged_ncp.csv \
        output/MyRegion_2022/ncp/ncp_uncertainty.csv \
        output/MyRegion_2022/ncp/ncp_uncertainty.png

If no arguments are provided the scripts fall back to default paths defined
near the top of each file, which can be edited for convenience.

---

## References

Bittig, H. C. et al. (2018). An alternative to static climatologies:
robust estimation of open ocean CO2 variables and nutrient concentrations
from T, S, and O2 data using Bayesian neural networks. Frontiers in Marine
Science, 5, 328.

Redfield, A. C. (1958). The biological control of chemical factors in the
environment. American Scientist, 46(3), 205-221.
