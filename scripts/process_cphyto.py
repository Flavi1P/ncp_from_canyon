"""
process_cphyto.py
-----------------
Per-float processing for phytoplankton carbon (Cphyto) from BGC-Argo Sprof files.

For each float with CHLA and BBP700 sensors, build 1 m-binned profiles from
surface to 200 m of {temp, sal, chla, bbp700}, fill NaN gaps <= 5 m by linear
interpolation, drop profiles with any larger gap, then compute Cphyto:

    bbp700_baseline = Briggs (2011) despike of bbp700  (minmax, window=11)
    bbp470          = bbp700_baseline * (700/470)^gamma,  gamma = 1.0
    cphyto          = 12128 * bbp470 + 0.59             (Graff et al. 2015)

Outputs one CSV per float + a manifest, mirroring the structure of
`process_sprof.R`.
"""


import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

# Local utils
sys.path.insert(0, str(Path(__file__).parent / "utils"))
from despike import despike                     # noqa: E402
from bbp_spectral import bbp_at_wavelength      # noqa: E402
from cphyto_graff2015 import cphyto_graff2015   # noqa: E402


DEPTH_MIN  = 0
DEPTH_MAX  = 200
DEPTH_STEP = 1
GAP_LIMIT  = 5      # m, max NaN gap tolerated after binning
BBP_WINDOW = 11     # samples (=metres here) for Briggs despike

# BGC-Argo ADJUSTED QC flags considered usable: 1=good, 2=probably good,
# 5=value changed (real-time correction applied — dominant for CHLA/BBP),
# 8=interpolated/estimated (e.g. CHLA NPQ correction).
QC_OK = {"1", "2", "3", "5", "8"}

TARGET_DEPTHS = np.arange(DEPTH_MIN, DEPTH_MAX + 1, DEPTH_STEP)


def _qc_mask(qc_arr) -> np.ndarray:
    """Boolean mask of samples with QC flag in QC_OK."""
    qc = np.asarray(qc_arr)
    if qc.dtype.kind in ("S", "O", "U"):
        flat = np.array([b.decode() if isinstance(b, bytes) else ("" if b is None else str(b))
                         for b in qc.ravel()]).reshape(qc.shape)
        return np.isin(flat, list(QC_OK))
    return np.isin(qc.astype(int), [int(f) for f in QC_OK])


# ── Per-profile cleaning ─────────────────────────────────────────────────────
def _bin_and_fill(df_raw: pd.DataFrame) -> pd.DataFrame | None:
    """Bin native samples to 1 m grid, interp gaps <= GAP_LIMIT, else return None.

    Input columns: depth, temp, sal, chla, bbp700 (non-NaN native samples).
    Output: a 201-row DataFrame on TARGET_DEPTHS with no NaN, or None if the
    profile cannot be filled within the gap tolerance.
    """
    df_raw = (
        df_raw.dropna(subset=["depth"])
              .assign(depth=lambda d: d["depth"].round().astype(int))
              .query(f"{DEPTH_MIN - GAP_LIMIT} <= depth <= {DEPTH_MAX + GAP_LIMIT}")
              .groupby("depth", as_index=False)
              .mean(numeric_only=True)
              .sort_values("depth")
    )
    if df_raw.empty:
        return None

    # Reindex onto a grid that extends GAP_LIMIT beyond [0, 200] so surface /
    # bottom edge gaps can be filled by linear interp instead of extrapolation.
    grid = np.arange(DEPTH_MIN - GAP_LIMIT, DEPTH_MAX + GAP_LIMIT + 1, DEPTH_STEP)
    df = df_raw.set_index("depth").reindex(grid)

    # Linear interp with a per-variable gap limit (pandas counts consecutive
    # NaNs, which equals metres on a 1 m grid). `limit_area='inside'` skips
    # leading/trailing NaNs at the extended-grid edges.
    df = df.interpolate(method="linear", limit=GAP_LIMIT,
                        limit_area="inside", axis=0)

    df = df.loc[DEPTH_MIN:DEPTH_MAX]
    if df.isna().any().any():
        return None
    return df.reset_index()


# ── Main per-float pipeline ──────────────────────────────────────────────────
def process_float(sprof_path: Path, out_dir: Path) -> dict | None:
    """Process one Sprof file; write CSV and return a manifest row, or None."""
    ds = xr.open_dataset(sprof_path)

    required = {"CHLA_ADJUSTED", "BBP700_ADJUSTED",
                "TEMP_ADJUSTED", "PSAL_ADJUSTED", "PRES_ADJUSTED"}
    if not required.issubset(ds.variables):
        ds.close()
        return None

    wmo_raw = ds["PLATFORM_NUMBER"].values[0]
    wmo = (wmo_raw.decode() if isinstance(wmo_raw, bytes) else str(wmo_raw)).strip()
    lon = ds["LONGITUDE"].values
    lat = ds["LATITUDE"].values
    date = pd.to_datetime(ds["JULD"].values).date

    pres = ds["PRES_ADJUSTED"].values
    temp = ds["TEMP_ADJUSTED"].values
    sal  = ds["PSAL_ADJUSTED"].values
    chla = ds["CHLA_ADJUSTED"].values
    bbp  = ds["BBP700_ADJUSTED"].values

    # QC masks — keep flags 1 & 2 per variable
    temp = np.where(_qc_mask(ds["TEMP_ADJUSTED_QC"].values), temp, np.nan)
    sal  = np.where(_qc_mask(ds["PSAL_ADJUSTED_QC"].values), sal,  np.nan)
    chla = np.where(_qc_mask(ds["CHLA_ADJUSTED_QC"].values), chla, np.nan)
    bbp  = np.where(_qc_mask(ds["BBP700_ADJUSTED_QC"].values), bbp, np.nan)
    ds.close()

    n_prof = pres.shape[0]
    kept = []
    for i in range(n_prof):
        if np.isnan(lon[i]) or np.isnan(lat[i]):
            continue
        prof = pd.DataFrame({
            "depth":   pres[i],
            "temp":    temp[i],
            "sal":     sal[i],
            "chla":    chla[i],
            "bbp700":  bbp[i],
        }).dropna(how="all", subset=["temp", "sal", "chla", "bbp700"])

        filled = _bin_and_fill(prof)
        if filled is None:
            continue

        # Despike bbp700 (Briggs 2011) and derive cphyto
        baseline, _ = despike(filled["bbp700"].to_numpy(),
                              window_size=BBP_WINDOW, spike_method="minmax")
        filled["bbp700_baseline"] = baseline
        filled["bbp470"]          = bbp_at_wavelength(baseline, 700.0, 470.0, 1.0)
        filled["cphyto"]          = cphyto_graff2015(filled["bbp470"])

        filled.insert(0, "lon", lon[i])
        filled.insert(1, "lat", lat[i])
        filled.insert(2, "date", date[i])
        filled.insert(3, "float_wmo", wmo)
        filled.insert(4, "prof_number", i + 1)
        kept.append(filled)

    if not kept:
        return None

    out = pd.concat(kept, ignore_index=True)
    csv_path = out_dir / f"argo_{wmo}_cphyto.csv"
    out.to_csv(csv_path, index=False)
    return {"wmo": wmo, "n_profiles": out["prof_number"].nunique(), "path": str(csv_path)}


# ── Entry point ──────────────────────────────────────────────────────────────
def main(sprof_dir: Path, out_dir: Path, manifest_path: Path) -> None:
    sprof_dir = Path(sprof_dir)
    out_dir   = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    files = sorted(sprof_dir.glob("*_Sprof.nc"))
    print(f"Found {len(files)} Sprof files in {sprof_dir}")

    rows = []
    for i, f in enumerate(files, 1):
        try:
            row = process_float(f, out_dir)
        except Exception as e:
            print(f"  [{i}/{len(files)}] {f.name}: FAILED ({e})")
            continue
        if row is None:
            print(f"  [{i}/{len(files)}] {f.name}: skipped (no CHLA/BBP or no valid profiles)")
        else:
            print(f"  [{i}/{len(files)}] {f.name}: {row['n_profiles']} profiles kept")
            rows.append(row)

    pd.DataFrame(rows).to_csv(manifest_path, index=False)
    print(f"Manifest -> {manifest_path}")


if __name__ == "__main__":
    if "snakemake" in globals():
        main(
            sprof_dir     = Path(snakemake.params["sprof_dir"]),                    # noqa: F821
            out_dir       = Path(snakemake.output["out_dir"]),                      # noqa: F821
            manifest_path = Path(snakemake.output["manifest"]),                     # noqa: F821
        )
    else:
        args = sys.argv[1:]
        sprof_dir = Path(args[0]) if len(args) > 0 else Path("data/raw")
        out_dir   = Path(args[1]) if len(args) > 1 else Path("data/NorthAtlantic_seas_comparison/intermediate/cphyto_profiles")
        manifest  = out_dir / "cphyto_manifest.csv"
        main(sprof_dir, out_dir, manifest)
