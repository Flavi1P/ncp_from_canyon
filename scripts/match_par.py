"""
match_par.py
------------
Match each Cphyto profile (by date + lon/lat) to the nearest pixel of the
corresponding daily MODIS-Aqua L3m PAR granule. Writes one row per profile.

Output columns: float_wmo, prof_number, date, lon, lat, par.
"""


import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr


PAR_VAR = "par"   # L3m variable name


def _nearest_idx(coord: np.ndarray, target: np.ndarray) -> np.ndarray:
    """Nearest-index lookup for a monotonic 1D coord (ascending or descending)."""
    if coord[0] > coord[-1]:
        return len(coord) - 1 - np.abs(coord[::-1][None, :] - target[:, None]).argmin(axis=1)
    return np.abs(coord[None, :] - target[:, None]).argmin(axis=1)


def match_par(
    cphyto_dir:   Path,
    par_manifest: Path,
    out_csv:      Path,
) -> None:
    cphyto_dir   = Path(cphyto_dir)
    par_manifest = Path(par_manifest)
    out_csv      = Path(out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    # Profile-level table (one row per profile)
    frames = []
    for csv in cphyto_dir.glob("argo_*_cphyto.csv"):
        df = pd.read_csv(csv, usecols=["float_wmo", "prof_number", "date", "lon", "lat"])
        frames.append(df.drop_duplicates(subset=["float_wmo", "prof_number"]))
    if not frames:
        print("No cphyto CSVs found.")
        pd.DataFrame(columns=["float_wmo", "prof_number", "date", "lon", "lat", "par"]).to_csv(out_csv, index=False)
        return
    profiles = pd.concat(frames, ignore_index=True)
    profiles["date"] = pd.to_datetime(profiles["date"]).dt.date.astype(str)
    print(f"{len(profiles)} profiles across {profiles['date'].nunique()} dates")

    par_idx = pd.read_csv(par_manifest)
    par_idx["date"] = par_idx["date"].astype(str)
    date_to_path = dict(zip(par_idx["date"], par_idx["path"]))

    # Group by date → open each PAR file once
    out_rows = []
    missing = 0
    for d, grp in profiles.groupby("date"):
        par_path = date_to_path.get(d)
        if par_path is None or not Path(par_path).exists():
            missing += len(grp)
            for _, r in grp.iterrows():
                out_rows.append({**r.to_dict(), "par": np.nan})
            continue

        ds = xr.open_dataset(par_path)
        lat = ds["lat"].values
        lon = ds["lon"].values
        par = ds[PAR_VAR].values   # shape (lat, lon)

        i_lat = _nearest_idx(lat, grp["lat"].to_numpy())
        i_lon = _nearest_idx(lon, grp["lon"].to_numpy())
        vals = par[i_lat, i_lon]
        # L3m flags fills with -32767 or NaN depending on encoding; ensure NaN
        vals = np.where(np.isfinite(vals), vals, np.nan)
        ds.close()

        for (_, r), v in zip(grp.iterrows(), vals):
            out_rows.append({**r.to_dict(), "par": float(v) if np.isfinite(v) else np.nan})

    out_df = pd.DataFrame(out_rows)
    out_df.to_csv(out_csv, index=False)
    n_ok = out_df["par"].notna().sum()
    print(f"Matched PAR: {n_ok}/{len(out_df)} profiles have a value  (missing files: {missing})")
    print(f"Output -> {out_csv}")


if __name__ == "__main__":
    if "snakemake" in globals():
        match_par(
            cphyto_dir   = Path(snakemake.input["cphyto_dir"]),     # noqa: F821
            par_manifest = Path(snakemake.input["par_manifest"]),   # noqa: F821
            out_csv      = Path(snakemake.output["matched_csv"]),   # noqa: F821
        )
    else:
        args = sys.argv[1:]
        cphyto_dir   = Path(args[0]) if len(args) > 0 else Path("data/NorthAtlantic_seas_comparison/intermediate/cphyto_profiles")
        par_manifest = Path(args[1]) if len(args) > 1 else Path("data/NorthAtlantic_seas_comparison/raw/par_download_manifest.csv")
        out_csv      = Path(args[2]) if len(args) > 2 else Path("data/NorthAtlantic_seas_comparison/intermediate/par_matched/par_matched.csv")
        match_par(cphyto_dir, par_manifest, out_csv)
