"""
download_par.py
---------------
Download bbox-subsetted daily MODIS-Aqua L3m PAR files (one per unique profile
date) using earthaccess. We open each granule lazily via `earthaccess.open()`,
slice to the union bbox of the configured basins, and write the resulting
subset locally. This keeps files ~100 kB each instead of ~8-25 MB.

Auth: user must run `earthaccess.login(persist=True)` once locally to create
~/.netrc. This script calls `earthaccess.login()` which reads that file.
"""


import sys
from pathlib import Path
from datetime import date

import pandas as pd
import xarray as xr
import earthaccess


SHORT_NAME = "MODISA_L3m_PAR"
VERSION    = "2022.0"
BBOX_PAD   = 0.5   # degrees of padding on each side of the basin bbox


def _bbox_from_config(cfg) -> tuple[float, float, float, float]:
    """Return (lon_min, lon_max, lat_min, lat_max) from config basins dict."""
    if "basins" in cfg:
        all_lons, all_lats = [], []
        for basin in cfg["basins"].values():
            for lon, lat in basin["polygon"]:
                all_lons.append(lon)
                all_lats.append(lat)
        return min(all_lons), max(all_lons), min(all_lats), max(all_lats)
    return cfg["lon_min"], cfg["lon_max"], cfg["lat_min"], cfg["lat_max"]


def _unique_dates_from_cphyto(cphyto_manifest: Path) -> list[str]:
    manifest = pd.read_csv(cphyto_manifest)
    dates: set[date] = set()
    for csv_path in manifest["path"]:
        dates.update(pd.read_csv(csv_path, usecols=["date"], parse_dates=["date"])["date"].dt.date.unique())
    return sorted(d.isoformat() for d in dates)


def _subset_and_save(src_obj, out_path: Path,
                     lon_min: float, lon_max: float,
                     lat_min: float, lat_max: float) -> bool:
    """Open a lazy granule handle, slice to bbox, write local NetCDF."""
    with xr.open_dataset(src_obj) as ds:
        # MODIS L3m lat is 90 -> -90 (descending); use slice accordingly.
        lat = ds["lat"].values
        lat_slice = slice(lat_max, lat_min) if lat[0] > lat[-1] else slice(lat_min, lat_max)
        sub = ds.sel(lat=lat_slice, lon=slice(lon_min, lon_max))
        if sub["lat"].size == 0 or sub["lon"].size == 0:
            return False
        # Load into memory then write — avoids dask issues with h5netcdf engine
        sub.load().to_netcdf(out_path)
    return True


def download_par(
    cphyto_manifest: Path,
    par_dir:         Path,
    manifest_path:   Path,
    bbox:            tuple[float, float, float, float],
    resolution_km:   int = 9,
) -> None:
    cphyto_manifest = Path(cphyto_manifest)
    par_dir         = Path(par_dir)
    manifest_path   = Path(manifest_path)
    par_dir.mkdir(parents=True, exist_ok=True)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)

    lon_min, lon_max, lat_min, lat_max = bbox
    lon_min -= BBOX_PAD; lon_max += BBOX_PAD
    lat_min -= BBOX_PAD; lat_max += BBOX_PAD
    print(f"Bbox subset: lon [{lon_min:.2f}, {lon_max:.2f}], lat [{lat_min:.2f}, {lat_max:.2f}]")

    dates = _unique_dates_from_cphyto(cphyto_manifest)
    print(f"Found {len(dates)} unique profile dates in {cphyto_manifest}")
    if not dates:
        pd.DataFrame(columns=["date", "filename", "path"]).to_csv(manifest_path, index=False)
        return

    res_str = f"{resolution_km}km"
    pattern = f"*.DAY.PAR.par.{res_str}.nc"

    # Pre-scan: split into already-cached and genuinely missing
    rows = []
    missing = []   # list of (date_str, expected_name, out_path)
    for d in dates:
        expected_name = f"AQUA_MODIS.{d.replace('-', '')}.L3m.DAY.PAR.par.{res_str}.subset.nc"
        out_path = par_dir / expected_name
        if out_path.exists():
            rows.append({"date": d, "filename": expected_name, "path": str(out_path)})
        else:
            missing.append((d, expected_name, out_path))

    print(f"{len(rows)}/{len(dates)} dates already cached, {len(missing)} to download")

    # Only authenticate + search when there are files to fetch
    if missing:
        # Force netrc strategy — the default 'all' strategy tries interactive prompts
        # first, which hangs/fails under Snakemake (no TTY).
        earthaccess.login(strategy="netrc")

        for i, (d, expected_name, out_path) in enumerate(missing, 1):
            try:
                granules = earthaccess.search_data(
                    short_name=SHORT_NAME, version=VERSION,
                    temporal=(d, d), granule_name=pattern,
                )
            except Exception as e:
                print(f"  [{i}/{len(missing)}] {d}: SEARCH FAILED ({e})")
                continue
            if not granules:
                print(f"  [{i}/{len(missing)}] {d}: no granule found")
                continue

            try:
                file_objs = earthaccess.open(granules)
                ok = _subset_and_save(file_objs[0], out_path,
                                      lon_min, lon_max, lat_min, lat_max)
            except Exception as e:
                print(f"  [{i}/{len(missing)}] {d}: OPEN/SUBSET FAILED ({e})")
                continue

            if ok and out_path.exists():
                print(f"  [{i}/{len(missing)}] {expected_name}: saved ({out_path.stat().st_size/1024:.0f} kB)")
                rows.append({"date": d, "filename": expected_name, "path": str(out_path)})
            else:
                print(f"  [{i}/{len(missing)}] {d}: bbox empty — skipped")

    pd.DataFrame(rows).to_csv(manifest_path, index=False)
    print(f"Manifest -> {manifest_path}  ({len(rows)}/{len(dates)} files available)")


def _snakemake_main():
    import traceback
    log_path = Path("data") / "_download_par_error.log"
    try:
        download_par(
            cphyto_manifest = Path(snakemake.input["cphyto_manifest"]),         # noqa: F821
            par_dir         = Path(snakemake.params["par_dir"]),                # noqa: F821
            manifest_path   = Path(snakemake.output["manifest"]),               # noqa: F821
            bbox            = _bbox_from_config(snakemake.config),              # noqa: F821
            resolution_km   = int(snakemake.config.get("par_resolution_km", 9)), # noqa: F821
        )
    except Exception:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        log_path.write_text(traceback.format_exc())
        raise


if __name__ == "__main__":
    if "snakemake" in globals():
        _snakemake_main()
    else:
        import yaml
        args = sys.argv[1:]
        cfg_path        = Path(args[0]) if len(args) > 0 else Path("config.yaml")
        cphyto_manifest = Path(args[1]) if len(args) > 1 else Path("data/NorthAtlantic_seas_comparison/intermediate/cphyto_profiles/cphyto_manifest.csv")
        par_dir         = Path(args[2]) if len(args) > 2 else Path("data/raw/modis_par")
        manifest        = Path(args[3]) if len(args) > 3 else Path("data/NorthAtlantic_seas_comparison/raw/par_download_manifest.csv")
        cfg = yaml.safe_load(cfg_path.read_text())
        download_par(cphyto_manifest, par_dir, manifest,
                     bbox=_bbox_from_config(cfg),
                     resolution_km=int(cfg.get("par_resolution_km", 9)))
