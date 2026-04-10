"""
download_sprof.py
-----------------
Download BGC-Argo Sprof files for floats with DOXY in a given
bounding box using the argo_gdac utility.
"""

from pathlib import Path
import sys
import pandas as pd
import os

# Make utils importable regardless of working directory
sys.path.insert(0, str(Path(__file__).parent / "utils"))
from float_download import argo_gdac


def download_sprof_files(
    lon_min:      float,
    lon_max:      float,
    lat_min:      float,
    lat_max:      float,
    sprof_dir:    Path,
    manifest_dir: Path,
    date_start = None,
    date_end   = None,
) -> None:
    """Download Sprof .nc files to a shared sprof_dir, write per-run
    manifest and wmo_list to manifest_dir."""

    from datetime import datetime

    sprof_dir    = Path(sprof_dir)
    manifest_dir = Path(manifest_dir)
    sprof_dir.mkdir(parents=True, exist_ok=True)
    manifest_dir.mkdir(parents=True, exist_ok=True)

    # Parse dates if provided as strings
    fmt = "%Y-%m-%d"
    start = datetime.strptime(date_start, fmt) if date_start else None
    end   = datetime.strptime(date_end,   fmt) if date_end   else None

    print(f"Querying GDAC index and downloading Sprof files → {sprof_dir}")
    wmoids, gdac_index, downloaded_filenames = argo_gdac(
        save_to            = str(sprof_dir) + "/",
        lat_range          = [lat_min, lat_max],
        lon_range          = [lon_min, lon_max],
        start_date         = start,
        end_date           = end,
        sensors            = "DOXY",
        overwrite_index    = False,   # keep index if already downloaded
        overwrite_profiles = False,   # skip floats already on disk
        verbose            = True,
    )

    print(f"\nFound {len(wmoids)} floats: {list(wmoids)}")
    print(f"Downloaded {len(downloaded_filenames)} Sprof files.")

    # Save WMO list and manifest to the per-run directory
    wmo_path = manifest_dir / "wmo_list.txt"
    wmo_path.write_text("\n".join(map(str, sorted(wmoids))))
    print(f"WMO list → {wmo_path}")

    manifest = pd.DataFrame({
        "wmo":      sorted(wmoids),
        "filename": [f"{w}_Sprof.nc" for w in sorted(wmoids)],
        "path":     [str(sprof_dir / f"{w}_Sprof.nc") for w in sorted(wmoids)],
    })
    manifest.to_csv(manifest_dir / "download_manifest.csv", index=False)
    print(f"Manifest → {manifest_dir / 'download_manifest.csv'}")


# ── Snakemake entry point ──────────────────────────────────────────────────────
if __name__ == "__main__":

    def main():
        """CLI entry point"""
        import yaml
        if len(sys.argv) < 2:
            print("Usage: python download_sprof.py config.yaml")
            sys.exit(1)
        cfg = yaml.safe_load(Path(sys.argv[1]).read_text())
        download_sprof_files(
            lon_min      = cfg["lon_min"],
            lon_max      = cfg["lon_max"],
            lat_min      = cfg["lat_min"],
            lat_max      = cfg["lat_max"],
            sprof_dir    = Path(cfg.get("raw_dir", "data/raw")),
            manifest_dir = Path(cfg["data_dir"]) / cfg["run_name"] / "raw",
            date_start   = cfg.get("date_start"),
            date_end     = cfg.get("date_end"),
        )

    if "snakemake" in globals():
        # Running inside Snakemake pipeline
        sys.path.insert(0, str(Path(__file__).parent / "utils"))
        download_sprof_files(
            lon_min      = snakemake.config["lon_min"],
            lon_max      = snakemake.config["lon_max"],
            lat_min      = snakemake.config["lat_min"],
            lat_max      = snakemake.config["lat_max"],
            sprof_dir    = Path(snakemake.config.get("raw_dir", "data/raw")),
            manifest_dir = Path(snakemake.output.manifest).parent,
            date_start   = snakemake.config.get("date_start"),
            date_end     = snakemake.config.get("date_end"),
        )
    else:
        # Running from CLI
        main()