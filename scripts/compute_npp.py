"""
compute_npp.py
--------------
Run the Carbon-based Productivity Model (CbPM, Westberry et al. 2008, as
implemented in scripts/utils/cbpm_argo.py) on every cphyto profile whose
surface PAR was matched in match_par.py.

Writes two CSVs:
  - npp_profiles.csv    depth-resolved (0..199 m), one row per (profile, depth)
  - npp_integrated.csv  one row per profile: surface PAR, zeu, and depth-integrated
                        NPP (mg C m-2 d-1) over 0-200 m and over 0-zeu.
"""


import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent / "utils"))
from cbpm_argo import cbpm_argo   # noqa: E402

# cbpm_argo internally sizes all output arrays to 200 — pass exactly 200 depths.
N_DEPTHS = 200


def _run_one(chl_z: np.ndarray, cphyto_z: np.ndarray,
             par_surface: float, d: pd.Timestamp, lat: float) -> dict | None:
    """Call cbpm_argo; return the depth-resolved arrays + integrated values."""
    if not np.isfinite(par_surface) or par_surface <= 0:
        return None
    try:
        pp_z, mu_z, par_z, prcnt_z, _, _, mzeu = cbpm_argo(
            chl_z.copy(), cphyto_z.copy(),
            float(par_surface),
            int(d.year), int(d.month), int(d.day),
            float(lat),
        )
    except Exception:
        return None

    # pp_z is volumetric NPP (mg C m-3 d-1). Integrate with a 1 m step.
    npp_0_200 = float(np.nansum(pp_z))
    if mzeu is None or not np.isfinite(mzeu):
        npp_0_zeu = np.nan
        zeu_out = np.nan
    else:
        zeu_out = int(mzeu)
        npp_0_zeu = float(np.nansum(pp_z[: zeu_out + 1]))

    return {
        "pp_z":     pp_z,
        "mu_z":     mu_z,
        "par_z":    par_z,
        "prcnt_z":  prcnt_z,
        "zeu_m":    zeu_out,
        "npp_int_0_200": npp_0_200,
        "npp_int_0_zeu": npp_0_zeu,
    }


def compute_npp(cphyto_dir: Path, par_matched_csv: Path,
                out_profiles_csv: Path, out_integrated_csv: Path) -> None:
    cphyto_dir         = Path(cphyto_dir)
    par_matched_csv    = Path(par_matched_csv)
    out_profiles_csv   = Path(out_profiles_csv)
    out_integrated_csv = Path(out_integrated_csv)
    out_profiles_csv.parent.mkdir(parents=True, exist_ok=True)

    par_df = pd.read_csv(par_matched_csv, parse_dates=["date"])
    par_df["float_wmo"] = par_df["float_wmo"].astype(str)
    par_key = par_df.set_index(["float_wmo", "prof_number"])["par"]

    profile_rows, integrated_rows = [], []
    skipped_no_par = skipped_bad_cbpm = 0

    for csv in sorted(cphyto_dir.glob("argo_*_cphyto.csv")):
        df = pd.read_csv(csv, parse_dates=["date"])
        df["float_wmo"] = df["float_wmo"].astype(str)

        for (wmo, prof), g in df.groupby(["float_wmo", "prof_number"]):
            g = g.sort_values("depth")
            # Take the first 200 rows (depths 0..199 m) — cbpm_argo hard-codes 200.
            if len(g) < N_DEPTHS:
                continue
            g = g.iloc[:N_DEPTHS]
            chl_z    = g["chla"].to_numpy(dtype=float)
            cphyto_z = g["cphyto"].to_numpy(dtype=float)
            lat = float(g["lat"].iloc[0])
            lon = float(g["lon"].iloc[0])
            d   = g["date"].iloc[0]

            par_surface = par_key.get((wmo, prof), np.nan)
            if not np.isfinite(par_surface):
                skipped_no_par += 1
                continue

            res = _run_one(chl_z, cphyto_z, par_surface, d, lat)
            if res is None:
                skipped_bad_cbpm += 1
                continue

            depth_arr = np.arange(N_DEPTHS)
            profile_rows.append(pd.DataFrame({
                "float_wmo":   wmo,
                "prof_number": prof,
                "date":        d,
                "lon":         lon,
                "lat":         lat,
                "depth":       depth_arr,
                "npp":         res["pp_z"],
                "mu":          res["mu_z"],
                "par_z":       res["par_z"],
                "prcnt_z":     res["prcnt_z"],
            }))
            integrated_rows.append({
                "float_wmo":      wmo,
                "prof_number":    prof,
                "date":           d,
                "lon":            lon,
                "lat":            lat,
                "par_surface":    par_surface,
                "zeu_m":          res["zeu_m"],
                "npp_int_0_200":  res["npp_int_0_200"],
                "npp_int_0_zeu":  res["npp_int_0_zeu"],
            })

    if profile_rows:
        pd.concat(profile_rows, ignore_index=True).to_csv(out_profiles_csv, index=False)
    else:
        pd.DataFrame(columns=["float_wmo", "prof_number", "date", "lon", "lat",
                              "depth", "npp", "mu", "par_z", "prcnt_z"]
                    ).to_csv(out_profiles_csv, index=False)

    pd.DataFrame(integrated_rows).to_csv(out_integrated_csv, index=False)
    print(f"NPP computed for {len(integrated_rows)} profiles  "
          f"(skipped: {skipped_no_par} no PAR, {skipped_bad_cbpm} cbpm error)")
    print(f"  profiles   -> {out_profiles_csv}")
    print(f"  integrated -> {out_integrated_csv}")


if __name__ == "__main__":
    if "snakemake" in globals():
        compute_npp(
            cphyto_dir         = Path(snakemake.input["cphyto_dir"]),         # noqa: F821
            par_matched_csv    = Path(snakemake.input["par_matched_csv"]),    # noqa: F821
            out_profiles_csv   = Path(snakemake.output["profiles_csv"]),      # noqa: F821
            out_integrated_csv = Path(snakemake.output["integrated_csv"]),    # noqa: F821
        )
    else:
        args = sys.argv[1:]
        base = Path("data/NorthAtlantic_seas_comparison")
        cphyto_dir      = Path(args[0]) if len(args) > 0 else base / "intermediate/cphyto_profiles"
        par_matched_csv = Path(args[1]) if len(args) > 1 else base / "intermediate/par_matched/par_matched.csv"
        out_profiles    = Path(args[2]) if len(args) > 2 else base / "intermediate/npp/npp_profiles.csv"
        out_integrated  = Path(args[3]) if len(args) > 3 else base / "intermediate/npp/npp_integrated.csv"
        compute_npp(cphyto_dir, par_matched_csv, out_profiles, out_integrated)
