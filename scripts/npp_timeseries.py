"""
npp_timeseries.py
-----------------
Bin NPP (depth-integrated, per profile) into time bins within a basin and
produce a mean/sd time series. Overlay against the corresponding NCP series.
"""


import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from shapely.geometry import Polygon, Point


INTEGRATION_COL_DEFAULT = "npp_int_0_zeu"   # override via CLI if needed


def _in_polygon(lon: pd.Series, lat: pd.Series, poly_verts) -> np.ndarray:
    coords = [(float(x), float(y)) for x, y in poly_verts]
    if coords[0] != coords[-1]:
        coords.append(coords[0])
    poly = Polygon(coords)
    return np.array([poly.contains(Point(x, y)) for x, y in zip(lon, lat)])


def _parse_time_step(ts: str) -> pd.Timedelta:
    # Accept '15 days', '10 days', '20 days', etc.
    return pd.Timedelta(ts.replace(" ", ""))


def _plot_basin_map(npp: pd.DataFrame, basin_polygon, basin_name: str,
                    integration_col: str, out_path: Path) -> None:
    """Scatter of profile locations within a basin, colored by NPP."""
    coords = [(float(x), float(y)) for x, y in basin_polygon]
    if coords[0] != coords[-1]:
        coords.append(coords[0])
    px = [c[0] for c in coords]
    py = [c[1] for c in coords]

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.plot(px, py, color="black", lw=1.2)
    sc = ax.scatter(npp["lon"], npp["lat"], c=npp[integration_col],
                    cmap="viridis", s=40, edgecolor="k", linewidth=0.3)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_aspect("equal")
    ax.set_title(f"{basin_name} — profile NPP ({integration_col})\n"
                 f"{len(npp)} profiles")
    cb = fig.colorbar(sc, ax=ax, shrink=0.8)
    cb.set_label("NPP  (mg C m⁻² d⁻¹)")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def npp_timeseries(
    integrated_csv: Path,
    ncp_csv:        Path,
    out_csv:        Path,
    out_png:        Path,
    out_map_png:    Path,
    basin_name:     str,
    basin_polygon,
    time_step:      str,
    integration_col: str = INTEGRATION_COL_DEFAULT,
) -> None:
    integrated_csv = Path(integrated_csv)
    out_csv        = Path(out_csv)
    out_png        = Path(out_png)
    out_map_png    = Path(out_map_png)
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    npp = pd.read_csv(integrated_csv, parse_dates=["date"])
    in_basin = _in_polygon(npp["lon"], npp["lat"], basin_polygon)
    npp = npp.loc[in_basin].copy()
    print(f"{basin_name}: {len(npp)} profiles after basin filter")

    if npp.empty:
        pd.DataFrame(columns=["date_bin", "npp_mean", "npp_sd", "n_profiles"]).to_csv(out_csv, index=False)
        for p in (out_png, out_map_png):
            fig, ax = plt.subplots(figsize=(6, 4))
            ax.set_title(f"{basin_name}: no profiles in basin")
            fig.savefig(p, dpi=150, bbox_inches="tight")
            plt.close(fig)
        return

    _plot_basin_map(npp, basin_polygon, basin_name, integration_col, out_map_png)

    # Time-bin NPP — floor dates to the start of each bin.
    dt = _parse_time_step(time_step)
    origin = npp["date"].min().normalize()
    bin_idx = ((npp["date"] - origin) // dt)
    npp["date_bin"] = origin + bin_idx * dt

    ts = (npp.groupby("date_bin")[integration_col]
              .agg(npp_mean="mean", npp_sd="std", n_profiles="count")
              .reset_index())
    ts.to_csv(out_csv, index=False)
    print(f"Time series written -> {out_csv}  ({len(ts)} bins)")

    # ── Plot — NPP mean ± sd, with NCP mean line overlaid if available ───────
    fig, ax1 = plt.subplots(figsize=(10, 5))
    ax1.fill_between(ts["date_bin"],
                     ts["npp_mean"] - ts["npp_sd"],
                     ts["npp_mean"] + ts["npp_sd"],
                     alpha=0.25, color="tab:green", label="NPP ±1σ")
    ax1.plot(ts["date_bin"], ts["npp_mean"], color="tab:green", lw=1.5, label="NPP mean")
    ax1.set_ylabel("NPP  (mg C m⁻² d⁻¹)", color="tab:green")
    ax1.tick_params(axis="y", labelcolor="tab:green")
    ax1.set_xlabel("Date")

    if Path(ncp_csv).exists():
        ncp = pd.read_csv(ncp_csv, parse_dates=["date_grid"])
        # NCP is in mmol C m-2 d-1; convert to mg C m-2 d-1 (× 12.011)
        ncp["NCP_mgC"] = ncp["NCP"] * 12.011
        ax2 = ax1.twinx()
        ax2.plot(ncp["date_grid"], ncp["NCP_mgC"], color="tab:blue", lw=1.2, label="NCP")
        ax2.set_ylabel("NCP  (mg C m⁻² d⁻¹)", color="tab:blue")
        ax2.tick_params(axis="y", labelcolor="tab:blue")

    ax1.set_title(f"{basin_name} — NPP ({integration_col}) vs NCP  |  {time_step} bins")
    fig.autofmt_xdate()
    fig.tight_layout()
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    print(f"Plot written -> {out_png}")


if __name__ == "__main__":
    if "snakemake" in globals():
        time_steps = snakemake.config.get("ncp_time_steps", ["15 days"])   # noqa: F821
        npp_timeseries(
            integrated_csv  = Path(snakemake.input["integrated_csv"]),     # noqa: F821
            ncp_csv         = Path(snakemake.input["ncp_csv"]),            # noqa: F821
            out_csv         = Path(snakemake.output["ts_csv"]),            # noqa: F821
            out_png         = Path(snakemake.output["plot_png"]),          # noqa: F821
            out_map_png     = Path(snakemake.output["map_png"]),           # noqa: F821
            basin_name      = snakemake.params["basin_name"],              # noqa: F821
            basin_polygon   = snakemake.params["basin_polygon"],           # noqa: F821
            time_step       = time_steps[0],
        )
    else:
        import yaml
        args = sys.argv[1:]
        cfg = yaml.safe_load(Path("config.yaml").read_text())
        basin = args[0] if len(args) > 0 else next(iter(cfg["basins"]))
        run   = cfg["run_name"]
        out_dir = Path(f"output/{run}/npp/{basin}")
        npp_timeseries(
            integrated_csv = Path(f"data/{run}/intermediate/npp/npp_integrated.csv"),
            ncp_csv        = Path(f"output/{run}/ncp/{basin}/ncp_results.csv"),
            out_csv        = out_dir / "npp_timeseries.csv",
            out_png        = out_dir / "npp_vs_ncp.png",
            out_map_png    = out_dir / "npp_map.png",
            basin_name     = basin,
            basin_polygon  = cfg["basins"][basin]["polygon"],
            time_step      = cfg["ncp_time_steps"][0],
        )
