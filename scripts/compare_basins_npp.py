"""
compare_basins_npp.py
---------------------
Cross-basin NPP visualisation:
  - basins_npp_comparison.png : mean NPP time series overlaid for all basins
  - basins_npp_map.png        : all profile locations on a single map,
                                colored by NPP, with every basin polygon drawn.
"""


import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from shapely.geometry import Polygon, Point


INTEGRATION_COL = "npp_int_0_zeu"


def _in_polygon(lon, lat, poly_verts) -> np.ndarray:
    coords = [(float(x), float(y)) for x, y in poly_verts]
    if coords[0] != coords[-1]:
        coords.append(coords[0])
    poly = Polygon(coords)
    return np.array([poly.contains(Point(x, y)) for x, y in zip(lon, lat)])


def compare_basins_npp(
    integrated_csv: Path,
    ts_csv_by_basin: dict[str, Path],
    basin_polygons:  dict[str, list],
    out_ts_png:      Path,
    out_map_png:     Path,
) -> None:
    out_ts_png  = Path(out_ts_png)
    out_map_png = Path(out_map_png)
    out_ts_png.parent.mkdir(parents=True, exist_ok=True)

    # ── Time-series overlay ──────────────────────────────────────────────────
    colors = plt.cm.Dark2.colors
    fig, ax = plt.subplots(figsize=(11, 5))
    for i, (basin, csv) in enumerate(ts_csv_by_basin.items()):
        csv = Path(csv)
        if not csv.exists():
            continue
        ts = pd.read_csv(csv, parse_dates=["date_bin"])
        if ts.empty:
            continue
        c = colors[i % len(colors)]
        ax.fill_between(ts["date_bin"],
                        ts["npp_mean"] - ts["npp_sd"].fillna(0),
                        ts["npp_mean"] + ts["npp_sd"].fillna(0),
                        alpha=0.18, color=c)
        ax.plot(ts["date_bin"], ts["npp_mean"], color=c, lw=1.5, label=basin)
    ax.set_ylabel("NPP  (mg C m⁻² d⁻¹)")
    ax.set_xlabel("Date")
    ax.set_title("NPP comparison across basins  (mean ± 1σ across profiles per time bin)")
    ax.axhline(0, color="grey", lw=0.5, ls="--")
    ax.legend(loc="upper right", frameon=False)
    fig.autofmt_xdate()
    fig.tight_layout()
    fig.savefig(out_ts_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Time-series comparison -> {out_ts_png}")

    # ── All-basin map ────────────────────────────────────────────────────────
    npp = pd.read_csv(integrated_csv, parse_dates=["date"])
    fig, ax = plt.subplots(figsize=(9, 7))

    # Draw each basin polygon + annotate
    for i, (basin, poly) in enumerate(basin_polygons.items()):
        coords = [(float(x), float(y)) for x, y in poly]
        if coords[0] != coords[-1]:
            coords.append(coords[0])
        px = [c[0] for c in coords]
        py = [c[1] for c in coords]
        c = colors[i % len(colors)]
        ax.plot(px, py, color=c, lw=1.5, label=basin)

    # All profiles (in any basin) colored by NPP
    mask_any = np.zeros(len(npp), dtype=bool)
    for poly in basin_polygons.values():
        mask_any |= _in_polygon(npp["lon"], npp["lat"], poly)
    pts = npp.loc[mask_any]
    sc = ax.scatter(pts["lon"], pts["lat"], c=pts[INTEGRATION_COL],
                    cmap="viridis", s=25, edgecolor="k", linewidth=0.2)
    cb = fig.colorbar(sc, ax=ax, shrink=0.75)
    cb.set_label("NPP  (mg C m⁻² d⁻¹)")

    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_aspect("equal")
    ax.set_title(f"Profile NPP across basins  ({INTEGRATION_COL}, N={len(pts)})")
    ax.legend(loc="lower right", frameon=False)
    fig.tight_layout()
    fig.savefig(out_map_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Basin map -> {out_map_png}")


if __name__ == "__main__":
    if "snakemake" in globals():
        ts_csvs = {b: p for b, p in
                   zip(snakemake.params["basin_names"],                # noqa: F821
                       snakemake.input["ts_csvs"])}                    # noqa: F821
        polys = {b: snakemake.config["basins"][b]["polygon"]           # noqa: F821
                 for b in snakemake.params["basin_names"]}             # noqa: F821
        compare_basins_npp(
            integrated_csv = Path(snakemake.input["integrated_csv"]),  # noqa: F821
            ts_csv_by_basin = ts_csvs,
            basin_polygons  = polys,
            out_ts_png  = Path(snakemake.output["ts_fig"]),            # noqa: F821
            out_map_png = Path(snakemake.output["map_fig"]),           # noqa: F821
        )
    else:
        import yaml
        cfg = yaml.safe_load(Path("config.yaml").read_text())
        run = cfg["run_name"]
        basins = list(cfg["basins"].keys())
        ts_csvs = {b: Path(f"output/{run}/npp/{b}/npp_timeseries.csv") for b in basins}
        polys = {b: cfg["basins"][b]["polygon"] for b in basins}
        compare_basins_npp(
            integrated_csv  = Path(f"data/{run}/intermediate/npp/npp_integrated.csv"),
            ts_csv_by_basin = ts_csvs,
            basin_polygons  = polys,
            out_ts_png  = Path(f"output/{run}/npp/basins_npp_comparison.png"),
            out_map_png = Path(f"output/{run}/npp/basins_npp_map.png"),
        )
