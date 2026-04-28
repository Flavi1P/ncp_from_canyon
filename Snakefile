configfile: "config.yaml"

RUN    = config["run_name"]
DATA   = f"{config['data_dir']}/{RUN}"
OUT    = f"{config['output_dir']}/{RUN}"
RAW    = config.get("raw_dir", "data/raw")
SHARED = config.get("shared_dir", "data/shared")
BASINS = list(config["basins"].keys())
FLOATS = [str(w) for w in config.get("floats", [])]

# ── Output targets ─────────────────────────────────────────────────────────────
_ncp_csv  = f"{OUT}/ncp/{{basin}}/ncp_results.csv"
_ncp_png  = f"{OUT}/ncp/{{basin}}/ncp_sensitivity.png"
_unc_csv  = f"{OUT}/ncp/{{basin}}/ncp_uncertainty.csv"
_unc_png  = f"{OUT}/ncp/{{basin}}/ncp_uncertainty.png"

_all_outputs = (
    expand(_ncp_csv, basin=BASINS) +
    expand(_ncp_png, basin=BASINS)
)
if config.get("uncertainty", False):
    _all_outputs += (
        expand(_unc_csv, basin=BASINS) +
        expand(_unc_png, basin=BASINS) +
        [f"{OUT}/ncp/basins_comparison.png",
         f"{OUT}/summary/_done.txt"]
    )
if FLOATS:
    _float_ncp_csv      = f"{OUT}/ncp_float/{{wmo}}/ncp_float.csv"
    _float_ncp_png      = f"{OUT}/ncp_float/{{wmo}}/ncp_float.png"
    _float_transect_png = f"{OUT}/ncp_float/{{wmo}}/nitrate_transect.png"
    _all_outputs += (
        expand(_float_ncp_csv,      wmo=FLOATS) +
        expand(_float_ncp_png,      wmo=FLOATS) +
        expand(_float_transect_png, wmo=FLOATS) +
        [f"{OUT}/ncp_float/floats_comparison.png"]
    )

if config.get("compute_npp", False):
    _npp_ts_csv  = f"{OUT}/npp/{{basin}}/npp_timeseries.csv"
    _npp_ts_png  = f"{OUT}/npp/{{basin}}/npp_vs_ncp.png"
    _npp_map_png = f"{OUT}/npp/{{basin}}/npp_map.png"
    _all_outputs += (
        expand(_npp_ts_csv,  basin=BASINS) +
        expand(_npp_ts_png,  basin=BASINS) +
        expand(_npp_map_png, basin=BASINS) +
        [f"{OUT}/npp/basins_npp_comparison.png",
         f"{OUT}/npp/basins_npp_map.png"]
    )

rule all:
    input: _all_outputs

# ── Upstream rules (basin-independent) ────────────────────────────────────────
rule download_sprof:
    output:
        manifest = f"{DATA}/raw/download_manifest.csv",
        wmo_list = f"{DATA}/raw/wmo_list.txt"
    script:
        "scripts/download_sprof.py"

rule process_sprof:
    input:
        manifest  = f"{DATA}/raw/download_manifest.csv",
        wmo_list  = f"{DATA}/raw/wmo_list.txt"
    params:
        sprof_dir  = RAW,
        shared_dir = SHARED
    output:
        manifest = f"{DATA}/intermediate/doxy_profiles/processing_manifest.csv"
    script:
        "scripts/process_sprof.R"

rule process_cphyto:
    input:
        manifest  = f"{DATA}/raw/download_manifest.csv",
        wmo_list  = f"{DATA}/raw/wmo_list.txt"
    params:
        sprof_dir  = RAW,
        shared_dir = SHARED
    output:
        manifest = f"{DATA}/intermediate/cphyto_profiles/cphyto_manifest.csv"
    script:
        "scripts/process_cphyto.py"

if config.get("compute_npp", False):
    rule download_par:
        input:
            cphyto_manifest = f"{DATA}/intermediate/cphyto_profiles/cphyto_manifest.csv"
        params:
            par_dir = f"{RAW}/modis_par"
        output:
            manifest = f"{DATA}/raw/par_download_manifest.csv"
        script:
            "scripts/download_par.py"

    rule match_par:
        input:
            cphyto_manifest = f"{DATA}/intermediate/cphyto_profiles/cphyto_manifest.csv",
            par_manifest    = f"{DATA}/raw/par_download_manifest.csv"
        output:
            matched_csv = f"{DATA}/intermediate/par_matched/par_matched.csv"
        script:
            "scripts/match_par.py"

    rule compute_npp:
        input:
            cphyto_manifest = f"{DATA}/intermediate/cphyto_profiles/cphyto_manifest.csv",
            par_matched_csv = f"{DATA}/intermediate/par_matched/par_matched.csv"
        output:
            profiles_csv   = f"{DATA}/intermediate/npp/npp_profiles.csv",
            integrated_csv = f"{DATA}/intermediate/npp/npp_integrated.csv"
        script:
            "scripts/compute_npp.py"

    rule npp_timeseries:
        input:
            integrated_csv = f"{DATA}/intermediate/npp/npp_integrated.csv",
            ncp_csv        = f"{OUT}/ncp/{{basin}}/ncp_results.csv"
        params:
            basin_name    = lambda wc: wc.basin,
            basin_polygon = lambda wc: config["basins"][wc.basin]["polygon"]
        output:
            ts_csv   = f"{OUT}/npp/{{basin}}/npp_timeseries.csv",
            plot_png = f"{OUT}/npp/{{basin}}/npp_vs_ncp.png",
            map_png  = f"{OUT}/npp/{{basin}}/npp_map.png"
        wildcard_constraints:
            basin = "|".join(BASINS)
        script:
            "scripts/npp_timeseries.py"

    rule compare_basins_npp:
        input:
            integrated_csv = f"{DATA}/intermediate/npp/npp_integrated.csv",
            ts_csvs        = expand(f"{OUT}/npp/{{basin}}/npp_timeseries.csv", basin=BASINS)
        params:
            basin_names = BASINS
        output:
            ts_fig  = f"{OUT}/npp/basins_npp_comparison.png",
            map_fig = f"{OUT}/npp/basins_npp_map.png"
        script:
            "scripts/compare_basins_npp.py"

rule predict_nitrate:
    input:
        manifest = f"{DATA}/intermediate/doxy_profiles/processing_manifest.csv"
    params:
        shared_dir = SHARED
    output:
        manifest = f"{DATA}/intermediate/nitrate_profiles/nitrate_manifest.csv"
    threads: workflow.cores
    script:
        "scripts/predict_nitrate.R"

rule merge_floats:
    input:
        manifest = f"{DATA}/intermediate/nitrate_profiles/nitrate_manifest.csv"
    output:
        merged_csv = f"{DATA}/intermediate/merged/merged_ncp.csv",
        map_png    = f"{OUT}/float_map.png"
    script:
        "scripts/merge_nitrate_files.R"

# ── Per-basin NCP rules ────────────────────────────────────────────────────────
rule compute_ncp:
    input:
        merged_csv = f"{DATA}/intermediate/merged/merged_ncp.csv"
    params:
        basin_name    = lambda wc: wc.basin,
        basin_polygon = lambda wc: config["basins"][wc.basin]["polygon"]
    output:
        out_dir         = directory(f"{OUT}/ncp/{{basin}}"),
        ncp_results     = f"{OUT}/ncp/{{basin}}/ncp_results.csv",
        sensitivity_png = f"{OUT}/ncp/{{basin}}/ncp_sensitivity.png"
    threads: workflow.cores
    resources:
        mem_mb = 16000
    wildcard_constraints:
        basin = "|".join(BASINS)
    script:
        "scripts/compute_ncp.R"

if config.get("uncertainty", False):
    rule compute_ncp_uncertainty:
        input:
            merged_csv = f"{DATA}/intermediate/merged/merged_ncp.csv",
            ncp_done   = f"{OUT}/ncp/{{basin}}/ncp_results.csv"
        params:
            basin_name    = lambda wc: wc.basin,
            basin_polygon = lambda wc: config["basins"][wc.basin]["polygon"]
        output:
            results_csv = f"{OUT}/ncp/{{basin}}/ncp_uncertainty.csv",
            plot_png    = f"{OUT}/ncp/{{basin}}/ncp_uncertainty.png"
        threads: workflow.cores
        resources:
            mem_mb = 16000
        wildcard_constraints:
            basin = "|".join(BASINS)
        script:
            "scripts/compute_ncp_uncertainty.R"

    rule compare_basins:
        input:
            expand(_unc_csv, basin=BASINS)
        params:
            basin_names = BASINS
        output:
            fig = f"{OUT}/ncp/basins_comparison.png"
        script:
            "scripts/compare_basins.R"

    rule summary_plots:
        input:
            basin_unc_csvs = expand(_unc_csv, basin=BASINS),
            basin_res_csvs = expand(_ncp_csv, basin=BASINS),
            float_csvs     = (expand(f"{OUT}/ncp_float/{{wmo}}/ncp_float.csv",
                                     wmo=FLOATS) if FLOATS else [])
        params:
            basin_names = BASINS,
            float_wmos  = FLOATS
        output:
            out_dir = directory(f"{OUT}/summary"),
            marker  = f"{OUT}/summary/_done.txt"
        script:
            "scripts/summary_plots.R"

# ── Per-float NCP rules ────────────────────────────────────────────────────────
if FLOATS:
    rule compute_ncp_float:
        input:
            merged_csv = f"{DATA}/intermediate/merged/merged_ncp.csv"
        params:
            wmo               = lambda wc: wc.wmo,
            zeu_default       = config.get("zeu_default", 40),
            float_time_window = config.get("float_time_window", None)
        output:
            results_csv  = f"{OUT}/ncp_float/{{wmo}}/ncp_float.csv",
            plot_png     = f"{OUT}/ncp_float/{{wmo}}/ncp_float.png",
            transect_png = f"{OUT}/ncp_float/{{wmo}}/nitrate_transect.png"
        threads: workflow.cores
        wildcard_constraints:
            wmo = "|".join(FLOATS)
        script:
            "scripts/compute_ncp_float.R"

    rule compare_floats:
        input:
            float_csvs = expand(f"{OUT}/ncp_float/{{wmo}}/ncp_float.csv", wmo=FLOATS)
        params:
            wmo_list = FLOATS
        output:
            fig = f"{OUT}/ncp_float/floats_comparison.png"
        script:
            "scripts/compare_floats_ncp.R"
