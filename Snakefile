configfile: "config.yaml"

RUN    = config["run_name"]
DATA   = f"{config['data_dir']}/{RUN}"
OUT    = f"{config['output_dir']}/{RUN}"
RAW    = config.get("raw_dir", "data/raw")
BASINS = list(config["basins"].keys())

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
        [f"{OUT}/ncp/basins_comparison.png"]
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
        sprof_dir = RAW
    output:
        out_dir  = directory(f"{DATA}/intermediate/doxy_profiles"),
        manifest = f"{DATA}/intermediate/doxy_profiles/processing_manifest.csv"
    script:
        "scripts/process_sprof.R"

rule predict_nitrate:
    input:
        in_dir   = f"{DATA}/intermediate/doxy_profiles",
        manifest = f"{DATA}/intermediate/doxy_profiles/processing_manifest.csv"
    output:
        out_dir  = directory(f"{DATA}/intermediate/nitrate_profiles"),
        manifest = f"{DATA}/intermediate/nitrate_profiles/nitrate_manifest.csv"
    threads: workflow.cores
    script:
        "scripts/predict_nitrate.R"

rule merge_floats:
    input:
        in_dir   = f"{DATA}/intermediate/nitrate_profiles",
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
    resources:
        mem_mb = 8000
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
        resources:
            mem_mb = 8000
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
