configfile: "config.yaml"

RUN  = config["run_name"]
DATA = f"{config['data_dir']}/{RUN}"
OUT  = f"{config['output_dir']}/{RUN}"
RAW  = config.get("raw_dir", "data/raw")

_all_outputs = [
    f"{OUT}/ncp/ncp_results.csv",
    f"{OUT}/ncp/ncp_sensitivity.png",
]
if config.get("uncertainty", False):
    _all_outputs += [
        f"{OUT}/ncp/ncp_uncertainty.csv",
        f"{OUT}/ncp/ncp_uncertainty.png",
    ]

rule all:
    input: _all_outputs

rule download_sprof:
    output:
        manifest = f"{DATA}/raw/download_manifest.csv",
        wmo_list = f"{DATA}/raw/wmo_list.txt"
    script:
        "scripts/download_sprof.py"

rule process_sprof:
    input:
        manifest = f"{DATA}/raw/download_manifest.csv",
        wmo_list = f"{DATA}/raw/wmo_list.txt"
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

rule compute_ncp:
    input:
        merged_csv = f"{DATA}/intermediate/merged/merged_ncp.csv"
    output:
        out_dir        = directory(f"{OUT}/ncp"),
        ncp_results    = f"{OUT}/ncp/ncp_results.csv",
        sensitivity_png= f"{OUT}/ncp/ncp_sensitivity.png"
    script:
        "scripts/compute_ncp.R"

if config.get("uncertainty", False):
    rule compute_ncp_uncertainty:
        input:
            merged_csv = f"{DATA}/intermediate/merged/merged_ncp.csv",
            ncp_done   = f"{OUT}/ncp/ncp_results.csv"
        output:
            results_csv = f"{OUT}/ncp/ncp_uncertainty.csv",
            plot_png    = f"{OUT}/ncp/ncp_uncertainty.png"
        script:
            "scripts/compute_ncp_uncertainty.R"