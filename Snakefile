configfile: "config.yaml"

RUN  = config["run_name"]
DATA = f"{config['data_dir']}/{RUN}"
OUT  = f"{config['output_dir']}/{RUN}"

rule all:
    input:
        f"{DATA}/intermediate/merged/merged_ncp.csv",
        f"{OUT}/float_map.png"

rule download_sprof:
    output:
        manifest = f"{DATA}/raw/download_manifest.csv",
        wmo_list = f"{DATA}/raw/wmo_list.txt"
    script:
        "scripts/download_sprof.py"

rule process_sprof:
    input:
        sprof_dir = f"{DATA}/raw",
        manifest = f"{DATA}/raw/download_manifest.csv",
        wmo_list = f"{DATA}/raw/wmo_list.txt"
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