configfile: "config.yaml"

RUN  = config["run_name"]
DATA = f"{config['data_dir']}/{RUN}"
OUT  = f"{config['output_dir']}/{RUN}"

rule all:
    input:
        f"{DATA}/intermediate/doxy_profiles/processing_manifest.csv"

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
    output:
        out_dir  = directory(f"{DATA}/intermediate/doxy_profiles"),
        manifest = f"{DATA}/intermediate/doxy_profiles/processing_manifest.csv"
    script:
        "scripts/process_sprof.R"