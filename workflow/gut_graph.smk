import os

if not os.path.exists('gut_graph_files'):
    os.mkdir('gut_graph_files')

configfile: "config/gut.yaml"

rule extract_metadata:
    input:
        config["metadata_file"]
    output:
        'gut_graph_files/metadata.pkl'
    conda:
        "envs/phage_clouds.yaml"
    script:
        "scripts/metadata_gut.py"
