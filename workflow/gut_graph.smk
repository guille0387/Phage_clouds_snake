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

rule split_dataset:
    input:
        fasta=config["sequence_file"],
        meta=config["metadata_file"]
    output:
        high=f'gut_graph_files/{config["sequence_file"].split(".")[0]}_high_quality.fa.gz',
        frag=f'gut_graph_files/{config["sequence_file"].split(".")[0]}_genome_fragment.fa.gz'
    conda:
        "envs/phage_clouds.yaml"
    script:
        "scripts/split_gut.py"
