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

rule ref_sketch:
    input:
        high=f'gut_graph_files/{config["sequence_file"].split(".")[0]}_high_quality.fa.gz',
        frag=f'gut_graph_files/{config["sequence_file"].split(".")[0]}_genome_fragment.fa.gz'
    output:
        high_sketch=f'gut_graph_files/{config["sequence_file"].split(".")[0]}_high_quality.msh',
        frag_sketch=f'gut_graph_files/{config["sequence_file"].split(".")[0]}_genome_fragment.msh'
    threads: 10
    run:
        high_prefix = os.path.splitext(output.high_sketch)[0]
        frag_prefix = os.path.splitext(output.frag_sketch)[0]
        shell(f'mash sketch -p 10 -i -k 15 -s 24000 -o {high_prefix} {input.high}')
        shell(f'mash sketch -p 10 -i -k 15 -s 24000 -o {frag_prefix} {input.frag}')

rule mash_dist:
    input:
        fasta='gut_graph_files/{subset}.fa.gz',
        sketch='gut_graph_files{subset}.msh'
    output:
        'gut_graph_files/{subset}.tbl'
    threads: 10
    run:
        print('Calculating all vs all mash distance tables (p-value reporting threshold 1e-10)...')
        shell("mash dist -p 10 -i -t -v 1e-10 {input.sketch} {input.fasta} > {output}")
        print('Done creating mash distance tables!')
