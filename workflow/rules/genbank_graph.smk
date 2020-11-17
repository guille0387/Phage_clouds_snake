import os

if not os.path.exists('genbank_graph_files'):
    os.mkdir('genbank_graph_files')

configfile: "config/genbank.yaml"

rule all:
    input:
        f'genbank_graph_files/{config["host"]}_{config["threshold"]}_clouds.html'

rule create_graph_files:
    input:
        plaingraph="genbank_graph_files/plain_graph.nxpkl.gz",
        annotgraph="genbank_graph_files/annotated_graph.nxpkl.gz"

rule get_fasta_metadata:
    input:
        config["input_dir"]
    output:
        fa="genbank_graph_files/all_records.fna",
        meta="genbank_graph_files/metadata.pkl"
    conda:
        "envs/phage_clouds.yaml"
    script:
        "scripts/extract_from_gbk.py"

rule ref_sketch:
    input:
        "genbank_graph_files/all_records.fna"
    output:
        "genbank_graph_files/all_records_k15_s24000.msh"
    threads: 10
    run:
        print('mash sketching parameters: kmer_size = 15, sketch_size = 24000')
        shell("mash sketch -p 10 -i -o {output} -k 15 -s 24000 {input}")

rule mash_dist:
    input:
        "genbank_graph_files/all_records_k15_s24000.msh",
        "genbank_graph_files/all_records.fna"
    output:
        "genbank_graph_files/mash_dist_k15_s24000.tbl"
    threads: 10
    run:
        print('Calculating all vs all mash distance table (p-value reporting threshold 1e-10)...')
        shell("mash dist -p 10 -i -t -v 1e-10 {input} > {output}")
        print(f'Done creating {output}')

rule get_network:
    input:
        mashtbl="genbank_graph_files/mash_dist_k15_s24000.tbl",
        annotdict="genbank_graph_files/metadata.pkl"
    output:
        plaingraph="genbank_graph_files/plain_graph.nxpkl.gz",
        annotgraph="genbank_graph_files/annotated_graph.nxpkl.gz"
    conda:
        "envs/phage_clouds.yaml"
    script:
        "scripts/generate_graph.py"

rule draw_host_centric:
    input:
        "genbank_graph_files/annotated_graph.nxpkl.gz"
    output:
        f'genbank_graph_files/{config["host"]}_{config["threshold"]}_clouds.html'
    conda:
        "envs/phage_clouds.yaml"
    params:
        host=config["host"],
        thres=config["threshold"]
    script:
        "scripts/host_centric_graph.py"
