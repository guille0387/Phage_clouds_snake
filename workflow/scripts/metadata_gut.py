#!/usr/bin/env python

import os, re, gzip, glob, pickle, ast, sys, csv
from collections import defaultdict
from ete3 import NCBITaxa
from tqdm import tqdm


def extract_meta(file_path, sep = '\t'):
    if os.path.splitext(file_path)[1] == '.gz':
        input_file = gzip.open(file_path, 'rt', newline = '')
    else:
        input_file = open(file_path, newline = '')
    annot_dir = defaultdict(dict)
    var_names = ['source', 'GPD_VC', 'genome_size', 'checkV_MIUViG', 'checkV_completion', 'checkV_host_cont']
    tax_ranks = ['genus', 'family', 'order', 'class', 'phylum']
    input_reader = csv.reader(input_file, delimiter = sep)
    next(input_reader)
    for line in tqdm(input_reader, desc = 'Extracting metadata from entries in GPD'):
        line_metadata = [*line[1:4], line[10], line[11], line[14]]
        annot_dir[line[0]].update(dict(zip(var_names, line_metadata)))
        host_lins = line[6].split(',')
        lowest_tax = set()
        tax_dir = defaultdict(list)
        for elem in host_lins:
            for taxon in reversed(elem.split('/')):
                if len(ncbi.get_name_translator([taxon])) != 0:
                    lowest_tax.add(ncbi.get_name_translator([taxon])[taxon][0])
                    break
        lowest_tax_lins = [ncbi.get_lineage(taxid) for taxid in lowest_tax]
        lowest_tax_lins_dict = [{ncbi.get_rank([x])[x]:ncbi.get_taxid_translator([x])[x] for x in lin} for lin in lowest_tax_lins]
        for rank in tax_ranks:
            tax_dir[f'host_{rank}'].extend(list({lin_dict.get(rank, 'not defined') for lin_dict in lowest_tax_lins_dict}))
        annot_dir[line[0]].update(tax_dir)
    return annot_dir

if __name__ == '__main__':

    metadata_file = snakemake.input[0]
    assert os.path.exists(metadata_file), "Metadata file not found"
    output_file = snakemake.output[0]

    if not os.path.exists('resources/ncbi_taxa.sqlite'):
        print("NCBI taxonomy file ncbi_taxa.sqlite not found, downloading NCBI's taxdump file...")
        urlretrieve('https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz', 'resources/taxdump.tar.gz')
        print('Done!')
        print('Generating ncbi_taxa.sqlite...')
        ncbi = NCBITaxa(dbfile = 'resources/ncbi_taxa.sqlite', taxdump_file = 'resources/taxdump.tar.gz')
        print('Done creating ncbi_taxa.sqlite!')
    else:
        print('Reading local instance of ncbi_taxa.sqlite...')
        ncbi = NCBITaxa(dbfile = 'resources/ncbi_taxa.sqlite')
        print('Done')

    metadata_dict = extract_meta(metadata_file)

    pickle.dump(metadata_dict, open(output_file, 'wb'))
    print('Done extracting metadata from entries in GPD')
