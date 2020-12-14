#!/usr/bin/env python

import os, gzip
import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser


frag_file = gzip.open(f'{file_prefix}_genome_fragment.fa.gz', 'wt')

def cat_seq(name_seq_pair, metadata_df = meta_df):
    name, seq = name_seq_pair
    name = name.split()[0]
    seq_index = metadata_df.index[metadata_df['GPD_id'] == name]
    category = list(metadata_df['checkV_MIUViG'].loc[seq_index])[0]
    return (name, seq, category)

if __name__ == '__main__':

    seq_file = snakemake.input.fasta
    meta_file = snakemake.input.meta
    high_file = snakemake.output.high
    frag_file = snakemake.output.frag

    file_prefix = seq_file.split('.')[0]

    if os.path.splitext(meta_file)[-1] == '.gz':
        with gzip.open(meta_file, 'rt') as input_file:
            meta_df = pd.read_csv(input_file, sep = '\t')
    else:
        with open(meta_file) as input_file:
            meta_df = pd.read_csv(input_file, sep = '\t')

    if os.path.splitext(seq_file)[-1] == '.gz':
        seq_map = map(cat_seq, SimpleFastaParser(gzip.open(seq_file, 'rt')))
    else:
        seq_map = map(cat_seq, SimpleFastaParser(open(seq_file)))

    high_handle = gzip.open(high_file, 'wt')
    frag_handle = gzip.open(frag_file, 'wt')

    for name, seq, cat in seq_map:
        if cat == 'High-quality':
            high_handle.write(f'>{name}\n{seq}\n')
        else:
            frag_handle.write(f'>{name}\n{seq}\n')

    high_handle.close()
    frag_handle.close()
