#!/usr/bin/env python

import os, glob, gzip
from Bio import SeqIO
from tqdm import tqdm

def get_fasta(seq_file_handle):
    seq_record = SeqIO.read(seq_file_handle, 'genbank')
    seq_record_acc = seq_record.id
    seq_record_nuc = seq_record.seq
    return (seq_record_acc, seq_record_nuc)

if __name__ == '__main__':

    gbk_dir = snakemake.input[0]

    assert os.path.exists(gbk_dir), "Path does not exist"

    gbk_files = glob.glob(os.path.join(gbk_dir, '*.gbk*'))

    assert len(gbk_files) > 0, 'No files with .gbk extension were found'

    fasta_file = snakemake.output[0]

    with open(fasta_file, 'w') as output_file:
        for item in tqdm(gbk_files, desc = 'Writing combined fasta file with all phage genome sequences'):
            if os.path.splitext(item)[1] == '.gbk':
                with open(item) as input_file:
                    accession, sequence = get_fasta(input_file)
                    output_file.write(f'>{accession}\n{sequence}\n')
            elif os.path.splitext(item)[1] == '.gz':
                with gzip.open(item, 'rt') as input_file:
                    accession, sequence = get_fasta(input_file)
                    output_file.write(f'>{accession}\n{sequence}\n')
    print(f'Done creating {fasta_file}')
