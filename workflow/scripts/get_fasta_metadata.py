#!/usr/bin/env python

import os, re, gzip, glob, pickle, ast, sys
import numpy as np
import pandas as pd
from Bio import SeqIO
from ete3 import NCBITaxa
from collections import Counter, defaultdict
from urllib.request import urlretrieve
from tqdm import tqdm

def get_fasta(seq_file):
    if os.path.splitext(seq_file)[1] == '.gbk':
        with open(seq_file) as input_file:
            seq_record = SeqIO.read(input_file, 'genbank')
    elif os.path.splitext(seq_file)[1] == '.gz':
        with gzip.open(seq_file, 'rt') as input_file:
            seq_record = SeqIO.read(input_file, 'genbank')
    seq_record_acc = seq_record.id
    seq_record_nuc = seq_record.seq
    return (seq_record_acc, seq_record_nuc)

def lineage_dict(taxid, rank_list):
    try:
        taxid_lineage = ncbi.get_lineage(taxid)
        taxid_lineage_ranks = ncbi.get_rank(taxid_lineage)
        taxid_lineage_names = ncbi.get_taxid_translator(taxid_lineage)
        taxid_lineage_dict = {item:[taxid_lineage_names[[x for x,y in taxid_lineage_ranks.items() if y == item][0]]] if item in taxid_lineage_ranks.values() else ['not defined'] for item in rank_list}
    except:
        taxid_lineage_dict = {item:['not defined'] for item in rank_list}
    return taxid_lineage_dict

def get_host_lineage(phage_record, tax_ranks):
    if all (x not in phage_record.features[0].qualifiers.keys() for x in ['host', 'lab_host']):
        host_taxon = phage_record.annotations['organism'].split()[0]
        if all (host_taxon != word for word in ['Enterobacteria', 'unidentified']):
            try:
                host_taxid = ncbi.get_name_translator([host_taxon])[host_taxon][0]
                host_lineage = ncbi.get_lineage(host_taxid)
                host_lineage_names = ncbi.get_taxid_translator(host_lineage)
                if any (x in host_lineage_names.values() for x in ['Bacteria', 'Archaea']):
                    phage_host_pair = (phage_record.id, host_taxid)
                else:
                    phage_host_pair = (phage_record.id, None)
            except:
                phage_host_pair = (phage_record.id, None)
        elif host_taxon == 'Enterobacteria':
            host_taxon = 'Enterobacteriaceae'
            host_taxid = ncbi.get_name_translator([host_taxon])[host_taxon][0]
            phage_host_pair = (phage_record.id, host_taxid)
        else:
            phage_host_pair = (phage_record.id, None)
    else:
        phage_host_taxid_list = [phage_record.id]
        hosts_list = [phage_record.features[0].qualifiers.get(y) for y in ['host', 'lab_host'] if y in phage_record.features[0].qualifiers.keys()]
        host_taxons = {x.split()[0] if 'Candidatus' not in x else x for x in [host for pair in hosts_list for host in pair]}
        for elem in host_taxons:
            try:
                host_taxid = ncbi.get_name_translator([elem])[elem][0]
                host_lineage = ncbi.get_lineage(host_taxid)
                host_lineage_names = ncbi.get_taxid_translator(host_lineage)
                if any (x in host_lineage_names.values() for x in ['Bacteria', 'Archaea']):
                    phage_host_taxid_list.append(host_taxid)
            except:
                continue
        if len(phage_host_taxid_list) > 1:
            phage_host_pair = tuple(phage_host_taxid_list)
        else:
            phage_host_pair = (phage_record.id, None)
    if len(phage_host_pair) > 2:
        counter = 1
        for item in phage_host_pair[1:]:
            if counter == 1:
                host_lineage_dict = lineage_dict(item, tax_ranks)
                counter += 1
            else:
                other_lineage_dict = lineage_dict(item, tax_ranks)
                for rank in tax_ranks:
                    host_lineage_dict[rank].extend(other_lineage_dict[rank])
    elif pd.isnull(phage_host_pair[1]):
        host_lineage_dict = {item:['not defined'] for item in tax_ranks}
    else:
        host_lineage_dict = lineage_dict(phage_host_pair[1], tax_ranks)
    return host_lineage_dict

def read_gbk(gbfile):
    d = {}
    phage_tax_ranks = ['genus', 'subfamily', 'family', 'order']
    host_tax_ranks = ['genus', 'family', 'order', 'class', 'phylum']
    if os.path.splitext(gbfile)[1] == '.gbk':
        with open(gbfile) as input_file:
            record = SeqIO.read(input_file, 'genbank')
    elif os.path.splitext(gbfile)[1] == '.gz':
        with gzip.open(gbfile, 'rt') as input_file:
            record = SeqIO.read(input_file, 'genbank')
    organism = record.annotations['organism']
    taxid = int({item.split(':')[0]:item.split(':')[1] for item in record.features[0].qualifiers['db_xref']}['taxon'])
    phage_taxonomy = lineage_dict(taxid, phage_tax_ranks)
    genome_size = len(record.seq)
    host_taxonomy = get_host_lineage(record, host_tax_ranks)
    d['organism'], d['taxid'], d['genome_size']= (organism, taxid, genome_size)
    d.update({f'phage_{key}':value for key,value in phage_taxonomy.items()})
    d.update({f'host_{key}':value for key,value in host_taxonomy.items()})
    return d

if __name__ == '__main__':

    gbk_dir = snakemake.input[0]
    assert os.path.exists(gbk_dir), "Path does not exist"
    gbk_files = glob.glob(os.path.join(gbk_dir, '*.gbk*'))
    assert len(gbk_files) > 0, 'No files with .gbk extension were found'

    fasta_file = snakemake.output.fa
    metadata_file = snakemake.output.meta
    annot_dir = defaultdict(dict)

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

    with open(fasta_file, 'w') as output_file:
        for item in tqdm(gbk_files, desc = 'Writing combined fasta file and extracting metadata'):
            if not os.path.splitext(item)[1] in ['.gbk', '.gz']:
                print(f'Extensions .gbk and .gz not detected in file {item}')
                continue
            else:
                accession, sequence = get_fasta(item)
                output_file.write(f'>{accession}\n{sequence}\n')
                annot_dir[accession].update(read_gbk(item))
    pickle.dump(annot_dir, open(metadata_file, 'wb'))
    print(f'Done creating {fasta_file} and extracting phage metadata')
