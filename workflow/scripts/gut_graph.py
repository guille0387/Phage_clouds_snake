#!/usr/bin/env python

import os, re, glob, sys, gzip, ast, pickle, csv
import networkx as nx
from tqdm import tqdm


#with open('gut_graph_files/metadata.pkl', 'rb') as metadata_file:
#    metadata_dict = pickle.load(metadata_file)

distance_files = glob.glob('gut_graph_files/*.tbl')

with open('gut_graph_files/edge_list.tsv', 'w', newline = '') as edge_file:
    dist_writer = csv.writer(edge_file, delimiter = '\t', quoting = csv.QUOTE_MINIMAL)
    for item in distance_files:
        with open(item, newline = '') as dist_file:
            dist_reader = csv.reader(dist_file, delimiter = '\t')
            for line in tqdm(dist_reader, desc = f'Writing edges from {item}'):
                if line[0] == line[1]:
                    continue
                else:
                    dist_writer.writerow(line[:3])
