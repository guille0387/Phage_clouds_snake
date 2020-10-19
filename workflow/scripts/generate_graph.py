#!/usr/bin/env python

import os, re, glob, sys, gzip, ast, pickle
import numpy as np
import pandas as pd
import networkx as nx
from collections import Counter
from tqdm import tqdm


mash_dist_file = snakemake.input.mashtbl
print('Creating phage intergenomic distance network...')
mash_dist_df = pd.read_csv(mash_dist_file, sep = '\t')
mash_dist_df = mash_dist_df.set_index('#query')

graph_file = snakemake.output.plaingraph
G = nx.from_numpy_matrix(mash_dist_df.values)
labels = mash_dist_df.columns.values
G = nx.relabel_nodes(G, dict(zip(range(len(labels)), labels)))
keep_edges = [(x,y) for x,y,z in G.edges(data = True) if pd.notnull(z['weight'])]
Gmod = nx.Graph(G.edge_subgraph(keep_edges))
singletons = []
for item in mash_dist_df.itertuples(index = True):
	non_null_values = [x for x in item[1:] if pd.notnull(x)]
	if len(non_null_values) == 1:
		singletons.append(item[0])
Gmod.add_nodes_from(singletons)
nx.write_gpickle(Gmod, graph_file)
print(f'Done! network saved as {graph_file}')

annot_dict_file = snakemake.input.annotdict
with open(annot_dict_file, 'rb') as metadata_file:
	annot_dir = pickle.load(metadata_file)

annot_graph_file = snakemake.output.annotgraph
for node in tqdm(Gmod.nodes, desc = 'Annotating phage intergenomic network'):
	Gmod.nodes[node].update(annot_dir[node])

nx.write_gpickle(Gmod, annot_graph_file)
print(f'Done! annotated network saved as {annot_graph_file}')
