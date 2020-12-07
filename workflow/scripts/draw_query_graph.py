#!/usr/bin/env python

import os, re, glob, sys, gzip, ast
import numpy as np
import pandas as pd
import networkx as nx
from networkx.drawing.nx_agraph import graphviz_layout
import matplotlib.pyplot as plt
from pyvis import network as net
from collections import Counter
from tqdm import tqdm

def draw_graph3(networkx_graph,threshold,notebook=False,show_buttons=True,only_physics_buttons=False):
    pyvis_graph = net.Network(notebook=notebook, height='1500px', width='1500px')
    pyvis_graph.force_atlas_2based()
    max_dist = max([z['weight'] for x,y,z in networkx_graph.edges(data = True)])
    for node,node_attrs in networkx_graph.nodes(data=True):
        pyvis_graph.add_node(node,**node_attrs)
    for source,target,edge_attrs in tqdm(networkx_graph.edges(data=True), desc="Drawing the graph"):
        if not 'value' in edge_attrs and not 'width' in edge_attrs and 'weight' in edge_attrs:
            edge_attrs['value']=max_dist-edge_attrs['weight']+0.1
        edge_attrs['color'] = 'lightgray'
        if edge_attrs['weight'] > threshold:
            edge_attrs['dashes'] = True
        pyvis_graph.add_edge(source,target,**edge_attrs)
    if show_buttons:
        if only_physics_buttons:
            pyvis_graph.show_buttons(filter_=['physics'])
        else:
            pyvis_graph.show_buttons()
    return pyvis_graph

# ### Resize the nodes to refelct the genome sizes
def set_node_sizes_for_plotting(graph, size_factor=3000):
    for node in graph.nodes():
        size = int(graph.nodes[node].get('genome_size', 40000) / size_factor)
        graph.nodes[node]['size'] = size
    return graph

if __name__ == '__main__':

    query_graph = snakemake.input[0]
    query_pyvis = snakemake.output[0]
    threshold = snakemake.params.thres
    control_board = snakemake.params.control

    print('Reading graph that includes query phages...')
    G = nx.read_gpickle(query_graph)
    print('Done reading graph!')

    print('Creating interactive graph file...')
    G_with_sizes = set_node_sizes_for_plotting(G)
    if control_board == 'off':
        pg = draw_graph3(G_with_sizes, threshold, show_buttons = False)
    else:
        pg = draw_graph3(G_with_sizes, threshold)
    pg.save_graph(query_pyvis)
    print('Done making interactive graph!')
