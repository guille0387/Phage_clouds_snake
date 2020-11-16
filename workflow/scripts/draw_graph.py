#!/usr/bin/env python
# coding: utf-8

# ## Import required libraries and functions

import os, re, glob, sys, gzip, ast
import numpy as np
import pandas as pd
import networkx as nx
from networkx.drawing.nx_agraph import graphviz_layout
import matplotlib.pyplot as plt
from pyvis.network import Network
from collections import Counter
from tqdm import tqdm

# ### Define function that returns a generator of subgraph from a user-defined graph

def connected_component_subgraphs(G):
    for node_list in nx.connected_components(G):
        yield G.subgraph(node_list)

# ### Define function for creating subgraph based on host genus and distance threshold

def mod_nx_graph(graph, thres, host):
    '''This function takes a phage distance graph previously annotated with their cognate host genera and generates a custom
    graph that only includes clouds whose edges weights are <= than thres, and that contain at least one phage targeting host'''
    retain_edges = [(x,y) for x,y,z in graph.edges(data = True) if z['weight'] <= thres]
    mod_graph = graph.edge_subgraph(retain_edges)
    keep_subgraphs = [item for item in connected_component_subgraphs(mod_graph) if any (host in y['host_genus'] for x,y in item.nodes(data = True))]
    joint_graph = nx.algorithms.operators.all.compose_all(keep_subgraphs)
    colorsd = {x:{'border':'#000000', 'background':'green'} if host in y['host_genus'] else {'border':'#000000', 'background':'red'} for x,y in joint_graph.nodes(data = True)}
    titlesd = {x: f'Target host genera:<br>{";".join(y["host_genus"])}<br>Genome size:<br>{y["genome_size"]} bp<br>Phage genus:<br>{y["phage_genus"]}<br>Accession:<br>{x}' for x,y in joint_graph.nodes(data = True)}
    nx.set_node_attributes(joint_graph, colorsd, 'color')
    nx.set_node_attributes(joint_graph, titlesd, 'title')
    nx.relabel_nodes(joint_graph, {x:y['organism'] for x,y in joint_graph.nodes(data = True)}, copy = False)
    #max_genome = max([y['genome_size'] for x,y in joint_graph.nodes(data = True)])
    #for node in joint_graph.nodes():
    #    joint_graph.nodes[node]['size'] = joint_graph.nodes[node]['genome_size'] / max_genome
    return joint_graph


# ### Define function to create pyvis visualization of graphs

def draw_graph3(networkx_graph,threshold,notebook=False,show_buttons=True,only_physics_buttons=False):
    from pyvis import network as net
    pyvis_graph = net.Network(notebook=notebook, height='1500px', width='1500px')
    pyvis_graph.force_atlas_2based()
    #pyvis_graph.inherit_edge_colors(False)

    for node,node_attrs in networkx_graph.nodes(data=True):
        pyvis_graph.add_node(node,**node_attrs)
    for source,target,edge_attrs in tqdm(networkx_graph.edges(data=True), desc="Drawing the graph"):
        if not 'value' in edge_attrs and not 'width' in edge_attrs and 'weight' in edge_attrs:
            edge_attrs['value']=threshold-edge_attrs['weight']+0.1
            #edge_attrs['color'] = '#000000'
            edge_attrs['color'] = 'lightgray'
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

    phage_host = snakemake.params.host
    dist_thres = snakemake.params.thres
    graph_file = snakemake.input[0]
    outfile = snakemake.output[0]

    assert isinstance(dist_thres, float), 'Entered threshold value is not float'

    print('Reading main full graph...')
    G = nx.read_gpickle(graph_file)
    print('Done!')

    assert dist_thres <= max([z['weight'] for x,y,z in G.edges(data = True)]), 'Entered threshold value is too high'

    print(f'Creating modified graph with {phage_host} as target host and {dist_thres} as mash distance threshold...')
    G_mod = mod_nx_graph(G, dist_thres, phage_host)
    print('Done!')
    print('Creating graph interactive file...')
    G_mod_with_sizes = set_node_sizes_for_plotting(G_mod)
    pg = draw_graph3(G_mod_with_sizes, dist_thres)
    pg.save_graph(outfile)
    print('Done!')
