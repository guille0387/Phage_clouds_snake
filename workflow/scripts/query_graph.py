import os, gzip
import pandas as pd
import networkx as nx
from Bio.SeqIO.FastaIO import SimpleFastaParser
from collections import defaultdict
from tqdm import tqdm

def connected_component_subgraphs(G):
    for node_list in nx.connected_components(G):
        yield G.subgraph(node_list)

def add_edges(refgraph, distfile, sizedir):
    if os.path.splitext(distfile)[1] == '.gz':
        with gzip.open(distfile, 'rt') as input_file:
            dist_gen = ((*line.strip().split('\t')[:2], float(line.strip().split('\t')[2])) for line in input_file if float(line.strip().split('\t')[2]) != 1)
    else:
        with open(distfile) as input_file:
            dist_gen = ((*line.strip().split('\t')[:2], float(line.strip().split('\t')[2])) for line in input_file if float(line.strip().split('\t')[2]) != 1)
    included_query = set()
    for target, query, dist in dist_gen:
        if target not in refgraph:
            continue
        else:
            refgraph.add_edge(target, query, weight = dist)
            included_query.add(query)
    nx.set_node_attributes(refgraph, {key:value for key,value in sizedir.items() if key in included_query}, 'genome_size')
    return (refgraph, included_query)

def retrieve_params(distfile, refgraph):
    if os.path.splitext(distfile)[1] == '.gz':
        with gzip.open(distfile, 'rt') as input_file:
            target_query_list = [(*line.strip().split('\t')[:2], float(line.strip().split('\t')[2])) for line in input_file if float(line.strip().split('\t')[2]) != 1]
    else:
        with open(distfile) as input_file:
            dist_gen = [(*line.strip().split('\t')[:2], float(line.strip().split('\t')[2])) for line in input_file if float(line.strip().split('\t')[2]) != 1]
    target_query_list = sorted(target_query_list, key = lambda x: x[2])
    hits_dir = defaultdict(list)
    for target, query, dist in target_query_list:
        if target not in refgraph:
            continue
        else:
            host = refgraph.nodes[target].get('host_genus')
            hits_dir[query].append((dist, target, host))
    best_hits = [hits_dir[x][0] for x in hits_dir.keys()]
    host_set = {y for x in best_hits for y in x[2]}
    host_set.discard('not defined')
    best_dist = [x[0] for x in best_hits]
    threshold = max(best_dist)
    return (threshold, host_set)

def reduce_graph(refgraph, threshold, query_names, upper_threshold):
    retain_edges = [(x,y) for x,y,z in refgraph.edges(data = True) if z['weight'] <= threshold]
    retain_edges_special = [(x,y) for x,y,z in refgraph.edges(data = True) if z['weight'] <= upper_threshold and (x in query_names or y in query_names)]
    retain_edges.extend(retain_edges_special)
    retain_edges = list(set(retain_edges))
    mod_graph = refgraph.edge_subgraph(retain_edges)
    return mod_graph

def get_host_subgraphs(modgraph, hosts, query_names, min_count = 1, keep_target_clouds_only = False):
    graphs = connected_component_subgraphs(modgraph)
    graphs_with_host = []
    for item in tqdm(graphs, desc = 'Searching subgraphs with matched hosts'):
        test = [x for x,y in item.nodes(data = True) if any(host in y.get('host_genus') if 'host_genus' in y.keys() else False for host in hosts)]
        if len(test) < min_count:
            continue
        else:
            graphs_with_host.append(item)
    if keep_target_clouds_only:
        graphs_with_host = [item for item in graphs_with_host if len(query_names.intersection(item.nodes)) > 0]
    G2 = nx.algorithms.operators.all.compose_all(graphs_with_host)
    return G2

def colour_graph(final_graph, hosts, query_names, query_colour = 'orange'):
    colorsd = {x:{'border':'#000000', 'background':'green'} if pd.notnull(y.get('host_genus', None)) and len(set(y.get['host_genus']).intersection(hosts)) > 0 elif x in query_names {'border':'#000000', 'background':query_colour} else {'border':'#000000', 'background':'red'} for x,y in final_graph.nodes(data = True)}
    nx.set_node_attributes(final_graph, colorsd, 'color')
    return final_graph

if __name__ == '__main__':
    fasta_file = snakemake.input.fasta
    dist_file = snakemake.input.dist
    ref_graph_file = snakemake.input.ref_graph
    out_graph_file = snakemake.output[0]
    main_threshold = snakemake.params.thres

    ref_graph = nx.read_gpickle(ref_graph_file)
    if os.path.splitext(fasta_file)[1] == '.gz':
        with gzip.open(fasta_file, 'rt') as input_file:
            query_sizes = {name:len(seq) for name,seq in SimpleFastaParser(input_file)}
    else:
        with open(fasta_file) as input_file:
            query_sizes = {name:len(seq) for name,seq in SimpleFastaParser(input_file)}

    query_graph, added_phages = add_edges(ref_graph, dist_file, query_sizes)
    upper_threshold, matched_host = retrieve_params(dist_file, ref_graph)
    subgraph = reduce_graph(query_graph, main_threshold, added_phages, upper_threshold)
    host_subgraph = get_host_subgraphs(subgraph, matched_host, added_phages, keep_target_clouds_only = True)
    host_subgraph_final = colour_graph(host_subgraph, matched_host, added_phages)
    nx.write_gpickle(host_subgraph_final, out_graph_file)
