#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html
#    python /debruijn/debruijn.py -i /data/eva71_hundred_reads.fq

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
import networkx as nx
import matplotlib
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics
import matplotlib.pyplot as plt
matplotlib.use("Agg")

__author__ = "Jane Schadtler-Law"
__copyright__ = "Universite Paris Cité"
__credits__ = ["Jane Schadtler-Law"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Jane Schadtler-Law"
__email__ = "jane.schadtler-law@etu.u-paris.fr"
__status__ = "Developpement"

def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=22, help="k-mer size (default 22)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    parser.add_argument('-f', dest='graphimg_file', type=str,
                        help="Save graph as image (png)")
    return parser.parse_args()


def read_fastq(fastq_file):
    with open(fastq_file, 'r') as f:
        for i in f:
            yield next(f).strip()
            next(f)
            next(f)

def cut_kmer(read, kmer_size):
    for i in range(len(read)):
        if i + kmer_size <= len(read): 
            yield read[i : i + kmer_size] 

def build_kmer_dict(fastq_file, kmer_size):
    kmer_occ = {}
    for seq in read_fastq(fastq_file):
        for kmer in cut_kmer(seq, kmer_size):
            if kmer in kmer_occ: 
                kmer_occ[kmer] += 1
            else:
                kmer_occ[kmer] = kmer_occ.get(kmer, 0) + 1
    return kmer_occ


def build_graph(kmer_dict):
    Graph = nx.DiGraph()
    for key in kmer_dict:
        Graph.add_edge(key[:-1], key[1:], weight = kmer_dict[key])
    return Graph


def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    for path in path_list:
        if delete_entry_node and not delete_sink_node:
            graph.remove_nodes_from(path[:-1])
        if delete_sink_node and not delete_entry_node:
            graph.remove_nodes_from(path[1:])
        if delete_entry_node and delete_sink_node: 
            graph.remove_nodes_from(path)
        if not delete_entry_node and not delete_sink_node:
            graph.remove_nodes_from(path[1:-1])  
    return graph

def std(data):
    return statistics.stdev([i for i in data])


def select_best_path(graph, path_list, path_length, weight_avg_list, 
                     delete_entry_node=False, delete_sink_node=False):
    if len(path_list) > 1 : 
        if std([i for i in weight_avg_list]) != 0 :
            path_to_keep = weight_avg_list.index(max(weight_avg_list))
    
        elif std([i for i in path_length]) != 0 :
            path_to_keep = path_length.index(max(path_length))
        
        else : 
            path_to_keep = randint(0, len(path_list))
        
        path_l = list(path_list)
        path_l.pop(path_to_keep)

        graph = remove_paths(graph, path_l, delete_entry_node, delete_sink_node)

    return graph

def path_average_weight(graph, path):
    return statistics.mean([d["weight"] for (u, v, d) in graph.subgraph(path).edges(data=True)])

def solve_bubble(graph, ancestor_node, descendant_node):
    paths = list(nx.all_simple_paths(graph, ancestor_node, descendant_node))
    weight_avg_list = [path_average_weight(graph, i) for i in paths]
    path_length = [len(i) for i in paths]
    return select_best_path(graph, paths, path_length, weight_avg_list)
    # map(lambda x : fun(x), paths)

from itertools import combinations

def simplify_bubbles(graph):
    bubble = False 
    for node in graph.nodes():
        list_predecessors = list(graph.predecessors(node))
        if len(list_predecessors) > 1:
            for node_anc_i,node_anc_j in combinations(list_predecessors, 2):
               node_ancestor = nx.lowest_common_ancestor(graph, node_anc_i, node_anc_j)
               if node_ancestor:
                bubble = True
                break
            if bubble:
                graph = simplify_bubbles(solve_bubble(graph, node_ancestor, node))
                break
    return graph


def solve_entry_tips(graph, starting_nodes):
    path_l = []
    path_len = []
    path_weight = []
    for node in starting_nodes:
        for des_node in nx.descendants(graph, node):
            pred_node = list(graph.predecessors(des_node))
            if len(pred_node) >= 2:
                for path in nx.all_simple_paths(graph, node, des_node):
                    path_l.append(path)
                    path_len.append(len(path))
                    path_weight.append(path_average_weight(graph, path))
    graphe = select_best_path(graph, path_l, path_len, path_weight, True, False)
    return graphe

def solve_out_tips(graph, ending_nodes):
    path_l = []
    path_len = []
    path_weight = []
    for node in ending_nodes:
        for pre_node in nx.ancestors(graph, node):
            des_node = list(graph.successors(pre_node))
            if len(des_node) >= 2:
                for path in nx.all_simple_paths(graph, pre_node, node):
                    path_l.append(path)
                    path_len.append(len(path))
                    path_weight.append(path_average_weight(graph, path))
    graphe = select_best_path(graph, path_l, path_len, path_weight, False, True)
    return graphe

def get_starting_nodes(graph):
    node_start = []
    for node in graph.nodes():
        if  not list(graph.predecessors(node)):
            node_start.append(node)
    return(node_start)

def get_sink_nodes(graph):
    node_end = []
    for node in graph.nodes():
        if not list(graph.successors(node)):
            node_end.append(node)
    return(node_end)

def get_contigs(graph, starting_nodes, ending_nodes):
    all_paths = []
    for start in starting_nodes:
        for end in ending_nodes:
            if list(nx.all_simple_paths(graph, start, end)): 
                for path in list(nx.all_simple_paths(graph, start, end)):
                    contig = path[0]
                    for node in range(1, len(path)):
                       contig = contig + path[node][1]
                    all_paths.append((contig, len(contig)))
    return(all_paths)


def save_contigs(contigs_list, output_file):
    with open(output_file, "w") as filout:
        for i, contig in enumerate(contigs_list):
            filout.write(f">contig_{i} len={contig[1]}\n")
            filout.write(f"{fill(contig[0])}\n")


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def draw_graph(graph, graphimg_file):
    """Draw the graph
    """                                    
    fig, ax = plt.subplots()
    elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]
    #print(elarge)
    esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]
    #print(elarge)
    # Draw the graph with networkx
    #pos=nx.spring_layout(graph)
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5, 
                           edge_color='b', style='dashed')
    #nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
    # save image
    plt.savefig(graphimg_file)


def save_graph(graph, graph_file):
    """Save the graph with pickle
    """
    with open(graph_file, "wt") as save:
            pickle.dump(graph, save)


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments

    args = get_arguments()
    # Lecture du fichier et construction du graphe
    dico = build_kmer_dict(fastq_file = args.fastq_file, kmer_size = args.kmer_size)
    graph = build_graph(dico)
    
    # Résolution des bulles
    graph = simplify_bubbles(graph)
    
    # Résolution des pointes d’entrée et de sortie
    starting_nodes = get_starting_nodes(graph)
    ending_nodes = get_sink_nodes(graph)
    graph = solve_entry_tips(graph, starting_nodes)
    graph = solve_out_tips(graph, ending_nodes)
    
    # Ecriture du/des contigs 
    contigs = get_contigs(graph, starting_nodes, ending_nodes)
    save_contigs(contigs, args.output_file)

    # Fonctions de dessin du graphe
    # A decommenter si vous souhaitez visualiser un petit 
    # graphe
    # Plot the graph
    if args.graphimg_file:
        draw_graph(graph, args.graphimg_file)
        save_graph(graph, args.graphimg_file)
    # Save the graph in file
    #if args.graph_file:
    #    


if __name__ == '__main__':
    main()
