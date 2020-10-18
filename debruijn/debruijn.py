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

#python debruijn.py -i ../data/eva71_two_reads.fq

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
import networkx as nx
import matplotlib as plt
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics

__author__ = ""
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Batista Reda"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
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
                        default=21, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    return parser.parse_args() 

def read_fastq(file):  
    """
    Read a fastq file .

        Parameters 
        ---------- 
        name : file
            contain de fastq file .

        Returns
        ------- 
        cursor  
            The sequence in the fastq file.
    """ 
    with open(file) as fastq_file : 
        for line in fastq_file: 
            yield(next(fastq_file)) 
            next(fastq_file) 
            next(fastq_file) 

def cut_kmer(sequence,kmer_size): 
    """
    Generator of kmer .

        Parameters 
        ---------- 
        name : sequence, kmer_size
            Contain the sequence. 
            the length of the kmer we want to produce

        Returns
        ------- 
          (sequence[seq:seq+kmer_size])
            A generator of kmer.
    """  
    for seq in range(len(sequence)-kmer_size+1): 
        yield(sequence[seq:seq+kmer_size]) 

def build_kmer_dict(fasta, kmer_size):  
    """ 
    Create a dictionnary with kmer as key and occurences as value .

        Parameters 
        ---------- 
        name : fasta, kmer_size
            the fastq file. 
            the length of the kmer we want to produce

        Returns
        ------- 
          k_mer_dict
            A dictionnary of kmer.
    """ 
    
    k_mer_dict={}
    for sequence in read_fastq(fasta): 
        for k_mer in cut_kmer(sequence, kmer_size): 
            if k_mer not in k_mer_dict.keys() : 
                k_mer_dict[k_mer]=0	
            k_mer_dict[k_mer]+=1   	
    return(k_mer_dict)  
    
def build_graph(k_mer_dict): 
    """
    Create a graph of de Bruijn .

        Parameters 
        ---------- 
        name : k_mer_dict
            A dictionnary of kmer

        Returns
        ------- 
          G
            A de Bruijn graph.
    """
    G = nx.DiGraph() 
    for kmer, poids in k_mer_dict.items(): 
        G.add_edge(kmer[:-1],kmer[1:], weight=poids) 
    #nx.draw(G,with_labels=True) 
    #plt.pyplot.show() 
    return(G)
    
def get_starting_node(graph): 
    """
    Get the starting node from a de Bruijn graph .

        Parameters 
        ---------- 
        name : graph 
        A de Bruijn graph

        Returns
        ------- 
          starting_node
            The starting node of a de Bruijn graph.
    """ 
    nodes = list(graph.nodes()) 
    starting_node = []
    for node in nodes:
        if not list(graph.predecessors(node)):
            starting_node.append(node) 
    return(starting_node) 

def get_sink_nodes(graph): 
    """
    Get the sink node from a de Bruijn graph .

        Parameters 
        ---------- 
        name : graph 
        A de Bruijn graph

        Returns
        ------- 
          sink_node
            The sink node of a de Bruijn graph.
    """ 
    nodes = list(graph.nodes()) 
    sink_node = []
    for node in nodes:
        if not list(graph.successors(node)):
            sink_node.append(node) 
    return(sink_node) 
  
def get_contigs(graph,st_node,sn_node): 
    """
    Create a liste of tuple of contigs .

        Parameters 
        ---------- 
        name : graph,st_node,sn_node 
        A de Bruijn graph 
        Starting node 
        Sink node

        Returns
        ------- 
          list_contigs. 
          A list of contigs
    """  
    list_contigs=[]
    for starting in st_node:
        for sink in sn_node: 
            path = list(nx.all_simple_paths(graph, starting, sink)) 
    contig=path[0][0]
    for i in range(1,len(path[0])): 
        contig = contig+path[0][i][-1] 
    list_contigs.append((contig, len(contig))) 
    return(list_contigs) 


def std(val):   
    """
    Give the standard deviation of a list .

        Parameters 
        ---------- 
        name : val 
        A list of values

        Returns
        ------- 
        statistics.stdev(val). 
        The standard deviation of the list
    """ 
    return (statistics.stdev(val))
    
def path_average_weight(graph,path): 
    """
    Give the average weight path .

        Parameters 
        ---------- 
        name : graph,path 
        A de Bruijn graph 
        A list of path

        Returns
        ------- 
        average_weight 
        The average weight path
    """
    poids = 0
    for w in path: 
        poids += G.out_degree(path, weight = "weight")
    average_weight = poids / (len(path) - 1)
    return(average_weight) 
    
def remove_path(graph,list_path,delete_entry_node,delete_sink_node): 
    """
    Give a graph with cleaned path .

        Parameters 
        ---------- 
        name : graph,path,delete_entry_node,delete_sink_node 
        A de Bruijn graph 
        A list of path 
        To indicate if a entry node is deleted 
        To indicate if a sink node is deleted
        

        Returns
        ------- 
        average_weight 
        A de Bruijn graph with cleaned path
    """
    for chemin in list_path: 
        if delete_entry_node: 
            graph.remove_node(chemin[0]) 
        if delete_sink_node: 
            graph.remove_node(chemin[-1]) 
    return(graph) 
    
def select_best_path(graph,list_path,list_lenght_path,list_average_weight,
    delete_entry_node=False,delete_sink_node=False): 
    """
    Give a graph with the best path .

        Parameters 
        ---------- 
        name : graph,list_path,list_lenght_path,list_average_weight,delete_entry_node,delete_sink_node 
        A de Bruijn graph 
        A list of path 
        A list of lenght path 
        A list of average weight
        To indicate if a entry node is deleted 
        To indicate if a sink node is deleted

        Returns
        ------- 
        average_weight 
        The average weight path
    """
    min_weight = [i for i, weight in enumerate(list_average_weight) if weight == min(list_average_weight)]
    if len(min_weight) == 1:
        return remove_paths(graph, [list_path[min_weight[0]]],
                            delete_entry_node, delete_sink_node)
    min_length = [i for i, length in enumerate(list_lenght_path) if length == min(list_lenght_path)]
    if len(min_length) == 1:
        return remove_paths(graph, [list_path[min_length[0]]],
                            delete_entry_node, delete_sink_node) 
    rand = randint(0, 1)
    return remove_paths(graph, [paths[rand]],
                        delete_entry_node, delete_sink_node) 
                        

def solve_bubble(graph,n_predecessor,n_successor): 
    """
    Delete the bubble in the graph .

        Parameters 
        ---------- 
        name : graph,n_predecessor,n_successor
        A predecessor node 
        A successor node

        Returns
        ------- 
        average_weight 
        Cleared graph from bubble
    """
    
    list_chemin = []
	list_lenght = []
	list_weight = []
	for chemin in nx.all_simple_paths(graph, ancestor, descendant): 
		list_chemin.append(chemin)
		list_lenght.append(len(chemin))
		list_weight.append(path_average_weight(graph, chemin))
	graph = select_best_path(graph, list_chemin, list_lenght, list_weight, delete_entry_node = False, delete_sink_node = False)
	return(graph)


    
    
    
#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
args = get_arguments()
kmer_dict=build_kmer_dict(args.fastq_file, args.kmer_size) 
dbj_graph=build_graph(kmer_dict) 
st_node=get_starting_node(dbj_graph) 
sn_node=get_sink_nodes(dbj_graph) 
contigs=get_contigs(dbj_graph,st_node,sn_node) 
if __name__ == '__main__':
	main() 	
	"""	
	for i in read_fastq("../data/eva71_two_reads.fq"): 
		print(i) 
		for j in cut_kmer(i,4): 
			print(j) 
	"""


