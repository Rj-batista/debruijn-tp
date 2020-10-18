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
	with open(file) as fastq_file : 
		for line in fastq_file: 
			yield(next(fastq_file)) 
			next(fastq_file) 
			next(fastq_file) 

def cut_kmer(sequence,kmer_size): 
	for seq in range(len(sequence)-kmer_size+1): 
		yield(sequence[seq:seq+kmer_size]) 

def build_kmer_dict(fasta, kmer_size): 
	k_mer_dict={} 
	for sequence in read_fastq(fasta): 
		for k_mer in cut_kmer(sequence, kmer_size): 
			if k_mer not in k_mer_dict.keys() : 
				k_mer_dict[k_mer]=0	
			k_mer_dict[k_mer]+=1   	
	return(k_mer_dict)  
    
def build_graph(k_mer_dict): 
    G = nx.DiGraph() 
    for kmer, poids in k_mer_dict.items(): 
        G.add_edge(kmer[:-1],kmer[1:], weight=poids) 
    #nx.draw(G,with_labels=True) 
    #plt.pyplot.show() 
    return(G)
    
def get_starting_node(graph): 
    nodes = list(graph.nodes()) 
    starting_node = []
    for node in nodes:
        if not list(graph.predecessors(node)):
            starting_node.append(node) 
    return(starting_node) 

def get_sink_nodes(graph): 
    nodes = list(graph.nodes()) 
    sink_node = []
    for node in nodes:
        if not list(graph.successors(node)):
            sink_node.append(node) 
    return(sink_node) 
  
def get_contigs(graph,st_node,sn_node): 
    list_contigs=[]
    for starting in st_node:
        for sink in sn_node: 
            path = list(nx.all_simple_paths(graph, starting, sink)) 
    contig=path[0][0]
    for i in range(1,len(path[0])): 
        contig = contig+path[0][i][-1] 
    list_contigs.append((contig, len(contig))) 
    return(list_contigs)

    
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


