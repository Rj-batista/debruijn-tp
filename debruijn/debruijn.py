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
    parser.add_argument('eva71_two_reads.fq', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('21', dest='kmer_size', type=int,
                        default=21, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    return parser.parse_args() 

def read_fastq(fastq_file) : 
    """
     Generate a sequence generator.

        Parameters 
        ---------- 
        name : args
            Contain the fastq file .

        Returns
        ------- 
        yield(next(fastq_file))  
            Sequence generator.
    """ 

with open(file) as fastq_file:
	for line in fastq_file:
		yield(next(fastq_file)) 
		next(fastq_file) 
		next(fastq_file) 

def cut_kmer(sequence,kmer_size) : 
	""" 
     Generate a k-mer generator.

        Parameters 
        ---------- 
        name : sequence,kmer-size
            The sequence to cut. 
			the size of the k-mer

        Returns
        ------- 
        k-mer_gen  
            K_mer generator. 
 
	"""   
	for seq in range(len(sequence)-kmer_size+1): 
		yield(sequence[i:i+kmer_size])


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments() 

	first_function=read_fastq('eva71_two_reads.fq')
if __name__ == '__main__':
    main()
