#!/usr/bin/env python3

doc = """
> tabulate_nnncgnnn_universe.py <

Script tabulates the number of distinct NNN CG NNN in the provided genome,
calculates a universal value for enrichment/depletion calculations.

Output automatically drops uninteresting scaffolds: lambda and pUC19
(control genomes).

Output format:
  1. NNNCGNNN sequence
  2. frequency

Autodetects whether input files are gzip-compressed from `.gz` suffix.
""".strip()

import argparse
import collections
from pathlib import Path
import re
import sys
import time

import parse_fasta

def benchmark_print(message):
    """
    Prints a timestamped message to stderr, to aid visual assessment of
    chokepoints in code.
    """
    print (f'[{time.asctime()}] {message}', file=sys.stderr)

parser = argparse.ArgumentParser(description=doc, 
                                 formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('genome_fasta', type=Path,
                    help='FASTA file of genome.')
parser.add_argument('-v', '--verbose', action='store_true',
                    help='prints diagnostic stuff to stderr.')

args = parser.parse_args()

# read sequences into memory
genome_seq = parse_fasta.get_all_sequences(args.genome_fasta, 'fasta')

# loop through sequences
nnncgnnn_tally = collections.Counter()
for seq in genome_seq:
    if args.verbose: benchmark_print(f'Parsing {seq}...')
    
    # ignore scaffolds from control genomes
    if seq in ['lambda', 'pUC19']: continue
    
    # use lookahead assertion so that string isn't consumed when matching
    # >>> x = 'AAACGCGAACGCGAAA'
    # >>> re.findall(r'\w{3}CG\w{3}', x)
    # ['AAACGCGA', 'ACGCGAAA']
    # >>> re.findall(r'(?=(\w{3}CG\w{3}))', x)
    # ['AAACGCGA', 'ACGCGAAC', 'GAACGCGA', 'ACGCGAAA']
    seq_nnncgnnn = re.findall(r'(?=(\w{3}CG\w{3}))', genome_seq[seq])
    
    # merge tallies together
    nnncgnnn_tally += collections.Counter(seq_nnncgnnn)

# print tally out
for nt in sorted(nnncgnnn_tally):
    print (nt, nnncgnnn_tally[nt], sep='\t')
