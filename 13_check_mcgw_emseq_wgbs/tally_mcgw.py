

#!/usr/bin/env python3

docstring = """
> tally_mcgw.py <

Given a FASTA file, break down sequences into 4-mers, then tally frequencies
of 4-mers that fit NCGN.

At the same time, check whether the 4-mers are MCGW or WCGK. As MCGW is not
palindromic, WCGK (revcomp of MCGW) implies that MCGW is on the non-reference
strand.
""".strip()

import argparse
import collections
from pathlib import Path

import parse_fasta

MCGW_4MERS = ['ACGA', 'ACGT', 'CCGA', 'CCGT']
WCGK_4MERS = ['ACGG', 'ACGT', 'TCGG', 'TCGT']
NCGN_4MERS = ['ACGA', 'ACGC', 'ACGG', 'ACGT',
              'CCGA', 'CCGC', 'CCGG', 'CCGT',
              'GCGA', 'GCGC', 'GCGG', 'GCGT',
              'TCGA', 'TCGC', 'TCGG', 'TCGT']

def count_kmer(seq, k=4):
    """
    Function parses sequence for 4-mers that fit NCGN.
    """
    # just in case, convert everything to uppercase
    seq = seq.upper()
    
    # start tallying
    kmer_count = {}
    for n in range(len(seq) - k + 1):
        subseq = seq[n:n+k]
        if subseq not in kmer_count:
            kmer_count[subseq] = 0
        
        kmer_count[subseq] += 1
    
    # convert dict into a Counter
    kmer_count = collections.Counter(kmer_count)
    
    return kmer_count

parser = argparse.ArgumentParser(
    description=docstring, formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument(
    'fasta_file', metavar='fasta_file', type=Path,
    help='FASTA file containing sequences for MCGW tallying.')
args = parser.parse_args()

 # parse the FASTA file
fasta_seqs = parse_fasta.get_all_sequences(args.fasta_file, 'fasta')

# break all sequences down into 4-mers, and tally them
fourmers = collections.Counter()
for seq in fasta_seqs:
    fourmers += count_kmer(fasta_seqs[seq])

# filter for NCGN
ncgn_dict = {mer:fourmers[mer] for mer in fourmers if mer in NCGN_4MERS}
mcgw_dict = {mer:fourmers[mer] for mer in fourmers if mer in MCGW_4MERS + WCGK_4MERS}

# print results table out
print ('motif', 'frequency', 'MCGW', sep='\t')
for mer in sorted(ncgn_dict):
    print (mer, f'{ncgn_dict[mer]:,}', '*' if mer in MCGW_4MERS + WCGK_4MERS else '', sep='\t')

print ('Total NCGN:', f'{sum(ncgn_dict.values()):,}', sep='\t')
print ('Total MCGW:', f'{sum(mcgw_dict.values()):,}', sep='\t')
print ('% MCGW:', f'{round(sum(mcgw_dict.values()) / sum(ncgn_dict.values()) * 100, 2)}%', sep='\t')
