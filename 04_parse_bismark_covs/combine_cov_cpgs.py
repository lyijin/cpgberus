#!/usr/bin/env python3

WINDOW_SIZE = 102
BP_UP_DOWN = int((WINDOW_SIZE - 2) / 2)
doc = f"""
> combine_cov_cpgs.py <

This script combines coverages from the C on Watson and C on Crick to produce a
per-CpG view of methylation.

Besides combining coverages and recalculating methylation levels, this script
also calculates coverage evenness, deviations in absolute methylation levels
from the original strands, and the GC% of the subsequence surrounding the
CpG dinucleotide.

Useful to check whether hemimethylation is present in the original cov file
(coverage approx even, but large delta absolute methylation levels).

Requires two (or more) input files:
  1. the genome fasta file (can be compressed) used in `bismark` mapping, to
     parse coords of CpGs in the genome
  2. one (or more) `bismark` .cov files from `bismark_methylation_extractor` 
     with per-position (not per-CpG) information

Produces a bismark-cov-esque file, with two extra columns and a header row.
Columns are
  1. scaffold
  2. dinuc start pos (1-based)
  3. dinuc end pos (1-based)
  4. combined meth %
  5. combined meth coverage
  6. combined unmeth coverage
  7. "coverage evenness"
     (0.5 means that the combined coverage was equally partitioned over both
      strands, 1 means coverage was completely dominated by a single strand)
  8. deviations in absolute methylation levels
     (absolute delta methylation % difference: abs(C-on-Watson - C-on-Crick))
  9. GC% of surrounding sequence
     (CpG dinucleotide centred in a window of {WINDOW_SIZE} bp)

There are no coverage filters with this script, do post-process the combined
tables in your $FAVE_SCRIPTING_LANGUAGE.

Table compilation automatically drops uninteresting scaffolds: lambda and pUC19
(control genomes).
""".strip()

import argparse
import csv
import gzip
from pathlib import Path
import sys
import time

import numpy as np

import parse_fasta

def benchmark_print(message):
    """
    Prints a timestamped message to stderr, to aid visual assessment of
    chokepoints in code.
    """
    print (f'[{time.asctime()}] {message}', file=sys.stderr)

def calc_gc_pct(scaf_nparray, dinuc_start):
    """
    Function calculates the GC% of the immediate context surrounding the CpG
    dinucleotide. Window size is controlled by WINDOW_SIZE constant defined
    at the beginning of script.
    """
    # "dinuc_start" too close to scaffold start
    if dinuc_start <= BP_UP_DOWN: return 'NA'
    
    # "dinuc_start" too close to scaffold end
    if dinuc_start + BP_UP_DOWN >= len(scaf_nparray): return 'NA'
    
    # all good!
    #
    # calculation check: window_start:window_end should equals WINDOW_SIZE
    #   (dinuc_start + 1 + BP_UP_DOWN) - (dinuc_start - 1 - BP_UP_DOWN)
    #   = 2 * BP_UP_DOWN + 2
    #   = WINDOW_SIZE
    window_start = dinuc_start - 1 - BP_UP_DOWN
    window_end = dinuc_start + 1 + BP_UP_DOWN
    acgt_bases = np.count_nonzero(scaf_nparray[window_start:window_end])
    cg_bases = np.count_nonzero(scaf_nparray[window_start:window_end] > 0)
    gc_pct = f'{cg_bases / acgt_bases * 100:.3f}'
    
    return gc_pct

def get_cpg_context(scaf_seq, dinuc_start):
    """
    Function returns the sequence context surrounding the CpG, to check
    whether there are methylation level biases wrt sequence context.
    
    Context is hardcoded as NNNCGNNN.
    """
    NUM_N = 3
    
    # "dinuc_start" too close to scaffold start
    if dinuc_start <= NUM_N: return 'NA'
    
    # "dinuc_start" too close to scaffold end
    if dinuc_start + NUM_N >= len(scaf_seq): return 'NA'
    
    # all good!
    #
    # calculation check: window_start:window_end should equals 8
    #   (dinuc_start + 1 + NUM_N) - (dinuc_start - 1 - NUM_N)
    #   = 2 * NUM_N + 2
    #   = 8
    window_start = dinuc_start - 1 - NUM_N
    window_end = dinuc_start + 1 + NUM_N
    
    return scaf_seq[window_start:window_end].upper()

parser = argparse.ArgumentParser(description=doc, 
                                 formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('genome_fasta', type=Path,
                    help='genome FASTA used in mapping (can be compressed).')
parser.add_argument('bismark_covs', type=Path, nargs='+',
                    help='one or more bismark cov files (can be compressed).')
parser.add_argument('-v', '--verbose', action='store_true',
                    help='prints diagnostic stuff to stderr.')

args = parser.parse_args()

# reads genome FASTA
if args.verbose: benchmark_print(f'Reading {args.genome_fasta.name}...')
genome_seq = parse_fasta.get_all_sequences(args.genome_fasta, 'fasta')

# to facilitate near O(n) lookups of whether a certain genomic pos is a CpG
# or not, and whether it's the C-on-Watson or C-on-Crick, and for GC%
# calculations, convert "genome_seq" to a NumPy dict with int8 values for
# memory-efficient lookups. values for "genome_int" are
#  -1: A/T
#   0: non-ACGT base
#   1: C/G not in CpG
#   2: C-on-Watson in CpG
#   3: C-on-Crick in CpG
trans_from = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
trans_to   = '/010001000000000000/000000'
trans_dict = str.maketrans(trans_from, trans_to)
genome_int = {}
for annot in genome_seq:
    genome_int[annot] = \
        genome_seq[annot].upper().replace('CG', '23').translate(trans_dict)
    
    # note that np.frombuffer converts characters to ASCII values; "- 48"
    # turns ASCII values back into digits. chr(47) == "/"; chr(48) == "0".
    genome_int[annot] = np.frombuffer(genome_int[annot].encode(), dtype=np.int8) - 48

# start parsing individual bismark covs
for bc in args.bismark_covs:
    if args.verbose: benchmark_print(f'Processing {bc.name}...')
    
    # have a buffer that is two-lines deep (i.e. "prev_line" and "curr_line")
    # if both lines belong to the same CpG dinuc, combine stats, flush it out
    # if both lines do not belong to the same dinuc, flush "prev_line" out
    # exception: when last line is encountered, flush both lines
    prev_line = []
    curr_line = []
    
    bc_file = gzip.open(bc, 'rt') if bc.suffix == '.gz' else open(bc)
    tsv_reader = csv.reader(bc_file, delimiter='\t')
    
    # print header
    if bc.suffix == '.gz':
        # strip last two suffixes (".cov.gz")
        outf = '.'.join(bc.name.split('.')[:-2]) + '.combined.tsv.gz'
    else:
        # strip last suffix (".cov")
        outf = '.'.join(bc.name.split('.')[:-1]) + '.combined.tsv.gz'
    outf = gzip.open(outf, 'wt')
    
    if args.verbose: benchmark_print(f'Outputting into {outf.name}...')
    
    print ('scaffold', 'dinuc_start', 'dinuc_end', 'meth_pct', 'meth_cov', 
           'unmeth_cov', 'evenness', 'abs_delta_meth_pct',
           'gc_pct', 'cpg_context_nnncgnnn', sep='\t', file=outf)
    
    for line in tsv_reader:
        # skip lines from scaffolds deemed uninteresting
        if line[0] in ['lambda', 'pUC19']: continue
        
        # at this point, "prev_line" is defo empty. conditional below checks
        # whether both lines are empty, or just "prev_line"
        if curr_line:   # only "prev_line" is empty
            prev_line = curr_line
            curr_line = line
        else:
            prev_line = line
            try:
                curr_line = next(tsv_reader)
            except StopIteration:
                # end-of-file reached. give dummy values to "curr_line" so
                # this loop can run one last time to process info in "prev_line"
                curr_line = ['dummy', 0, 0, 0, 0, 0]
        
        # check whether both lines belong to the same CpG dinuc
        prev_scaf, prev_pos, _, _, prev_meth, prev_unmeth = prev_line
        prev_pos, prev_meth, prev_unmeth = int(prev_pos), int(prev_meth), int(prev_unmeth)
        
        curr_scaf, curr_pos, _, _, curr_meth, curr_unmeth = curr_line
        curr_pos, curr_meth, curr_unmeth = int(curr_pos), int(curr_meth), int(curr_unmeth)
        
        # check whether prev_scaf:prev_pos is NOT in a CpG dinuc
        if genome_int[prev_scaf][prev_pos - 1] < 2:    # convert from 1-based to 0-based
            # flush "prev_line", compute nothing
            prev_line = []
            continue
        
        # past this stage, prev_scaf:prev_pos IS in a CpG
        #        "prev_pos" and "curr_pos" are adjacent                   "prev_pos" is C-on-Watson
        #  ---------------------------------------------------     ----------------------------------------
        if prev_scaf == curr_scaf and prev_pos + 1 == curr_pos and genome_int[prev_scaf][prev_pos - 1] == 2:
            # both lines in buffer are from same dinuc!
            combined_meth = prev_meth + curr_meth
            combined_unmeth = prev_unmeth + curr_unmeth
            combined_beta = combined_meth / (combined_meth + combined_unmeth) * 100
            
            prev_cov = prev_meth + prev_unmeth
            curr_cov = curr_meth + curr_unmeth
            combined_evenness = max(prev_cov, curr_cov) / (prev_cov + curr_cov)
            delta_meth_pct = abs(prev_meth / (prev_meth + prev_unmeth) - \
                                 curr_meth / (curr_meth + curr_unmeth)) * 100
            gc_pct = calc_gc_pct(genome_int[prev_scaf], prev_pos)
            seq_context =  get_cpg_context(genome_seq[prev_scaf], prev_pos)
            
            print (prev_scaf, prev_pos, curr_pos,
                   f'{combined_beta:.3f}',
                   combined_meth, combined_unmeth,
                   f'{combined_evenness:.6f}', f'{delta_meth_pct:.3f}',
                   gc_pct, seq_context, sep='\t', file=outf)
            
            prev_line = []
            curr_line = []
        else:
            # "prev_line" and "curr_line" are from a different dinuc
            if genome_int[prev_scaf][prev_pos - 1] == 2:
                dinuc_start, dinuc_end = prev_pos, prev_pos + 1
            else:
                dinuc_start, dinuc_end = prev_pos - 1, prev_pos
            
            # in this case, coverage is totally dominated by a strand (1),
            # and delta_meth_pct == prev_meth_pct
            prev_meth_pct = prev_meth / (prev_meth + prev_unmeth) * 100
            gc_pct = calc_gc_pct(genome_int[prev_scaf], dinuc_start)
            seq_context =  get_cpg_context(genome_seq[prev_scaf], dinuc_start)
            
            print (prev_scaf, dinuc_start, dinuc_end,
                   f'{prev_meth_pct:.3f}',
                   prev_meth, prev_unmeth,
                   1, f'{prev_meth_pct:.3f}', 
                   gc_pct, seq_context, sep='\t', file=outf)
            
            prev_line = []
    
    bc_file.close()
