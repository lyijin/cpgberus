#!/usr/bin/env python3

doc = """
> parse_mpileup_tsvs.py <

Script parses mpileup tsvs (or compressed tsvs) and outputs, in tabular form:

  0. scaffold name
  1. coord
  2. total per-base coverage
  3. reads mapping to Watson strand (+) of reference
  4. reads mapping to Crick strand (-) of reference
  5. number of A bases at coordinate
  6. number of C bases at coordinate
  7. number of G bases at coordinate
  8. number of T bases at coordinate
  (note: Ns are not considered)

Input tables are assumed to be produced with the command
  `samtools mpileup -A -B -d 0 -Q 0 -x ${file} | gzip > ${output_file}`

The reason behind the flags is that `samtools mpileup` tries to filter for
"good quality" bases, but `samtools depth` just calculates the depth, quality
be damned. Using both tools with default settings will result in different
numbers, see https://github.com/samtools/samtools/issues/1129.

`samtools mpileup` settings:
  -A, --count-orphans     do not discard anomalous read pairs
  -B, --no-BAQ            disable BAQ (per-Base Alignment Quality)
  -d, --max-depth INT     max per-file depth; avoids excessive memory usage [8000]
  -Q, --min-BQ INT        skip bases with baseQ/BAQ smaller than INT [13]
  -x, --ignore-overlaps   disable read-pair overlap detection
""".strip()

import argparse
import csv
import gzip
from pathlib import Path
import re
import sys
import time

def benchmark_print(message):
    """
    Prints a timestamped message to stderr, to aid visual assessment of
    chokepoints in code.
    """
    print (f'[{time.asctime()}] {message}', file=sys.stderr)

parser = argparse.ArgumentParser(description=doc,
                                 formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('mpileup_tsv', metavar='mpileup_tsv', type=Path,
                    nargs='?', default=sys.stdin,
                    help='tsv/gzip-compressed tsv file.')
parser.add_argument('-v', '--verbose', action='store_true',
                    help='prints diagnostic stuff to stderr.')

args = parser.parse_args()

# check whether input file has ".gz" at the end; if it has, then file is 
# assumed to be gzip-compressed
if args.mpileup_tsv.name == '<stdin>':
    tsv_reader = csv.reader(args.mpileup_tsv, delimiter='\t', quoting=csv.QUOTE_NONE)
elif args.mpileup_tsv.suffix == '.gz':
    tsv_reader = csv.reader(gzip.open(args.mpileup_tsv, 'rt'), delimiter='\t',
                            quoting=csv.QUOTE_NONE)
else:
    tsv_reader = csv.reader(open(args.mpileup_tsv), delimiter='\t',
                            quoting=csv.QUOTE_NONE)

if args.verbose:
    benchmark_print('Parsing {args.mpileup_tsv.name}...')

lineno = 0
for row in tsv_reader:
    # indicate progress
    lineno += 1
    if args.verbose and lineno % 1_000_000 == 0:
        benchmark_print('Reading line {lineno:,}...')
    
    scaf, pos, refbase, depth, bases, *_ = row
    pos = int(pos)
    depth = int(depth)
    
    # https://en.wikipedia.org/wiki/Pileup_format
    # `mpileup` adds additional text to denote "start of read" (^.) and
    # "end of read" ($). for "start of read", the ASCII character following
    # caret is the Phred+33 quality of the base. strip these text to retain
    # sequenced bases
    bases = re.sub(r'\^.', '', bases)
    bases = bases.replace('$', '')
    
    # "*" means base was deleted, but mentioned in the previous lines. both
    # `samtools depth` and `samtools mpileup` considers "*" as +1 coverage for
    # said base, but as I can't assign "*" to Watson/Crick/ACGT buckets, I'm
    # more inclined to ignore them. i.e., this can be tested and should be True
    #   n_star = bases.count('*')
    #   assert depth == n_total + n_star (after n_total is calculated)
    bases = bases.replace('*', '')
    
    # code to deal with indels, which complicate per-base coverage calcs
    # greatly. indels are coded e.g. '+2AG' (insertion of 'AG' in Watson)
    # or '-2ct' (deletion of 'CT' on Crick). in both cases, it has the pattern
    #   [- or +] [numbers] [text]
    # and [numbers] is the # of characters stripped out at the start of [text]
    bases = bases.replace('-', '+')
    bases = bases.split('+')
    
    temp = ''
    for b in bases:
        match = re.match(r'(\d+)(\w+)', b)
        if match:
            n = int(match[1])
            temp += match[2][n:]
        else:
            temp += b
    
    bases = temp
    
    n_watson = bases.count('A') + bases.count('C') + bases.count('G') + bases.count('T') + bases.count('.')
    n_crick = bases.count('a') + bases.count('c') + bases.count('g') + bases.count('t') + bases.count(',')
    n_total = n_watson + n_crick
    
    n_refbase = bases.count('.') + bases.count(',')
    
    n_a = bases.count('A') + bases.count('a') + (n_refbase if refbase in 'Aa' else 0)
    n_c = bases.count('C') + bases.count('c') + (n_refbase if refbase in 'Cc' else 0)
    n_g = bases.count('G') + bases.count('g') + (n_refbase if refbase in 'Gg' else 0)
    n_t = bases.count('T') + bases.count('t') + (n_refbase if refbase in 'Tt' else 0)
    
    print (scaf, pos, n_total, n_watson, n_crick, n_a, n_c, n_g, n_t, sep='\t')
