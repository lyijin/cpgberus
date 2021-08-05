#!/usr/bin/env python3

doc = """
> tabulate_cov_coverages_methpct.grch38.py <

Script uses pandas to outer-join `bismark` cov files provided as input, and
produces a compressed tab-separated file as output.

Input files should not have a header, output files WILL have headers to make
downstream table reading easier.

Output format:
  1. scaffold name
  2. starting pos (1-based, like the input cov files)
  3. ending pos (same as starting pos)
  4. file1
  5. file2
  ...

Autodetects whether input files are gzip-compressed from `.gz` suffix.
""".strip()

import argparse
import gzip
from pathlib import Path
import sys
import time

import pandas as pd

def benchmark_print(message):
    """
    Prints a timestamped message to stderr, to aid visual assessment of
    chokepoints in code.
    """
    print (f'[{time.asctime()}] {message}', file=sys.stderr)

parser = argparse.ArgumentParser(description=doc, 
                                 formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('bismark_covs', type=Path, nargs='+',
                    help='one or more bismark cov files (can be compressed).')
parser.add_argument('-v', '--verbose', action='store_true',
                    help='prints diagnostic stuff to stderr.')

args = parser.parse_args()

# start parsing individual bismark covs
cov_df = pd.DataFrame()
methpct_df = pd.DataFrame()
for n, bc in enumerate(args.bismark_covs):
    if args.verbose: benchmark_print(f'Processing {bc.name}...')
    
    bc_id = bc.name.split('.')[0]
    
    bc_df = pd.read_table(
        bc, header=None,
        names=['scaffold', 'start_pos', 'end_pos', 
               f'{bc_id}_methpct', 'meth_reads', 'unmeth_reads'])
    
    # calculate coverage
    bc_df[f'{bc_id}_cov'] = bc_df['meth_reads'] + bc_df['unmeth_reads']
    
    # no need to merge when dealing with first file
    if not n:
        cov_df = bc_df[['scaffold', 'start_pos', 'end_pos', f'{bc_id}_cov']]
        methpct_df = bc_df[['scaffold', 'start_pos', 'end_pos', f'{bc_id}_methpct']]
    else:
        cov_df = pd.merge(
            cov_df,
            bc_df[['scaffold', 'start_pos', 'end_pos', f'{bc_id}_cov']],
            how='outer')
        methpct_df = pd.merge(
            methpct_df,
            bc_df[['scaffold', 'start_pos', 'end_pos', f'{bc_id}_methpct']],
            how='outer')

# irritatingly, anytime NaN appears, the column will be casted to np.float64
# tidy up tables after the merging
cov_df = cov_df.fillna(0)
cov_cols = [x for x in cov_df.columns if '_cov' in x]
cov_df[cov_cols] = cov_df[cov_cols].astype(int)
cov_df.columns = [x.replace('_cov', '') for x in cov_df.columns]

# for meth %ages, reduce the number of decimal points, save some storage space;
# then coerce into string to ensure saving floats won't result in stuff like
# 3.33000000001 later on
methpct_cols = [x for x in methpct_df.columns if '_methpct' in x]
methpct_df[methpct_cols] = methpct_df[methpct_cols].round(3).astype(str)
methpct_df.columns = [x.replace('_methpct', '') for x in methpct_df.columns]

if args.verbose:
    benchmark_print(f'Union of all files produces {len(cov_df)} rows.')

# print to designated output files
if args.verbose: benchmark_print(f'Writing coverages out...')
cov_df.to_csv('all.coverages.tsv.gz', sep='\t', index=False)

if args.verbose: benchmark_print(f'Writing methylation percentages out...')
methpct_df.to_csv('all.methpcts.tsv.gz', sep='\t', index=False)
