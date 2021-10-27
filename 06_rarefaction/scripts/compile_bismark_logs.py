#!/usr/bin/env python3

"""
> compile_bismark_logs.py <

Compile a table from all of the logs produced by Bismark during the mapping-
deduplication-methylation extractor steps.

Script does not take any arguments, expects folders/filenames to follow what
was produced by the `bismsmark` pipeline.
"""
import csv
import gzip
from pathlib import Path
import re

# create dict to store all stats
bismark_stats = {}

# stop script if there aren't valid files in  `05_renamed_covs/` subfolder
renamed_covs = list(Path('05_renamed_covs/').glob('*.cov.gz'))
assert renamed_covs, 'Could not find "*.cov.gz" files in the `05_renamed_covs` subfolder!'

# work backwards--check what the final names are in `05_renamed_cov/`, then
# track their provenances...
for rc in renamed_covs:
    # parse the final, renamed version of the covs
    sample_id = rc.name.split('.')[0]
    mapped_genome = rc.name.split('.')[1]
    unique_id = f'{sample_id}|{mapped_genome}'
    
    bismark_stats[unique_id] = {}
    bismark_stats[unique_id]['00|sample_id'] = sample_id
    bismark_stats[unique_id]['01|mapped_genome'] = mapped_genome
    
    # .resolve() follows the symbolic link to the actual .cov.gz. use
    # spaghetti code to discover original ID associated with the file
    orig_id = rc.resolve().name.replace('.deduplicated.bismark.cov.gz', '')
    bismark_stats[unique_id]['02|orig_id'] = orig_id
    
    # parse `bismark_methylation_extraction` log
    meth_extract_log = Path(f'04_meth_extract_{mapped_genome}/{orig_id}.deduplicated_splitting_report.txt')
    log_content = open(meth_extract_log).read()
    bismark_stats[unique_id]['bme_read_counts'] = \
        int(re.search(r'Processed (\d+) lines in total', log_content)[1])
    bismark_stats[unique_id]['09|CpG_meth'] = \
        re.search(r'C methylated in CpG context:\s+(.+%)', log_content)[1]
    bismark_stats[unique_id]['10|CHG_meth'] = \
        re.search(r'C methylated in CHG context:\s+(.+%)', log_content)[1]
    bismark_stats[unique_id]['11|CHH_meth'] = \
        re.search(r'C methylated in CHH context:\s+(.+%)', log_content)[1]
    
    # calculate lambda/pUC19 stats from the cov, which are expected to be
    # non-methylated and fully methylated respectively. this only applies
    # if mapping was done vs. "hg38_lambda_puc"
    bismark_stats[unique_id]['12|lambda_reads'] = 'NA'
    bismark_stats[unique_id]['13|lambda_meth'] = 'NA'
    bismark_stats[unique_id]['14|pUC19_reads'] = 'NA'
    bismark_stats[unique_id]['15|pUC19_meth'] = 'NA'
    if mapped_genome == 'grch38p13_lambda_puc':
        n_reads = {'lambda': 0, 'pUC19': 0} 
        n_meth = {'lambda': 0, 'pUC19': 0}
        n_coverage = {'lambda': 0, 'pUC19': 0}
        
        tsv_reader = csv.reader(gzip.open(rc, 'rt'), delimiter='\t')
        for line in tsv_reader:
            if line[0] not in ['lambda', 'pUC19']: continue
            
            n_reads[line[0]] += 1
            n_meth[line[0]] += int(line[4])
            n_coverage[line[0]] += int(line[4]) + int(line[5])
        
        if n_reads['lambda']:
            bismark_stats[unique_id]['12|lambda_reads'] = n_reads['lambda']
            bismark_stats[unique_id]['13|lambda_meth'] = \
                f"{n_meth['lambda'] / n_coverage['lambda'] * 100:.2f}%"
        
        if n_reads['pUC19']:
            bismark_stats[unique_id]['14|pUC19_reads'] = n_reads['pUC19']
            bismark_stats[unique_id]['15|pUC19_meth'] = \
                f"{n_meth['pUC19'] / n_coverage['pUC19'] * 100:.2f}%"
    
    # parse `deduplicate_bismark` logs. deduplication is SKIPPED for RRBS data,
    # hence the if/else to handle that
    dedup_log = Path(f'03_dedup_{mapped_genome}/{orig_id}.deduplication_report.txt')
    if dedup_log.is_file():
        log_content = open(dedup_log).read()
        
        pre_dedup = re.search(r'Total number of alignments.*?:\s+(\d+)', log_content)[1]
        post_dedup = re.search(r'leftover sequences:\s+(\d+)', log_content)[1]
        bismark_stats[unique_id]['06|pre-dedup_reads'] = f'{int(pre_dedup):,}'
        bismark_stats[unique_id]['07|post-dedup_reads'] = f'{int(post_dedup):,}'
        bismark_stats[unique_id]['08|post-dedup_pct_kept'] = \
            f'{int(post_dedup) / int(pre_dedup) * 100:.2f}%'
    else:
        bismark_stats[unique_id]['06|pre-dedup_reads'] = \
            f"{bismark_stats[unique_id]['bme_read_counts']:,}"
        bismark_stats[unique_id]['07|post-dedup_reads'] = \
            f"{bismark_stats[unique_id]['bme_read_counts']:,}"
        bismark_stats[unique_id]['08|post-dedup_pct_kept'] = 'NA'
    
    # parse `bismark` log
    lib_type_upper = orig_id[:-3] + orig_id[-3:].upper()
    mapping_log = Path(f'02_map_vs_{mapped_genome}/{lib_type_upper}_report.txt')
    
    log_content = open(mapping_log).read()
    n_trim = int(re.search(r'Sequence.*?analysed in total:\s+(\d+)', log_content)[1])
    n_align = int(re.search(r'Number of.*?alignments with.*?:\s+(\d+)', log_content)[1])
    
    bismark_stats[unique_id]['03|trimmed_reads'] = f'{n_trim:,}'
    bismark_stats[unique_id]['04|aligned_reads'] = f'{n_align:,}'
    bismark_stats[unique_id]['05|mapping_rate'] = f'{n_align / n_trim * 100:.2f}%'

# print parsed information into a table
impt_columns = sorted(x for x in bismark_stats[unique_id].keys() if '|' in x)

# remove "nn|..." from every element in `impt_columns`
print (*[x[3:] for x in impt_columns], sep='\t')

for uid in sorted(bismark_stats):
    print (*[bismark_stats[uid][x] for x in impt_columns], sep='\t')
