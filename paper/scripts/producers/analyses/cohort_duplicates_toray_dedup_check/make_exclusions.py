#!/usr/bin/env python3
"""Build evaluation_exclusions.tsv from the near-duplicate audit outputs (coordinator, 2026-10-12).
Inputs (same folder): toray_duplicate_rows.tsv (from nn_pairs.R output via the R one-liner recorded in README.md:
rows of the 197 Toray pairs) and all_units_near_duplicate_pairs.tsv (from all_cohorts.R)."""
from run import analysis_directory as _public_analysis_directory, public_table as _public_table
import csv, sys
from pathlib import Path
d = Path(__file__).resolve().parent
rows = list(csv.DictReader(open(d / 'toray_duplicate_rows.tsv'), delimiter='\t'))
pr = list(csv.DictReader(open(d / 'all_units_near_duplicate_pairs.tsv'), delimiter='\t'))
out = d / (_public_analysis_directory(sys.argv[1] if len(sys.argv) > 1 else 'evaluation_exclusions.tsv'))
with open(out, 'w') as f:
    f.write('unit\trow_name\tsample_id\tgroup_id\tpair_id\tcentred_log2_distance\treason\n')
    for r in rows:
        f.write(f"TORAY_PUBLIC_CANONICAL\t{r['row']}\t{r['row'].split('::')[-1]}\t{r['group']}\tT{r['pair']}\t{r['nn_dist']}\tnear-identical profile deposited under another accession (residual duplicate after identifier and exact-profile deduplication)\n")
    for i, p in enumerate(pr):
        for k, g in (('a', 'ga'), ('b', 'gb')):
            f.write(f"{p['unit']}\t{p[k]}\t{p[k]}\t{p[g]}\tP{i+1}\t{p['d']}\texact duplicate profile under two sample accessions\n")
