#!/usr/bin/env python3
"""Combine exploratory pooled sensitivity analyses into one table (Additional file 8)."""
import csv
from pathlib import Path
A = Path('./analysis')
cols = ['analysis','method','k','estimate','ci_lower','ci_upper','pi_lower','pi_upper','i2','p_above0','p_above005','by_above0','by_above005']
labels = {'primary': 'Primary prespecified synthesis', 'hemolysis_associated': 'Excluding four hemolysis-associated units',
          'repeated_profile_units': 'Excluding the two units with repeated profiles per person',
          'unresolved_or_normalized_input': 'Excluding units with unresolved scale or timing or normalized input',
          'cancer_status_unconfirmed': 'Excluding contrasts without confirmed cancer status', 'toray': 'Excluding the Toray unit',
          'gse85589_log_scale': 'Excluding GSE85589', 'linear_gse85589': 'GSE85589 re-run on linear scale'}
rows = []
for r in csv.DictReader(open(A/'sensitivity_exclusions/sensitivity_exclusions.tsv'), delimiter='\t'):
    r['analysis'] = labels[r['analysis']]; rows.append({c: r.get(c, '') for c in cols})
for r in csv.DictReader(open(A/'gse85589_linear_compute/pooled_with_linear_gse85589.tsv'), delimiter='\t'):
    key = 'linear_gse85589'
    rows.append({c: (labels[key] if c == 'analysis' else r.get(c, '')) for c in cols})
for r in rows:
    for c in cols[2:]:
        if r[c] in ('NA', 'nan'): r[c] = ''
out = A/'sensitivity_exclusions/pooled_sensitivity_combined.tsv'
with open(out, 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=cols, delimiter='\t'); w.writeheader(); w.writerows(rows)
print(len(rows), 'rows ->', out)
