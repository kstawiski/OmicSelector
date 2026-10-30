"""Stream public inputs; keep only named features and the authors' first filter."""
from run import analysis_directory as _public_analysis_directory, public_table as _public_table
import csv, gzip, json, pathlib, hashlib, socket, platform, datetime, os
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
P = pathlib.Path('.')
O = ROOT / 'inputs'
O.mkdir(exist_ok=True)
rules = json.loads((P/'analysis/opportunity/manuscript_drive/mirbit_fidelity/rules.json').read_text())
features = [r[k] for r in rules['rules'] for k in ('feature_a','feature_b')]
rows, names, chosen = [], [], {}
source = P/'data_bucket/GSE211692/raw/GSE211692_processed_data_matrix.txt.gz'
with gzip.open(source, 'rt') as f:
    header = next(csv.reader([next(f)], delimiter='\t'))[1:]
    for line in f:
        name, _, data = line.rstrip('\n').partition('\t')
        name = name.strip('"')
        a = np.fromstring(data, sep='\t')
        assert len(a) == len(header), (name, len(a))
        if np.mean(a) > 5:
            rows.append(a); names.append(name)
        if name in features:
            assert name not in chosen
            chosen[name] = a
assert len(chosen) == 8
with (O/'train_eight.tsv').open('w') as f:
    w=csv.writer(f,delimiter='\t'); w.writerow(['title']+features)
    w.writerows(zip(header, *[chosen[n] for n in features]))
with (O/'train_filter.tsv').open('w') as f:
    w=csv.writer(f,delimiter='\t'); w.writerow(['feature']+header)
    w.writerows([n]+list(a) for n,a in zip(names,rows))
(O/'rules.json').write_text(json.dumps(rules,indent=2))
for src,name in [(P/'data_bucket/GSE211692/metadata.tsv','train_metadata.tsv'),
                 (ROOT/'out/evidence/ST2.csv','ST2.csv'),
                 (ROOT/'out/evidence/ST4.csv','ST4.csv'),
                 (ROOT/'out/evidence/meta_train_full.csv','meta_train_full.csv'),
                 (ROOT/'out/evidence/metadata85589test.csv','metadata85589test.csv')]:
    (O/_public_analysis_directory(name)).write_bytes(src.read_bytes())
print(json.dumps({'features_first_filter':len(names),'samples':len(header),'feature_min':{n:float(chosen[n].min()) for n in features}}))
(ROOT/'out/evidence/preparation_identity.json').write_text(json.dumps({'hostname':socket.gethostname(),'python':platform.python_version(),'numpy':np.__version__,'cpu_threads':1,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source':str(source),'purpose':'streaming input preparation'},indent=2))
