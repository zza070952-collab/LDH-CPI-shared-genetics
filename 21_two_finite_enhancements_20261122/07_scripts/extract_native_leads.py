"""Extract exactly the locked lead rows from the public original release on GRCh38."""
import csv,gzip,json
from pathlib import Path
n=Path(__file__).resolve().parents[1]
leads=set(json.loads((n/'00_admin/ANALYSIS_PLAN_LOCK.json').read_text())['leads'])
source=n.parents[2]/'shared_resources/raw/finngen_r12/finngen_R12_M13_LUMBAR_PROLAPSE.gz'
rows=[]
with gzip.open(source,'rt') as f:
    reader=csv.DictReader(f,delimiter='\t')
    for row in reader:
        if row['rsids'] in leads:rows.append(row)
    fields=reader.fieldnames
assert len(rows)==len(leads)==18 and {r['rsids'] for r in rows}==leads
with (n/'02_external/LDH_RAW_GRCH38_LEADS.tsv').open('w',newline='') as f:
    writer=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(rows)

