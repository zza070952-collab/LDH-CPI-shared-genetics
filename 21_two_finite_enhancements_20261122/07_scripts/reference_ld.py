"""Extract pairwise LD for fixed-region published sentinels using the frozen EUR503 reference."""
from pathlib import Path
import csv,hashlib,zipfile,json
import numpy as np
N=Path(__file__).resolve().parents[1];A=next((N/'00_preanalysis').glob('LDH_*'));D=N/'02_external'
read=lambda p:list(csv.DictReader(p.open(),delimiter='\t'))
pairs=set()
for r in read(A/'external/MVP_ST4_fixed18_region_overlaps.tsv'):pairs.add((r['primary_lead'],r['rsID']))
for r in read(A/'external/ESTBB_published64_per_variant_effects.tsv'):
    for x in read(A/'frozen_20261113/01_evidence_tables/LOCUS_MASTER_18.tsv'):
        if x['locus_id'] in r['current_regions']:pairs.add((x['lead_snp'],r['rsid']))
for rs in ['rs4456576','rs12185377','rs4384683','rs17487130','rs12956563']:pairs.add(('rs62100775',rs))
arc=Path('/analysis/projects/sciatica_multomics_bakeoff/raw/phase3_final/magma/g1000_eur.zip')
h=hashlib.sha256()
with arc.open('rb') as fp:
    for part in iter(lambda:fp.read(4<<20),b''):h.update(part)
assert h.hexdigest()=='83a48fd9dcaa0b9a874b18c63143a4ede93f05505b215b0bd8790130a0d7a954'
ids=set.union(*(set(p) for p in pairs));idx={};dose={}
with zipfile.ZipFile(arc) as z:
    with z.open('g1000_eur.bim') as f:
        for i,line in enumerate(f):
            r=line.decode().split()
            if r[1] in ids:idx[i]=r
    with z.open('g1000_eur.bed') as f:
        assert f.read(3)==bytes([108,27,1])
        for i,r in sorted(idx.items()):
            f.seek(3+126*i);buf=np.frombuffer(f.read(126),dtype=np.uint8)
            raw=((buf[:,None]>>np.array([0,2,4,6]))&3).ravel()[:503]
            dose[r[1]]=(np.array([2,np.nan,1,0])[raw],r)
rows=[]
for a,b in sorted(pairs):
    v=np.nan
    if a in dose and b in dose:
        x,ra=dose[a];y,rb=dose[b];ok=np.isfinite(x+y);v=np.corrcoef(x[ok],y[ok])[0,1]
    rows.append([a,b,v,v*v,dose[a][1][4] if a in dose else 'NA',dose[b][1][4] if b in dose else 'NA','AVAILABLE' if np.isfinite(v) else 'NOT_IN_REFERENCE'])
with (D/'FROZEN_EUR503_PAIRWISE_LD.tsv').open('w') as f:
    w=csv.writer(f,delimiter='\t');w.writerow(['lead','other','signed_r_A1_A1','r2','lead_A1','other_A1','status']);w.writerows(rows)
print(len(rows),'pairs')
