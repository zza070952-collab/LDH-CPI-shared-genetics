from pathlib import Path
import csv,json,hashlib
N=Path(__file__).resolve().parents[1]
A=N/'00_preanalysis/LDH_CPI_Preanalysis_20261122'
read=lambda p:list(csv.DictReader(p.open(encoding='utf-8-sig'),delimiter='\t'))
for r in read(A/'FILE_MANIFEST_SHA256.tsv'):
    p=A/r['path']
    assert p.stat().st_size==int(r['bytes']) and hashlib.sha256(p.read_bytes()).hexdigest()==r['sha256'],p
m=read(A/'frozen_20261113/06_supplement/GWAS_SOURCE_MANIFEST.tsv')
lock=dict(version='20261122; after exploratory orientation-free preanalysis',GWAS=[{k:r[k] for k in ['resource_id','N','Neff','full_QC_path','full_QC_build','full_QC_sha256']} for r in m if r['resource_id'] in ['FINNGEN_R12_LDH','UKB_CPI_WANG2025']],tissues=['Nerve_Tibial','Brain_Cortex','Brain_Spinal_cord_cervical_c-1','Whole_Blood'],primary='MASHR',sensitivity='EN standard; CPI phi conditional on documented quantitative N and observed h2; LDH phi withheld',coverage='>=80% nonzero weights, complete if n<=2; compatible unique rsID and allele; covariance positive',families='per trait and model family BH pooled across four tissues; common-predictor intersection-union max P, BH across eligible shared predictors; compare cross-tissue gene overlap separately',leads=[r['lead_snp'] for r in read(A/'frozen_20261113/01_evidence_tables/LOCUS_MASTER_18.tsv')],external='MVP EUR pha013136; EstBB only published sentinel ancillary; no overlapping meta; exact SNP preferred, high LD r2>=0.8, phase needed for direction; Bonferroni denominator 18',download='one repair and one official alternate; do not retrieve unattributed multi-GB archives')
p=N/'00_admin/ANALYSIS_PLAN_LOCK.json';data=json.dumps(lock,indent=2).encode()
if p.exists():assert p.read_bytes()==data
else:p.write_bytes(data)
print('Evidence and lock verified',hashlib.sha256(data).hexdigest())
