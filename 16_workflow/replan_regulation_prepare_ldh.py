#!/usr/bin/env python3
"""Prepare FinnGen-only LDH SMR input on the full eight-cell QTL SNP universe.

Uses audited rsID+allele mapping to GRCh37 EUR reference; does not relabel the
original GRCh38 coordinates as hg19 and is not a full-GWAS liftover product.
"""
from pathlib import Path
import argparse,csv,gzip,hashlib,json,math,re
from collections import Counter
P2=Path('/analysis/projects/JTM_recent_highIF_direct_clone/P2_LDH_shared_architecture')
R=P2/'11_NatCommun_strict_clone_Pro';M=R/'08_smr_heidi_coloc';Q=M/'bryois_full_smr_universe'
B=Path('/analysis/projects/hipOA_spatial_multiomics_clean/resources/active/extracted/magma_g1000_eur/g1000_eur')
G=P2/'01_data_raw/finngen_R12_M13_LUMBAR_PROLAPSE.gz'
EXPECTED='7f2d32eaf6040914bea07c4193775923a20a5cffc5ef6f00e004ca4d1e278392'
CELLS=['Astrocytes','Endothelial.cells','Excitatory.neurons','Inhibitory.neurons','Microglia','Oligodendrocytes','OPCs...COPs','Pericytes']
p=argparse.ArgumentParser();p.add_argument('--ld-checksum-manifest');p.add_argument('--normalized-gwas');p.add_argument('--normalized-sha256');a=p.parse_args()
assert not(R/'00_admin/PAUSE_REQUESTED').exists(),'PAUSE_REQUESTED: no new stage'
for c in CELLS:assert Path(str(Q/c)+'.conversion.done').exists(),f'Full eight-panel QTL universe not complete: {c}'
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()
checks=[]
if a.normalized_gwas:
    SRC=Path(a.normalized_gwas);h=sha(SRC);assert a.normalized_sha256 and h==a.normalized_sha256,'Unapproved normalized GWAS hash';checks.append(dict(file=str(SRC),sha256=h))
else:
    SRC=G;h=sha(G);assert h==EXPECTED,'Unapproved GWAS hash';checks.append(dict(file=str(G),sha256=h))
needed=set()
for c in CELLS:
    fn=Path(str(Q/c)+'.esi');checks.append(dict(file=str(fn),sha256=sha(fn)))
    with fn.open() as f:
        for ln in f:needed.add(ln.split()[1])
ref={};amb=set()
with Path(str(B)+'.bim').open() as f:
    for ln in f:
        x=ln.split()
        if x[1] not in needed:continue
        rr=(int(x[0]),int(x[3]),x[4],x[5])
        if x[1] in ref and ref[x[1]]!=rr:amb.add(x[1])
        ref[x[1]]=rr
for s in amb:ref.pop(s,None)
qc=Counter();kept={};duplicates=set();locations={};comp=str.maketrans('ACGT','TGCA')
with gzip.open(SRC,'rt') as f:
    for r in csv.DictReader(f,delimiter='\t'):
        qc['input_rows']+=1
        rsids=r['SNP'] if a.normalized_gwas else r['rsids']
        candidates=[s for s in rsids.replace(';',',').split(',') if s in ref]
        if not candidates:continue
        qc['QTL_union_rsID_match']+=1
        try:
            if a.normalized_gwas:c=int(r['CHR']);bp38='NA';b=float(r['BETA']);se=float(r['SE']);pv=float(r['P']);af=float(r['EAF'])
            else:c=int(r['#chrom']);bp38=int(r['pos']);b=float(r['beta']);se=float(r['sebeta']);pv=float(r['pval']);af=float(r['af_alt'])
        except ValueError:qc['invalid_numeric']+=1;continue
        a1=(r['A1'] if a.normalized_gwas else r['alt']).upper();a2=(r['A2'] if a.normalized_gwas else r['ref']).upper()
        if c not in range(1,23) or len(a1)!=1 or len(a2)!=1 or a1 not in 'ACGT' or a2 not in 'ACGT' or a1==a2:qc['invalid_allele_or_chr']+=1;continue
        if not(math.isfinite(b) and math.isfinite(se) and se>0 and 0<pv<=1 and .01<af<.99):qc['invalid_beta_SE_P_or_MAF']+=1;continue
        if a1+a2 in ('AT','TA','CG','GC'):qc['palindrome_excluded']+=1;continue
        snp=next((s for s in candidates if c==ref[s][0] and (set([a1,a2])==set(ref[s][2:]) or set([a1.translate(comp),a2.translate(comp)])==set(ref[s][2:]))),None)
        if snp is None:qc['allele_or_chromosome_mismatch']+=1;continue
        rr=ref[snp]
        if a.normalized_gwas and int(r['BP'])!=rr[1]:qc['normalized_hg19_vs_LD_position_mismatch']+=1;continue
        if c==6 and 25e6<=rr[1]<=35e6:qc['MHC_snp_excluded']+=1;continue
        if snp in kept:duplicates.add(snp);continue
        loc=(c,rr[1])
        if loc in locations:duplicates.add(snp);duplicates.add(locations[loc])
        locations[loc]=snp
        kept[snp]=dict(SNP=snp,A1=a1,A2=a2,freq=af,b=b,se=se,p=pv,N=388514,chr=c,bp_hg19=rr[1],bp_original_hg38=bp38)
for snp in duplicates:kept.pop(snp,None)
qc['duplicate_snp_or_position_excluded']=len(duplicates);qc['final_SMR_rows']=len(kept);qc['reference_ambiguous_rsID']=len(amb)
out=M/'gwas';out.mkdir(exist_ok=True)
with (out/'LDH_FinnGen_R12.QTL_universe.smr.ma').open('w') as f:
    w=csv.DictWriter(f,fieldnames=['SNP','A1','A2','freq','b','se','p','N'],delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(kept.values())
with gzip.open(out/'LDH_FinnGen_R12.QTL_universe.coordinate_audit.tsv.gz','wt') as f:
    w=csv.DictWriter(f,fieldnames=list(next(iter(kept.values()))),delimiter='\t');w.writeheader();w.writerows(kept.values())
cached={}
if a.ld_checksum_manifest:
    for r in csv.DictReader(open(a.ld_checksum_manifest),delimiter='\t'):
        st=Path(r['file']).stat()
        if st.st_size==int(r['size']) and st.st_mtime_ns==int(r['mtime_ns']):cached[r['file']]=r['sha256']
for ext in ['.fam','.bim','.bed']:
    f=str(B)+ext;checks.append(dict(file=f,sha256=cached[f] if f in cached else sha(f)))
with (out/'LDH_SMR_INPUT_CHECKSUMS.tsv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=['file','sha256'],delimiter='\t');w.writeheader();w.writerows(checks)
meta=dict(trait='LDH',cohorts='FinnGen_R12_only',N=388514,cases=35290,controls=353224,original_build='GRCh38',module_reference_build='GRCh37',coordinate_method='Approved full-GWAS liftOver/QC followed by rsID+allele+BP match' if a.normalized_gwas else 'rsID_plus_alleles_to_audited_reference; NOT_genomewide_liftOver',source=str(SRC),reference=str(B),reference_sample_count=sum(1 for l in Path(str(B)+'.fam').open()),effect_allele='ALT harmonized to target forward strand',effect_frequency='af_alt; EAF not MAF',INFO='NOT_REPORTED_IN_SOURCE',qc=dict(qc),coloc_eligible='NO_QTL_HEIDI_UNION_ONLY; use full-window normalized GWAS extraction for coloc')
(out/'LDH_SMR_INPUT_QC.json').write_text(json.dumps(meta,indent=2));(out/'LDH_SMR_INPUT.done').write_text('SUCCESS\n');print(json.dumps(meta,indent=2))
