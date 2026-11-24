"""Prepare the unchanged fixed chr18 region using existing GWAS and EUR archive.

No analysis threshold depends on an observed SuSiE or colocalization result.
No reference download or historical-file mutation occurs.
"""
import csv
import gzip
import hashlib
import json
import shutil
import zipfile
from pathlib import Path
import numpy as np
import pandas as pd

R=Path('/analysis/projects/JTM_recent_highIF_direct_clone/P2_LDH_shared_architecture/11_NatCommun_strict_clone_Pro')
NEW=R/'19_computational_finalize_JTM_20261103'
O=NEW/'03_signal_resolution'
O.mkdir(parents=True,exist_ok=True)
SRC=R/'17_evidence_tier_revision_20260912/01_locus_evidence/LDH_CPI_p05_chr18_001.harmonized.tsv.gz'
ARC=Path('/analysis/projects/sciatica_multomics_bakeoff/raw/phase3_final/magma/g1000_eur.zip')
ARC_SHA='83a48fd9dcaa0b9a874b18c63143a4ede93f05505b215b0bd8790130a0d7a954'

def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for chunk in iter(lambda:f.read(8*1024**2),b''):h.update(chunk)
    return h.hexdigest()

sig={'source':str(SRC),'source_sha256':sha(SRC),'reference_archive':str(ARC),'reference_archive_sha256':ARC_SHA,'window':[18,50358109,51055069],'build':'GRCh37','LDH_N':388514,'CPI_N_scheme':'source min/median/max; no best-N selection'}
done=O/'CHR18_PREPARED.json'
if done.exists():
    old=json.loads(done.read_text())
    assert old['signature']==sig
    assert all(sha(Path(p))==v for p,v in old['outputs'].items())
    print('REUSE_VERIFIED_CHR18_PREPARATION',flush=True)
    raise SystemExit(0)
assert sha(ARC)==ARC_SHA
d=pd.read_csv(SRC,sep='\t')
assert len(d)==1357 and d.SNP.is_unique
assert d.CHR_LDH.eq(18).all() and d.CHR_CPI.eq(18).all()
assert d.BP_LDH.between(50358109,51055069).all()
assert np.array_equal(d.BP_LDH,d.BP_CPI)
# Historical revision_locus_coloc.R already flips BETA_CPI to LDH alleles but
# deliberately retains original CPI allele-label columns. Do not flip twice.
assert d.alignment.isin(['SAME','FLIP','COMPLEMENT','COMPLEMENT_FLIP']).all()
assert d.N_LDH.eq(388514).all()
assert np.all(np.isfinite(d[['BETA_LDH','SE_LDH','BETA_CPI','SE_CPI']]))
assert (d[['SE_LDH','SE_CPI']]>0).all().all()
ids=set(d.SNP);bims=[];indices=[];seen=set()
with zipfile.ZipFile(ARC) as z:
    fam=z.read('g1000_eur.fam').decode().splitlines()
    assert len(fam)==503
    with z.open('g1000_eur.bim') as h:
        for i,line in enumerate(h):
            row=line.decode().split()
            if row[1] in ids:
                assert row[1] not in seen
                bims.append(row);indices.append(i);seen.add(row[1])
    bm=pd.DataFrame(bims,columns=['CHR','SNP','CM','BP','A1','A2'])
    bm.CHR=bm.CHR.astype(int);bm.BP=bm.BP.astype(int)
    a=d.set_index('SNP').loc[bm.SNP].reset_index()
    assert np.array_equal(a.BP_LDH,bm.BP) and bm.CHR.eq(18).all()
    stride=(len(fam)+3)//4
    bed=z.open('g1000_eur.bed')
    assert bed.read(3)==bytes([108,27,1])
    offset=3; packed=[]
    for idx in indices:
        target=3+idx*stride
        skip=target-offset
        while skip:
            consume=min(skip,8*1024**2)
            v=bed.read(consume);assert len(v)==consume
            skip-=consume;offset+=consume
        v=bed.read(stride);assert len(v)==stride
        packed.append(np.frombuffer(v,dtype=np.uint8));offset+=stride
    bed.close()
enc=np.array(packed)
codes=np.stack([(enc>>k)&3 for k in [0,2,4,6]],axis=2).reshape(len(bm),-1)[:,:len(fam)]
g=np.array([2.,np.nan,1.,0.])[codes]
missing=np.isnan(g)
sign=[];align=[];comp=str.maketrans('ACGT','TGCA')
for ar,br in zip(a.itertuples(),bm.itertuples()):
    assert {ar.A1_LDH,ar.A2_LDH} not in [{'A','T'},{'C','G'}], 'Historical harmonization unexpectedly retained ambiguous SNP'
    if (ar.A1_LDH,ar.A2_LDH)==(br.A1,br.A2):s,tag=1,'SAME'
    elif (ar.A1_LDH,ar.A2_LDH)==(br.A2,br.A1):s,tag=-1,'SWAPPED'
    elif (ar.A1_LDH,ar.A2_LDH)==(br.A1.translate(comp),br.A2.translate(comp)):s,tag=1,'STRAND'
    elif (ar.A1_LDH,ar.A2_LDH)==(br.A2.translate(comp),br.A1.translate(comp)):s,tag=-1,'STRAND_SWAPPED'
    else:raise ValueError('Allele mismatch: '+ar.SNP)
    sign.append(s);align.append(tag)
g=np.where(missing,np.nanmean(g,axis=1)[:,None],g)
g=np.where(np.array(sign)[:,None]==1,g,2-g)
assert np.isfinite(g).all()
nonmono=np.std(g,axis=1)>0
ex=[]
for snp in sorted(ids-set(bm.SNP)):ex.append({'SNP':snp,'reason':'ABSENT_FROM_EXISTING_503_EUR_REFERENCE_BY_SNP_ID'})
for snp in a.SNP[~nonmono]:ex.append({'SNP':snp,'reason':'MONOMORPHIC_EXISTING_REFERENCE'})
pd.DataFrame(ex,columns=['SNP','reason']).to_csv(O/'CHR18_REFERENCE_EXCLUSIONS.tsv',sep='\t',index=False)
a=a.loc[nonmono].reset_index(drop=True);g=g[nonmono];bm=bm.loc[nonmono].reset_index(drop=True)
out=pd.DataFrame({'SNP':a.SNP,'CHR':18,'BP':a.BP_LDH,'A1':a.A1_LDH,'A2':a.A2_LDH,'beta_LDH':a.BETA_LDH,'SE_LDH':a.SE_LDH,'P_LDH':a.P_LDH,'N_LDH':a.N_LDH,'MAF_LDH':np.minimum(a.EAF_LDH,1-a.EAF_LDH),'beta_CPI':a.BETA_CPI,'SE_CPI':a.SE_CPI,'P_CPI':a.P_CPI,'N_CPI':a.N_CPI,'MAF_CPI':np.minimum(a.EAF_CPI,1-a.EAF_CPI),'reference_EAF_effect_allele':g.mean(axis=1)/2,'reference_alignment':np.array(align)[nonmono]})
assert len(out)>100 and out.SNP.is_unique
out.to_csv(O/'CHR18_TWO_TRAIT_COMMON_FIXED_WINDOW.tsv',sep='\t',index=False)
with gzip.open(O/'CHR18_effect_allele_genotypes.tsv.gz','wt') as h:
    np.savetxt(h,g.T,fmt='%.9g',delimiter='\t',header='\t'.join(out.SNP),comments='')
ld=np.corrcoef(g)
assert np.max(np.abs(ld-ld.T))<1e-10 and np.max(np.abs(np.diag(ld)-1))<1e-10
np.save(O/'CHR18_signed_r.npy',ld)
q={'n_input_fixed_region':len(d),'n_reference_matched':len(bm)+int((~nonmono).sum()),'n_missing_reference':len(ids-set(seen)),'n_monomorphic_removed':int((~nonmono).sum()),'n_final':len(out),'n_reference_samples':503,'N':{t:{'min':float(out['N_'+t].min()),'median':float(out['N_'+t].median()),'max':float(out['N_'+t].max()),'unique':int(out['N_'+t].nunique())} for t in ['LDH','CPI']},'LD_symmetry_max_abs':float(np.max(np.abs(ld-ld.T))),'LD_diagonal_max_abs_error':float(np.max(np.abs(np.diag(ld)-1))),'dense_LD_bytes':int(ld.nbytes),'missing_reference_calls_mean_imputed':int(missing.sum()),'source_SNP_selection':'all 1357 frozen-region aligned SNPs; no P-value pruning','CPI_beta_encoding':'already signed to LDH alleles in original revision_locus_coloc.R; original CPI label columns retained; no second beta flip','source_LDH_reference_EAF_abs_diff_max':float(np.max(np.abs(out.reference_EAF_effect_allele-a.EAF_LDH))),'reference_scope':'existing official MAGMA 1000G Phase3 EUR503 archive; separate from UKB LAVA reference'}
(O/'CHR18_PREPARATION_QC.json').write_text(json.dumps(q,indent=2)+'\n')
cfg=f'''# Fixed before any RSS diagnostic or SuSiE model fit
fixed_locus_id: LDH_CPI_p05_chr18_001
build: GRCh37
chromosome: 18
start: 50358109
end: 51055069
boundary_sensitivity: none
LDH_N: 388514
CPI_N_values: [{q['N']['CPI']['min']}, {q['N']['CPI']['median']}, {q['N']['CPI']['max']}]
CPI_N_scope: scalar min median max sensitivity only
signed_LD_source: existing MAGMA 1000G Phase3 European 503 reference
archive_sha256: {ARC_SHA}
RSS_s_max: 0.1
RSS_s_threshold_origin: pre-existing project operational QC not an official universal threshold
SuSiE_L: 10
credible_set_coverage: 0.95
min_abs_corr: 0.5
max_iterations: 200
estimate_residual_variance: false
repeat_until_convergence: false
priors_p1_p2: [0.0001, 0.0001]
primary_p12: 0.00001
sensitivity_p12: [0.000001, 0.00005]
coloc_version: 5.2.2
susieR_version: 0.12.35
threads: 4
memory_limit_gib: 32
new_QTL_panels: 0
'''
(O/'CHR18_RUN_CONFIG.yaml').write_text(cfg)
lib=O/'R_library';lib.mkdir(exist_ok=True)
oldlib=R/'16_workflow/envs/coloc_exact_zero/library'
for package in ['coloc','susieR']:
    if not (lib/package).exists():shutil.copytree(oldlib/package,lib/package)
outs=[O/'CHR18_TWO_TRAIT_COMMON_FIXED_WINDOW.tsv',O/'CHR18_effect_allele_genotypes.tsv.gz',O/'CHR18_signed_r.npy',O/'CHR18_REFERENCE_EXCLUSIONS.tsv',O/'CHR18_PREPARATION_QC.json',O/'CHR18_RUN_CONFIG.yaml']
done.write_text(json.dumps({'signature':sig,'qc':q,'outputs':{str(p):sha(p) for p in outs}},indent=2)+'\n')
print(json.dumps(q,indent=2),flush=True)
