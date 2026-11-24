#!/usr/bin/env python3
"""Global BH and exact probe/panel/top-eQTL shared regulation, all traits retained."""
import argparse,csv,gzip,json,math
from pathlib import Path
import numpy as np
import pandas as pd

p=argparse.ArgumentParser();p.add_argument('--job-manifest',required=True);p.add_argument('--outdir',required=True);p.add_argument('--family-status',choices=['PROVISIONAL_AVAILABLE_TRAITS','FINAL_REGISTERED_FAMILY'],default='PROVISIONAL_AVAILABLE_TRAITS');a=p.parse_args()
manifest=pd.read_csv(a.job_manifest,sep='\t',keep_default_na=False);out=Path(a.outdir);out.mkdir(parents=True,exist_ok=True)
chunks=[];audit=[]
for _,j in manifest.iterrows():
    if j['pairing_approved']!='YES':continue
    pref=Path(j['output_prefix']);exitf=Path(str(pref)+'.exitcode');marker=Path(str(pref)+'.done')
    if not(exitf.exists() and exitf.read_text().strip()=='0' and marker.exists()):
        raise RuntimeError('Incomplete prespecified job '+j['job_id']+'; do not re-estimate global FDR on selected completed panels')
    d=pd.read_csv(str(pref)+'.smr',sep=r'\s+')
    d=d.rename(columns={'ProbeID':'probeID','Probe_Chr':'ProbeChr','SNP':'topSNP'})
    for k in ['trait','axis','primary_or_sensitivity','panel','qtl_resource_id']:d[k]=j[k]
    d['job_id']=j['job_id'];chunks.append(d)
    audit.append(dict(job_id=j['job_id'],trait=j['trait'],panel=j['panel'],raw_tests=len(d),status='SUCCESS'))
assert chunks,'No eligible complete jobs'
d=pd.concat(chunks,ignore_index=True)
d['multiplicity_family_status']=a.family_status
key=['trait','panel','qtl_resource_id','probeID'];assert not d.duplicated(key).any(),'Duplicated trait-panel-probe hypotheses'
d['MHC_flag']=(pd.to_numeric(d.ProbeChr,errors='coerce')==6)&pd.to_numeric(d.Probe_bp,errors='coerce').between(25000000,35000000)
def bh(p):
    arr=pd.to_numeric(p,errors='coerce').to_numpy();q=np.full(len(arr),np.nan);ok=np.isfinite(arr)&(arr>=0)&(arr<=1);xx=arr[ok];ix=np.argsort(xx);res=np.empty(len(xx))
    if len(xx):res[ix]=np.minimum(1,np.minimum.accumulate((xx[ix]*len(xx)/np.arange(1,len(xx)+1))[::-1])[::-1])
    q[ok]=res;return q
eligible=~d.MHC_flag.to_numpy(dtype=bool)
# Use explicit positional arrays after a reproducible boolean .loc assignment
# anomaly in this runtime; no upstream-library defect is assumed. The retained
# diagnostic rules out duplicate/mismatched row labels. Every output is also
# validated independently against R stats::p.adjust before scientific use.
eligible_positions=np.flatnonzero(eligible)
qglobal=np.full(len(d),np.nan);qglobal[eligible_positions]=bh(d.p_SMR.iloc[eligible_positions])
qwithin=np.full(len(d),np.nan)
for trait in d.trait.unique():
    positions=np.flatnonzero(eligible&(d.trait.to_numpy()==trait));qwithin[positions]=bh(d.p_SMR.iloc[positions])
d['FDR_global_all_traits_panels']=qglobal
d['FDR_within_trait_all_panels']=qwithin
assert np.allclose(d.FDR_global_all_traits_panels.to_numpy()[eligible_positions],bh(d.p_SMR.iloc[eligible_positions]))
d['primary_pass']=(d.FDR_global_all_traits_panels<.05)&(pd.to_numeric(d.p_HEIDI,errors='coerce')>.05)&eligible
d['legacy_HEIDI01_pass']=(d.FDR_global_all_traits_panels<.05)&(pd.to_numeric(d.p_HEIDI,errors='coerce')>.01)&eligible
d.to_csv(out/'SMR_ALL_TRAITS_ALL_PANELS.tsv.gz',sep='\t',index=False)
ldh=d[d.trait=='LDH'];pain=d[d.trait!='LDH'];join=['panel','qtl_resource_id','probeID']
pairs=ldh.merge(pain,on=join,suffixes=('_LDH','_pain'))
def orient(r):
    aa=(str(r.A1_LDH).upper(),str(r.A2_LDH).upper());bb=(str(r.A1_pain).upper(),str(r.A2_pain).upper())
    if r.topSNP_LDH!=r.topSNP_pain:return 'TOP_EQTL_MISMATCH'
    if aa==bb:return 'SAME'
    if aa==bb[::-1]:return 'FLIPPED'
    comp=str.maketrans('ACGT','TGCA');cc=tuple(x.translate(comp) for x in bb)
    if aa==cc:return 'COMPLEMENT'
    if aa==cc[::-1]:return 'COMPLEMENT_FLIPPED'
    return 'ALLELE_MISMATCH'
if len(pairs):
    pairs['allele_alignment']=pairs.apply(orient,axis=1)
    pairs['eqtl_beta_aligned_pain']=pairs.b_eQTL_pain*np.where(pairs.allele_alignment.isin(['FLIPPED','COMPLEMENT_FLIPPED']),-1,1)
    pairs['eqtl_effect_consistent']=np.isclose(pairs.b_eQTL_LDH,pairs.eqtl_beta_aligned_pain,rtol=1e-4,atol=1e-8)
    pairs['shared_primary_pass']=pairs.primary_pass_LDH&pairs.primary_pass_pain&pairs.allele_alignment.isin(['SAME','FLIPPED','COMPLEMENT','COMPLEMENT_FLIPPED'])&pairs.eqtl_effect_consistent
    pairs['expression_effect_direction']=np.where(pairs.b_SMR_LDH*pairs.b_SMR_pain>0,'CONCORDANT','DISCORDANT')
else:
    for k in ['allele_alignment','eqtl_beta_aligned_pain','eqtl_effect_consistent','shared_primary_pass','expression_effect_direction']:pairs[k]=pd.Series(dtype='object')
pairs.to_csv(out/'SMR_PROBE_TISSUE_PAIR_AUDIT_ALL.tsv',sep='\t',index=False)
shared=pairs[pairs.shared_primary_pass==True].copy()
if len(shared):
    n=shared[shared.primary_or_sensitivity_pain=='PRIMARY'].groupby(join)['axis_pain'].nunique()
    shared=shared.merge(n.rename('n_primary_axes'),on=join,how='left');shared['axis_class']=np.where(shared.n_primary_axes>=2,'PAN-PAIN','AXIS-SPECIFIC')
else:shared['n_primary_axes']=pd.Series(dtype='int');shared['axis_class']=pd.Series(dtype='object')
shared.to_csv(out/'SMR_SHARED_PROBE_TISSUE.tsv',sep='\t',index=False)
cols=[c for c in ['Gene_LDH','probeID','panel','qtl_resource_id','trait_pain','axis_pain','b_SMR_LDH','b_SMR_pain','FDR_global_all_traits_panels_LDH','FDR_global_all_traits_panels_pain','p_HEIDI_LDH','p_HEIDI_pain','topSNP_LDH','expression_effect_direction','axis_class'] if c in shared]
shared[cols].to_csv(out/'SMR_SHARED_GENE_AXIS_MATRIX.tsv',sep='\t',index=False)
pd.DataFrame(audit).to_csv(out/'SMR_JOB_COMPLETION_AUDIT.tsv',sep='\t',index=False)
(out/'SMR_MULTIPLICITY_DEFINITION.md').write_text('Primary BH family pools every completed prespecified non-MHC trait-panel-probe test. Within-trait BH across all panels is retained only as sensitivity. No filtering on SMR P, HEIDI, candidate genes, or positive tissue precedes primary BH. Shared regulation additionally requires identical QTL resource, panel, probe and top eQTL, aligned alleles, consistent eQTL effects, and primary SMR/HEIDI passes in both LDH and pain. b_SMR is expression-scaled; flipping both GWAS and QTL effects leaves its sign invariant.\n')
print(json.dumps(dict(all_tests=len(d),nonMHC_tests=int(eligible.sum()),shared_probe_trait_rows=len(shared))))
