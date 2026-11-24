"""All-test multiplicity summaries; positional convergence is not causal mapping."""
import json
from pathlib import Path
import numpy as np,pandas as pd
from replan_run_task import digest
R=Path(__file__).resolve().parents[1];D=R/'06_fuma_magma_otg/magma_original_gwas'
def bh(p):
 p=np.asarray(p,float);o=np.argsort(p);v=np.empty(len(p));v[o]=np.minimum(1,np.minimum.accumulate((p[o]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]);return v
g=pd.read_csv(R/'06_fuma_magma_otg/magma_method_audit/official_resources/MAGMAgenev92.resource',sep=r'\s+',names=['GENE','chr','start','end','strand','symbol'])
allg=[];allt=[];audit=[]
for trait in ['LDH','CPI','SCIATICA']:
 if not (D/(trait+'.COMPLETE.json')).exists():continue
 for kind,suf in [('gene','.genes.out'),('GTEx54','_GTEx54.gsa.out'),('GTEx30','_GTEx30.gsa.out'),('MSigDB7','_MSigDB7.gsa.out')]:
  p=D/(trait+suf);done=json.loads((D/(trait+'_'+kind+'.done.json')).read_text());assert done['actual_exitcode']==0
  assert done['outputs'][str(p)]==digest(p)
  z=pd.read_csv(p,sep=r'\s+',comment='#');assert np.isfinite(z.P).all() and z.P.between(0,1).all()
  z['trait']=trait;z['family']=kind;z['BH_FDR_within_trait_family']=bh(z.P);z['Bonferroni_P_within_trait_family']=np.minimum(1,z.P*len(z));z['scope']='UNIVARIATE_ORIGINAL_GWAS_NOT_SHARED_PROOF'
  if kind=='gene':z=z.merge(g[['GENE','symbol']],on='GENE',how='left',validate='one_to_one');allg.append(z)
  else:allt.append(z)
  audit.append(dict(trait=trait,family=kind,n_tests=len(z),n_BH05=int((z.BH_FDR_within_trait_family<.05).sum()),n_Bonf05=int((z.Bonferroni_P_within_trait_family<.05).sum()),input=str(p),sha256=digest(p)))
if allg:pd.concat(allg,ignore_index=True).to_csv(D/'MAGMA_ALL_GENES.tsv',sep='\t',index=False)
if allt:pd.concat(allt,ignore_index=True).to_csv(D/'MAGMA_ALL_GENESETS_TISSUES.tsv',sep='\t',index=False)
pd.DataFrame(audit).to_csv(D/'MAGMA_MULTIPLICITY_SUMMARY.tsv',sep='\t',index=False)
if allg:
 gg=pd.concat(allg);conv=[]
 for pair,trait in [('LDH_CPI_low_overlap','CPI'),('LDH_SCIATICA_PANUKB_low_overlap','SCIATICA')]:
  p=R/'06_fuma_magma_otg/new_pair_interval_annotation'/pair/'LOCUS_POSITIONAL_AND_NEAREST_GENES.tsv'
  if not p.exists():continue
  loc=pd.read_csv(p,sep='\t');loc=loc[(loc.threshold==.05)&(loc.gene_body_overlaps_locus==True)]
  for t in ['LDH',trait]:
   zz=gg[gg.trait==t][['GENE','P','BH_FDR_within_trait_family','Bonferroni_P_within_trait_family']].rename(columns=lambda x:t+'_'+x if x!='GENE' else 'ensembl')
   loc=loc.merge(zz,on='ensembl',how='left',validate='many_to_one')
  loc['both_original_GWAS_gene_BH05']=(loc.LDH_BH_FDR_within_trait_family<.05)&(loc[trait+'_BH_FDR_within_trait_family']<.05)
  loc['interpretation']='Shared locus positional overlap plus original GWAS MAGMA; correlated evidence, not causal gene or independent replication'
  conv.append(loc)
 if conv:pd.concat(conv,ignore_index=True).to_csv(D/'MAGMA_SHARED_LOCUS_CONVERGENCE.tsv',sep='\t',index=False)
print(pd.DataFrame(audit).to_string(index=False))
