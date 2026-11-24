"""Correct symbol normalization in the newly prepared NeuN background only.

Preserve every historical output and the initial new H5AD. No scDRS score from
this atlas existed before the correction.
"""
from pathlib import Path
import hashlib,json,os
import anndata as ad,numpy as np,pandas as pd
from scipy import sparse
P=Path('/analysis/projects/JTM_recent_highIF_direct_clone/P2_LDH_shared_architecture')
N=P/'11_NatCommun_strict_clone_Pro/19_computational_finalize_JTM_20261103';S=N/'06_scDRS';D=N/'05_DRG'
mp=pd.read_csv(S/'EXISTING_H5_GENE_ID_SYMBOL_MAP.tsv',sep='\t');symbol_ids=mp.groupby('gene').ensembl.agg(list).to_dict()
mapping={1:'donor_1',2:'donor_2',3:'donor_3',4:'donor_4',5:'donor_3',6:'donor_5'}
original=S/'GSE168243.neun.raw_counts.h5ad';backup=S/'GSE168243.neun.raw_counts.pre_symbol_correction.h5ad'
if (S/'NEUN_SYMBOL_CORRECTION_QC.json').exists():raise SystemExit('Correction already completed; refusing rerun')
assert not list((S/'GSE168243').glob('*/**/*.full_score.gz')),'Cannot replace input after scoring'
old=ad.read_h5ad(backup if backup.exists() else original);mats=[];ob=[];names=None;records=[]
for prep in range(1,7):
 p=next((P/'05_results_intermediate/mda/drg/GSE168243_RAW').glob('*prep'+str(prep)+'.csv.gz'))
 d=pd.read_csv(p,index_col=0);d.index=d.index.astype(str);n_before=len(d);raw_library=d.sum(axis=0).to_numpy();d=d.groupby(level=0).sum()
 if names is None:names=d.index.to_numpy()
 else:assert np.array_equal(names,d.index)
 old_ix=old.obs.preparation.eq('prep'+str(prep)).to_numpy();oldlib=np.asarray(old.X[old_ix].sum(axis=1)).ravel();assert np.array_equal(oldlib,raw_library)
 assert np.array_equal(d.sum(axis=0).to_numpy(),raw_library)
 x=sparse.csr_matrix(d.values.T,dtype='float32');mats.append(x)
 ob.append(pd.DataFrame({'donor':mapping[prep],'preparation':f'prep{prep}','cell_class':'NeuN-enriched neuronal nuclei','GSM':f'GSM{5134536+prep}'},index=[f'prep{prep}_'+str(c) for c in d.columns]))
 records.append({'preparation':f'prep{prep}','n_raw_symbols':n_before,'n_unique_exact_symbols':len(d),'n_valid_dotted_symbols':int(d.index.str.contains('.',regex=False).sum()),'old_vs_new_library_totals_identical':True})
ids=[symbol_ids[g][0] if len(symbol_ids.get(g,[]))==1 else 'UNMAPPED_SYMBOL:'+g for g in names];assert len(ids)==len(set(ids))
a=ad.AnnData(X=sparse.vstack(mats),obs=pd.concat(ob),var=pd.DataFrame({'gene_symbol':names},index=ids));assert a.obs_names.equals(old.obs_names)
if not backup.exists():
 temp=S/'GSE168243.neun.raw_counts.corrected.partial.h5ad';a.write_h5ad(temp,compression='gzip');os.rename(original,backup);os.replace(temp,original)
targets=pd.read_csv(D/'DRG_CANDIDATE_SCOPE.tsv',sep='\t');exprpath=D/'DRG_CANDIDATE_DONOR_EXPRESSION.tsv';expr=pd.read_csv(exprpath,sep='\t');neu=expr[expr.dataset.eq('GSE168243')]
lib=np.asarray(a.X.sum(axis=1)).ravel();genes=a.var.gene_symbol.to_numpy();checks=[];revised=[];differences=[]
for _,r in neu.iterrows():
 ix=np.flatnonzero(a.obs.donor.eq(r.donor).to_numpy());jx=np.flatnonzero(genes==r.gene)
 new=r.to_dict()
 if len(jx)==0:
  new.update(n_expressing=np.nan,percent_expressing=np.nan,mean_log1p_CP10K=np.nan,total_UMI=np.nan,status='GENE_NOT_IN_DEPOSITED_MATRIX');revised.append(new);continue
 v=np.asarray(a.X[ix][:,jx].sum(axis=1)).ravel().astype(float)
 det=int((v>0).sum());mean=float(np.log1p(v/lib[ix]*1e4).mean());n=len(ix)
 new.update(n_cells=n,n_expressing=det,percent_expressing=100*det/n,mean_log1p_CP10K=mean,total_UMI=float(v.sum()),status='MEASURED' if n>=20 else 'MEASURED_LOW_CELL_COVERAGE',primary_coverage=n>=20)
 revised.append(new)
 if not np.isclose(det,r.n_expressing,equal_nan=True) or not np.isclose(mean,r.mean_log1p_CP10K,rtol=0,atol=1e-10,equal_nan=True):
  differences.append({'donor':r.donor,'gene':r.gene,'old_n_expressing':r.n_expressing,'corrected_n_expressing':det,'old_mean_log1p_CP10K':r.mean_log1p_CP10K,'corrected_mean_log1p_CP10K':mean,'reason':'Preserve exact deposited gene symbols; no unverified gene-symbol suffix merging'})
 checks.append(r.gene)
prev=D/'DRG_CANDIDATE_DONOR_EXPRESSION.pre_symbol_correction.tsv'
if not prev.exists():expr.to_csv(prev,sep='\t',index=False,na_rep='NA')
expr=pd.concat([expr[~expr.dataset.eq('GSE168243')],pd.DataFrame(revised)],ignore_index=True);expr.to_csv(exprpath,sep='\t',index=False,na_rep='NA')
pd.DataFrame(differences).to_csv(D/'NEUN_CANDIDATE_SYMBOL_CORRECTIONS.tsv',sep='\t',index=False,na_rep='NA')
summary=[]
for (ds,cl,g),z in expr.groupby(['dataset','cell_class','gene']):
 good=z[z.primary_coverage & z.status.eq('MEASURED')]
 summary.append({'dataset':ds,'cell_class':cl,'gene':g,'n_donors_total':len(z),'n_donors_eligible':len(good),'median_percent_expressing':good.percent_expressing.median(),'min_percent_expressing':good.percent_expressing.min(),'max_percent_expressing':good.percent_expressing.max(),'median_mean_log1p_CP10K':good.mean_log1p_CP10K.median(),'nuclei_in_eligible_donors':good.n_cells.sum(),'scope':z.scope.iloc[0]})
pd.DataFrame(summary).to_csv(D/'DRG_CANDIDATE_CLASS_SUMMARY.tsv',sep='\t',index=False,na_rep='NA')
result={'status':'CORRECTED_BEFORE_SCORING','cause':'Legacy code stripped dot suffixes from gene symbols, merging valid distinct symbols or unverified duplicate-name aliases','initial_background_genes':old.n_vars,'corrected_background_genes':a.n_vars,'n_nuclei':a.n_obs,'candidate_donor_rows_rechecked':len(checks),'candidate_donor_rows_changed':len(differences),'genes_with_changed_candidate_expression':sorted({x['gene'] for x in differences}),'nucleus_library_totals_changed':False,'historical_results_modified':False,'original_new_h5ad_preserved_at':str(backup),'replicate_data':records,'scDRS_input':str(original),'scope':'All original gene symbols preserved; no phenotype or score-dependent selection'}
(S/'NEUN_SYMBOL_CORRECTION_QC.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
