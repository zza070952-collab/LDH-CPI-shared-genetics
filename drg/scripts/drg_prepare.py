"""Prepare adult author-labelled DRG references and donor-level candidate expression.

Only existing P2 raw counts/metadata and frozen positional candidate mapping are read.
No clustering or disease-state inference is performed.
"""
import csv, gzip, hashlib, io, json, os, re, tarfile
from pathlib import Path
import h5py, anndata as ad, numpy as np, pandas as pd
from scipy import sparse

P=Path('/analysis/projects/JTM_recent_highIF_direct_clone/P2_LDH_shared_architecture')
R=P/'11_NatCommun_strict_clone_Pro'
N=R/'19_computational_finalize_JTM_20261103'
D=N/'05_DRG'; S=N/'06_scDRS'; D.mkdir(exist_ok=True,parents=True); S.mkdir(exist_ok=True,parents=True)

def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
def tsv(df,p):df.to_csv(p,sep='\t',index=False,na_rep='NA')
def dec(a):return np.array([x.decode() if isinstance(x,bytes) else str(x) for x in a])
def soft(g):
 p=P/f'01_data_raw/human_drg_final/{g}/{g}_family.full.soft.gz'; result=[]; row=None; design=''
 with gzip.open(p,'rt') as h:
  for line in h:
   line=line.rstrip()
   if line.startswith('!Series_overall_design = '):design=line.split(' = ',1)[1]
   if line.startswith('^SAMPLE = '):
    if row:result.append(row)
    row={'dataset':g,'GSM':line.split(' = ',1)[1]}
   elif row and line.startswith('!Sample_title = '):row['title']=line.split(' = ',1)[1]
   elif row and line.startswith('!Sample_characteristics_ch1 = '):
    x=line.split(' = ',1)[1]
    if ': ' in x:
     k,v=x.split(': ',1);row[k.lower()]=v
  if row:result.append(row)
 return result, design, p

pospath=R/'06_fuma_magma_otg/new_pair_interval_annotation/LDH_CPI_low_overlap/LOCUS_POSITIONAL_AND_NEAREST_GENES.tsv'
pos=pd.read_csv(pospath,sep='\t');pos=pos[pos.threshold.eq(.05)]
cmatpath=R/'17_evidence_tier_revision_20260912/04_gene_evidence/CANDIDATE_CLAIM_EVIDENCE_MATRIX.tsv'
cmat=pd.read_csv(cmatpath,sep='\t')
targets=sorted(set(pos.gene)|set(cmat.gene)); primary=set(pos.gene)
targetdf=pd.DataFrame({'gene':targets,'scope':['PRIMARY_FIXED_18' if x in primary else 'HISTORICAL_CONTEXT_ONLY' for x in targets]})
tsv(targetdf,D/'DRG_CANDIDATE_SCOPE.tsv')
ex=pos[pos.locus_id.eq('LDH_CPI_p05_chr18_001')][['gene','ensembl','locus_id']].drop_duplicates();ex['source_file']=str(pospath)
l2gpath=R/'06_fuma_magma_otg/otg_new_loci_zero/L2G_NEW_INTERVAL_MATCHES.tsv'
l2g=pd.read_csv(l2gpath,sep='\t');lx=l2g[l2g.locus_id.eq('LDH_CPI_p05_chr18_001')][['gene','ensembl','locus_id']].drop_duplicates();lx['source_file']=str(l2gpath)
ex=pd.concat([ex,lx],ignore_index=True).drop_duplicates(['gene','ensembl','locus_id']);assert 'DCC' in set(ex.gene)
tsv(ex,S/'CHR18_EXCLUDED_GENES_FROZEN_BEFORE_SCORE.tsv')
provenance=[]; metadata=[]; all_expr=[]; datasetqc=[]

# Dataset 1: retain exact existing author-metadata/raw-count intersection and adults.
raw=P/'01_data_raw/human_drg_final/GSE201586/GSE201586_RAW.tar'
cache=R/'00_admin/zero_process_evidence_20260909/DRG_DUPLICATE_H5_CACHE.json'
cachedata=json.loads(cache.read_text());assert sha(raw)==cachedata['archive_sha256']
meta_path=P/'06_results_tables/P2_GSE201586_CELL_METADATA.tsv.gz'
meta=pd.read_csv(meta_path,sep='\t');assert meta.cell_id.is_unique
adults_path=R/'10_drg_spinal_validation/adult_metadata_20260909/GSE201586_ADULT_DONOR_ELIGIBILITY.tsv'
adult_table=pd.read_csv(adults_path,sep='\t');adultset=set(adult_table.loc[adult_table.adult_status.eq('INCLUDE_ADULT'),'donor'])
assert len(adultset)==5
softrows,design,softpath=soft('GSE201586');soft_by={x['GSM']:x for x in softrows}
for x in softrows:
 donor=x['donor #']; expected=adult_table.set_index('donor').loc[donor,'age_reported'];assert expected==x['age']
 adultstatus=adult_table.set_index('donor').loc[donor,'adult_status']
 metadata.append({'dataset':'GSE201586','GSM':x['GSM'],'preparation':x['title'],'donor':donor,'age':x['age'],'sex':x['sex'],'tissue':'adult human DRG' if donor in adultset else 'human DRG','platform':'10x snRNA-seq','selection':'not NeuN selected; FACS or density-gradient nuclei isolation','adult_primary_status':adultstatus,'donor_mapping_source':str(softpath),'cross_dataset_donor_overlap':'NOT_ESTABLISHED'})

counts_path=S/'GSE201586.adult_author.raw_counts.h5ad'
if counts_path.exists():
 a=ad.read_h5ad(counts_path)
else:
 mats=[];ob=[];q=[];ens=None;names=None
 with tarfile.open(raw) as tf:
  for rec in cachedata['files']:
   member=rec['archive_member'];sam=Path(member).name.split('_')[1]
   sub=meta[meta.donor.isin(adultset)&meta['orig.ident'].eq(sam)].set_index('cell_id')
   if sub.empty:continue
   b=tf.extractfile(member).read();assert hashlib.sha256(b).hexdigest()==rec['sha256']
   with h5py.File(io.BytesIO(b),'r') as h:
    m=h['matrix'];en=dec(m['features']['id'][:]);nm=dec(m['features']['name'][:]);bar=dec(m['barcodes'][:]);x=sparse.csc_matrix((m['data'][:],m['indices'][:],m['indptr'][:]),shape=tuple(m['shape'][:]))
   if ens is None:ens=en;names=nm
   else:assert np.array_equal(ens,en) and np.array_equal(names,nm)
   ids=np.array([v+'_'+sam for v in bar]);keep=np.isin(ids,sub.index);c=ids[keep];sub=sub.loc[c].copy();x=x[:,keep].T.tocsr();assert x.shape[0]==len(sub)
   assert np.isfinite(x.data).all() and (x.data>=0).all() and np.all(x.data==np.floor(x.data))
   ob.append(pd.DataFrame({'donor':sub.donor.astype(str),'preparation':sam,'cell_class':sub.celltype.astype(str),'GSM':sub.GEO_sample.astype(str)},index=c));mats.append(x)
   q.append({'sample':sam,'donor':sub.donor.iloc[0],'author_cells':len(meta[meta['orig.ident'].eq(sam)&meta.donor.isin(adultset)]),'matched_count_cells':len(c),'excluded_author_cells_without_raw_count':len(meta[meta['orig.ident'].eq(sam)&meta.donor.isin(adultset)])-len(c),'raw_member_sha256':rec['sha256']});print('raw matched',sam,len(c),flush=True)
 a=ad.AnnData(X=sparse.vstack(mats).astype('float32'),obs=pd.concat(ob),var=pd.DataFrame({'gene_symbol':names},index=ens))
 assert a.obs_names.is_unique and a.var_names.is_unique
 a.write_h5ad(counts_path,compression='gzip');tsv(pd.DataFrame(q),D/'GSE201586_RAW_METADATA_INTERSECTION.tsv')

gene_map=pd.DataFrame({'ensembl':a.var_names,'gene':a.var.gene_symbol.values});tsv(gene_map,S/'EXISTING_H5_GENE_ID_SYMBOL_MAP.tsv')
symbol_ids=gene_map.groupby('gene').ensembl.agg(list).to_dict()

def expression(a,dataset):
 lib=np.asarray(a.X.sum(axis=1)).ravel();assert (lib>0).all()
 genes=a.var.gene_symbol.astype(str).to_numpy();target_ix={g:np.flatnonzero(genes==g) for g in targets}
 rows=[];groupqc=[]
 for (donor,cl),ix in a.obs.groupby(['donor','cell_class'],observed=True).indices.items():
  x=a.X[ix].tocsc();n=len(ix);groupqc.append({'dataset':dataset,'donor':str(donor),'cell_class':str(cl),'n_cells':n,'primary_coverage':n>=20})
  for g,jx in target_ix.items():
   if len(jx)==0:
    rows.append({'dataset':dataset,'donor':str(donor),'cell_class':str(cl),'gene':g,'n_cells':n,'n_expressing':np.nan,'percent_expressing':np.nan,'mean_log1p_CP10K':np.nan,'total_UMI':np.nan,'status':'GENE_NOT_IN_DEPOSITED_MATRIX','primary_coverage':n>=20});continue
   v=np.asarray(x[:,jx].sum(axis=1)).ravel().astype(float);det=int((v>0).sum());mean=float(np.log1p(v/lib[ix]*1e4).mean())
   rows.append({'dataset':dataset,'donor':str(donor),'cell_class':str(cl),'gene':g,'n_cells':n,'n_expressing':det,'percent_expressing':100*det/n,'mean_log1p_CP10K':mean,'total_UMI':v.sum(),'status':'MEASURED' if n>=20 else 'MEASURED_LOW_CELL_COVERAGE','primary_coverage':n>=20})
 return rows,groupqc

e,q=expression(a,'GSE201586');all_expr+=e;datasetqc+=q
tsv(a.obs.reset_index(names='cell_id'),D/'GSE201586_AUTHOR_CELL_METADATA_USED.tsv')
provenance += [{'file':str(p),'sha256':sha(p)} for p in [raw,meta_path,adults_path,softpath,cache]]
del a

# Dataset 2: author-deposited neuronal nuclei; no inferred neuronal subtype labels.
softrows,design2,softpath=soft('GSE168243');assert 'hDRG3 and hDRG5 are from the same donor' in design2
mapping={1:'donor_1',2:'donor_2',3:'donor_3',4:'donor_4',5:'donor_3',6:'donor_5'}
counts_path=S/'GSE168243.neun.raw_counts.h5ad'
if counts_path.exists():a=ad.read_h5ad(counts_path)
else:
 mats=[];obs=[];original_names=None
 for prep in range(1,7):
  files=list((P/'05_results_intermediate/mda/drg/GSE168243_RAW').glob('*prep'+str(prep)+'.csv.gz'));assert len(files)==1;p=files[0]
  # These are gene symbols, not versioned ENSG IDs. Dots are valid parts of
  # symbols such as RP11-34P13.7 and must be preserved. Only identical labels
  # may be collapsed by summing counts.
  d=pd.read_csv(p,index_col=0);d.index=d.index.astype(str);d=d.groupby(level=0).sum()
  assert np.isfinite(d.values).all() and (d.values>=0).all() and np.all(d.values==np.floor(d.values))
  if original_names is None:original_names=d.index.to_numpy()
  else:assert np.array_equal(original_names,d.index)
  x=sparse.csr_matrix(d.values.T,dtype='float32');ids=[f'prep{prep}_'+str(c) for c in d.columns]
  obs.append(pd.DataFrame({'donor':mapping[prep],'preparation':f'prep{prep}','cell_class':'NeuN-enriched neuronal nuclei','GSM':f'GSM{5134536+prep}'},index=ids));mats.append(x)
  provenance.append({'file':str(p),'sha256':sha(p)})
 var_ids=[]
 for g in original_names:
  hits=symbol_ids.get(g,[]);var_ids.append(hits[0] if len(hits)==1 else 'UNMAPPED_SYMBOL:'+g)
 assert len(var_ids)==len(set(var_ids))
 a=ad.AnnData(X=sparse.vstack(mats),obs=pd.concat(obs),var=pd.DataFrame({'gene_symbol':original_names},index=var_ids))
 a.write_h5ad(counts_path,compression='gzip')
for prep,x in enumerate(softrows,1):metadata.append({'dataset':'GSE168243','GSM':x['GSM'],'preparation':f'prep{prep}','donor':mapping[prep],'age':x['age'],'sex':x['gender'],'tissue':'adult human DRG','platform':'10x snRNA-seq v2' if prep<=3 else '10x snRNA-seq v3','selection':'NeuN-enriched neuronal nuclei','adult_primary_status':'INCLUDE_ADULT','donor_mapping_source':str(softpath)+'; Series_overall_design: hDRG3 and hDRG5 same donor','cross_dataset_donor_overlap':'NOT_ESTABLISHED'})
e,q=expression(a,'GSE168243');all_expr+=e;datasetqc+=q
tsv(a.obs.reset_index(names='cell_id'),D/'GSE168243_CELL_METADATA_USED.tsv');provenance.append({'file':str(softpath),'sha256':sha(softpath)})

expr=pd.DataFrame(all_expr).merge(targetdf,on='gene',how='left');tsv(expr,D/'DRG_CANDIDATE_DONOR_EXPRESSION.tsv')
groups=pd.DataFrame(datasetqc);tsv(groups,D/'DRG_DONOR_CLASS_CELL_COUNTS.tsv')
metadata=pd.DataFrame(metadata)
cellcounts=pd.concat([pd.read_csv(D/'GSE201586_AUTHOR_CELL_METADATA_USED.tsv',sep='\t'),pd.read_csv(D/'GSE168243_CELL_METADATA_USED.tsv',sep='\t')])
sizes=cellcounts.groupby('GSM').size();metadata['n_cells_used']=metadata.GSM.map(sizes).fillna(0).astype(int)
tsv(metadata,D/'DRG_DONOR_METADATA.tsv')
summary=[]
for (ds,cl,g),z in expr.groupby(['dataset','cell_class','gene']):
 good=z[z.primary_coverage & z.status.eq('MEASURED')]
 summary.append({'dataset':ds,'cell_class':cl,'gene':g,'n_donors_total':len(z),'n_donors_eligible':len(good),'median_percent_expressing':good.percent_expressing.median(),'min_percent_expressing':good.percent_expressing.min(),'max_percent_expressing':good.percent_expressing.max(),'median_mean_log1p_CP10K':good.mean_log1p_CP10K.median(),'nuclei_in_eligible_donors':good.n_cells.sum(),'scope':'PRIMARY_FIXED_18' if g in primary else 'HISTORICAL_CONTEXT_ONLY'})
tsv(pd.DataFrame(summary),D/'DRG_CANDIDATE_CLASS_SUMMARY.tsv')
for p in [pospath,cmatpath,l2gpath]:provenance.append({'file':str(p),'sha256':sha(p)})
tsv(pd.DataFrame(provenance),D/'DRG_INPUT_PROVENANCE.tsv')
qc={'status':'COMPLETE','n_primary_candidates':len(primary),'n_total_candidates':len(targets),'expression_rows':len(expr),'metadata_rows':len(metadata),'datasets':{ds:{'n_donors':int(g.donor.nunique()),'n_nuclei':int(g.n_cells.sum()),'n_donor_class_groups':len(g),'n_eligible_donor_class_groups':int(g.primary_coverage.sum())} for ds,g in groups.groupby('dataset')},'GSE168243_design_verified':design2,'cross_dataset_donor_independence':'NOT_ESTABLISHED','gene_symbol_ambiguity_rule':'GSE168243 symbols map to existing GSE201586 H5 ENSG only if uniquely mapped; others retained in full background with UNMAPPED_SYMBOL prefix','candidate_expression_scale':'log1p(CP10K), donor-class means; no cell-level disease-effect tests','minimum_donor_class_nuclei':20}
(D/'DRG_MODULE_STATUS.json').write_text(json.dumps(qc,indent=2));print(json.dumps(qc,indent=2),flush=True)
