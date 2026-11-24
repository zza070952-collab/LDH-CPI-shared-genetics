"""Integrate full published-model families. Keep unavailable tests and the candidate universe explicit."""
from pathlib import Path
import hashlib,json,sqlite3
import numpy as np
import pandas as pd
N=Path(__file__).resolve().parents[1]
O=N.parent/'20_JTM_literature_benchmarked_finalize_20261113'
D=N/'01_twas'; I=N/'03_integrated_evidence';I.mkdir(exist_ok=True)
def read(p):return pd.read_csv(p,sep='\t',low_memory=False)
def save(d,p):d.to_csv(p,sep='\t',na_rep='NA',index=False)
def bh(p):
    a=np.array(p); idx=np.argsort(a,kind='stable'); q=np.empty(len(a));q[idx]=np.minimum(1,np.minimum.accumulate((a[idx]*len(a)/np.arange(1,len(a)+1))[::-1])[::-1]);return q
def adjust(df,col,mask,group,out):
    df[out]=np.nan;df[out+'_N']=pd.Series(pd.NA,index=df.index,dtype='Int64')
    for _,sub in df.loc[mask & df[col].notna()].groupby(group):
        df.loc[sub.index,out]=bh(sub[col]);df.loc[sub.index,out+'_N']=len(sub)
lock=json.loads((N/'00_admin/ANALYSIS_PLAN_LOCK.json').read_text())
paths=sorted(D.glob('*.tests.tsv.gz'));assert len(paths)==8,paths
raw=pd.concat([read(p) for p in paths],ignore_index=True)
catalog=[];zero=[]
for tissue in lock['tissues']:
    for family in ['mashr','en']:
        filename=family+'_'+tissue+'.db'
        db=next(iter([p for root in [N/'model_cache',D/'models'] for p in root.rglob(filename)]))
        con=sqlite3.connect('file:'+db.as_posix()+'?mode=ro',uri=True)
        meta=pd.read_sql('select * from extra',con)
        nonzero=set(r[0] for r in con.execute('select distinct gene from weights where weight!=0'));con.close()
        for _,r in meta.iterrows():
            catalog.append([family,tissue,r.gene,r['n.snps.in.model'],r['pred.perf.R2'],r['pred.perf.pval'],r['pred.perf.qval']])
            if r.gene not in nonzero:
                for tr in ['LDH','CPI']:
                    for scope in ['trait','common']:
                        zero.append(dict(family=family,tissue=tissue,ensembl=r.gene,symbol=r.genename,trait=tr,predictor=scope,n_weight_snps=0,n_in_covariance=0,n_in_GWAS=0,n_allele_aligned=0,status='NO_NONZERO_WEIGHTS',subset_hash='EMPTY'))
raw=pd.concat([raw,pd.DataFrame(zero)],ignore_index=True).merge(pd.DataFrame(catalog,columns=['family','tissue','ensembl','db_model_snps','model_performance_R2','model_performance_P','model_performance_q']),on=['family','tissue','ensembl'],how='left',validate='many_to_one')
raw['stable_id']=raw.ensembl.str.split('.').str[0]
keys=['family','tissue','stable_id'];assert not raw.duplicated(keys+['trait','predictor']).any()
s=raw[raw.predictor=='trait'].copy(); eligible=s.status=='ELIGIBLE'
adjust(s,'P',eligible,['family','trait'],'q_trait_global4tissue')
adjust(s,'CPI_phi_P',eligible & s.family.eq('en') & s.trait.eq('CPI'),['family','trait'],'q_phi_global4tissue')
c=raw[raw.predictor=='common'].copy()
def wide(d,cols):
    a=d[d.trait=='LDH'][keys+cols].set_index(keys).add_suffix('_LDH')
    b=d[d.trait=='CPI'][keys+cols].set_index(keys).add_suffix('_CPI')
    assert a.index.equals(b.index)
    return a.join(b).reset_index()
j=wide(c,['P','Z','effect','status','subset_hash','n_allele_aligned','coverage_fraction'])
assert (j.subset_hash_LDH==j.subset_hash_CPI).all()
ok=(j.status_LDH=='ELIGIBLE')&(j.status_CPI=='ELIGIBLE')
j['P_conjunction']=j[['P_LDH','P_CPI']].max(axis=1).where(ok)
j['same_tissue_Z_direction']=np.where(~ok,'INELIGIBLE',np.where(j.Z_LDH*j.Z_CPI>0,'SAME_SIGN',np.where(j.Z_LDH*j.Z_CPI<0,'OPPOSITE_SIGN','ZERO')))
adjust(j,'P_conjunction',ok,['family'],'q_conjunction')
w=wide(s,['Z','P','effect','status','n_allele_aligned','coverage_fraction','q_trait_global4tissue','q_trait_global4tissue_N','CPI_phi_P','q_phi_global4tissue'])
w=w.rename(columns={col:'trait_'+col for col in w if col not in keys})
design=s[s.trait=='LDH'][keys+['ensembl','symbol','n_weight_snps','n_in_covariance','phi','model_performance_R2','model_performance_P','model_performance_q']]
j=j.merge(w,on=keys,validate='one_to_one').merge(design,on=keys,validate='one_to_one')
save(j,D/'TWAS_SAME_TISSUE_CONJUNCTION.tsv')
raw=raw.merge(s[keys+['trait','predictor','q_trait_global4tissue','q_trait_global4tissue_N','q_phi_global4tissue','q_phi_global4tissue_N']],on=keys+['trait','predictor'],how='left')
raw=raw.merge(j[keys+['q_conjunction','P_conjunction','q_conjunction_N']].assign(predictor='common'),on=keys+['predictor'],how='left')
save(raw,D/'TWAS_ALL_TESTS.tsv.gz')
save(raw[keys+'ensembl symbol trait predictor n_weight_snps n_in_covariance n_in_GWAS n_allele_aligned coverage_fraction var_g subset_hash status'.split()],D/'MODEL_COVERAGE_AND_HARMONIZATION.tsv')
frozen=read(O/'03_gene_and_tissue/GENE_EVIDENCE_MATRIX.tsv')
pos=read(O/'03_gene_and_tissue/POSITIONAL_MAPPING_ALL.tsv')
stable=pos[['gene','ensembl']].drop_duplicates()
symbolids=j[['symbol','stable_id']].drop_duplicates()
rows=[]
for row in frozen.itertuples():
    fixed=str(row.in_fixed_primary_18).lower()=='true'
    p=set(stable.loc[stable.gene==row.gene,'ensembl'].dropna())
    q=set(symbolids.loc[symbolids.symbol==row.gene,'stable_id'].dropna())
    ids=p if len(p)==1 else q
    assert len(ids)<=1,(row.gene,ids)
    sid=next(iter(ids),None)
    for fam in ['mashr','en']:
        for tissue in lock['tissues']:
            rows.append(dict(gene=row.gene,scope='CURRENT101' if fixed else 'HISTORICAL6',stable_id=sid,family=fam,tissue=tissue,ID_source='POSITIONAL_ENSEMBL' if len(p)==1 else 'MODEL_SYMBOL' if q else 'UNRESOLVED_SYMBOL',locus_id=row.locus_id))
focus=pd.DataFrame(rows).merge(j,on=keys,how='left',validate='many_to_one')
for a in ['LDH','CPI']:
    for k in ['status_','trait_status_']:focus[k+a]=focus[k+a].fillna('NO_MODEL')
assert focus.groupby('scope').gene.nunique().to_dict()=={'CURRENT101':101,'HISTORICAL6':6}
save(focus,D/'TWAS_CURRENT101_AND_HISTORICAL6.tsv')
cols=['ensembl','n_weight_snps','trait_status_LDH','trait_status_CPI','trait_Z_LDH','trait_Z_CPI','trait_q_trait_global4tissue_LDH','trait_q_trait_global4tissue_CPI','trait_CPI_phi_P_CPI','trait_q_phi_global4tissue_CPI','P_conjunction','q_conjunction','same_tissue_Z_direction']
sens=focus.pivot(index=['scope','gene','tissue'],columns='family',values=cols)
sens.columns=[b+'_'+a for a,b in sens.columns];save(sens.reset_index(),D/'TWAS_MODEL_SENSITIVITY.tsv')
out=[]
for scope,frame in [('ALL',j),('CURRENT101',focus[focus.scope=='CURRENT101']),('HISTORICAL6',focus[focus.scope=='HISTORICAL6'])]:
    for fam,group in frame.groupby('family'):
        gene='stable_id' if scope=='ALL' else 'gene'
        l=group.trait_q_trait_global4tissue_LDH<.05;cpi=group.trait_q_trait_global4tissue_CPI<.05
        both=set(group.loc[l,gene])&set(group.loc[cpi,gene]);same=set(group.loc[l&cpi,gene]);conj=group.q_conjunction<.05
        record=dict(scope=scope,family=fam,model_rows=int(group.ensembl.notna().sum()),model_genes=group.loc[group.ensembl.notna(),gene].nunique(),
          LDH_eligible_rows=(group.trait_status_LDH=='ELIGIBLE').sum(),CPI_eligible_rows=(group.trait_status_CPI=='ELIGIBLE').sum(),common_eligible_rows=group.q_conjunction.notna().sum(),
          LDH_sig_rows=l.sum(),LDH_sig_genes=group.loc[l,gene].nunique(),CPI_sig_rows=cpi.sum(),CPI_sig_genes=group.loc[cpi,gene].nunique(),
          any_tissue_dual_genes=len(both),dual_same_traitq_rows=(l&cpi).sum(),dual_same_traitq_genes=len(same),cross_tissue_only_genes=len(both-same),conjunction_sig_rows=conj.sum(),conjunction_sig_genes=group.loc[conj,gene].nunique(),
          conjunction_opposite_rows=(conj & group.same_tissue_Z_direction.eq('OPPOSITE_SIGN')).sum(),CPI_phi_tested=group.trait_q_phi_global4tissue_CPI.notna().sum(), CPI_phi_sig_rows=(group.trait_q_phi_global4tissue_CPI<.05).sum())
        out.append(record)
counts=pd.DataFrame(out);save(counts,D/'TWAS_FAMILY_COUNTS.tsv')
phi=pd.DataFrame([['mashr','both',None,None,'NO_MATCHING_PHI','Ordinary S-PrediXcan'],['en','LDH',None,None,'NOT_ESTIMABLE_UNDER_CURRENT_ASSUMPTIONS','Binary trait; N effective vs total and observed vs liability scale unresolved'],['en','CPI',134627,.1012,'SENSITIVITY_FULL_PREDICTORS_ONLY','Observed quantitative LDSC h2; maximum per-SNP N approximation; variance inflation 1+N*h2*phi; own BH family']],columns=['family','trait','N','h2','status','notes'])
save(phi,D/'TWAS_PHI_PARAMETER_AUDIT.tsv')
summary=[]
for gene,group in focus.groupby('gene'):
    d={'gene':gene}
    for fam,f in group.groupby('family'):
        d[fam+'_eligible_common_tissues']=int(f.q_conjunction.notna().sum())
        d[fam+'_conjunction_min_q']=f.q_conjunction.min()
        d[fam+'_conjunction_tissues']='|'.join(f.loc[f.q_conjunction<.05,'tissue'])
        d[fam+'_LDH_sig_tissues']='|'.join(f.loc[f.trait_q_trait_global4tissue_LDH<.05,'tissue'])
        d[fam+'_CPI_sig_tissues']='|'.join(f.loc[f.trait_q_trait_global4tissue_CPI<.05,'tissue'])
    summary.append(d)
save(frozen.merge(pd.DataFrame(summary),on='gene',validate='one_to_one'),I/'UPDATED_GENE_EVIDENCE_MATRIX.tsv')
ext=read(N/'02_external/EXTERNAL_FIXED18_ALL_STATUS.tsv')
locus=read(O/'01_evidence_tables/LOCUS_MASTER_18.tsv')
links=read(N/'02_external/EXTERNAL_EXACT_AND_PROXY_SUPPORT.tsv')
ext['selected_published_proxies']=ext.lead_snp.map(links[links.evidence_type=='PUBLISHED_VARIANT_SUPPORT'].groupby('lead').variant.agg('|'.join))
save(ext,N/'02_external/EXTERNAL_FIXED18_ALL_STATUS.tsv')
assert list(locus.lead_snp)==lock['leads']
save(locus.merge(ext[[x for x in ext if x not in locus or x=='lead_snp']],on='lead_snp',validate='one_to_one'),I/'UPDATED_LOCUS_MASTER_18.tsv')
print(counts.to_string(index=False))
print(focus[(focus.gene=='DCC')&focus.ensembl.notna()][['family','gene','n_weight_snps','trait_n_allele_aligned_LDH','trait_n_allele_aligned_CPI','trait_P_LDH','trait_P_CPI','q_conjunction','trait_CPI_phi_P_CPI']].to_string(index=False))

