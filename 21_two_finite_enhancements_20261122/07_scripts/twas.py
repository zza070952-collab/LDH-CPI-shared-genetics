"""Bounded S-PrediXcan using the official association kernel with allele-aligned per-gene inputs."""
import os,sys,gzip,sqlite3,hashlib,json,logging,csv
from pathlib import Path
from collections import defaultdict,Counter
import numpy as np
import pandas as pd
from scipy.stats import norm
N=Path(__file__).resolve().parents[1];B=N.parent;A=next((N/'00_preanalysis').glob('LDH_*'));D=N/'01_twas';M=D/'models'
commit='fc735bf01a1cf43f25cef466439dcfa04712dd50'
software=B.parent.parent/'P4_LSS_drug_repurposing/10_workflow/vendor'/('MetaXcan-'+commit)/'software'
sys.path.insert(0,str(software))
from metax.metaxcan import AssociationCalculation as official
from metax import Constants
from metax.PredictionModel import WDBQF
logging.getLogger().setLevel(logging.WARNING)
family=sys.argv[1]
lock=json.loads((N/'00_admin/ANALYSIS_PLAN_LOCK.json').read_text())
files=pd.read_csv(A/'twas/MODEL_DOWNLOAD_CHECKSUMS.tsv',sep='\t').set_index('filename')
def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for buf in iter(lambda:f.read(4<<20),b''):h.update(buf)
    return h.hexdigest()
def write(df,p):df.to_csv(p,index=False,sep='\t',na_rep='NA')
models={};target=set()
for t in lock['tissues']:
    stem=family+'_'+t;db=M/(stem+'.db');cov=M/(stem+'.txt.gz')
    for p in [db,cov]:assert sha(p)==files.loc[p.name,'sha256'],p
    con=sqlite3.connect('file:'+str(db)+'?mode=ro',uri=True)
    w=pd.read_sql('select * from weights where weight != 0',con);info=pd.read_sql('select * from extra',con);con.close()
    models[t]=(w,info,cov);target.update(w.rsid.unique())
# Stream full frozen QC inputs, preserving signed beta/SE and actual sample sizes.
inputs={};qc=[];p_diffs={}
for trait,gw in zip(['LDH','CPI'],lock['GWAS']):
    path=Path(gw['full_QC_path']);assert sha(path)==gw['full_QC_sha256']
    data={};duplicate=set();total=0
    with gzip.open(path,'rt') as f:
        cols=f.readline().strip().split('\t');assert cols[0]=='SNP'
        for line in f:
            total+=1;rs=line.split('\t',1)[0]
            if rs not in target:continue
            if rs in data:duplicate.add(rs);continue
            r=dict(zip(cols,line.rstrip().split('\t')))
            data[rs]=(r['A1'],r['A2'],float(r['BETA']),float(r['SE']),float(r['P']),float(r['N']))
    for rs in duplicate:del data[rs]
    vals=np.array([(v[2]/v[3],v[4]) for v in data.values() if v[3]>0 and v[4]>0])
    diff=np.abs(np.log10(vals[:,1])-(np.log(2)+norm.logsf(abs(vals[:,0])))/np.log(10))
    p_diffs[trait]={'rows_scanned':total,'unique_model_rsids_present':len(data),'duplicated_dropped':len(duplicate),'median_abs_log10P_difference':float(np.median(diff)),'max_abs_log10P_difference':float(diff.max()),'abs_log10P_difference_gt1':int(sum(diff>1)),'P_zero':sum(v[4]==0 for v in data.values())}
    inputs[trait]=data
    print('GWAS',trait,total,len(data),flush=True)
(N/'00_admin'/('GWAS_INPUT_QC_'+family+'.json')).write_text(json.dumps(p_diffs,indent=2))
class Context:
    def provide_calculation(self,gene):return self.n,self.i,self.cov,None
    def get_n_in_covariance(self,gene):return self.nc
ctx=Context()
def state(v,ea,oa,ambiguous):
    if ambiguous:return 'NONUNIQUE_MODEL_RSID',None
    if v is None:return 'MISSING_GWAS',None
    a,b,beta,se,p,n=v
    if len(ea)!=1 or len(oa)!=1 or (set(ea+oa)-set('ACGT')):return 'NON_BIALLELIC_SNV',None
    if set(ea+oa) in [{'A','T'},{'C','G'}]:return 'PALINDROMIC',None
    if not np.isfinite(beta) or not np.isfinite(se) or se<=0:return 'INVALID_ESTIMATE',None
    comp=str.maketrans('ATCG','TAGC')
    for pair,mult,op in [((ea,oa),1,'MATCH'),((oa,ea),-1,'SWAP'),((ea.translate(comp),oa.translate(comp)),1,'COMPLEMENT'),((oa.translate(comp),ea.translate(comp)),-1,'COMPLEMENT_SWAP')]:
        if (a,b)==pair:return op,(beta*mult,beta/se*mult)
    return 'ALLELE_MISMATCH',None
for t,(w,extra,covfile) in models.items():
    out=D/(family+'_'+t+'.tests.tsv.gz');done=out.with_suffix('.sha256')
    if done.is_file() and out.is_file() and sha(out)==done.read_text().strip():continue
    info=extra.set_index('gene')
    key='varID' if family=='mashr' else 'rsid'
    bad=set(w.groupby('rsid').varID.nunique().loc[lambda x:x>1].index)
    variations=w[['rsid','varID','ref_allele','eff_allele']].drop_duplicates()
    align={};audit=[]
    for v in variations.itertuples():
        for tr in inputs:
            op,pair=state(inputs[tr].get(v.rsid),v.eff_allele,v.ref_allele,v.rsid in bad)
            align[(v.varID,tr)]=pair
            audit.append((v.rsid,v.varID,tr,op))
    write(pd.DataFrame(audit,columns=['source_rsid','model_varID','trait','alignment']),D/(family+'_'+t+'.alignment.tsv.gz'))
    blocks=defaultdict(list)
    with gzip.open(covfile,'rt') as f:
        assert f.readline().split()[:3]==['GENE','RSID1','RSID2']
        for line in f:
            g,a,b,val=line.split();blocks[g].append((a,b,float(val)))
    output=[];checks=[]
    for g,weight in w.groupby('gene',sort=False):
        gene=info.loc[g,'genename'] if g in info.index else None
        phi=info.loc[g,'phi'] if g in info.index and 'phi' in info else np.nan
        ids=weight[key].tolist();vs=weight.varID.tolist(); rs=weight.rsid.tolist()
        assert len(ids)==len(set(ids))
        positions={a:i for i,a in enumerate(ids)}
        C=np.full((len(ids),len(ids)),np.nan)
        for a,b,value in blocks.get(g,[]):
            if a in positions and b in positions:
                i,j=positions[a],positions[b];C[i,j]=C[j,i]=value
        diag=np.isfinite(C.diagonal()) & (C.diagonal()>0)
        mask={tr:np.array([align[(v,tr)] is not None for v in vs]) & diag for tr in inputs}
        weights=weight.weight.to_numpy();both=mask['LDH'] & mask['CPI']
        for tr in inputs:
            for scope,select in [('trait',mask[tr]),('common',both)]:
                ii=np.flatnonzero(select);used=len(ii);nw=len(ids)
                status='NO_ALIGNED_SNPS';z=p=var=cal=effect=np.nan
                if used:
                    matrix=C[np.ix_(ii,ii)]
                    if np.isfinite(matrix).all():
                        betas,zs=np.array([align[(vs[i],tr)] for i in ii]).T
                        ctx.n=nw;ctx.nc=int(diag.sum());ctx.cov=matrix
                        ctx.i={Constants.SNP:ii,Constants.ZSCORE:zs,Constants.BETA:betas,WDBQF.K_WEIGHT:weights[ii]}
                        result=official.association(g,ctx)
                        z,effect,var=result[1:4];p=2*norm.sf(abs(z))
                        status='ELIGIBLE' if np.isfinite(z) and var>0 else 'INVALID_VARIANCE'
                        if status=='ELIGIBLE' and (used/nw<.8 or (nw<=2 and used<nw)):status='LOW_COVERAGE'
                        if family=='en' and tr=='CPI' and used==nw and np.isfinite(phi):cal=2*norm.sf(abs(z)/np.sqrt(1+134627*.1012*phi))
                        if gene in ['DCC','APEH','BSN']:
                            ref_z=sum(weights[ii]*np.sqrt(matrix.diagonal())*zs)/np.sqrt(weights[ii]@matrix@weights[ii])
                            assert np.isclose(ref_z,z);checks.append([g,gene,tr,scope,z,ref_z])
                    else:status='INCOMPLETE_COVARIANCE'
                if not used and not diag.any():status='COVARIANCE_UNAVAILABLE'
                output.append([family,t,g,gene,tr,scope,nw,int(diag.sum()),sum(s in inputs[tr] for s in rs),used,used/nw,var,z,p,effect,status,phi,cal,hashlib.sha256('|'.join(vs[i] for i in ii).encode()).hexdigest()[:20]])
    columns='family tissue ensembl symbol trait predictor n_weight_snps n_in_covariance n_in_GWAS n_allele_aligned coverage_fraction var_g Z P effect status phi CPI_phi_P subset_hash'.split()
    write(pd.DataFrame(output,columns=columns),out)
    write(pd.DataFrame(checks,columns=['ensembl','symbol','trait','scope','official_Z','checked_Z']),D/(family+'_'+t+'.checks.tsv'))
    done.write_text(sha(out));print('DONE',t,len(output),flush=True)

