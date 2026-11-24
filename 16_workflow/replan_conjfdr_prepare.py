#!/usr/bin/env python3
"""Full-GWAS adapter: exact GRCh37 positions, all duplicate IDs removed, real N."""
import gc,gzip,json,subprocess,sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import norm
from scipy.io import loadmat
from replan_run_task import digest
ROOT=Path(__file__).resolve().parents[1];P2=ROOT.parent
OUT=ROOT/'05_conjfdr/LDH_CPI_low_overlap';OUT.mkdir(parents=True,exist_ok=True)
REF=P2/'01_data_raw/conjfdr_reference/9545380.ref'
CONVERTER=P2/'10_workflow/vendor/python_convert-e46ebdfafd495c1420c7f8a4740a0da75c94d84d/sumstats.py'
for label,source,ncol in [('LDH','FINNGEN_R12_LDH','NEFF'),('CPI','UKB_CPI_WANG2025','N')]:
 inp=ROOT/'02_qc/gwas'/f'{source}.GRCh37.QC.tsv.gz';done=OUT/f'{label}.input_qc.json';out=OUT/f'{label}.tsv.gz';mat=OUT/f'{label}.mat'
 signature={'gwas_sha256':digest(inp),'reference_sha256':digest(REF),'adapter_sha256':digest(__file__),'converter_sha256':digest(CONVERTER),'N_definition':ncol}
 if done.exists():
  old=json.loads(done.read_text());assert old['signature']==signature,'Input changed; new stage required'
  assert digest(out)==old['tsv_sha256'] and digest(mat)==old['mat_sha256'],'Successful input changed'
  print('REUSE_EXACT_INPUT',label,flush=True);continue
 ids=pd.read_csv(inp,sep='\t',usecols=['SNP'])['SNP'];dups=set(ids[ids.duplicated(keep=False)]);del ids;gc.collect()
 ref=pd.read_csv(REF,sep='\t',usecols=['SNP','CHR','BP'],dtype={'CHR':'int8','BP':'int32'})
 rd=ref.SNP.duplicated(keep=False);nrd=int(rd.sum());ref=ref.loc[~rd].set_index('SNP').rename(columns={'CHR':'REF_CHR','BP':'REF_BP'})
 counts={'rows_read':0,'duplicate_input_ids':len(dups),'duplicate_reference_rows':nrd,'matched_id':0,'coordinate_mismatch':0,'invalid':0,'MHC_removed':0,'written':0}
 nmin=float('inf');nmax=0.;first=True
 with gzip.open(out,'wt',compresslevel=6) as f:
  for x in pd.read_csv(inp,sep='\t',chunksize=250000):
   counts['rows_read']+=len(x);x=x[~x.SNP.isin(dups)].join(ref,on='SNP',how='inner');counts['matched_id']+=len(x)
   same=(x.CHR==x.REF_CHR)&(x.BP==x.REF_BP);counts['coordinate_mismatch']+=int((~same).sum());x=x.loc[same].copy()
   good=np.isfinite(x.BETA)&np.isfinite(x.P)&(x.P>0)&(x.P<=1)&np.isfinite(x[ncol])&(x[ncol]>0)&(x.EAF>.01)&(x.EAF<.99)&x.A1.isin(list('ACGT'))&x.A2.isin(list('ACGT'))&(x.A1!=x.A2)
   counts['invalid']+=int((~good).sum());x=x.loc[good].copy();mhc=(x.CHR==6)&x.BP.between(25000000,34000000);counts['MHC_removed']+=int(mhc.sum());x=x.loc[~mhc].copy()
   x['PVAL']=x.P;x['Z']=np.sign(x.BETA)*norm.isf(x.P/2);x['N']=x[ncol]
   assert np.isfinite(x.Z).all(),'Nonfinite Z; no arbitrary P floor'
   if len(x):nmin=min(nmin,float(x.N.min()));nmax=max(nmax,float(x.N.max()))
   x[['SNP','CHR','BP','A1','A2','PVAL','Z','N']].to_csv(f,sep='\t',index=False,header=first);first=False;counts['written']+=len(x)
   print(label,counts,flush=True)
 del ref;gc.collect();assert counts['written']>1000000
 cmd=[sys.executable,str(CONVERTER),'mat','--sumstats',str(out),'--ref',str(REF),'--out',str(mat),'--force']
 (OUT/f'{label}.convert_command.json').write_text(json.dumps(cmd,indent=2))
 with (OUT/f'{label}.converter.log').open('a') as log:rc=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT).returncode
 (OUT/f'{label}.converter.exitcode').write_text(str(rc)+'\n');assert rc==0,'Official converter failed'
 data=loadmat(mat);assert all(k in data for k in ['zvec','logpvec','nvec']);assert all(data[k].shape==(9545380,1) for k in ['zvec','logpvec','nvec'])
 valid=np.isfinite(data['zvec'])&np.isfinite(data['logpvec'])&np.isfinite(data['nvec']);assert valid.sum()>1000000
 assert np.all(data['nvec'][valid]>0) and np.all(data['logpvec'][valid]>=0)
 assert float(np.nanmin(data['nvec']))==nmin and float(np.nanmax(data['nvec']))==nmax
 done.write_text(json.dumps({'status':'SUCCESS','signature':signature,'counts':counts,'N_min':nmin,'N_max':nmax,'mat_finite_all3':int(valid.sum()),'tsv_sha256':digest(out),'mat_sha256':digest(mat),'P_handling':'Require 0<P<=1; no replacement/floor. Ambiguous Z set NA by official converter; reference allele alignment official.'},indent=2))
 del data;gc.collect()
