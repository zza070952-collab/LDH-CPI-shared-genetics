#!/usr/bin/env python3
"""Run each new LOW-pair MiXeR stage once; preserve failures and quantitative QC."""
import csv,gzip,json,math,os,subprocess,sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import norm
from replan_run_task import digest
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'04_ldsc_mixer/LDH_CPI_protocol13';OUT.mkdir(parents=True,exist_ok=True)
V=ROOT/'16_workflow/vendor/mixer_protocol13_3d17e8a';LIB=V/'src/build/lib/libbgmg.so';LD=ROOT/'04_ldsc_mixer/reference_protocol13'
REF=Path('/analysis/projects/hipOA_spatial_multiomics_clean/resources/active/extracted/ldsc_phase3_plink/1000G_EUR_Phase3_plink')
records=[]
for ch in range(1,23):
 d=json.loads((LD/f'chr{ch:02}.done.json').read_text());p=LD/f'1000G.EUR.QC.{ch}.run4.ld'
 assert d['actual_exitcode']==0 and d['output_sha256']==digest(p) and d['fingerprint']['lib_sha256']==digest(LIB)
 for ext,h in d['fingerprint']['inputs'].items():assert digest(REF/f'1000G.EUR.QC.{ch}.{ext}')==h
 commit=d['fingerprint']['source_commit'];assert commit=='3d17e8a72e91982b186bc069a0506ef35d873b9e' or (ch==22 and commit=='3d17e8a72e91982b186bc069a0506ef35d873b9e4' and (ROOT/'00_admin/MIXER_REFERENCE_PROVENANCE_CORRECTION.md').exists())
 records.append({'chr':ch,'ld_sha256':digest(p),'bim_sha256':d['fingerprint']['inputs']['bim']})
(OUT/'REFERENCE22_VALIDATED.json').write_text(json.dumps(records,indent=2))
for label,source,ncol in [('LDH','FINNGEN_R12_LDH','NEFF'),('CPI','UKB_CPI_WANG2025','N')]:
 inp=ROOT/'02_qc/gwas'/f'{source}.GRCh37.QC.tsv.gz';out=OUT/f'{label}.full_nonMHC.tsv.gz';done=OUT/f'{label}.input.done.json'
 signature={'gwas':digest(inp),'code':digest(__file__),'ncol':ncol}
 if done.exists():
  old=json.loads(done.read_text());assert old['signature']==signature and digest(out)==old['sha256'];continue
 counts={'read':0,'written':0,'MHC':0,'invalid':0};seen=set();first=True;nmin=float('inf');nmax=0
 with gzip.open(out,'wt',compresslevel=6) as f:
  for x in pd.read_csv(inp,sep='\t',chunksize=250000):
   counts['read']+=len(x);assert not x.SNP.duplicated().any() and not any(s in seen for s in x.SNP),'Duplicate rsID in normalized QC input';seen.update(x.SNP)
   mask=np.isfinite(x.BETA)&np.isfinite(x.P)&(x.P>0)&(x.P<=1)&np.isfinite(x[ncol])&(x[ncol]>0)&(x.EAF>.01)&(x.EAF<.99)
   counts['invalid']+=int((~mask).sum());x=x.loc[mask].copy();mhc=(x.CHR==6)&x.BP.between(25000000,34000000);counts['MHC']+=int(mhc.sum());x=x.loc[~mhc].copy()
   x['Z']=np.sign(x.BETA)*norm.isf(x.P/2);x['N']=x[ncol];assert np.isfinite(x.Z).all()
   x[['SNP','N','Z','A1','A2']].to_csv(f,sep='\t',index=False,header=first);first=False;counts['written']+=len(x)
   if len(x):nmin=min(nmin,float(x.N.min()));nmax=max(nmax,float(x.N.max()))
 assert counts['written']>6000000
 done.write_text(json.dumps({'status':'SUCCESS','signature':signature,'sha256':digest(out),'counts':counts,'N_min':nmin,'N_max':nmax,'scope':'Full imputed GWAS, NOT HapMap3 or pleioFDR-reference-restricted; fitting uses official MAF>.05/r2.8/2M extracts; all-SNP test uses no extract.'},indent=2));del seen
COMMON=['--lib',str(LIB),'--bim-file',str(REF/'1000G.EUR.QC.@.bim'),'--ld-file',str(LD/'1000G.EUR.QC.@.run4.ld')]
stages=[]
def finite_params(x):
 if isinstance(x,dict):return all(finite_params(v) for v in x.values())
 if isinstance(x,list):return all(finite_params(v) for v in x)
 return math.isfinite(x) if isinstance(x,(float,int)) else True
def runstage(name,kind,args,inputs,output):
 checkpoint=OUT/(name+'.stage.json');cmd=[sys.executable,str(V/'precimed/mixer.py'),kind]+COMMON+args
 sig={'cmd':cmd,'inputs':{str(p):digest(p) for p in inputs},'LD22':records,'library':digest(LIB),'controller':digest(__file__)}
 if checkpoint.exists():
  prior=json.loads(checkpoint.read_text());assert prior['signature']==sig,'Changed stage inputs: explicit new stage needed'
  if prior['status']=='SUCCESS':assert output.is_file() and digest(output)==prior['output_sha256']
  stages.append(prior);return prior['status']=='SUCCESS'
 if (ROOT/'00_admin/PAUSE_REQUESTED').exists():raise SystemExit(75)
 with (OUT/(name+'.stdout.log')).open('a') as log:rc=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT).returncode
 (OUT/(name+'.exitcode')).write_text(str(rc)+'\n');row={'name':name,'kind':kind,'signature':sig,'actual_exitcode':rc,'status':'FAILED','reason':'Nonzero exit or missing output'}
 if rc==0 and output.is_file():
  try:
   if kind=='snps':
    with output.open() as f:ids=[s.strip() for s in f if s.strip()]
    assert len(ids)>100000 and len(ids)==len(set(ids));row['n_snps']=len(ids)
   else:
    data=json.loads(output.read_text());assert data['analysis']==('univariate' if kind.endswith('1') else 'bivariate') and data.get('params') and finite_params(data['params']) and data.get('ci')
    row['params']=data['params'];row['point_estimates']={k:v.get('point_estimate') for k,v in data['ci'].items()};row['optimization_diagnostics']=data.get('optimize')
    row['interpretation']='Valid numerical output, NOT proof of stability. Boundary/AIC/replicate consistency must be assessed across all3 seeds.'
   row.update(status='SUCCESS',reason='Exit0 and schema/core-numeric QC',output_sha256=digest(output))
  except Exception as e:row['reason']='OUTPUT_QC: '+str(e)
 checkpoint.write_text(json.dumps(row,indent=2));stages.append(row)
 print(name,row['status'],row['reason'],flush=True)
 (OUT/'STAGE_STATUS.json').write_text(json.dumps([{k:v for k,v in r.items() if k not in ['signature','optimization_diagnostics']} for r in stages],indent=2))
 return row['status']=='SUCCESS'
for seed in [9527,20260904,314159]:
 snps=OUT/f'maf05_rand2M_r2p8.seed{seed}.snps';name=f'snps_{seed}'
 if not runstage(name,'snps',['--out',str(snps),'--maf','.05','--subset','2000000','--r2','.8','--seed',str(seed)],[],snps):continue
 fits={}
 for trait in ['LDH','CPI']:
  full=OUT/f'{trait}.full_nonMHC.tsv.gz';base=OUT/f'{trait}.fit.seed{seed}';out=Path(str(base)+'.json')
  fits[trait]=out if runstage(f'{trait}_fit1_{seed}','fit1',['--trait1-file',str(full),'--extract',str(snps),'--out',str(base),'--seed',str(seed),'--threads','4'],[full,snps],out) else None
  if fits[trait]:
   base=OUT/f'{trait}.test.seed{seed}';runstage(f'{trait}_test1_{seed}','test1',['--trait1-file',str(full),'--load-params-file',str(fits[trait]),'--out',str(base),'--seed',str(seed),'--threads','4'],[full,fits[trait]],Path(str(base)+'.json'))
 if not all(fits.values()):continue
 fulls=[OUT/f'{t}.full_nonMHC.tsv.gz' for t in ['LDH','CPI']];base=OUT/f'LDH_CPI.fit.seed{seed}';out=Path(str(base)+'.json')
 ok=runstage(f'LDH_CPI_fit2_{seed}','fit2',['--trait1-file',str(fulls[0]),'--trait2-file',str(fulls[1]),'--trait1-params-file',str(fits['LDH']),'--trait2-params-file',str(fits['CPI']),'--extract',str(snps),'--out',str(base),'--seed',str(seed),'--threads','4'],fulls+list(fits.values())+[snps],out)
 if ok:
  base=OUT/f'LDH_CPI.test.seed{seed}';runstage(f'LDH_CPI_test2_{seed}','test2',['--trait1-file',str(fulls[0]),'--trait2-file',str(fulls[1]),'--load-params-file',str(out),'--out',str(base),'--seed',str(seed),'--threads','4'],fulls+[out],Path(str(base)+'.json'))
(OUT/'MODULE_STATUS.json').write_text(json.dumps({'controller':'COMPLETE','expected_stages':21,'terminal_stages':len(stages),'success':sum(r['status']=='SUCCESS' for r in stages),'failed':sum(r['status']=='FAILED' for r in stages),'quantitative_interpretation':'REQUIRES_REPLICATE_AIC_BOUNDARY_REVIEW','old11258_reused':False},indent=2))
