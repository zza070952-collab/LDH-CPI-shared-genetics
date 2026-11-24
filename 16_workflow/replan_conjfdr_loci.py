#!/usr/bin/env python3
"""New CPI core only: audited SNP export and explicit local FUMA-equivalent LD loci."""
import csv,gzip,hashlib,json,math,os,resource,shutil,subprocess,sys,time,traceback
from pathlib import Path
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='2',MKL_NUM_THREADS='1')
import numpy as np
import pandas as pd
from scipy.io import loadmat
ROOT=Path(sys.argv[1]); PAIR=ROOT/'05_conjfdr/LDH_CPI_low_overlap'; OUT=PAIR/'ld_loci'; OUT.mkdir(exist_ok=True)
LOG=ROOT/'logs/conjfdr_cpi_loci';LOG.mkdir(exist_ok=True)
REF=ROOT.parent/'01_data_raw/conjfdr_reference/9545380.ref'
GENO=Path('/analysis/projects/hipOA_spatial_multiomics_clean/resources/active/extracted/magma_g1000_eur/g1000_eur')
PLINK=Path('/analysis/miniconda3/envs/pops_env/bin/plink')
SHA_EXPECT={'bed':'fadce7321e622ed20e3cb7798b3e54699c7f342752c642a14ddb982d93d8e08c','bim':'b7dfeb30cef0573059ae43dd85a6924e46e418d6fe03df6b0ff1166c79cfa131','fam':'79dee21226fd04b5a889cbbeacd08079e1740e7d63b0a4dea115527e8fea4dae'}
PARAMETERS={'r2_definition':'squared dosage correlation in every stage, FUMA precomputed-LD definition','indep_r2':.6,'lead_r2':.1,'merge_bp':250000,'candidate_r2':.6,'candidate_window':'chromosome-wide','MAF':'>0.01','exclude':'GRCh37 6:25000000-34000000'}
RESULT_COLS=['SNP','CHR','BP','A1_REF_EFFECT','A2_REF_OTHER','conjFDR','condFDR_LDH_given_CPI','condFDR_CPI_given_LDH','LDH_BETA_RAW_ALIGNED','LDH_SE_RAW','LDH_P_RAW','CPI_BETA_RAW_ALIGNED','CPI_SE_RAW','CPI_P_RAW','MHC_flag','status']
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(8*1024**2),b''):h.update(b)
 return h.hexdigest()
def atom(p,text):
 p=Path(p);t=p.with_name(p.name+'.tmp');t.write_text(text);os.replace(t,p)
def js(p,x):atom(p,json.dumps(x,indent=2,allow_nan=False)+'\n')
def tab(p,rows,cols):
 with open(str(p)+'.tmp','w',newline='') as f:
  w=csv.DictWriter(f,cols,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rows)
 os.replace(str(p)+'.tmp',p)
def guard():
 if (ROOT/'00_admin/PAUSE_REQUESTED').exists():raise RuntimeError('PAUSE_REQUESTED')
 assert shutil.disk_usage(ROOT).free>=80*1024**3,'80 GiB free reserve'
 assert sum(p.stat().st_size for p in OUT.rglob('*') if p.is_file())<950*1024**2,'1 GiB new output budget (950 MiB safety stop)'
def state(status,**kwargs):js(LOG/'STATUS.json',dict(status=status,updated=time.strftime('%Y-%m-%dT%H:%M:%S%z'),**kwargs))
def matchsign(a,b,x,y):
 if a==x and b==y:return 1
 if a==y and b==x:return -1
 tr=str.maketrans('ACGT','TGCA')
 if len(a)==len(b)==len(x)==len(y)==1:
  if a.translate(tr)==x and b.translate(tr)==y:return 1
  if a.translate(tr)==y and b.translate(tr)==x:return -1
 return 0
def fmt(v):return format(float(v),'.17g') if np.isfinite(v) else 'NA'
def mergeblocks(blocks):
 merged=[]
 for g in sorted(blocks,key=lambda x:(x['start'],x['end'],x['lead'])):
  if merged and g['start']-merged[-1]['end']<=250000:
   q=merged[-1];q['end']=max(q['end'],g['end']);q['leads'].append(g['lead']);q['indep']+=g['indep'];q['candidates'].update(g['candidates'])
  else:merged.append(dict(start=g['start'],end=g['end'],leads=[g['lead']],indep=g['indep'][:],candidates=g['candidates'].copy()))
 return merged
def software_tests():
 assert matchsign('A','G','A','G')==1 and matchsign('G','A','A','G')==-1
 assert matchsign('T','C','A','G')==1 and matchsign('C','T','A','G')==-1
 assert matchsign('AC','A','A','AC')==-1 and matchsign('AC','A','G','C')==0
 def b(s,e,n):return dict(start=s,end=e,lead=n,indep=[n],candidates={n:s})
 assert len(mergeblocks([b(1,100,'a'),b(250100,250200,'b')]))==1
 assert len(mergeblocks([b(1,100,'a'),b(250101,250200,'b')]))==2
 assert len(mergeblocks([b(1,1000000,'a'),b(100,200,'b'),b(1200000,1500000,'c')]))==1
 assert len(mergeblocks([b(1,10,'a'),b(250010,250020,'b'),b(500020,500030,'c')]))==1
 ds={s:dict(conjFDR=p,BP=i) for i,(s,p) in enumerate([('a',.001),('b',.002),('c',.003)])}
 chosen,mapping=greedy(set(ds),{'a':{'b'},'b':{'a','c'},'c':{'b'}},ds)
 assert chosen==['a','c'] and mapping=={'a':'a','b':'a','c':'c'}
 lead,mapping=greedy(set(chosen),{'a':{'c'},'c':{'a'}},ds)
 assert lead==['a'] and mapping=={'a':'a','c':'a'}
 js(OUT/'SOFTWARE_UNIT_TESTS.json',dict(status='PASS',scope='Allele sign and inclusive/transitive interval-merging software tests; toy intervals are NOT biological data or results.'))
def export(core,qc,lock):
 done=OUT/'EXPORT_QC.json'
 if done.exists():
  d=json.loads(done.read_text());assert d['lock']==lock and all(sha(OUT/p)==h for p,h in d['hashes'].items());return
 state('RUNNING_EXPORT');guard()
 ref=pd.read_csv(REF,sep='\t',usecols=['SNP','CHR','BP','A1','A2'],dtype={'SNP':str,'CHR':'int8','BP':'int32','A1':'category','A2':'category'})
 assert len(ref)==9545380 and not ref.SNP.duplicated().any()
 ix=pd.Index(ref.SNP);n=len(ref);cc=ref.CHR.to_numpy();bp=ref.BP.to_numpy();ra=ref.A1.to_numpy();rb=ref.A2.to_numpy()
 m=loadmat(core,variable_names=['fdrmat','fdrmat12','fdrmat21','chrnumvec','posvec'])
 conj,co12,co21=(m[k].ravel() for k in ['fdrmat','fdrmat12','fdrmat21'])
 assert np.array_equal(m['chrnumvec'].ravel(),cc) and np.array_equal(m['posvec'].ravel(),bp),'reference row order mismatch'
 valid=np.isfinite(conj)&np.isfinite(co12)&np.isfinite(co21)
 assert int(valid.sum())==qc['finite_both_cond_and_conj']
 assert np.allclose(conj[valid],np.maximum(co12[valid],co21[valid]),atol=1e-12,rtol=0)
 arrays={};counts={}
 unmatched=OUT/'FULLQC_UNMATCHED_REFERENCE.tsv.gz'
 with gzip.open(str(unmatched)+'.tmp','wt',compresslevel=6) as uf:
  uf.write('trait\tSNP\tCHR\tBP\treason\n')
  for trait,file in [('LDH','FINNGEN_R12_LDH'),('CPI','UKB_CPI_WANG2025')]:
   g=ROOT/f'02_qc/gwas/{file}.GRCh37.QC.tsv.gz';meta=json.loads((PAIR/f'{trait}.input_qc.json').read_text())
   assert sha(g)==meta['signature']['gwas_sha256']
   assert meta['signature']['reference_sha256']==lock['reference_sha256']
   vals=np.full((n,3),np.nan);seen=np.zeros(n,dtype='uint8');cnt={'read':0,'aligned':0,'unmatched':0,'coordinate_mismatch':0,'allele_mismatch':0,'duplicate':0}
   for chunk in pd.read_csv(g,sep='\t',usecols=['SNP','CHR','BP','A1','A2','BETA','SE','P'],chunksize=150000,keep_default_na=False):
    chunk=chunk[['SNP','CHR','BP','A1','A2','BETA','SE','P']]
    inds=ix.get_indexer(chunk.SNP);cnt['read']+=len(chunk)
    for row,i in zip(chunk.itertuples(index=False,name=None),inds):
     snp,c,p,a,b,beta,se,pval=row;reason=''
     if i<0:reason='UNMATCHED_REFERENCE_ID';cnt['unmatched']+=1
     elif seen[i]:vals[i,:]=np.nan;seen[i]=2;reason='DUPLICATE_ALL_EXCLUDED';cnt['duplicate']+=1
     else:
      seen[i]=1
      if int(c)!=cc[i] or int(p)!=bp[i]:reason='COORDINATE_MISMATCH';cnt['coordinate_mismatch']+=1
      else:
       sg=matchsign(a,b,ra[i],rb[i])
       if not sg:reason='ALLELE_MISMATCH';cnt['allele_mismatch']+=1
       else:
        try:vals[i]=[sg*float(beta),float(se),float(pval)];cnt['aligned']+=1
        except ValueError:reason='INVALID_NUMERIC'
     if reason:uf.write(f'{trait}\t{snp}\t{c}\t{p}\t{reason}\n')
    guard();print('ALIGN',trait,cnt,flush=True)
   assert sha(PAIR/f'{trait}.mat')==meta['mat_sha256'],'pre-GC input MAT checksum changed'
   pre=loadmat(PAIR/f'{trait}.mat',variable_names=['zvec','logpvec']);z=pre['zvec'].ravel();lp=pre['logpvec'].ravel()
   assert len(z)==n==len(lp) and np.all(np.isfinite(z[valid])&np.isfinite(lp[valid])),'tested core not closed to pre-GC inputs'
   check=np.isfinite(z)&np.isfinite(lp)&(vals[:,2]>0)&(vals[:,2]<=1)
   assert check.sum()>1000000
   assert np.allclose(-np.log10(vals[check,2]),lp[check],rtol=1e-8,atol=1e-8),'raw P differs from locked official-converter input'
   nz=check&(z!=0)&(vals[:,0]!=0);assert np.all(np.sign(z[nz])==np.sign(vals[nz,0])),'raw BETA effect allele inconsistent with converter Z'
   cnt['converter_P_and_BETA_sign_checks']=int(check.sum());arrays[trait]=vals;counts[trait]=cnt;del pre,z,lp,seen
   cnt['aligned_rows_before_duplicate_exclusion']=cnt.pop('aligned');cnt['final_finite_aligned_rows']=int(np.isfinite(vals).all(axis=1).sum())
   assert np.isfinite(vals[valid,:]).all() and np.all((vals[valid,2]>0)&(vals[valid,2]<=1)) and np.all(vals[valid,1]>0), 'Tested core SNP missing original aligned BETA/SE/P'
 os.replace(str(unmatched)+'.tmp',unmatched)
 cols=RESULT_COLS
 dest=OUT/'CONJFDR_SNP_LEVEL_ALL.tsv.gz';sig=[]
 with gzip.open(str(dest)+'.tmp','wt',compresslevel=6) as f:
  f.write('\t'.join(cols)+'\n')
  for i in range(n):
   mhc=cc[i]==6 and 25000000<=bp[i]<=34000000
   row=[ref.SNP.iat[i],str(cc[i]),str(bp[i]),ra[i],rb[i],fmt(conj[i]),fmt(co12[i]),fmt(co21[i])]+[fmt(x) for t in ['LDH','CPI'] for x in arrays[t][i]]+[str(int(mhc)),'TESTED' if valid[i] else 'NOT_TESTABLE']
   f.write('\t'.join(row)+'\n')
   if valid[i] and conj[i]<.05 and not mhc:sig.append(dict(zip(cols,row)))
   if i%500000==0:guard();print('EXPORT',i,n,flush=True)
 os.replace(str(dest)+'.tmp',dest);tab(OUT/'SIGNIFICANT_SNP_INPUT.tsv',sig,cols)
 assert len(sig)==qc['conj_lt_0_05_SNPs_NOT_LOCI'],'Unexpected finite significant MHC or core count drift'
 d=dict(status='SUCCESS',lock=lock,n_reference_rows=n,n_tested=int(valid.sum()),n_primary_nonMHC=len(sig),n_sensitivity_nonMHC=sum(float(x['conjFDR'])<.01 for x in sig),counts=counts,raw_effect_definition='Original full-QC BETA sign-aligned to reference A1; SE/P unchanged. Neither GC Z nor derived Z is used as BETA.',hashes={p.name:sha(p) for p in [dest,OUT/'SIGNIFICANT_SNP_INPUT.tsv',unmatched]})
 js(done,d);guard()
def plinkrun(args,prefix,expects,allow_empty_maf=False):
 done=Path(str(prefix)+'.done.json');cmd=[str(PLINK),*map(str,args),'--threads','2','--memory','6000','--out',str(prefix)]
 if done.exists():
  d=json.loads(done.read_text());assert d['cmd']==cmd and all(sha(p)==h for p,h in d['hashes'].items());return d
 guard();js(str(prefix)+'.command.json',cmd)
 def rss(pid):
  try:
   for line in Path(f'/proc/{pid}/status').read_text().splitlines():
    if line.startswith('VmRSS:'):return int(line.split()[1])*1024
  except FileNotFoundError:return 0
  return 0
 child_budget=8*1024**3-rss(os.getpid())-256*1024**2
 assert child_budget>1024**3,'Insufficient memory headroom for LD process'
 with open(str(prefix)+'.stdout.log','a') as f:
  proc=subprocess.Popen(cmd,stdout=f,stderr=subprocess.STDOUT,preexec_fn=lambda:resource.setrlimit(resource.RLIMIT_AS,(child_budget,child_budget)))
  while proc.poll() is None:
   try:
    guard();assert rss(os.getpid())+rss(proc.pid)<8*1024**3-128*1024**2,'8 GiB task-total RSS safety stop'
   except Exception:
    proc.terminate()
    try:proc.wait(timeout=15)
    except subprocess.TimeoutExpired:proc.kill();proc.wait()
    atom(str(prefix)+'.exitcode',str(proc.returncode)+'\n');raise
   time.sleep(5)
  rc=proc.returncode
 atom(str(prefix)+'.exitcode',str(rc)+'\n')
 if rc!=0 and allow_empty_maf:
  logpath=Path(str(prefix)+'.log');txt=logpath.read_text() if logpath.exists() else ''
  if ('All variants removed due to minor allele threshold(s)' in txt) or ('No variants remaining' in txt and ('minor allele' in txt or 'after --maf' in txt)):
   d=dict(cmd=cmd,exit_status=rc,status='ZERO_ELIGIBLE_LD_REFERENCE_NOT_A_SUCCESSFUL_LD_CALCULATION',hashes={str(logpath):sha(logpath)})
   js(done,d);return d
 assert rc==0,f'PLINK failed {prefix}'
 for p in expects:assert Path(p).is_file() and Path(p).stat().st_size>0,f'Missing expected {p}'
 if '--clump' in cmd:
  cl=Path(str(prefix)+'.clumped')
  if cl.exists():expects=list(expects)+[cl]
  else:assert 'No significant --clump results' in Path(str(prefix)+'.log').read_text(), 'Missing clumped output without explicit zero-result log'
 d=dict(cmd=cmd,exit_status=rc,status='SUCCESS',hashes={str(p):sha(p) for p in expects});js(done,d);guard();return d
def clumpread(p):
 if not Path(p).exists():return []
 with open(p) as f:
  lines=[x.split() for x in f if x.strip()]
 return [dict(zip(lines[0],r)) for r in lines[1:]]
def greedy(nodes,adj,sd):
 """FUMA-style P-ranked independence on one explicitly defined LD relation."""
 mapping={};chosen=[]
 for snp in sorted(nodes,key=lambda x:(float(sd[x]['conjFDR']),int(sd[x]['BP']),x)):
  if snp in mapping:continue
  chosen.append(snp);mapping[snp]=snp
  for other in adj.get(snp,set()):
   if other in nodes and other not in mapping:mapping[other]=snp
 assert set(mapping)==set(nodes)
 for snp in chosen:assert not (set(chosen)-{snp})&adj.get(snp,set()),'retained variants are not LD independent'
 return chosen,mapping
def process_ld(lock):
 state('RUNNING_LD');guard();sig=list(csv.DictReader(open(OUT/'SIGNIFICANT_SNP_INPUT.tsv'),delimiter='\t'));sd={r['SNP']:r for r in sig}
 matched={};dups=set()
 with open(str(GENO)+'.bim') as f:
  for line in f:
   c,s,g,p,a,b=line.split()
   if s not in sd:continue
   if s in matched:dups.add(s)
   matched[s]=(c,p,a,b)
 missing=[];eligible=[]
 for s,r in sd.items():
  reason=None
  if s not in matched:reason='NO_LD_REFERENCE_ID'
  elif s in dups:reason='DUPLICATE_LD_REFERENCE_ID'
  else:
   c,p,a,b=matched[s]
   if int(c)!=int(r['CHR']) or int(p)!=int(r['BP']):reason='LD_REFERENCE_COORDINATE_MISMATCH'
   elif not matchsign(a,b,r['A1_REF_EFFECT'],r['A2_REF_OTHER']):reason='LD_REFERENCE_ALLELE_MISMATCH'
  if reason:missing.append(dict(r,reason=reason))
  else:eligible.append(r)
 tab(OUT/'SIGNIFICANT_UNMATCHED_LD_REFERENCE.tsv',missing,RESULT_COLS+['reason'])
 tab(OUT/'SIGNIFICANT_LD_MATCHED.tsv',eligible,RESULT_COLS)
 excludes=OUT/'MHC.exclude.range';atom(excludes,'6\t25000000\t34000000\tMHC_NON_PRIMARY\n')
 summaries=[];allind=[];alllead=[];allloc=[];allmaf=[]
 canddest=OUT/'CONJFDR_CANDIDATE_LD_PROXIES.tsv.gz'
 candfile=gzip.open(str(canddest)+'.tmp','wt')
 cw=csv.DictWriter(candfile,['threshold','locus_id','chr','independent_SNP','candidate_SNP','candidate_BP','r2','proxy_in_significant_set'],delimiter='\t');cw.writeheader()
 for chrom in range(1,23):
  ss=[x for x in eligible if int(x['CHR'])==chrom]
  if not ss:continue
  cd=OUT/f'chr{chrom}';cd.mkdir(exist_ok=True);sf=cd/'sig.tsv';tab(sf,ss,['SNP','conjFDR']);ids=cd/'sig.snps';atom(ids,''.join(r['SNP']+'\n' for r in ss))
  base=['--bfile',GENO,'--chr',chrom,'--maf','0.010000000001','--exclude','range',excludes]
  fp=cd/'maf_qc';mafstate=plinkrun(base+['--extract',ids,'--freq','--write-snplist'],fp,[Path(str(fp)+'.log'),Path(str(fp)+'.frq'),Path(str(fp)+'.snplist')],allow_empty_maf=True)
  freqid=set(Path(str(fp)+'.snplist').read_text().split()) if mafstate['exit_status']==0 else set();assert freqid<={r['SNP'] for r in ss}
  mafex=[dict(r,reason='GENOTYPE_MAF_NOT_GT_0_01_OR_REFERENCE_FILTER') for r in ss if r['SNP'] not in freqid];allmaf+=mafex
  tab(cd/'SIGNIFICANT_MAF_FILTERED.tsv',mafex,list(ss[0])+['reason'])
  adjacency={}
  for code,cut in [('r06',.6),('r01',.1)]:
   adj={s:set() for s in freqid};ap=cd/('significant_pairs_'+code)
   if len(freqid)>1:
    agz=Path(str(ap)+'.ld.gz')
    plinkrun(base+['--extract',ids,'--r2','gz','--ld-window','99999999','--ld-window-kb','1000000','--ld-window-r2',cut],ap,[Path(str(ap)+'.log'),agz])
    with gzip.open(agz,'rt') as af:
     ah=next(af).split()
     for line in af:
      r=dict(zip(ah,line.split()));a,b=r['SNP_A'],r['SNP_B'];assert a in freqid and b in freqid
      # Exact threshold membership is computed by PLINK, not inferred from rounded displayed r2.
      if a!=b:adj[a].add(b);adj[b].add(a)
   adjacency[code]=adj
  assert all(adjacency['r06'][s]<=adjacency['r01'][s] for s in freqid),'LD threshold relation is not nested'
  sets={}
  for tag,thr in [('p05',.05),('p01',.01)]:
   nodes={s for s in freqid if float(sd[s]['conjFDR'])<thr}
   indep,sigmap=greedy(nodes,adjacency['r06'],sd);leads,mapping=greedy(set(indep),adjacency['r01'],sd)
   tab(cd/(tag+'.independent_mapping.tsv'),[dict(SNP=s,independent_SNP=v) for s,v in sigmap.items()],['SNP','independent_SNP'])
   tab(cd/(tag+'.lead_mapping.tsv'),[dict(independent_SNP=s,lead_SNP=v) for s,v in mapping.items()],['independent_SNP','lead_SNP'])
   sets[tag]=([sd[s] for s in indep],[sd[s] for s in leads],mapping)
  indep05={r['SNP'] for r in sets['p05'][0]}
  assert {r['SNP'] for r in sets['p01'][0]}<=indep05,'Nested significance did not preserve greedy independent subset'
  if not indep05:
   js(cd/'CHROMOSOME_QC.json',dict(status='NO_LD_ELIGIBLE_REFERENCE',chr=chrom,n_sig=len(ss),n_MAF_excluded=len(mafex),n_independent_loci_resolved=0,interpretation='Reference coverage gap; NOT evidence of no shared association'))
   continue
  iq=cd/'indep_union.snps';atom(iq,''.join(s+'\n' for s in sorted(indep05)))
  pp=cd/'proxies';gz=Path(str(pp)+'.ld.gz')
  plinkrun(base+['--r2','gz','--ld-snp-list',iq,'--ld-window','99999999','--ld-window-kb','1000000','--ld-window-r2','.6'],pp,[Path(str(pp)+'.log'),gz])
  proxies={s:{s:(int(sd[s]['BP']),1.)} for s in indep05}
  with gzip.open(gz,'rt') as f:
   fields=next(f).split()
   for line in f:
    r=dict(zip(fields,line.split()));a,b=r['SNP_A'],r['SNP_B'];v=float(r['R2'])
    if a not in proxies:continue
    assert int(r['CHR_A'])==chrom==int(r['CHR_B']);assert v>=.6-1e-6
    proxies[a][b]=(int(r['BP_B']),v)
  for tag,thr in [('p05',.05),('p01',.01)]:
   inds,leads,mapping=sets[tag];groups={}
   for s,lead in mapping.items():
    d=groups.setdefault(lead,{'lead':lead,'indep':[],'candidates':{}});d['indep'].append(s);d['candidates'].update({p:b for p,(b,r) in proxies[s].items()})
   blocks=[]
   for lead,g in groups.items():g['start']=min(g['candidates'].values());g['end']=max(g['candidates'].values());blocks.append(g)
   merged=mergeblocks(blocks)
   for k,g in enumerate(merged,1):
    locus=f'LDH_CPI_{tag}_chr{chrom}_{k:03d}';best=min(g['leads'],key=lambda x:(float(sd[x]['conjFDR']),int(sd[x]['BP']),x));r=sd[best]
    allloc.append(dict(trait_pair='LDH_CPI_LOW',threshold=thr,locus_id=locus,chr=chrom,start=g['start'],end=g['end'],lead_SNP=best,lead_position=r['BP'],conjFDR=r['conjFDR'],n_independent_significant=len(g['indep']),n_lead_SNPs=len(g['leads']),n_candidate_LD_proxies=len(g['candidates']),lead_SNPs=';'.join(g['leads']),independent_significant_SNPs=';'.join(g['indep']),MHC_flag=0,status='LOCAL_FUMA_EQUIVALENT'))
    for s in g['indep']:
     allind.append(dict(sd[s],threshold=thr,locus_id=locus,lead_SNP=mapping[s]))
     for p,(bp1,r2) in proxies[s].items():cw.writerow(dict(threshold=thr,locus_id=locus,chr=chrom,independent_SNP=s,candidate_SNP=p,candidate_BP=bp1,r2=r2,proxy_in_significant_set=int(p in sd and float(sd[p]['conjFDR'])<thr)))
    for s in g['leads']:alllead.append(dict(sd[s],threshold=thr,locus_id=locus,is_locus_representative=int(s==best)))
   assert all(b['start']-a['end']>250000 for a,b in zip(merged,merged[1:])), 'unmerged loci remain'
  js(cd/'CHROMOSOME_QC.json',dict(status='SUCCESS',chr=chrom,n_sig=len(ss),n_indep05=len(indep05),n_indep01=len(sets['p01'][0]),proxy_queries_chromosome_wide=True));guard()
 for thr in [.05,.01]:
  summaries.append(dict(trait_pair='LDH_CPI_LOW',conjFDR_threshold=thr,n_significant_nonMHC_SNPs=sum(float(r['conjFDR'])<thr for r in sig),n_unmatched_reference=sum(float(r['conjFDR'])<thr for r in missing),n_genotype_MAF_excluded=sum(float(r['conjFDR'])<thr for r in allmaf),n_independent_significant_SNPs=sum(r['threshold']==thr for r in allind),n_lead_SNPs=sum(r['threshold']==thr for r in alllead),n_independent_genomic_loci=sum(r['threshold']==thr for r in allloc),status='SUCCESS_LOCAL_FUMA_EQUIVALENT'))
  d=summaries[-1];d['reference_coverage_status']='COMPLETE' if d['n_unmatched_reference']+d['n_genotype_MAF_excluded']==0 else ('NO_LD_ELIGIBLE_REFERENCE' if d['n_independent_significant_SNPs']==0 else 'PARTIAL_REFERENCE_COVERAGE')
 tab(OUT/'SIGNIFICANT_GENOTYPE_MAF_EXCLUDED.tsv',allmaf,RESULT_COLS+['reason'])
 tab(OUT/'CONJFDR_INDEPENDENT_SIG_SNPS.tsv',allind,RESULT_COLS+['threshold','locus_id','lead_SNP'])
 tab(OUT/'CONJFDR_LEAD_SNPS.tsv',alllead,RESULT_COLS+['threshold','locus_id','is_locus_representative'])
 tab(OUT/'CONJFDR_GENOMIC_LOCI.tsv',allloc,['trait_pair','threshold','locus_id','chr','start','end','lead_SNP','lead_position','conjFDR','n_independent_significant','n_lead_SNPs','n_candidate_LD_proxies','lead_SNPs','independent_significant_SNPs','MHC_flag','status'])
 candfile.close();os.replace(str(canddest)+'.tmp',canddest)
 tab(OUT/'CONJFDR_PAIR_SUMMARY.tsv',summaries,list(summaries[0]));guard()
 outputs=[p for p in OUT.rglob('*') if p.is_file() and p.name not in ['FINAL_QC.json'] and not p.name.endswith('.tmp')]
 js(OUT/'FINAL_QC.json',dict(status='SUCCESS',lock=lock,summary=summaries,other_axes='QUEUED_OR_BLOCKED_NO_NEW_CORE_INPUT',hashes={str(p.relative_to(OUT)):sha(p) for p in outputs}));atom(LOG/'.done','SUCCESS\n');atom(LOG/'.exitcode','0\n');state('SUCCESS',summaries=summaries)
def main():
 software_tests()
 prelock={'code_sha256':sha(__file__),'parameters':PARAMETERS,'LD_sha256':SHA_EXPECT}
 prepath=OUT/'PRE_RESULT_PARAMETER_LOCK.json'
 if prepath.exists():assert json.loads(prepath.read_text())==prelock,'Pre-result LD parameters changed; explicit new module required'
 else:js(prepath,prelock)
 state('WAITING_FOR_NEW_CORE')
 while not ((ROOT/'logs/conjfdr_cpi/.done').exists() and (PAIR/'CORE_QC.json').exists()):
  if (ROOT/'00_admin/PAUSE_REQUESTED').exists():state('PAUSED');return
  time.sleep(30)
 qc=json.loads((PAIR/'CORE_QC.json').read_text());assert qc['status']=='SUCCESS' and qc['new_pair_500_prunings'] is True
 core=PAIR/'results/result_core.mat';assert sha(core)==qc['result_sha256']
 lock={'code_sha256':sha(__file__),'core_sha256':qc['result_sha256'],'reference_sha256':sha(REF),'LD_sha256':SHA_EXPECT,'parameters':PARAMETERS}
 lock['plink_binary_sha256']=sha(PLINK)
 lock['preGC_inputs']={}
 for trait in ['LDH','CPI']:
  meta=json.loads((PAIR/f'{trait}.input_qc.json').read_text());actual=sha(PAIR/f'{trait}.mat')
  assert meta['mat_sha256']==actual and meta['signature']['reference_sha256']==lock['reference_sha256']
  lock['preGC_inputs'][trait]={'mat_sha256':actual,'input_qc_sha256':sha(PAIR/f'{trait}.input_qc.json'),'gwas_sha256':meta['signature']['gwas_sha256']}
 lp=OUT/'INPUT_PARAMETER_LOCK.json'
 if lp.exists():assert json.loads(lp.read_text())==lock,'Changed postprocessing input/code; refuse checkpoint reuse'
 else:js(lp,lock)
 final=OUT/'FINAL_QC.json'
 if final.exists():
  d=json.loads(final.read_text());assert d['lock']==lock and all(sha(OUT/p)==h for p,h in d['hashes'].items());atom(LOG/'.done','SUCCESS\n');atom(LOG/'.exitcode','0\n');state('SUCCESS_ALREADY_VERIFIED');return
 audit=[]
 for ext,expected in SHA_EXPECT.items():
  p=Path(str(GENO)+'.'+ext);actual=sha(p);assert actual==expected,f'LD reference checksum mismatch {ext}';audit.append(dict(file=str(p),size=p.stat().st_size,sha256=actual))
 fam=[x.split() for x in open(str(GENO)+'.fam')];assert len(fam)==503 and all(x[2]=='0' and x[3]=='0' for x in fam), '503 founder LD reference required'
 tab(OUT/'LD_REFERENCE_CHECKSUMS.tsv',audit,['file','size','sha256'])
 atom(OUT/'PLINK_VERSION.txt',subprocess.check_output([str(PLINK),'--version'],text=True));atom(OUT/'PYTHON_VERSIONS.json',json.dumps(dict(python=sys.version,numpy=np.__version__,pandas=pd.__version__)))
 export(core,qc,lock);process_ld(lock)
if __name__=='__main__':
 resource.setrlimit(resource.RLIMIT_AS,(8*1024**3,8*1024**3))
 resource.setrlimit(resource.RLIMIT_FSIZE,(950*1024**2,950*1024**2))
 try:main()
 except Exception as e:
  state('FAILED',reason=str(e),traceback=traceback.format_exc());atom(LOG/'.exitcode','1\n');traceback.print_exc();raise
