#!/usr/bin/env python3
"""Trait-independent seven-column LDSC annotations, checkpointed by chromosome."""
import bisect, csv, gzip, hashlib, json, math, os, resource, subprocess, sys, time
from pathlib import Path
root=Path(sys.argv[1]); od=root/'03_sldsc_scavenge/ldscores';od.mkdir(exist_ok=True)
ann=root/'03_sldsc_scavenge/annotations'; ref=root/'03_sldsc_scavenge/ld_reference_audit'
base=Path('/analysis/projects/hipOA_spatial_multiomics_clean/resources/active/extracted')
py='/analysis/miniconda3/envs/ldsc/bin/python'; code='/analysis/projects/sciatica_multomics_bakeoff/work/software/src/ldsc/ldsc.py'
labels=['AST','ExNeu','IN','MG','OL','OPC','VASC']
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='1'
resource.setrlimit(resource.RLIMIT_AS,(8*1024**3,8*1024**3))
def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
def gate():
    if (root/'00_admin/PAUSE_REQUESTED').exists():raise SystemExit(75)
    n=sum(p.stat().st_size for p in od.rglob('*') if p.is_file())
    if n>300*1024**2:raise RuntimeError('LDscore output exceeds allocated 300 MiB')
    st=os.statvfs(str(od))
    if st.f_bavail*st.f_frsize<80*1024**3:raise RuntimeError('Free disk below 80 GiB safety floor')
def write_tsv(path,rows):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
while not (ref/'reference_audit.done').exists():gate();time.sleep(15)
assert (ann/'annotations.done').exists()
manifest=list(csv.DictReader((ann/'SLDSC_ANNOTATION_MANIFEST.tsv').open(),delimiter='\t'))
eligible=[r['cell_type'] for r in manifest if r['eligible_S_LDSC']=='TRUE' and int(r['n_peaks_hg19'])>0]
assert eligible==labels,eligible
lock={'labels':labels,'LDSC':'1.0.1','python':py,'code':code,'code_sha256':sha(code),'LD_window_cM':1,'threads':1,'memory_GiB':8,'maximum_output_MiB':300,'coordinate_overlap':'BED start <= BIM BP - 1 < BED end','annotation_type':'binary sum matrix counts >= 10 in author Adult label; cells library >1000','print_snps':'exact chromosome HM3 non-MHC weight SNP identifiers','construction_scope':'explicit input reconstruction supplement; NOT byte-identical mother annotations','source_BED_sha256':{x:sha(ann/(x+'.hg19.bed')) for x in labels}}
(od/'LD_SCORE_PARAMETERS_LOCK.json').write_text(json.dumps(lock,indent=2)+'\n')
for name in [py,code]:assert Path(name).is_file()
with (od/'software_versions.txt').open('w') as f:
    for cmd in [[py,'--version'],[py,'-c','import numpy,pandas,scipy,bitarray; print(numpy.__version__,pandas.__version__,scipy.__version__,bitarray.__version__)']]:subprocess.run(cmd,stdout=f,stderr=f,check=True)
with (od/'snATAC_Adult_final.ldcts').open('w') as f:
    for ct in labels:f.write(ct+'\t'+str(od/ct/(ct+'.'))+'\n')
allqc=[]
for c in range(1,23):
    gate();cp=od/f'chr{c:02d}';cp.mkdir(exist_ok=True);done=cp/'ldscore.done'
    if done.exists():
        prior=json.loads((cp/'qc.json').read_text())
        for file,h in prior['output_sha256'].items():assert sha(file)==h,('completed output changed',file)
        allqc.extend(prior['annotation_counts']);print('SKIP VERIFIED CHR',c,flush=True);continue
    bedpre=base/'ldsc_phase3_plink/1000G_EUR_Phase3_plink'/f'1000G.EUR.QC.{c}'
    bim=Path(str(bedpre)+'.bim'); intervals=[]
    for ct in labels:
        rows=[]
        with (ann/(ct+'.hg19.bed')).open() as f:
            for line in f:
                t=line.split()
                if t[0]=='chr'+str(c):rows.append((int(t[1]),int(t[2])))
        merged=[]
        for s,e in sorted(rows):
            if merged and s<=merged[-1][1]:merged[-1]=(merged[-1][0],max(e,merged[-1][1]))
            else:merged.append((s,e))
        intervals.append(([t[0] for t in merged],[t[1] for t in merged]))
    ap=cp/'Adult.annot.gz';counts=[0]*7;n=0
    with gzip.open(ap,'wt') as out,bim.open() as f:
        out.write('\t'.join(labels)+'\n')
        for line in f:
            bp=int(line.split()[3])-1;values=[]
            for i,(starts,ends) in enumerate(intervals):
                j=bisect.bisect_right(starts,bp)-1;v=int(j>=0 and bp<ends[j]);counts[i]+=v;values.append(str(v))
            out.write('\t'.join(values)+'\n');n+=1
    weights=base/'ldsc_phase3_weights/1000G_Phase3_weights_hm3_no_MHC'/f'weights.hm3_noMHC.{c}.l2.ldscore.gz'
    expected=[]
    with gzip.open(weights,'rt') as f,(cp/'print_snps.txt').open('w') as out:
        header=next(f).split();j=header.index('SNP')
        for line in f:snp=line.split()[j];expected.append(snp);out.write(snp+'\n')
    prefix=cp/'Adult';cmd=[py,code,'--l2','--bfile',str(bedpre),'--ld-wind-cm','1','--annot',str(ap),'--thin-annot','--print-snps',str(cp/'print_snps.txt'),'--out',str(prefix)]
    (cp/'command.json').write_text(json.dumps(cmd,indent=2)+'\n');print('START CHR',c,'reference_snps',n,'annotation_counts',counts,flush=True)
    with (cp/'task.stdout.log').open('a') as log:res=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT)
    (cp/'task.exitcode').write_text(str(res.returncode)+'\n')
    if res.returncode:raise RuntimeError('LDSC process failed chromosome '+str(c))
    logtext=(cp/'Adult.log').read_text()
    if 'Traceback' in logtext or 'Error' in logtext or 'Total time elapsed' not in logtext:raise RuntimeError('LDSC log not successful chromosome '+str(c))
    score=Path(str(prefix)+'.l2.ldscore.gz');mf=Path(str(prefix)+'.l2.M');m50f=Path(str(prefix)+'.l2.M_5_50')
    m=[float(x) for x in mf.read_text().split()];m50=[float(x) for x in m50f.read_text().split()];assert len(m)==len(m50)==7 and all(v>=0 for v in m50)
    assert all(abs(x-y)<1e-4 for x,y in zip(m,counts)),('M differs from binary counts',c)
    split=[]
    for ct in labels:
        dd=od/ct;dd.mkdir(exist_ok=True);split.append(gzip.open(dd/f'{ct}.{c}.l2.ldscore.gz','wt'))
    seen=[]
    with gzip.open(score,'rt') as f:
        h=next(f).split();assert h[:3]==['CHR','SNP','BP'] and h[3:]==[ct+'L2' for ct in labels],h
        for i,o in enumerate(split):o.write('\t'.join(h[:3]+[h[i+3]])+'\n')
        for line in f:
            t=line.split();assert len(t)==10 and all(math.isfinite(float(v)) for v in t[3:]);seen.append(t[1])
            for i,o in enumerate(split):o.write('\t'.join(t[:3]+[t[i+3]])+'\n')
    for o in split:o.close()
    assert seen==expected,('SNP output mismatch with weights',c,len(seen),len(expected))
    chrqc=[];outputs=[ap,score,mf,m50f]
    for i,ct in enumerate(labels):
        dd=od/ct
        for suffix,val in [('M',m[i]),('M_5_50',m50[i])]:p=dd/f'{ct}.{c}.l2.{suffix}';p.write_text(str(val)+'\n');outputs.append(p)
        outputs.append(dd/f'{ct}.{c}.l2.ldscore.gz')
        chrqc.append(dict(chr=c,cell_type=ct,n_reference_snps=n,n_annotated_snps=counts[i],annotation_fraction=counts[i]/n,n_common_annotated_snps=m50[i],n_hm3_nonmhc_ldscores=len(seen),status='SUCCESS'))
    record={'chr':c,'annotation_counts':chrqc,'input_sha256':{'bim':sha(bim),'weights':sha(weights),'annotations':sha(ap)},'output_sha256':{str(p):sha(p) for p in outputs},'exit_status':0,'finite_ldscores':True,'exact_HM3_nonMHC_order':True}
    (cp/'qc.json').write_text(json.dumps(record,indent=2)+'\n');done.write_text(time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())+'\n');allqc.extend(chrqc)
    write_tsv(od/'LD_SCORE_CHROMOSOME_QC.tsv',allqc);gate();print('SUCCESS CHR',c,flush=True)
write_tsv(od/'LD_SCORE_CHROMOSOME_QC.tsv',allqc)
(od/'ldscores.done').write_text(time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())+'\n')
print('SUCCESS ALL 22 CHROMOSOMES / 7 ANNOTATIONS',flush=True)
