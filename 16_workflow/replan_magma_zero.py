"""User-authorized separate original-P MAGMA, reproducible local FUMA context."""
import csv,gzip,hashlib,json,os,shutil,subprocess,sys
from pathlib import Path
import numpy as np
import pandas as pd
from replan_run_task import digest
R=Path(__file__).resolve().parents[1];D=R/'06_fuma_magma_otg/magma_original_gwas';D.mkdir(exist_ok=True)
A=R/'06_fuma_magma_otg/magma_method_audit';resources=A/'official_resources';lock=json.loads((A/'MAGMA_METHOD_LOCK.json').read_text())
exe=Path(lock['magma_binary']);assert digest(exe)==lock['magma_binary_sha256'];ref=Path(lock['ld_reference']['prefix'])
for name in ['MAGMAgenev92','GTEx54v8','GTEx30v8','MSigDB7']:
 q=json.loads((resources/(name+'.provenance.json')).read_text());assert digest(resources/(name+'.resource'))==q['sha256']
def guard(extra=0):assert shutil.disk_usage(R).free>=80*1024**3+extra,'DISK_RESERVE'
def stage(name,args,inputs,outputs):
 p=D/(name+'.done.json');sig={'cmd':[str(exe)]+list(map(str,args)),'inputs':{str(x):digest(x) for x in inputs},'binary':digest(exe)}
 if p.exists():
  q=json.loads(p.read_text());assert q['signature']==sig and all(digest(x)==h for x,h in q['outputs'].items());return
 guard(32*1024**2)
 with (D/(name+'.stdout.log')).open('a') as f:rc=subprocess.run(sig['cmd'],cwd=D,stdout=f,stderr=subprocess.STDOUT).returncode
 (D/(name+'.exitcode')).write_text(str(rc)+'\n');assert rc==0 and all(x.is_file() and x.stat().st_size>100 for x in outputs),name
 for x in outputs:
  if x.name.endswith(('.genes.out','.gsa.out')):
   z=pd.read_csv(x,sep=r'\s+',comment='#');assert 'P' in z and len(z)>0 and np.isfinite(z.P).all() and z.P.between(0,1).all()
 p.write_text(json.dumps({'status':'SUCCESS','signature':sig,'actual_exitcode':rc,'outputs':{str(x):digest(x) for x in outputs}},indent=2));print('SUCCESS',name,flush=True)
gene=D/'Ensembl92.protein_coding.loc'
if not gene.exists():
 with gene.open('w') as f:
  for line in (resources/'MAGMAgenev92.resource').open():f.write('\t'.join(line.split()[:4])+'\n')
annot=D/'Ensembl92_window0';annotfile=Path(str(annot)+'.genes.annot')
stage('annotation',['--annotate','window=0,0','--snp-loc',str(ref)+'.bim','--gene-loc',gene,'--out',annot],[gene,Path(str(ref)+'.bim')],[annotfile])
trait=sys.argv[1];source={'LDH':'FINNGEN_R12_LDH','CPI':'UKB_CPI_WANG2025','SCIATICA':'PANUKB_SCIATICA_1476_EUR'}[trait]
gwas=R/f'02_qc/gwas/{source}.GRCh37.QC.tsv.gz';zipped=D/(trait+'.original_P_N.tsv.gz');prep=D/(trait+'.input.json')
sig={'gwas_sha256':digest(gwas),'method':'All valid original P and total/perSNP N, exclude GRCh37 chr6:25-34Mb; no conjFDR-as-P; Ensembl92 gene body window0','N':'N_not_NEFF'}
if prep.exists():q=json.loads(prep.read_text());assert q['signature']==sig and q['sha256']==digest(zipped)
else:
 n=0;excluded=0
 with gzip.open(str(zipped)+'.tmp','wt',compresslevel=6) as f:
  f.write('SNP\tP\tN\n')
  with gzip.open(gwas,'rt') as g:
   for i,row in enumerate(csv.DictReader(g,delimiter='\t')):
    if i%250000==0:guard(32*1024**2)
    if int(row['CHR'])==6 and 25000000<=int(row['BP'])<=34000000:excluded+=1;continue
    assert 0<float(row['P'])<=1 and float(row['N'])>0
    f.write(row['SNP']+'\t'+row['P']+'\t'+row['N']+'\n');n+=1
 os.replace(str(zipped)+'.tmp',zipped);assert n>6000000
 prep.write_text(json.dumps({'signature':sig,'sha256':digest(zipped),'n':n,'MHC_excluded':excluded},indent=2))
# Plain input is an explicitly disposable duplicate cache of the retained gzip input.
plain=D/(trait+'.duplicate_input_cache.tsv');gzsha=hashlib.sha256()
if not plain.exists():
 with gzip.open(zipped,'rb') as f,plain.open('wb') as out:
  for chunk in iter(lambda:f.read(8*1024**2),b''):guard(16*1024**2);out.write(chunk);gzsha.update(chunk)
else:
 with gzip.open(zipped,'rb') as f:
  for chunk in iter(lambda:f.read(8*1024**2),b''):gzsha.update(chunk)
assert digest(plain)==gzsha.hexdigest()
prefix=D/trait;genes=Path(str(prefix)+'.genes.out');raw=Path(str(prefix)+'.genes.raw')
stage(trait+'_gene',['--bfile',ref,'--pval',plain,'use=SNP,P','ncol=N','--gene-annot',annotfile,'--out',prefix],[zipped,annotfile,Path(str(ref)+'.bed'),Path(str(ref)+'.bim'),Path(str(ref)+'.fam')],[genes,raw])
assert plain.resolve().parent==D.resolve() and plain.name==trait+'.duplicate_input_cache.tsv' and digest(plain)==gzsha.hexdigest()
(D/(trait+'.cache_cleanup.json')).write_text(json.dumps({'path':str(plain),'retained_compressed_input':str(zipped),'decompressed_sha256':gzsha.hexdigest(),'gzip_sha256':digest(zipped),'reason':'Byte-identical regenerable duplicate input cache only; original GWAS and all results/logs retained'},indent=2));plain.unlink()
for label,flag,file,extra in [('GTEx54','--gene-covar','GTEx54v8',['--model','direction-covar=greater','condition-hide=Average']),('GTEx30','--gene-covar','GTEx30v8',['--model','direction-covar=greater','condition-hide=Average']),('MSigDB7','--set-annot','MSigDB7',[])]:
 out=D/(trait+'_'+label);stage(trait+'_'+label,['--gene-results',raw,flag,resources/(file+'.resource')]+extra+['--out',out],[raw,resources/(file+'.resource')],[Path(str(out)+'.gsa.out')])
(D/(trait+'.COMPLETE.json')).write_text(json.dumps({'status':'SUCCESS','trait':trait,'scope':'Original GWAS univariate MAGMA; not itself shared-axis proof; shared-locus convergence follows','gene_rows':len(pd.read_csv(genes,sep=r'\s+',comment='#'))},indent=2))
