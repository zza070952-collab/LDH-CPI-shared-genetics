#!/usr/bin/env python3
"""Run one fully specified SMR job; exact input/version/checkpoint gate."""
import argparse,csv,hashlib,json,re,subprocess,sys,fcntl
from pathlib import Path

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()
p=argparse.ArgumentParser();p.add_argument('--job-manifest',required=True);p.add_argument('--job-id',required=True);a=p.parse_args()
jobs=list(csv.DictReader(open(a.job_manifest),delimiter='\t')); jj=[j for j in jobs if j['job_id']==a.job_id]
assert len(jj)==1,'Job must be uniquely preregistered'
j=jj[0]
replan=Path('/analysis/projects/JTM_recent_highIF_direct_clone/P2_LDH_shared_architecture/11_NatCommun_strict_clone_Pro')
assert not(replan/'00_admin/PAUSE_REQUESTED').exists(),'PAUSE_REQUESTED: do not begin a new SMR job'
assert j['pairing_approved']=='YES','Await audited GWAS pairing'
assert j['gwas_build']=='GRCh37' and j['qtl_build']=='GRCh37','Build mismatch'
assert j['full_heidi_coverage']=='YES','P<1e-5 lite QTLs cannot provide complete HEIDI input'
assert j['full_probe_universe']=='YES','Candidate-selected QTL universe cannot yield formal global FDR'
out=Path(j['output_prefix']);out.parent.mkdir(parents=True,exist_ok=True)
job_lock=open(str(out)+'.job.lock','a')
fcntl.flock(job_lock.fileno(),fcntl.LOCK_EX)
assert not(replan/'00_admin/PAUSE_REQUESTED').exists(),'PAUSE_REQUESTED while waiting for job lock'
version=subprocess.run([j['smr_binary'],'--help'],capture_output=True,text=True).stdout
match=re.search(r'Version\s+([\d.]+)',version,re.I)
assert match,'Cannot verify SMR version'
actual=match.group(1); expected=j['smr_version']
assert actual==expected,'Binary version mismatch'
if actual not in ('1.02','1.0.2'):
    deviation=Path(j.get('version_deviation_file',''))
    assert deviation.is_file() and deviation.stat().st_size>0,'Nonmother version requires explicit recorded deviation; no silent replacement'
required=[j['smr_binary'],j['gwas_path'],j['qtl_prefix']+'.besd',j['qtl_prefix']+'.esi',j['qtl_prefix']+'.epi',j['bfile']+'.bed',j['bfile']+'.bim',j['bfile']+'.fam',j['pairing_matrix']]
cached={}
if j.get('ld_checksum_manifest'):
    for r in csv.DictReader(open(j['ld_checksum_manifest']),delimiter='\t'):
        stt=Path(r['file']).stat()
        if stt.st_size==int(r['size']) and stt.st_mtime_ns==int(r['mtime_ns']):cached[r['file']]=r['sha256']
checks={f:cached[f] if f in cached else sha(f) for f in required}
if j.get('gwas_sha256'):assert checks[j['gwas_path']]==j['gwas_sha256'],'GWAS changed since approval'
cmd=[j['smr_binary'],'--bfile',j['bfile'],'--gwas-summary',j['gwas_path'],'--beqtl-summary',j['qtl_prefix'],'--peqtl-smr','5e-8','--peqtl-heidi','1.57e-3','--ld-lower-limit','0.05','--ld-upper-limit','0.9','--heidi-min-m','3','--heidi-max-m','20','--cis-wind','2000','--maf','0.01','--thread-num','2','--out',str(out)]
if j.get('exclude_probe'):cmd+=['--exclude-probe',j['exclude_probe']];checks[j['exclude_probe']]=sha(j['exclude_probe'])
stamp=dict(command=cmd,input_checksums={k:v for k,v in checks.items() if k!=j['pairing_matrix']},version=actual,approved_trait=j['trait'],qtl_resource_id=j['qtl_resource_id'])
st=Path(str(out)+'.provenance.json');done=Path(str(out)+'.done');exitf=Path(str(out)+'.exitcode')
if done.exists() and exitf.exists() and exitf.read_text().strip()=='0' and st.exists() and json.loads(st.read_text())==stamp:
    res=Path(str(out)+'.smr');assert res.is_file() and res.stat().st_size>0
    print('REUSE_EXACT',j['job_id']);sys.exit(0)
assert not done.exists(),'Existing completed job inputs differ; choose a new output prefix'
st.write_text(json.dumps(stamp,indent=2));Path(str(out)+'.software_versions.tsv').write_text('software\tversion\nSMR\t'+actual+'\n')
with open(str(out)+'.input_checksums.tsv','w') as f:
    f.write('file\tsha256\n');f.writelines(k+'\t'+v+'\n' for k,v in checks.items())
with open(str(out)+'.run.log','w') as log:
    result=subprocess.run(['/usr/bin/prlimit','--as=17179869184','--']+cmd,stdout=log,stderr=subprocess.STDOUT)
exitf.write_text(str(result.returncode)+'\n')
assert result.returncode==0,'SMR failed; see run log, preserve checkpoint'
res=Path(str(out)+'.smr');assert res.exists() and res.stat().st_size>0,'Missing result even though process exited zero'
header=res.open().readline().split()
assert 'p_SMR' in header and 'p_HEIDI' in header,'Invalid output schema'
done.write_text('SUCCESS\n')
