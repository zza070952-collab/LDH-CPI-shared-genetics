#!/usr/bin/env python3
"""Isolated LOW-overlap CPI run; official statistics, checkpoint-only additions."""
import difflib,json,os,shutil,subprocess,sys,urllib.request
from pathlib import Path
from replan_run_task import digest
ROOT=Path(__file__).resolve().parents[1];P2=ROOT.parent;W=ROOT/'16_workflow'
OUT=ROOT/'05_conjfdr/LDH_CPI_low_overlap';OUT.mkdir(parents=True,exist_ok=True)
OLD=P2/'10_workflow/vendor/pleiofdr-0da963cac22fe8de9030166d7aea974bcbb0a367';V=OUT/'vendor';AUDIT=OUT/'software_audit';AUDIT.mkdir(exist_ok=True)
OCT=P2/'envs/conjfdr_octave/bin/octave';PKG=P2/'envs/conjfdr_octave/octave_packages'
REF=P2/'01_data_raw/conjfdr_reference/ref9545380_1kgPhase3eur_LDr2p1.mat'
fingerprint={str(p):digest(p) for p in [Path(__file__),W/'replan_conjfdr_prepare.py',W/'replan_conjfdr_checkpoint.m',W/'replan_conjfdr_restore.m',REF,ROOT/'02_qc/gwas/FINNGEN_R12_LDH.GRCh37.QC.tsv.gz',ROOT/'02_qc/gwas/UKB_CPI_WANG2025.GRCh37.QC.tsv.gz']}
lock=OUT/'RUN_INPUT_LOCK.json'
if lock.exists():assert json.loads(lock.read_text())==fingerprint,'Checkpoint inputs/software changed; refuse reuse'
else:lock.write_text(json.dumps(fingerprint,indent=2))
if not V.exists():shutil.copytree(OLD,V)
# Preserve and compare inherited compatibility code to the exact upstream commit.
for name in ['FastPrune.m','random_prune_idx_amd_fb.m','pleiotropy_analysis.m']:
 original=AUDIT/(name+'.upstream')
 if not original.exists():
  with urllib.request.urlopen('https://raw.githubusercontent.com/precimed/pleiofdr/0da963cac22fe8de9030166d7aea974bcbb0a367/'+name,timeout=90) as r:body=r.read()
  assert b'<html' not in body.lower() and len(body)>100;original.write_bytes(body)
 inherited=(OLD/name).read_text();(AUDIT/(name+'.inherited')).write_text(inherited)
 (AUDIT/(name+'.inherited.diff')).write_text(''.join(difflib.unified_diff(original.read_text().splitlines(True),inherited.splitlines(True),fromfile='upstream/'+name,tofile='inherited/'+name)))
source=(OLD/'random_prune_idx_amd_fb.m').read_text();restore=(W/'replan_conjfdr_restore.m').read_text();checkpoint=(W/'replan_conjfdr_checkpoint.m').read_text()
assert source.count('for k=1:niter')==1 and source.count('prunemask(:, k) = isfinite(tmplogp_p);')==1
patched=source.replace('for k=1:niter',restore+'\nfor k=first_iteration:niter').replace('prunemask(:, k) = isfinite(tmplogp_p);','prunemask(:, k) = isfinite(tmplogp_p);\n'+checkpoint)
(V/'random_prune_idx_amd_fb.m').write_text(patched)
analysis=(OLD/'pleiotropy_analysis.m').read_text();assert analysis.count('%% PLOT QQ AND ENRICHMENT')==1
analysis=analysis.split('%% PLOT QQ AND ENRICHMENT')[0]+"\nsave('-v7',fullfile(outputdir,'result_core.mat.tmp'),'fdrmat','fdrmat12','fdrmat21','fdrvec0','logpvec1','logpmat2','zvec1','zmat2','excludevec','mafvec','chrnumvec','posvec','lookup12','lookup21');\nmovefile(fullfile(outputdir,'result_core.mat.tmp'),fullfile(outputdir,'result_core.mat'));\ndiary off;\n"
(V/'pleiotropy_analysis.m').write_text(analysis)
for name in ['random_prune_idx_amd_fb.m','pleiotropy_analysis.m']:
 (AUDIT/(name+'.checkpoint.diff')).write_text(''.join(difflib.unified_diff((OLD/name).read_text().splitlines(True),(V/name).read_text().splitlines(True),fromfile='inherited/'+name,tofile='checkpointed/'+name)))
(AUDIT/'CHECKSUMS.json').write_text(json.dumps({str(p.relative_to(V)):digest(p) for p in sorted(V.glob('*.m'))},indent=2))
(AUDIT/'IMPLEMENTATION.md').write_text('''# LOW-overlap LDH–UKB CPI conjFDR implementation
Official pleioFDR commit 0da963c; inherited Octave sparse-CSC compatibility is explicitly diffed, not an algorithm replacement. Core cond_FDR_amd, pleioFDR_amd, GCcorrect_logpvec and lookup code are unmodified.
The new pair uses 500 NEW prunings, seed9527; old 500-result files are never read or written. Every50 columns are compressed and persisted with exact RNG state. Resume restores contiguous validated input-locked chunks. A numerical toy regression test requires identical masks and final RNG for uninterrupted and resumed runs (test data never become biological results).
MAF>.01 and MHC6:25–34Mb excluded by normalized input, following new user QC. Mother fit exclusions6:25119106–33854733 and8:7.2–12.5Mb retained; exclude_from_discovery=false as mother. Ambiguous SNPs excluded; intergenic random-pruning GC, not standardGC; no sample-overlap correction for this LOW pair. Both conditional directions and conjunction saved before plot/locus reporting. Plot generation and true LD genomic-locus merging remain separate incomplete modules; SNP count is not locus count.
Official Python converter aligns alleles; adapter first requires exact GRCh37 coordinate and removes ALL duplicate rsIDs, preserves P and CPI per-SNP N / LDH effectiveN. P<=0 is excluded and counted rather than assigned a made-up floor.
''')
prunes=OUT/'prune_checkpoints';prunes.mkdir(exist_ok=True)
toy=OUT/'checkpoint_test';toy.mkdir(exist_ok=True)
octhead=f"pkg prefix '{PKG}'; pkg local_list '{PKG}/octave_packages'; pkg load statistics; pkg load nan;"
testcode=octhead+f"addpath('{V}'); setenv('P2_PRUNE_CHECKPOINT_DIR',''); rng(9527); refmask=random_prune_idx_amd_fb(100,sparse(eye(20)),true(20,1),'default'); refstate=rng(); setenv('P2_PRUNE_CHECKPOINT_DIR','{toy}'); rng(9527); first=random_prune_idx_amd_fb(50,sparse(eye(20)),true(20,1),'default'); clear first; resumed=random_prune_idx_amd_fb(100,sparse(eye(20)),true(20,1),'default'); assert(isequal(refmask,resumed)); assert(isequal(refstate,rng())); fprintf('CHECKPOINT_EXACT_RNG_MASK_PASS\\n');"
testdone=AUDIT/'checkpoint_test.done'
if not testdone.exists():
 with (AUDIT/'checkpoint_test.log').open('a') as f:rc=subprocess.run([str(OCT),'--no-gui','--quiet','--eval',testcode],stdout=f,stderr=subprocess.STDOUT).returncode
 assert rc==0,'Checkpoint numerical test failed';testdone.write_text('PASS exact mask and RNG, synthetic software test only\n')
subprocess.run([sys.executable,str(W/'replan_conjfdr_prepare.py')],check=True)
config=f'''reffile={REF}
traitfolder={OUT}
traitfile1=LDH.mat
traitname1=LDH
traitfiles={{'CPI.mat'}}
traitnames={{'CPI'}}
mlibrary={V}
outputdir={OUT}/results
stattype=conjfdr
fdrthresh=0.05
randprune=true
randprune_n=500
exclude_chr_pos=[6 25119106 33854733; 8 7200000 12500000]
exclude_from_discovery=false
mafthresh=0.01
exclude_ambiguous_snps=true
onscreen=false
manh_plot=false
refinfo=
reset_pruneidx=true
randprune_repeats=default
pthresh=1
perform_gc=true
use_standard_gc=false
randprune_gc=true
dummy_zscore=false
exit_matlab_upon_completion=false
'''
(V/'config.txt').write_text(config);(OUT/'results').mkdir(exist_ok=True)
result=OUT/'results/result_core.mat'
if not result.exists():
 env=os.environ.copy();env['P2_PRUNE_CHECKPOINT_DIR']=str(prunes)
 cmd=[str(OCT),'--no-gui','--quiet','--eval',octhead+"rng(9527); run('runme.m');"]
 (OUT/'octave_command.json').write_text(json.dumps(cmd,indent=2))
 with (OUT/'octave.log').open('a') as log:rc=subprocess.run(cmd,cwd=V,env=env,stdout=log,stderr=subprocess.STDOUT).returncode
 (OUT/'octave.exitcode').write_text(str(rc)+'\n');assert rc==0,'Octave failed; preserve pruning checkpoints'
assert result.stat().st_size>1000000 and (OUT/'octave.exitcode').read_text().strip()=='0'
from scipy.io import loadmat
import numpy as np
m=loadmat(result,variable_names=['fdrmat','fdrmat12','fdrmat21']);a,b,c=(m[k].ravel() for k in ['fdrmat','fdrmat12','fdrmat21']);valid=np.isfinite(a)&np.isfinite(b)&np.isfinite(c)
assert len(a)==9545380 and valid.sum()>1000000
assert np.allclose(a[valid],np.maximum(b[valid],c[valid]),rtol=0,atol=1e-12)
assert np.all((a[valid]>=0)&(a[valid]<=1))
(OUT/'CORE_QC.json').write_text(json.dumps({'status':'SUCCESS','finite_both_cond_and_conj':int(valid.sum()),'conj_lt_0_05_SNPs_NOT_LOCI':int(np.sum(valid&(a<.05))),'conj_lt_0_01_SNPs_NOT_LOCI':int(np.sum(valid&(a<.01))),'result_sha256':digest(result),'locus_status':'NOT_YET_LD_MERGED','new_pair_500_prunings':True},indent=2))
