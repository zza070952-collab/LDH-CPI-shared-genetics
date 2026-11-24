"""Resume completed pruning, with memory-bounded diagnostic column sums only."""
import difflib,json,os,subprocess,sys
from pathlib import Path
import numpy as np
from scipy.io import loadmat
from replan_run_task import digest
R=Path(__file__).resolve().parents[1];pair=sys.argv[1];assert pair in ['LDH_CPI_low_overlap','LDH_SCIATICA_PANUKB_low_overlap']
D=R/'05_conjfdr'/pair;V=D/'vendor';A=D/'software_audit';result=D/'results/result_core.mat'
assert not result.exists(),'Core exists: verify, do not repeat'
lock=json.loads((D/'RUN_INPUT_LOCK.json').read_text())
assert all(digest(p)==h for p,h in lock.items()),'Input fingerprint changed'
before={str(p):digest(p) for p in sorted((D/'prune_checkpoints').glob('prune_*.mat'))};assert len(before)==10
for p in before:
 x=loadmat(p,variable_names=['k']);assert int(x['k'].item())==int(Path(p).stem.split('_')[1])
for trait in ['LDH','CPI' if 'CPI' in pair else 'SCIATICA']:
 q=json.loads((D/(trait+'.input_qc.json')).read_text());assert q['status']=='SUCCESS' and digest(D/(trait+'.mat'))==q['mat_sha256']
p=V/'pleiotropy_analysis.m';backup=A/'pleiotropy_analysis.m.pre_zero_repair'
if not backup.exists():backup.write_bytes(p.read_bytes())
src=backup.read_text();old="fprintf('mean (std) variants per random pruning iteration = %.2f (%.2f)\\n', mean(sum(pruneidx)), ...\n    std(sum(pruneidx)))"
assert src.count(old)==1
new="prune_counts=zeros(1,size(pruneidx,2));\nfor prune_col=1:size(pruneidx,2), prune_counts(prune_col)=sum(pruneidx(:,prune_col)); end\nfprintf('mean (std) variants per random pruning iteration = %.2f (%.2f)\\n',mean(prune_counts),std(prune_counts));"
patched=src.replace(old,new);p.write_text(patched)
(A/'POST500_MEMORY_DIAGNOSTIC_ONLY.diff').write_text(''.join(difflib.unified_diff(src.splitlines(True),patched.splitlines(True))))
# Restore hook verifies all10 chunks and makes first_iteration501: no pruning body executes.
random_code=(V/'random_prune_idx_amd_fb.m').read_text();assert 'for k=first_iteration:niter' in random_code and 'Checkpoint gap' in random_code
octave=R.parent/'envs/conjfdr_octave/bin/octave';pkg=R.parent/'envs/conjfdr_octave/octave_packages'
head=f"pkg prefix '{pkg}'; pkg local_list '{pkg}/octave_packages'; pkg load statistics; pkg load nan;"
test="x=reshape(mod(1:2000,3)==0,100,20);y=zeros(1,20);for j=1:20,y(j)=sum(x(:,j));end;assert(isequal(y,sum(x)));assert(mean(y)==mean(sum(x)));assert(std(y)==std(sum(x)));"
env=os.environ.copy();env['P2_PRUNE_CHECKPOINT_DIR']=str(D/'prune_checkpoints')
env['OCTAVE_HOME']=str(R.parent/'envs/conjfdr_octave')
cmd=[str(octave),'--no-gui','--quiet','--eval',head+test+"rng(9527);run('runme.m');"]
(A/'POST500_RESUME_COMMAND.json').write_text(json.dumps({'cmd':cmd,'before_checkpoints':before,'change':'Diagnostic column sums only; no FDR/GC/statistical parameter changes; original64GiB limit increased to96GiB','pruning_iterations_to_execute':0},indent=2))
with (D/'post500_resume.log').open('a') as f:rc=subprocess.run(cmd,cwd=V,env=env,stdout=f,stderr=subprocess.STDOUT).returncode
(D/'post500_resume.exitcode').write_text(str(rc)+'\n')
assert all(digest(p)==h for p,h in before.items()),'Pruning checkpoint unexpectedly modified'
assert rc==0 and result.is_file()
m=loadmat(result,variable_names=['fdrmat','fdrmat12','fdrmat21']);a,b,c=(m[k].ravel() for k in ['fdrmat','fdrmat12','fdrmat21']);v=np.isfinite(a)&np.isfinite(b)&np.isfinite(c)
assert len(a)==9545380 and v.sum()>1000000 and np.allclose(a[v],np.maximum(b[v],c[v]),atol=1e-12,rtol=0) and np.all((a[v]>=0)&(a[v]<=1))
(D/'CORE_QC.json').write_text(json.dumps({'status':'SUCCESS','finite_both_cond_and_conj':int(v.sum()),'conj_lt_0_05_SNPs_NOT_LOCI':int(np.sum(v&(a<.05))),'conj_lt_0_01_SNPs_NOT_LOCI':int(np.sum(v&(a<.01))),'result_sha256':digest(result),'locus_status':'NOT_YET_LD_MERGED','pruning_checkpoint_hashes_unchanged':True,'pruning_rerun':False,'actual_postprocess_exit':rc},indent=2))
print('POST500_QC_PASS',pair)
