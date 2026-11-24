#!/usr/bin/env python3
"""Generate the pinned protocol-1.3 LD reference, checkpointed chromosome-wise."""
import json,os,shutil,subprocess,sys
from pathlib import Path
from replan_run_task import digest
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'04_ldsc_mixer/reference_protocol13';OUT.mkdir(parents=True,exist_ok=True)
V=ROOT/'16_workflow/vendor/mixer_protocol13_3d17e8a';LIB=V/'src/build/lib/libbgmg.so'
REF=Path('/analysis/projects/hipOA_spatial_multiomics_clean/resources/active/extracted/ldsc_phase3_plink/1000G_EUR_Phase3_plink')
chromosomes=[int(c) for c in sys.argv[1:]];assert chromosomes and all(1<=c<=22 for c in chromosomes)
for ch in chromosomes:
 if (ROOT/'00_admin/PAUSE_REQUESTED').exists():raise SystemExit(75)
 prefix=REF/f'1000G.EUR.QC.{ch}';out=OUT/f'1000G.EUR.QC.{ch}.run4.ld';done=OUT/f'chr{ch:02}.done.json'
 cmd=[sys.executable,str(V/'precimed/mixer.py'),'ld','--lib',str(LIB),'--bfile',str(prefix),'--out',str(out),'--r2min','.05','--ldscore-r2min','.05','--ld-window-kb','30000']
 fingerprint={'command':cmd,'source_commit':'3d17e8a72e91982b186bc069a0506ef35d873b9e','lib_sha256':digest(LIB),'inputs':{ext:digest(str(prefix)+'.'+ext) for ext in ['bed','bim','fam']}}
 if done.exists():
  old=json.loads(done.read_text());assert old['fingerprint']==fingerprint and digest(out)==old['output_sha256'];print('REUSE_EXACT_LD',ch,flush=True);continue
 if shutil.disk_usage(ROOT).free<80*1024**3+512*1024**2:raise SystemExit('DISK_RESERVE_NO_NEW_CHROMOSOME')
 (OUT/f'chr{ch:02}.command.json').write_text(json.dumps(cmd,indent=2))
 print('START_LD',ch,flush=True)
 with (OUT/f'chr{ch:02}.stdout.log').open('a') as log:rc=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT).returncode
 (OUT/f'chr{ch:02}.exitcode').write_text(str(rc)+'\n');assert rc==0,'LD command failed'
 assert out.is_file() and out.stat().st_size>1000000
 logpath=Path(str(out)+'.log');log=logpath.read_text();assert '<save_ld_matrix(' in log,'Missing successful serialization log'
 done.write_text(json.dumps({'status':'SUCCESS','fingerprint':fingerprint,'bytes':out.stat().st_size,'output_sha256':digest(out),'log_sha256':digest(logpath),'actual_exitcode':rc},indent=2))
 print('SUCCESS_LD',ch,out.stat().st_size,flush=True)
