#!/usr/bin/env python3
"""Finish only the locked LDH+CPI atomic job family; do not claim full preregistered completion."""
from pathlib import Path
import csv,json,subprocess,sys,time
R=Path('/analysis/projects/JTM_recent_highIF_direct_clone/P2_LDH_shared_architecture/11_NatCommun_strict_clone_Pro');M=R/'08_smr_heidi_coloc'
manifest=M/'SMR_JOBS_LDH_CPI.tsv';rows=list(csv.DictReader(manifest.open(),delimiter='\t'))
assert len(rows)==16 and {x['trait'] for x in rows}=={'LDH','CPI'}
while True:
    states=[]
    for j in rows:
        f=Path(j['output_prefix']);e=Path(str(f)+'.exitcode');d=Path(str(f)+'.done')
        if e.exists() and e.read_text().strip()!='0':raise SystemExit('FAILED_JOB_PRESERVED '+j['job_id'])
        states.append(d.exists() and e.exists() and e.read_text().strip()=='0')
    if all(states):break
    if (R/'00_admin/PAUSE_REQUESTED').exists():raise SystemExit('PAUSED_NO_NEW_STAGE')
    disp=(M/'CPI_DISPATCH.log').read_text() if (M/'CPI_DISPATCH.log').exists() else ''
    if 'DISPATCH_COMPLETE' in disp:raise SystemExit('Dispatcher complete but output family incomplete: inspect markers')
    time.sleep(30)
out=M/'shared_LDH_CPI_provisional'
cmd=[sys.executable,str(R/'16_workflow/replan_regulation_shared_smr.py'),'--job-manifest',str(manifest),'--outdir',str(out),'--family-status','PROVISIONAL_AVAILABLE_TRAITS']
code=subprocess.run(cmd).returncode
(M/'LDH_CPI_SHARED_FAMILY.exitcode').write_text(str(code)+'\n')
if code:raise SystemExit(code)
subprocess.run(['/usr/bin/Rscript',str(R/'16_workflow/replan_regulation_validate_bh.R'),str(out)],check=True)
subprocess.run([sys.executable,str(R/'16_workflow/replan_regulation_result_report.py')],check=True)
(M/'LDH_CPI_SHARED_FAMILY.done').write_text('SUCCESS_AVAILABLE_LDH_CPI_ONLY_NOT_ALL_PRIMARY_AXES\n')
subprocess.run([sys.executable,str(R/'16_workflow/replan_regulation_status.py')],check=True)
print('AVAILABLE_FAMILY_SUCCESS: global BH over all LDH+CPI panels; full registered family remains incomplete',flush=True)
