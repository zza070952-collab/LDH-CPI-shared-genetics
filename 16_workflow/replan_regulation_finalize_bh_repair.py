from pathlib import Path
import difflib,hashlib,json,subprocess
R=Path('/analysis/projects/JTM_recent_highIF_direct_clone/P2_LDH_shared_architecture/11_NatCommun_strict_clone_Pro');M=R/'08_smr_heidi_coloc';D=M/'shared_LDH_CPI_provisional';S=M/'audit_invalid_BH_assignment_20260908'
before=json.loads((S/'RAW_SMR_PRE_REPAIR_SHA256.json').read_text());checks=[]
for x in before:
    h=hashlib.sha256(Path(x['file']).read_bytes()).hexdigest();assert h==x['sha256'];checks.append(dict(x,after_sha256=h,status='UNCHANGED'))
(S/'RAW_SMR_AFTER_REPAIR_SHA256.json').write_text(json.dumps(checks,indent=2))
code=R/'16_workflow/replan_regulation_shared_smr.py'
(S/'BH_MINIMAL_REPAIR.diff').write_text(''.join(difflib.unified_diff((S/'shared_smr_BEFORE_PATCH.py').read_text().splitlines(True),code.read_text().splitlines(True),fromfile='before_shared_smr.py',tofile='after_shared_smr.py')))
with (S/'BOOLEAN_ASSIGNMENT_DIAGNOSTIC.log').open('w') as f:subprocess.run(['/analysis/projects/JTM_recent_highIF_direct_clone/P2_LDH_shared_architecture/envs/mixer_v2.2.1/bin/python',str(R/'16_workflow/replan_regulation_debug_bh.py')],stdout=f,stderr=subprocess.STDOUT,check=True)
assert (D/'BH_INDEPENDENT_R_VALIDATION.done').exists()
report='''# BH postprocessing incident and repair — 2026-09-08

## Original error and detection

The first combined LDH+CPI summary wrote incorrect row-specific global BH values. A biological sanity check on TGFA found raw P=2.071207e-6 but written global FDR=0.64291817. Independent R stats::p.adjust returned 0.007549549515. The whole first joint summary was suspended and retained in this directory as invalid; raw SMR outputs were never changed.

## Reproduction and scope

The standalone Python BH function matches R. In this installed runtime (pandas1.5.0, NumPy1.23.3), the retained Boolean .loc assignment path reproduces a maximum row discrepancy of0.9219016 even though input P arrays match exactly, all row labels are unique and exactly0..7289, columns are unique, and Boolean mask values are all True. Fresh mask construction and whole-column assignment return correct values. A four-row synthetic unit test also passes. This evidence identifies the problematic assignment path but does not establish an upstream pandas defect; the lower-level root cause remains unassigned. The self-contained reproducer, original data snapshot, script and diagnostic log are retained.

## Minimal repair and validation

The postprocessor now computes BH using explicit integer row positions and assigns the complete output columns. No raw P, HEIDI, family membership, threshold, algorithm or SMR job was changed. All7,290 global and within-trait BH values were independently checked against R stats::p.adjust at absolute tolerance1e-12; every primary-pass and shared-pass Boolean was also recomputed in R. All16 original SMR output SHA256 values are unchanged. Only the inexpensive summary was rerun.

## Corrected result

LDH:3,638 tests,3 global FDR<0.05 plus HEIDI>0.05 signals. CPI:3,652 tests,2 such signals. The erroneous snapshot had misassigned3 to CPI and2 to LDH. The total strict shared signal count remains0 after independent verification. TGFA OPC/COP: LDH FDR0.007549549515, CPI FDR0.219769551; not a shared pass. This is the completed available LDH+CPI family, not all prespecified pain axes or all GTEx resources.

No completed long-running SMR calculation was repeated. The repair does not loosen the preregistration or promote any candidate. Future family summaries must pass the independent R validation before use.
'''
(M/'SMR_BH_REPAIR_LOG.md').write_text(report)
(M/'LDH_CPI_SHARED_FAMILY.done').write_text('SUCCESS_AVAILABLE_LDH_CPI_ONLY_NOT_ALL_PRIMARY_AXES; INDEPENDENT_R_BH_VALIDATION_PASS\n')
print('REPAIR_CLOSED: 16 raw SHA unchanged; 7290 full-row R validation passed')
