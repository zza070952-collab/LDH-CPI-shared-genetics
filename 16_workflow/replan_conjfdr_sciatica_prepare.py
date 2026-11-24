"""Same audited full-GWAS converter; only new LOW sciatica input is converted."""
import difflib,json
from pathlib import Path
from replan_run_task import digest
R=Path(__file__).resolve().parents[1];W=R/'16_workflow';D=R/'05_conjfdr/LDH_SCIATICA_PANUKB_low_overlap';D.mkdir(exist_ok=True)
base=W/'replan_conjfdr_prepare.py';assert digest(base)=='b607c22bd831fe933406d9beda216c7bfefc3d3286d058df8add11c6b2c69514'
old=R/'05_conjfdr/LDH_CPI_low_overlap';q=json.loads((old/'LDH.input_qc.json').read_text())
assert q['status']=='SUCCESS' and q['signature']['adapter_sha256']==digest(base)
assert q['signature']['gwas_sha256']==digest(R/'02_qc/gwas/FINNGEN_R12_LDH.GRCh37.QC.tsv.gz')
assert q['signature']['N_definition']=='NEFF'
assert q['signature']['reference_sha256']==digest(R.parent/'01_data_raw/conjfdr_reference/9545380.ref')
assert q['signature']['converter_sha256']==digest(R.parent/'10_workflow/vendor/python_convert-e46ebdfafd495c1420c7f8a4740a0da75c94d84d/sumstats.py')
for name,key in [('LDH.mat','mat_sha256'),('LDH.tsv.gz','tsv_sha256')]:
 p=old/name;assert digest(p)==q[key];link=D/name
 if link.is_symlink():assert link.resolve()==p
 else:assert not link.exists();link.symlink_to(p)
(D/'LDH_EXACT_INPUT_REUSE.json').write_text(json.dumps({'source_qc':str(old/'LDH.input_qc.json'),'verified':q,'conversion_rerun':False,'scope':'Same LDH full GWAS, software, reference, MAF, MHC and effective N; consumer never edits these links.'},indent=2)+'\n')
qc_link=D/'LDH.input_qc.json'
if qc_link.is_symlink():assert qc_link.resolve()==old/'LDH.input_qc.json'
else:assert not qc_link.exists();qc_link.symlink_to(old/'LDH.input_qc.json')
src=base.read_text();needle="[('LDH','FINNGEN_R12_LDH','NEFF'),('CPI','UKB_CPI_WANG2025','N')]";assert src.count(needle)==1
new=src.replace('LDH_CPI_low_overlap','LDH_SCIATICA_PANUKB_low_overlap').replace(needle,"[('SCIATICA','PANUKB_SCIATICA_1476_EUR','NEFF')]")
(D/'PREPARE_SCOPE_ONLY.diff').write_text(''.join(difflib.unified_diff(src.splitlines(True),new.splitlines(True),fromfile=str(base),tofile='SCIATICA_adapter')))
exec(compile(new,str(base), 'exec'),{'__file__':__file__,'__name__':'__main__'})
