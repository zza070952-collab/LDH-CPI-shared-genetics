"""New sciatica fits with exact, input-verified reuse of shared LDH and SNP-selection stages."""
import difflib,json
from pathlib import Path
from replan_run_task import digest
R=Path(__file__).resolve().parents[1];W=R/'16_workflow';base=W/'replan_mixer_cpi_chain.py'
BASESHA='3da6a414bd622cdbf375516e49b3a3290b0ca68c56d2b6206b6bc55d3d4f10d4';assert digest(base)==BASESHA
old=R/'04_ldsc_mixer/LDH_CPI_protocol13';out=R/'04_ldsc_mixer/LDH_SCIATICA_protocol13';out.mkdir(exist_ok=True)
assert (R/'logs/mixer_cpi_models/.done').exists(),'Wait for existing CPI chain; never repeat LDH model work'
q=json.loads((old/'LDH.input.done.json').read_text());assert q['status']=='SUCCESS' and q['signature']['code']==BASESHA and q['signature']['ncol']=='NEFF'
assert q['signature']['gwas']==digest(R/'02_qc/gwas/FINNGEN_R12_LDH.GRCh37.QC.tsv.gz') and q['sha256']==digest(old/'LDH.full_nonMHC.tsv.gz')
p=out/'LDH.full_nonMHC.tsv.gz'
if p.is_symlink():assert p.resolve()==old/p.name
else:assert not p.exists();p.symlink_to(old/p.name)
(out/'LDH_FULL_INPUT_REUSE.json').write_text(json.dumps(q,indent=2)+'\n')
src=base.read_text();needle="[('LDH','FINNGEN_R12_LDH','NEFF'),('CPI','UKB_CPI_WANG2025','N')]";assert src.count(needle)==1
new=src.replace(needle,"[('SCIATICA','PANUKB_SCIATICA_1476_EUR','NEFF')]").replace('LDH_CPI_protocol13','LDH_SCIATICA_protocol13').replace("'CPI'","'SCIATICA'").replace('LDH_CPI.','LDH_SCIATICA.').replace('LDH_CPI_','LDH_SCIATICA_')
hook=" if checkpoint.exists():\n"
assert new.count(hook)==1
reuse=''' if not checkpoint.exists() and (name.startswith('snps_') or name.startswith('LDH_fit1_') or name.startswith('LDH_test1_')):
  oldout=ROOT/'04_ldsc_mixer/LDH_CPI_protocol13';origin=oldout/(name+'.stage.json')
  assert origin.exists(),'Missing existing LDH/SNP stage; no automatic repeat'
  prior=json.loads(origin.read_text());ps=prior['signature']
  def canon(v):return str(v).replace(str(OUT),str(oldout))
  assert ps['cmd']==[canon(v) for v in cmd] and ps['LD22']==records and ps['library']==digest(LIB)
  assert ps['controller']=='3da6a414bd622cdbf375516e49b3a3290b0ca68c56d2b6206b6bc55d3d4f10d4'
  assert ps['inputs']=={canon(p):digest(p) for p in inputs}
  if prior['status']=='SUCCESS':
   oldfile=Path(canon(output));assert digest(oldfile)==prior['output_sha256']
   if output.is_symlink():assert output.resolve()==oldfile
   else:assert not output.exists();output.symlink_to(oldfile)
  row=dict(prior);row['signature']=sig;row['reused_from']=str(origin);row['reuse_reason']='Identical LD22, library, full LDH GWAS, effective N, extract, seed, parameters and upstream checksums; only output paths differ'
  checkpoint.write_text(json.dumps(row,indent=2));stages.append(row);return row['status']=='SUCCESS'
'''
new=new.replace(hook,reuse+hook)
(out/'PAIR_AND_REUSE_ONLY.diff').write_text(''.join(difflib.unified_diff(src.splitlines(True),new.splitlines(True),fromfile=str(base),tofile='SCIATICA_chain')))
exec(compile(new,str(base),'exec'),{'__file__':__file__,'__name__':'__main__'})
