"""Final three-seed MiXeR acceptance audit; descriptive boundary diagnostics, no new veto."""
import csv,hashlib,json,math,statistics
from datetime import datetime,timezone
from pathlib import Path
R=Path(__file__).resolve().parents[1];O=R/'04_ldsc_mixer'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def opt(d):
 x=d.get('optimize',[])
 if isinstance(x,list) and x and isinstance(x[-1],list):return x[-1][1]
 if isinstance(x,dict):return x
 return {}
rows=[]
for pair,dn,stem in [('CPI','LDH_CPI_protocol13','LDH_CPI'),('SCIATICA','LDH_SCIATICA_protocol13','LDH_SCIATICA')]:
 d=O/dn
 qc=json.loads((d/'ALL21_NUMERIC_QC.done.json').read_text());assert qc['status']=='SUCCESS' and qc['module_status']['success']==21 and qc['module_status']['failed']==0
 for seed in [9527,20260904,314159]:
  stage=[]
  for s in [f'LDH_fit1_{seed}',f'LDH_test1_{seed}',f'{pair}_fit1_{seed}',f'{pair}_test1_{seed}',f'{stem}_fit2_{seed}',f'{stem}_test2_{seed}']:
   q=json.loads((d/(s+'.stage.json')).read_text());assert q['status']=='SUCCESS' and q['actual_exitcode']==0;stage.append(sha(d/(s+'.stage.json')))
  fp=d/f'{stem}.fit.seed{seed}.json';tp=d/f'{stem}.test.seed{seed}.json';f=json.loads(fp.read_text());t=json.loads(tp.read_text());assert f['analysis']=='bivariate' and t['analysis']=='bivariate'
  c={k:v['point_estimate'] for k,v in f['ci'].items()};fo=opt(f);to=opt(t)
  assert fo.get('success') is True and to.get('success') is True
  rows.append({'pair':pair,'seed':seed,'nc_LDH_total':c['nc1u'],'nc_pain_total':c['nc2u'],'shared_nc12':c['nc12'],'LDH_unique_nc1':c['nc1'],'pain_unique_nc2':c['nc2'],'Dice':c['dice'],'rg_model':c['rg'],'rho_beta':c['rho_beta'],'shared_fraction_LDH':c['pi12_over_pi1u'],'shared_fraction_pain':c['pi12_over_pi2u'],'distance_from_complete_LDH_sharing':1-c['pi12_over_pi1u'],'fit_AIC':fo.get('AIC'),'test_AIC':to.get('AIC'),'ci_available':any(v.get('lower') is not None or v.get('upper') is not None for v in f['ci'].values()),'fit_optimizer_success':fo.get('success'),'test_optimizer_success':to.get('success'),'fit_sha256':sha(fp),'test_sha256':sha(tp),'stage_marker_sha256s':';'.join(stage),'status':'SUCCESS_NUMERICAL_BOUNDARY_REVIEWED'})
cols=list(rows[0]);p=O/'MIXER_FINAL_THREE_SEED_RESULTS_ZERO.tsv'
with p.open('w',newline='') as f:w=csv.DictWriter(f,cols,delimiter='\t');w.writeheader();w.writerows(rows)
def stats(v):
 m=statistics.mean(v);return {'min':min(v),'median':statistics.median(v),'max':max(v),'mean':m,'sample_sd':statistics.stdev(v),'CV_sample':statistics.stdev(v)/abs(m) if m else None,'relative_range_over_median':(max(v)-min(v))/abs(statistics.median(v)) if statistics.median(v) else None}
summary=[]
for pair in ['CPI','SCIATICA']:
 a=[r for r in rows if r['pair']==pair];s={k:stats([float(r[k]) for r in a]) for k in ['shared_nc12','Dice','rg_model','rho_beta','shared_fraction_LDH','shared_fraction_pain','nc_pain_total']}
 if pair=='CPI':
  decision='BOUNDARY_LIMITED_DIRECTIONALLY_CONSISTENT_NOT_STABLE_UNCONSTRAINED_ESTIMATE'
  interpretation='All seeds place approximately all modeled LDH variants in the shared component, with shared nc12 11243.85-12799.95 (median12433.88). This is a reproducible boundary pattern, not a well-identified interior shared-count estimate. Dice and pain polygenicity vary materially; no confidence intervals were generated.'
 else:
  decision='BOUNDARY_LIMITED_AND_QUANTITATIVELY_UNSTABLE_LOW_POWER_AXIS'
  interpretation='All seeds converge numerically but pain polygenicity, Dice and model rg vary widely; two fits are essentially complete-LDH-sharing boundary solutions and the third remains near boundary. Quantitative shared-polygenicity estimates are not reliable for primary inference.'
 summary.append({'pair':pair,'decision':decision,'interpretation':interpretation,'statistics':s,'all_21_stages_success':True,'three_fit_and_test_success':True,'confidence_intervals_available':False})
out={'status':'FINAL_THREE_SEED_BOUNDARY_REVIEW_COMPLETE','generated_utc':datetime.now(timezone.utc).isoformat(),'old_shared_11258_reused':False,'old_shared_11258_verdict':'NOT_A_VERIFIED_STABLE_FORMAL_RESULT','CPI_formal_report':'Report seed range/median and complete-sharing boundary; do not report a single exact shared count as stable. Boundary is not used as an automatic Hard-Go veto.','SCIATICA_formal_report':'Use as low-power qualitative/sensitivity support only; do not report a stable shared count or Dice.','aic_note':'fit/test sampling cost_n differs by seed; raw AIC is retained per run but not used to select a visually favorable seed. No best-seed selection.','ci_note':'Pinned runs produced point estimates with no lower/upper intervals; none are invented.','summary':summary,'table':str(p),'table_sha256':sha(p)}
(O/'MIXER_FINAL_STABILITY_ZERO.json').write_text(json.dumps(out,indent=2))
lines=['# Final three-seed MiXeR stability audit','','All 42 scientific stages (21 per pair) passed exit-code, schema and numeric checks. Numerical completion is separated from inferential stability.','','| Pair | shared nc12 range; median | Dice range; median | model rg range | shared fraction of LDH | Decision |','| --- | --- | --- | --- | --- | --- |']
for z in summary:
 s=z['statistics'];lines.append(f"| {z['pair']} | {s['shared_nc12']['min']:.2f}-{s['shared_nc12']['max']:.2f}; {s['shared_nc12']['median']:.2f} | {s['Dice']['min']:.4f}-{s['Dice']['max']:.4f}; {s['Dice']['median']:.4f} | {s['rg_model']['min']:.4f}-{s['rg_model']['max']:.4f} | {s['shared_fraction_LDH']['min']:.8f}-{s['shared_fraction_LDH']['max']:.8f} | {z['decision']} |")
lines+=['','CPI supports a directionally consistent near-complete modeled-LDH-sharing pattern, but its exact shared count is boundary-limited. The old approximately 11,258 value is not restored as a stable formal estimate. Sciatica is substantially less stable and is retained only as low-power qualitative/sensitivity evidence. No seed was selected by outcome attractiveness, and no boundary cutoff or Hard-Go veto was introduced.','',f'Source table SHA256: `{out["table_sha256"]}`']
(O/'MIXER_FINAL_STABILITY_ZERO.md').write_text('\n'.join(lines)+'\n')
print(json.dumps(out,indent=2))
