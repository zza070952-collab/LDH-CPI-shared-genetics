#!/usr/bin/env python3
"""Mother-code LDSC h2 and h2-cts with checkpointed true outputs."""
import csv,gzip,hashlib,json,math,os,subprocess,sys,time
from pathlib import Path
r=Path(sys.argv[1]);od=r/'03_sldsc_scavenge/sldsc_salo';ld=r/'03_sldsc_scavenge/ldscores';ss=r/'03_sldsc_scavenge/salo_hm3/SALO2024_LDH_META.strict_hm3.sumstats.gz'
base='/analysis/projects/hipOA_spatial_multiomics_clean/resources/active/extracted/'
py='/analysis/miniconda3/envs/ldsc/bin/python';code='/analysis/projects/sciatica_multomics_bakeoff/work/software/src/ldsc/ldsc.py'
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
assert (ld/'ldscores.done').exists() and (r/'03_sldsc_scavenge/salo_hm3/salo_hm3.done').exists()
for c in range(1,23):
    cp=ld/f'chr{c:02d}';assert (cp/'ldscore.done').exists();q=json.loads((cp/'qc.json').read_text())
    for p,h in q['output_sha256'].items():assert sha(p)==h
    for suffix in ['annot.gz','l2.ldscore.gz','l2.M','l2.M_5_50']:assert (ld/f'Adult.{c}.{suffix}').is_file()
prints=od/'hm3_non_MHC.list.txt'
with prints.open('w') as out:
    for c in range(1,23):out.write((ld/f'chr{c:02d}/print_snps.txt').read_text())
baseline=str(od/'baselineLD_hm3_nonMHC/baselineLD.');weights=base+'ldsc_phase3_weights/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.';freq=base+'ldsc_phase3_frq/1000G_Phase3_frq/1000G.EUR.QC.'
assert (od/'S_LDSC_JOIN_REPAIR.json').exists()
lock={'trait':'SALO2024_LDH_META Track A','N':829699,'inference':'observed-scale LDH SNP heritability and genetic annotation enrichment; no liability transform','benchmark':'STRICT_CLONE_BENCHMARK_CORTEX','donors':2,'source':'author Adult ages25,40, PMID37824647','LDSC_version':'1.0.1','primary_family':'seven author-lineage h2-cts one-sided positive-coefficient tests; BH across all seven','sensitivity_family':'seven custom h2 enrichment P values BH across all seven; all 97 baseline estimates also retained','global_family':'seven CTS + seven h2 custom-enrichment P values; BH across 14; explicitly reported separately','negative_results':'retain all coefficients, SE, P, FDR, negative tau and lack of enrichment','source_GWAS_sha256':sha(ss),'LDSC_code_sha256':sha(code),'software_interpreter':py,'baseline':baseline,'weights':weights,'frequencies':freq,'threads':1,'memory_GiB':8,'no_force_inversion':True}
(od/'SLDSC_ANALYSIS_LOCK.json').write_text(json.dumps(lock,indent=2)+'\n')
common=['--w-ld-chr',weights,'--frqfile-chr',freq,'--print-snps',str(prints)]
jobs=[('h2',[py,code,'--h2',str(ss),'--ref-ld-chr',str(ld/'Adult.')+','+baseline]+common+['--overlap-annot','--print-coefficients','--out',str(od/'SALO_Adult_joint')],od/'SALO_Adult_joint.results'),('cts',[py,code,'--h2-cts',str(ss),'--ref-ld-chr',baseline,'--ref-ld-chr-cts',str(ld/'snATAC_Adult_final.ldcts')]+common+['--out',str(od/'SALO_Adult_CTS')],od/'SALO_Adult_CTS.cell_type_results.txt')]
for name,cmd,out in jobs:
    if (r/'00_admin/PAUSE_REQUESTED').exists():raise SystemExit(75)
    done=od/(name+'.done');qcp=od/(name+'.qc.json')
    if done.exists():q=json.loads(qcp.read_text());assert out.exists() and sha(out)==q['output_sha256'];print('SKIP VERIFIED',name,flush=True);continue
    previous_exit=od/(name+'.exitcode')
    if previous_exit.exists():previous_exit.rename(od/(name+'.exitcode.previous.'+time.strftime('%Y%m%dT%H%M%S')))
    (od/(name+'.command.json')).write_text(json.dumps(cmd,indent=2)+'\n');print('START',name,time.strftime('%FT%T'),flush=True)
    with (od/(name+'.stdout.log')).open('a') as log:res=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT)
    (od/(name+'.exitcode')).write_text(str(res.returncode)+'\n')
    if res.returncode:raise RuntimeError('LDSC nonzero exit: '+name)
    logfile=Path(str(out).replace('.results','.log').replace('.cell_type_results.txt','.log'));text=logfile.read_text()
    if 'Traceback' in text or 'Total time elapsed' not in text:raise RuntimeError('LDSC incomplete/exception log: '+name)
    assert out.is_file() and out.stat().st_size>50
    with out.open() as f:header=next(f).strip().split();rows=[line.strip().split() for line in f if line.strip()]
    expect=104 if name=='h2' else 7;assert len(rows)==expect,(name,len(rows),expect)
    expected_universal_annotation_NA=0
    for row in rows:
        assert len(row)==len(header),(name,row)
        for j,v in enumerate(row[1:],1):
            # Official LDSC makes Enrichment_p undefined for the universal baseline
            # (M_annotation == M_total, identically zero enrichment contrast SE).
            if name=='h2' and header[j]=='Enrichment_p' and v=='NA' and float(row[header.index('Prop._SNPs')])==1.0:
                expected_universal_annotation_NA+=1;continue
            assert math.isfinite(float(v)),(name,'unexpected nonfinite',header[j],v)
    q={'status':'SUCCESS','rows':len(rows),'columns':header,'output':str(out),'output_sha256':sha(out),'log_sha256':sha(logfile),'exitcode':0,'all_estimable_numeric_values_finite':True,'expected_universal_annotation_undefined_enrichment_P':expected_universal_annotation_NA}
    qcp.write_text(json.dumps(q,indent=2)+'\n');done.write_text(time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())+'\n');print('SUCCESS',name,flush=True)
