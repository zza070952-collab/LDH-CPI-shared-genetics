#!/usr/bin/env python3
"""Prepare all eight Bryois SMR universes without GWAS/candidate selection.

Retain all probes having a cis-eQTL P<5e-8 and all their eQTL rows satisfying
the preregistered HEIDI inclusion P<1.57e-3. This loses no SNP eligible for
SMR/HEIDI, but is NOT a full-window coloc input. Full nominal source remains
read-only for targeted coloc extraction. BESD conversion uses official SMR
data-management only, with version separately recorded from inference.
"""
from pathlib import Path
import argparse,csv,gzip,hashlib,json,math,os,re,shutil,subprocess
from collections import defaultdict
from scipy.stats import norm

P2=Path('/analysis/projects/JTM_recent_highIF_direct_clone/P2_LDH_shared_architecture')
R=P2/'11_NatCommun_strict_clone_Pro';OUT=R/'08_smr_heidi_coloc/bryois_full_smr_universe'
BR=P2.parent/'P1_OA_csMR/01_data_raw/brain_eqtl_bryois192_zenodo7276971'
GTF=Path('/analysis/projects/sciatica_multomics_bakeoff/work/phase3_final/gsmap_resource/genome_annotation/gtf/gencode.v46lift37.basic.annotation.gtf')
SMR=str(OUT.parent/'software/smr_Linux')
CELLS=['Astrocytes','Endothelial.cells','Excitatory.neurons','Inhibitory.neurons','Microglia','Oligodendrocytes','OPCs...COPs','Pericytes']
p=argparse.ArgumentParser();p.add_argument('--prepare-only',action='store_true');p.add_argument('--cell',choices=CELLS);p.add_argument('--max-new-gb',type=float,default=2);a=p.parse_args()
OUT.mkdir(parents=True,exist_ok=True)
def safe_space():
    used=sum(x.stat().st_size for x in OUT.rglob('*') if x.is_file())
    if shutil.disk_usage(OUT).free<80*1024**3:raise RuntimeError('DISK_FLOOR_80_GIB')
    if used>a.max_new_gb*1024**3:raise RuntimeError('MODULE_DISK_BUDGET_EXCEEDED')
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()
safe_space()
annotation={};gchecksum=sha(GTF)
with GTF.open() as f:
    for line in f:
        if line.startswith('#'):continue
        z=line.rstrip().split('\t')
        if len(z)<9 or z[2]!='gene' or not z[0].removeprefix('chr').isdigit():continue
        gid=re.search(r'gene_id "([^"]+)"',z[8]);gn=re.search(r'gene_name "([^"]+)"',z[8])
        if not gid:continue
        ensg=gid.group(1).split('.')[0];chrom=int(z[0].removeprefix('chr'))
        bp=int(z[3]) if z[6]=='+' else int(z[4])
        annotation[ensg]=(chrom,bp,gn.group(1) if gn else ensg,z[6])
manifest=[]
for cell in ([a.cell] if a.cell else CELLS):
    safe_space();pref=OUT/cell;marker=Path(str(pref)+'.conversion.done')
    if marker.exists():
        print('PRESERVE_ALREADY_CONVERTED',cell,flush=True);continue
    retained=defaultdict(list);top=defaultdict(lambda:1.0);source_n=0;checks=[]
    for chrom in range(1,23):
        src=BR/f'{cell}.{chrom}.gz';checks.append(dict(file=str(src),sha256=sha(src)))
        with gzip.open(src,'rt') as f:
            for line in f:
                x=line.split();source_n+=1
                if len(x)<5:continue
                try:pv=float(x[3]);b=float(x[4])
                except ValueError:continue
                if not(0<pv<=1 and math.isfinite(b)):continue
                if pv<top[x[0]]:top[x[0]]=pv
                if pv<1.57e-3 and b!=0:retained[x[0]].append((x[1],b,pv))
        print(cell,'chr',chrom,'source_rows',source_n,'candidate_heidi_rows',sum(map(len,retained.values())),flush=True)
    eligible={g for g,v in top.items() if v<5e-8}
    needed={s for g in eligible for s,b,pv in retained[g]}
    metadata={}
    with gzip.open(BR/'snp_pos.txt.gz','rt') as f:
        for m in csv.DictReader(f,delimiter='\t'):
            if m['SNP'] in needed:metadata[m['SNP']]=m
    qfile=Path(str(pref)+'.qfile');nrows=0;qc=defaultdict(int);gene_qc=[]
    with qfile.open('w') as f:
        f.write('SNP\tChr\tBP\tA1\tA2\tFreq\tProbe\tProbe_Chr\tProbe_bp\tGene\tOrientation\tb\tse\tp\n')
        for gene in sorted(eligible):
            ensg='ENSG'+gene.split('_ENSG')[-1] if '_ENSG' in gene else gene
            ann=annotation.get(ensg)
            if ann is None:qc['missing_hg19_probe_annotation']+=1;continue
            if ann[0]==6 and 25e6<=ann[1]<=35e6:qc['MHC_probe_excluded']+=1;continue
            written=0;seen=set()
            for snp,b,pv in retained[gene]:
                m=metadata.get(snp)
                if m is None:qc['missing_snp_metadata']+=1;continue
                if snp in seen:qc['duplicate_snp_probe']+=1;continue
                try:cc,bb=m['SNP_id_hg19'].split(':');cc=int(cc.removeprefix('chr'));bb=int(bb)
                except ValueError:qc['invalid_hg19_snp']+=1;continue
                a1,a2=m['effect_allele'],m['other_allele']
                if a1 not in ('A','C','G','T') or a2 not in ('A','C','G','T') or a1==a2:qc['invalid_allele']+=1;continue
                if cc!=ann[0] or abs(bb-ann[1])>2e6:qc['outside_prespecified_cis_window']+=1;continue
                se=abs(b)/norm.isf(pv/2)
                if not(math.isfinite(se) and se>0):qc['invalid_SE']+=1;continue
                # Official SNP metadata reports MAF only, not EAF. NA is intentional.
                f.write('\t'.join(map(str,[snp,cc,bb,a1,a2,'NA',gene,ann[0],ann[1],ann[2],ann[3],b,se,pv]))+'\n')
                seen.add(snp);nrows+=1;written+=1
            gene_qc.append(dict(panel=cell,probe=gene,top_P_source=top[gene],included_rows=written,gene_annotation='GENCODEv46lift37',probe_hg19=ann[1]))
    safe_space()
    checks.append(dict(file=str(GTF),sha256=gchecksum));checks.append(dict(file=str(BR/'snp_pos.txt.gz'),sha256=sha(BR/'snp_pos.txt.gz')))
    for suffix,rows in [('input_checksums.tsv',checks),('probe_qc.tsv',gene_qc)]:
        with open(str(pref)+'.'+suffix,'w') as f:
            ww=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');ww.writeheader();ww.writerows(rows)
    Path(str(pref)+'.qc.json').write_text(json.dumps(dict(panel=cell,source_rows=source_n,source_probes=len(top),eligible_probes=len(eligible),rows=nrows,qc=dict(qc),frequency='NA_NOT_INFERRED_FROM_MAF',coverage='ALL_SMR_AND_HEIDI_ELIGIBLE_SNPS; NOT_COLOC_FULL_CIS',N_study=192),indent=2))
    if not a.prepare_only:
        cmd=[SMR,'--qfile',str(qfile),'--make-besd','--out',str(pref)]
        with open(str(pref)+'.conversion.log','w') as log:rr=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT)
        Path(str(pref)+'.conversion.exitcode').write_text(str(rr.returncode)+'\n')
        if rr.returncode!=0:raise RuntimeError('BESD data conversion failed; source preserved')
        assert Path(str(pref)+'.besd').stat().st_size>0
        Path(str(pref)+'.software_versions.tsv').write_text('step\tsoftware\tversion\nBESD_data_format_conversion_only\tSMR\t1.03\nSMR_inference\tSMR\tNOT_RUN_BY_CONVERSION_SCRIPT\n')
        marker.write_text('SUCCESS_DATA_FORMAT_ONLY\n')
    manifest.append(dict(panel=cell,qfile=str(qfile),qfile_bytes=qfile.stat().st_size,full_heidi_coverage='YES',full_probe_universe='YES_AFTER_DOCUMENTED_ANNOTATION_QC',coloc_input='NO_USE_ORIGINAL_FULL_NOMINAL'))
    print('PREPARED',cell,nrows,flush=True)
if manifest:
    with (OUT/'BRYOIS_SMR_UNIVERSE_BUILD_MANIFEST.tsv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(manifest[0]),delimiter='\t');w.writeheader();w.writerows(manifest)
