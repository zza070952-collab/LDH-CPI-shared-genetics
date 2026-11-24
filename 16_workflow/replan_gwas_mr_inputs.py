"""Independent exposure-only IV selection and 503-EUR PLINK clumping, both ways.

This prepares real inputs for the mother's TwoSampleMR algorithms. It does not
substitute custom MR estimators or select instruments using outcome significance.
"""
import collections,csv,gzip,hashlib,json,pathlib,subprocess
R=pathlib.Path('/analysis/projects/JTM_recent_highIF_direct_clone/P2_LDH_shared_architecture/11_NatCommun_strict_clone_Pro');O=R/'08_mr/LDH_CPI';O.mkdir(parents=True,exist_ok=True)
REF=pathlib.Path('/analysis/projects/hipOA_spatial_multiomics_clean/resources/active/extracted/magma_g1000_eur/g1000_eur');PLINK='/analysis/miniconda3/envs/pops_env/bin/plink'
SPECS={'LDH':{'id':'FINNGEN_R12_LDH','sha256':'aef6a58605765206dbfa30eac5b1998e1ac68ddd8da875bc3fe9d4bbfbafdbc4'},'CPI':{'id':'UKB_CPI_WANG2025','sha256':'1edb5aab9445960c4915dd393eb580e103e142fcf6cd07acee1d00a1585541f3'}}
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(8*1024**2),b''):h.update(b)
    return h.hexdigest()
def table(p,rows,fields):
    with p.open('w') as f:w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)
if (R/'00_admin/PAUSE_REQUESTED').exists():raise SystemExit('PAUSE_REQUESTED')
nref=sum(1 for _ in pathlib.Path(str(REF)+'.fam').open())
if nref!=503:raise RuntimeError('Reference is not exactly 503 samples')
allgws={};counts={};audit=[]
for trait,spec in SPECS.items():
    path=R/'02_qc/gwas'/(spec['id']+'.GRCh37.QC.tsv.gz');assert sha(path)==spec['sha256'],'GWAS checksum mismatch'
    spec['path']=path;d={};c=collections.Counter()
    with gzip.open(path,'rt') as f:
        for x in csv.DictReader(f,delimiter='\t'):
            c['QC_rows']+=1
            if float(x['P'])>=5e-8:continue
            c['GWS_rows']+=1;x['F_stat']=(float(x['BETA'])/float(x['SE']))**2
            if int(x['CHR'])==6 and 25000000<=int(x['BP'])<=34000000:c['GWS_MHC_excluded']+=1;continue
            if x['F_stat']<10:c['GWS_F_below10']+=1;continue
            if x['SNP'] in d:raise RuntimeError('Duplicate SNP in validated QC file')
            d[x['SNP']]=x
    allgws[trait]=d;counts[trait]=c
union=set().union(*(set(d) for d in allgws.values()));matches={};duplicates=set()
with pathlib.Path(str(REF)+'.bim').open() as f:
    for line in f:
        ch,snp,cm,bp,a1,a2=line.split()
        if snp not in union:continue
        if snp in matches:duplicates.add(snp)
        matches[snp]=(int(ch),int(bp),a1,a2)
leads={}
for trait,d in allgws.items():
    eligible={}
    for snp,x in d.items():
        reason='REFERENCE_ELIGIBLE'
        if snp not in matches:reason='NOT_IN_REFERENCE'
        elif snp in duplicates:reason='DUPLICATE_REFERENCE_RSID'
        else:
            ch,bp,a1,a2=matches[snp]
            if (ch,bp)!=(int(x['CHR']),int(x['BP'])):reason='REFERENCE_POSITION_MISMATCH'
            elif {a1,a2}!={x['A1'],x['A2']}:reason='REFERENCE_ALLELE_MISMATCH'
        counts[trait][reason]+=1
        audit.append({'exposure':trait,'SNP':snp,'CHR':x['CHR'],'BP':x['BP'],'P':x['P'],'F_stat':x['F_stat'],'reference_QC':reason,'clump_status':'NOT_ELIGIBLE','sample_overlap':'LOW'})
        if reason=='REFERENCE_ELIGIBLE':eligible[snp]=x
    prefix=O/(trait+'_independent_IVs')
    table(O/(trait+'_GWS_for_clump.tsv'),[{'SNP':s,'P':x['P']} for s,x in eligible.items()],['SNP','P'])
    (O/(trait+'_GWS_for_clump.snps')).write_text(''.join(s+'\n' for s in eligible))
    chosen=[]
    if eligible:
        cmd=[PLINK,'--bfile',str(REF),'--extract',str(O/(trait+'_GWS_for_clump.snps')),'--clump',str(O/(trait+'_GWS_for_clump.tsv')),'--clump-p1','5e-8','--clump-p2','5e-8','--clump-r2','0.001','--clump-kb','10000','--threads','2','--memory','8192','--out',str(prefix)]
        prefix.with_suffix('.command.json').write_text(json.dumps(cmd,indent=2));p=subprocess.run(cmd);prefix.with_suffix('.exitcode').write_text(str(p.returncode)+'\n')
        if p.returncode:raise RuntimeError('PLINK clump failed for '+trait)
        path=prefix.with_suffix('.clumped')
        if not path.exists():raise RuntimeError('Missing clumped output with nonzero eligible IVs')
        with path.open() as f:
            header=None
            for line in f:
                bits=line.split()
                if not bits:continue
                if header is None:header=bits;continue
                row=dict(zip(header,bits))
                if row['SNP'] not in eligible:raise RuntimeError('Unrecognized clump lead')
                chosen.append(row['SNP'])
    leads[trait]=chosen;counts[trait]['independent_IVs']=len(chosen)
    for x in audit:
        if x['exposure']==trait and x['reference_QC']=='REFERENCE_ELIGIBLE':x['clump_status']='INDEPENDENT_LEAD' if x['SNP'] in chosen else 'LD_PRUNED'
    table(O/(trait+'_INDEPENDENT_IVS.tsv'),[eligible[s] for s in chosen],list(next(iter(d.values())).keys()) if d else ['SNP','CHR','BP','A1','A2','BETA','SE','P','EAF','N','NEFF','INFO','F_stat'])
lookup={}
for trait,spec in SPECS.items():
    other='CPI' if trait=='LDH' else 'LDH';needed=set(leads[other]);found={}
    with gzip.open(spec['path'],'rt') as f:
        for x in csv.DictReader(f,delimiter='\t'):
            if x['SNP'] in needed:found[x['SNP']]=x
    lookup[trait]=found
pair_fields=['SNP','exposure','outcome','CHR','BP','A1_exposure','A2_exposure','BETA_exposure','SE_exposure','P_exposure','EAF_exposure','N_exposure','F_stat','A1_outcome','A2_outcome','BETA_outcome','SE_outcome','P_outcome','EAF_outcome','N_outcome','status','sample_overlap']
for exposure,outcome in [('LDH','CPI'),('CPI','LDH')]:
    rows=[]
    for snp in leads[exposure]:
        x=allgws[exposure][snp];y=lookup[outcome].get(snp);row=dict(SNP=snp,exposure=exposure,outcome=outcome,CHR=x['CHR'],BP=x['BP'],F_stat=x['F_stat'],status='READY_FOR_TWOSAMPLEMR_HARMONISE' if y else 'OUTCOME_SNP_NOT_FOUND',sample_overlap='LOW')
        for k in ['A1','A2','BETA','SE','P','EAF','N']:
            row[k+'_exposure']=x[k];row[k+'_outcome']=y[k] if y else 'NA'
        rows.append(row)
    counts[exposure]['outcome_SNP_found']=sum(x['status']=='READY_FOR_TWOSAMPLEMR_HARMONISE' for x in rows)
    table(O/(exposure+'_to_'+outcome+'.instrument_pairs.tsv'),rows,pair_fields)
table(O/'MR_INSTRUMENT_QC.tsv',audit,['exposure','SNP','CHR','BP','P','F_stat','reference_QC','clump_status','sample_overlap'])
result={'status':'SUCCESS_INPUT_PREPARATION_NOT_MR_ESTIMATION','reference_n':nref,'reference_build':'GRCh37','counts':{k:dict(v) for k,v in counts.items()},'thresholds':{'instrument_p_exclusive':5e-8,'F_inclusive':10,'clump_r2':.001,'clump_kb':10000,'MHC':'chr6:25-34Mb'},'selection':'Exposure-only significance and F; clump each exposure independently; no outcome-P instrument selection','software':'PLINK1.90b6.21','GWAS_sha256':{k:v['sha256'] for k,v in SPECS.items()},'MR_estimation':'PENDING_FIXED_TWOSAMPLEMR','Steiger':'BLOCKED_UNSPECIFIED_POPULATION_PREVALENCE_FOR_BINARY_LDH'}
(O/'MR_INPUTS_QC.json').write_text(json.dumps(result,indent=2));print(json.dumps(result),flush=True)
