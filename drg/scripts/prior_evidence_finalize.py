"""Bounded reconciliation of published locus evidence and preserved L2G records.

No genetic model is fitted and no historical result is modified.
"""
from pathlib import Path
import csv,gzip,hashlib,json,re
import openpyxl

R=Path('/analysis/projects/JTM_recent_highIF_direct_clone/P2_LDH_shared_architecture/11_NatCommun_strict_clone_Pro')
N=R/'19_computational_finalize_JTM_20261103'
D=N/'02_locus_evidence';D.mkdir(parents=True,exist_ok=True)
G=N/'04_gene_evidence';G.mkdir(parents=True,exist_ok=True)
S=D/'novelty_sources'
OT=R/'06_fuma_magma_otg/otg_new_loci_zero'

def read(p):
 with Path(p).open() as h:return list(csv.DictReader(h,delimiter='\t'))
def write(p,rr):
 fields=list(dict.fromkeys(k for r in rr for k in r))
 with Path(p).open('w') as h:
  w=csv.DictWriter(h,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rr)
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(2**20),b''):h.update(b)
 return h.hexdigest()

fixed=read(R/'17_evidence_tier_revision_20260912/01_locus_evidence/LOCUS_EVIDENCE_18.tsv')
assert len(fixed)==18
fby={x['locus_id']:x for x in fixed}
orig={x['locus_id']:x for x in read(R/'06_fuma_magma_otg/local_fuma_final_zero/FUMA_LOCI_ALL.tsv') if x['locus_id'] in fby}
all_l2g=read(OT/'OTG_L2G_CANDIDATE_ALL.tsv')
matches=[x for x in read(OT/'L2G_NEW_INTERVAL_MATCHES.tsv') if x['locus_id'] in fby]
trace=[]
for m in matches:
 sid=m['dependence_group'].split('|')[1]
 rr=[x for x in all_l2g if x['gene']==m['gene'] and x['study_id']==m['study'] and x['study_locus_id']==sid]
 assert len(rr)==1,(m,rr)
 a=rr[0];assert abs(float(a['L2G_score'])-float(m['L2G_score']))<1e-12
 z=dict(a);z.update(locus_id=m['locus_id'],GRCh37_chr=m['GRCh37_chr'],GRCh37_bp=m['GRCh37_bp'],rsid_mapped=m['rsid'],mapping_rule='Exact rsID occurs once in previously audited GRCh37 EUR reference and falls within unchanged fixed locus',same_primary_LDH_endpoint=a['study_id']=='FINNGEN_R12_M13_LUMBAR_PROLAPSE',same_primary_Wang_CPI='NOT_ESTABLISHED',signal_relationship='SAME_FIXED_REGION_SPECIFIC_SIGNAL_UNRESOLVED',support_category='OTHER_TRAIT_REGION_BACKGROUND',LD_with_current_lead='UNKNOWN_NOT_ESTIMATED',primary_lead=fby[m['locus_id']]['lead_SNP'],primary_lead_GRCh37_bp=fby[m['locus_id']]['lead_position'],source_file=str(OT/'OTG_L2G_CANDIDATE_ALL.tsv'),coordinate_source_file=str(OT/'L2G_NEW_INTERVAL_MATCHES.tsv'),accessed='2026-11-13')
 trace.append(z)
write(G/'L2G_TRACEABILITY.tsv',trace)

# Original CPI paper: six primary GWAS intervals in GRCh37 (exclude no-MRI GWAS).
wb=openpyxl.load_workbook(S/'Wang2025_SupplementaryData2.xlsx',read_only=True,data_only=True)
wang=[]
for rowno,v in enumerate(wb.worksheets[0].iter_rows(min_row=4,max_row=9,max_col=26,values_only=True),4):
 assert isinstance(v[0],str) and v[0].startswith('rs')
 wang.append(dict(study='Wang2025_CPI',doi='10.1038/s41467-025-64904-y',rsid=v[0],chr=int(v[1]),position=int(v[2]),start=int(v[3]),end=int(v[4]),build='GRCh37',gene_label=v[12],P=v[9],source_file=str(S/'Wang2025_SupplementaryData2.xlsx'),source_row=rowno))

# Salo supplementary Table 2 is GRCh38. Obtain exact GRCh37 mappings from the
# already harmonized Salo GWAS, requiring agreement with published source build.
wb=openpyxl.load_workbook(S/'Salo2024_SourceData.xlsx',read_only=True,data_only=True)
salo=[]
for rowno,v in enumerate(wb['Supplementary Table 2'].iter_rows(min_row=6,max_col=13,values_only=True),6):
 rs=str(v[4] or '').strip().replace('*','')
 if not re.fullmatch('rs[0-9]+',rs):continue
 coord=str(v[3] or '')
 if not re.fullmatch('[0-9]+:[0-9]+',coord):continue
 ch,bp=coord.split(':')
 salo.append(dict(study='Salo2024_LDH',doi='10.1038/s41467-024-53467-z',rsid=rs,signal_type='CONDITIONAL_SECONDARY' if '*' in str(v[2]) else 'PRIMARY',source_chr=int(ch),source_position=int(bp),source_build='GRCh38',gene_label=v[2] or '',P=v[9],source_file=str(S/'Salo2024_SourceData.xlsx'),source_sheet='Supplementary Table 2',source_row=rowno))
want={x['rsid'] for x in salo};mapped={}
qc=R/'02_qc/gwas/SALO2024_LDH_META.GRCh37.QC.tsv.gz'
prior_map=D/'novelty_PUBLISHED_PRIMARY_LEADS.tsv'
prior={v['rsid']:v for v in read(prior_map) if v['study']=='Salo2024_LDH'} if prior_map.exists() else {}
if set(prior)==want:
 for x in salo:
  old=prior[x['rsid']]
  assert int(old['source_chr'])==x['source_chr'] and int(old['source_position'])==x['source_position']
  x.update(chr=int(old['chr']) if old['chr']!='NA' else 'NA',position=int(old['position']) if old['position']!='NA' else 'NA',build='GRCh37',coordinate_status=old['coordinate_status'],coordinate_source_file=old['coordinate_source_file'])
else:
 with gzip.open(qc,'rt') as h:
  for row in csv.DictReader(h,delimiter='\t'):
   if row['SNP'] in want:mapped.setdefault(row['SNP'],[]).append(row)
 for x in salo:
  rr=[v for v in mapped.get(x['rsid'],[]) if int(v['SOURCE_CHR'])==x['source_chr'] and int(v['SOURCE_BP'])==x['source_position']]
  if len(rr)==1:
   x.update(chr=int(rr[0]['CHR']),position=int(rr[0]['BP']),build='GRCh37',coordinate_status='UNIQUE_EXACT_RSID_AND_SOURCE_COORDINATE_MATCH_IN_FROZEN_QC',coordinate_source_file=str(qc))
  else:x.update(chr='NA',position='NA',build='GRCh37',coordinate_status='UNRESOLVED_COORDINATE_OR_NOT_RETAINED_IN_QC',coordinate_source_file=str(qc))
published=wang+salo
write(D/'novelty_PUBLISHED_PRIMARY_LEADS.tsv',published)

rows=[];summary=[]
for f in fixed:
 lid=f['locus_id'];ch=int(f['chr']);start=int(f['start']);end=int(f['end'])
 wr=[x for x in wang if x['chr']==ch and x['start']<=end and start<=x['end']]
 sr=[x for x in salo if x['chr']!='NA' and x['chr']==ch and start<=x['position']<=end]
 tx=[x for x in trace if x['locus_id']==lid]
 for study,records,rule in [('Wang2025_CPI',wr,'Published primary GRCh37 region overlaps unchanged fixed region'),('Salo2024_LDH',sr,'Published lead, exactly mapped from GRCh38 to GRCh37 using frozen study QC, lies in fixed region')]:
  if records:
   for a in records:
    rows.append(dict(locus_id=lid,lead_snp=f['lead_SNP'],evidence_type='PUBLISHED_SINGLE_TRAIT_REGION',study=study,doi=a['doi'],study_lead=a['rsid'],study_gene_label=a['gene_label'],study_P=a['P'],build='GRCh37',position=a['position'],match_rule=rule,match_status='MATCH',same_signal_with_current_lead='EXACT_SENTINEL' if a['rsid']==f['lead_SNP'] else 'UNKNOWN_WITHOUT_LD_OR_COLOC',evidence_scope='Region or exact variant association; gene label is not causal gene proof',cohort_overlap='FinnGen+EstBB+UKB meta includes both source biobanks' if study=='Salo2024_LDH' else 'Same UKB primary CPI GWAS',independent_replication='NO',source_file=a['source_file'],source_row=a['source_row'],search_cutoff='2026-11-13'))
  else:
   rows.append(dict(locus_id=lid,lead_snp=f['lead_SNP'],evidence_type='PUBLISHED_SINGLE_TRAIT_REGION',study=study,doi='10.1038/s41467-025-64904-y' if study=='Wang2025_CPI' else '10.1038/s41467-024-53467-z',match_rule=rule,match_status='NO_MATCH_IN_THIS_PUBLISHED_LEAD_TABLE',same_signal_with_current_lead='UNKNOWN',evidence_scope='No matched lead under stated rule is not proof of a novel region',independent_replication='NO',search_cutoff='2026-11-13'))
 # Curated relevant-study sentinels already retrieved from official OTG/GWAS
 # metadata. Preserve exact provenance; do not substitute for main GWAS.
 seen=set()
 for a in tx:
  key=a['study_id'],a['study_locus_id']
  if key in seen:continue
  seen.add(key)
  rows.append(dict(locus_id=lid,lead_snp=f['lead_SNP'],evidence_type='CURATED_PUBLISHED_RELATED_TRAIT_REGION',study=a['study_id'],doi='',pmid=a['pubmed_id'],study_trait=a['trait'],study_lead=a['lead_rsids'],build='GRCh37',position=a['GRCh37_bp'],match_rule=a['mapping_rule'],match_status='MATCH',same_signal_with_current_lead='EXACT_SENTINEL' if f['lead_SNP'] in a['lead_rsids'].split(';') else 'UNKNOWN_WITHOUT_LD_OR_COLOC',evidence_scope='Related-phenotype regional annotation; publication and CS recorded; no established LDH-CPI shared signal inferred',cohort_overlap='Potentially shared UKB or FinnGen samples, or different phenotype in external cohort',independent_replication='NOT_TESTED',source_file=a['source_file'],database_release=a['release'],source_study_locus_id=a['study_locus_id'],search_cutoff='2026-11-13'))
 category='PUBLISHED_MATCH_BOTH_PRIMARY_STUDIES' if wr and sr else 'PUBLISHED_PRIMARY_CPI_REGION' if wr else 'PUBLISHED_PRIMARY_LDH_LEAD' if sr else 'PUBLISHED_RELATED_TRAIT_REGION' if tx else 'NO_MATCH_IN_BOUNDED_PUBLISHED_TABLES'
 summary.append(dict(locus_id=lid,lead_snp=f['lead_SNP'],published_LDH_match=len(sr)>0,published_CPI_match=len(wr)>0,prior_LDH_association='MATCHED_PUBLISHED_LDH_LEAD_IN_FIXED_REGION' if sr else 'NO_LEAD_MATCH_UNDER_FIXED_RULE_NOT_PROOF_OF_NOVELTY',prior_pain_association='MATCHED_PRIMARY_CPI_REGION' if wr else 'MATCHED_RELATED_PAIN_OR_DISC_STUDY' if tx else 'NOT_IDENTIFIED_IN_BOUNDED_SOURCES',prior_shared_signal_evidence='NO_PREVIOUS_EXACT_LDH_CPI_SHARED_SIGNAL_ESTABLISHED_IN_BOUNDED_REVIEW',novelty_category=category,current_GWAS_interval_category=orig[lid]['original_GWAS_interval_category'],n_L2G_related_study_loci=len(seen),published_LDH_leads=';'.join(x['rsid'] for x in sr),published_CPI_leads=';'.join(x['rsid'] for x in wr),same_signal_claim='No signal equivalence inferred from coordinate overlap or gene label',search_cutoff='2026-11-13',source_file='NOVELTY_AND_PRIOR_EVIDENCE.tsv'))

# Prior DCC mechanism is literature context only, separate from our Results.
rows.append(dict(locus_id='LDH_CPI_p05_chr18_001',lead_snp='rs62100775',evidence_type='PRIOR_EXPERIMENTAL_BIOLOGY',study='Pan2026_DCC_male_mouse',doi='10.1038/s41467-026-75423-9',pmid='42448693',match_status='GENE_BIOLOGY_ONLY',same_signal_with_current_lead='NOT_TESTED',evidence_scope='Published male-mouse endplate/DRG/dorsal-horn mechanism, not a result of the current human genetic study',independent_replication='NO_HUMAN_SHARED_SIGNAL_TEST',source_file='https://pubmed.ncbi.nlm.nih.gov/42448693/',search_cutoff='2026-11-13'))
for a in salo:
 if a['gene_label']=='DCC':
  rows.append(dict(locus_id='LDH_CPI_p05_chr18_001',lead_snp='rs62100775',evidence_type='PUBLISHED_DCC_GENE_LABEL_LDH',study='Salo2024_LDH',doi=a['doi'],pmid='39511132',study_lead=a['rsid'],study_gene_label='DCC',study_P=a['P'],build=a['source_build'],position=a['source_position'],match_rule='Primary paper explicitly names DCC; source lead recorded; exact current-signal relationship unresolved',match_status='GENE_LABEL_PRIOR_ASSOCIATION',same_signal_with_current_lead='UNKNOWN',evidence_scope='Published LDH locus labelled DCC; not proof of DCC gene causality or identity with current shared signal',independent_replication='NO',source_file=a['source_file'],source_row=a['source_row'],search_cutoff='2026-11-13'))
for a in summary:
 if a['locus_id']=='LDH_CPI_p05_chr18_001':a['prior_LDH_association']='DCC_LOCUS_PREVIOUSLY_REPORTED_IN_SALO2024_EXACT_CURRENT_SIGNAL_RELATIONSHIP_UNRESOLVED'
cat_path=D/'novelty_CATALOG_EXACT_SENTINELS.tsv'
if cat_path.exists():
 for a in read(cat_path):
  rows.append(dict(locus_id=a['locus_id'],lead_snp=a['rsid'],evidence_type='GWAS_CATALOG_EXACT_SENTINEL',study=a['accession_id'],pmid=a['pubmed_id'],study_trait=a['trait'],study_lead=a['rsid'],study_P=a['p_value'],match_rule='Exact rsID in GWAS Catalog v2; queried current official API; no new LD inference',match_status='MATCH',same_signal_with_current_lead='EXACT_SENTINEL',evidence_scope=a['support_scope'],independent_replication='NOT_TESTED',source_file=a['source_url'],search_cutoff='2026-11-13'))
write(D/'NOVELTY_AND_PRIOR_EVIDENCE.tsv',rows)
write(D/'novelty_LOCUS_PRIOR_SUMMARY.tsv',summary)
counts={k:sum(x['novelty_category']==k for x in summary) for k in sorted({x['novelty_category'] for x in summary})}
status=dict(search_cutoff='2026-11-13',fixed_loci=18,L2G_traceable_rows=len(trace),published_Salo_rows=len(salo),published_Salo_unique_rsid=len({x['rsid'] for x in salo}),published_Salo_mapped=sum(x['chr']!='NA' for x in salo),published_Wang_primary_regions=6,novelty_counts=counts,global_first_claim_allowed=False,genetic_models_rerun=False)
(D/'novelty_AUDIT_SUMMARY.json').write_text(json.dumps(status,indent=2))
print(json.dumps(status,indent=2))
print('DCC_MAX_L2G',[(x['study_id'],x['L2G_score'],x['lead_rsids'],x['GRCh37_bp']) for x in trace if x['gene']=='DCC'])
