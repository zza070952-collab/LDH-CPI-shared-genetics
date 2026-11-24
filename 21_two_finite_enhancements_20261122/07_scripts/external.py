"""Validate the original public rows and compare published sentinels with fixed leads."""
from pathlib import Path
import csv,io,json,math
import numpy as np,pandas as pd
from openpyxl import load_workbook
N=Path(__file__).resolve().parents[1];O=N/'02_external';A=next((N/'00_preanalysis').glob('LDH_*'));E=A/'external'
read=lambda p:pd.read_csv(p,sep='\t')
write=lambda d,p:d.to_csv(O/p,sep='\t',index=False,na_rep='NA')
f=read(A/'frozen_20261113/01_evidence_tables/LOCUS_MASTER_18.tsv')
ld=read(O/'FROZEN_EUR503_PAIRWISE_LD.tsv')
native=read(O/'LDH_RAW_GRCH38_LEADS.tsv').set_index('rsids')
raw=(E/'MVP_dbGaP_EUR_public.txt').read_text(); pub=read(io.StringIO(raw[raw.index('ID\tAnalysis ID'):]))
assert len(pub)==17857 and pub['P-value'].max()<1e-5
exact=f.merge(pub[pub['SNP ID'].isin(f.lead_snp)],left_on='lead_snp',right_on='SNP ID',how='left',validate='1:1')
known=exact['P-value'].notna()
for _,r in exact[known].iterrows():
    row=native.loc[r.lead_snp]; assert row.pos==r['Chr Position']+1 and row['#chrom']==r['Chr ID']
    pair=set([r.Allele1,r.Allele2]);a=set([row.ref,row.alt]);assert pair==a or pair=={x.translate(str.maketrans('ACGT','TGCA')) for x in a}
exact=exact[['locus_id','lead_snp','chromosome','lead_position','harmonized_effect_allele','other_allele','LDH_beta','CPI_beta','Chr Position','Allele1','Allele2','Coded Allele','P-value','|&beta;|','SE']].rename(columns={'P-value':'MVP_P','|&beta;|':'MVP_abs_beta','Chr Position':'MVP_GRCh38_pos'})
exact=exact.rename(columns={'MVP_GRCh38_pos':'MVP_Chr_Position_as_deposited'})
exact['GRCh38_one_based_verified_against_FinnGen']=exact.MVP_Chr_Position_as_deposited+1
exact['MVP_status']=np.where(known,'EXACT_ASSOCIATION_UNSIGNED','NOT_AVAILABLE_IN_TRUNCATED_PUBLIC_FILE');exact['MVP_adjusted_P']=(18*exact.MVP_P).clip(upper=1);exact['exact_direction']='UNASSESSABLE';exact['MVP_N']=436683
est=read(E/'ESTBB_published64_per_variant_effects.tsv')
assert len(est)==64 and est.beta.notna().sum()==58 and not set(est.rsid)&set(f.lead_snp)
for x in est.itertuples():
    if np.isfinite(x.beta):assert np.isclose(x.P_derived_Wald,math.erfc(abs(x.beta/x.se)/math.sqrt(2)))
write(est,'EstBB_64_selected.tsv')
exact['EstBB_exact_status']='UNAVAILABLE_IN_SELECTED_SENTINEL_TABLE'
write(exact,'EXTERNAL_FIXED18_ALL_STATUS.tsv')
wb=load_workbook(E/'41591_2024_2839_MOESM3_ESM.xlsx',read_only=True,data_only=True)
all_st=[]
for sheet in ['ST2','ST4']:
    rows=list(wb[sheet].values);d=pd.DataFrame(rows[3:],columns=rows[2]);d=d.loc[pd.to_numeric(d.chr,errors='coerce').notna(),['rsID','SNP_GRCh37','chr','GRCh37_bp','beta','se','p']].copy()
    d['wald_P']=[math.erfc(abs(float(b)/float(s))/math.sqrt(2)) for b,s in zip(d.beta,d.se)]
    d['abs_log10P_delta']=abs(np.log10(d.wald_P.clip(lower=1e-300))-np.log10(d.p.astype(float)))
    d['QC']=np.where(d.abs_log10P_delta<=.1,'APPROX_WALD_COMPATIBLE','WALD_DISCREPANT');d['effect_allele_definition']='Allele pair recorded; no explicit effect-allele column';d['sheet']=sheet
    all_st.append(d)
write(pd.concat(all_st),'MVP_PUBLISHED_ROW_QC.tsv')
evid=[]
for r in exact[known].itertuples():
    evid.append(dict(lead=r.lead_snp,variant=r.lead_snp,resource='MVP EUR public',evidence_type='EXACT_ASSOCIATION_UNSIGNED',P=r.MVP_P,adjusted_P=r.MVP_adjusted_P,direction='UNASSESSABLE_UNSIGNED'))
for r in ld.itertuples():
    if not r.r2>=.8:continue
    pair=est[est.rsid==r.other]
    if len(pair) and np.isfinite(pair.iloc[0].beta):
        s=pair.iloc[0];d=f[f.lead_snp==r.lead].iloc[0]
        assert r.other_A1 in (s.EA,s.OA) and r.lead_A1 in (d.harmonized_effect_allele,d.other_allele)
        beta=s.beta*(1 if s.EA==r.other_A1 else -1)*(1 if d.harmonized_effect_allele==r.lead_A1 else -1)*np.sign(r.signed_r_A1_A1)
        evid.append(dict(lead=r.lead,variant=r.other,resource='EstBB selected Salo sentinel',evidence_type='PUBLISHED_VARIANT_SUPPORT',r2=r.r2,P=s.P_derived_Wald,adjusted_P=min(1,18*s.P_derived_Wald),direction='CONCORDANT_LDH' if beta*d.LDH_beta>0 else 'DISCORDANT_LDH',aligned_beta=beta))
    s=all_st[1].loc[all_st[1].rsID==r.other]
    if len(s):evid.append(dict(lead=r.lead,variant=r.other,resource='MVP EUR ST4 selected',evidence_type='PUBLISHED_VARIANT_SUPPORT',r2=r.r2,P=s.iloc[0].p,direction='EFFECT_ALLELE_UNRESOLVED',QC=s.iloc[0].QC))
write(pd.DataFrame(evid),'EXTERNAL_EXACT_AND_PROXY_SUPPORT.tsv')
audit=[['pha013136.1','MVP EUR',436683,'Independent related-pain cohort','Truncated at P<1e-5; absolute beta only','5 exact / 13 unavailable'],['FullData index','MVP EUR',436683,'Same cohort','Set archives; no verified target member in available listing','No signed exact result'],['ST4','MVP EUR',436683,'Same cohort','Published select signals; no explicit effect allele','Separate selected-sentinel evidence'],['ST2','MVP across ancestries',598339,'Overlaps EUR','158 numeric rows; correction refers to 162 index variants; beta/SE/P discrepancies','Source QC only'],['Source Data Fig 42','EstBB','Study-specific','Separate LDH cohort, meta-selected sentinels','64 selected; 58 with estimates; no exact frozen lead','Prior sentinel context'],['Salo meta','FinnGen/UKB/EstBB','Meta','Overlaps discovery cohorts','Full summary covers combined cohorts','Excluded from independent support']]
write(pd.DataFrame(audit,columns=['resource','sample','N','overlap','available_scope','interpretation']),'EXTERNAL_RESOURCE_AND_OVERLAP_AUDIT.tsv')
print('MVP exact',known.sum(),'; selected LD-supported comparisons',len(evid)-known.sum(),'; ST QC',[(x.QC.value_counts().to_dict()) for x in all_st])

