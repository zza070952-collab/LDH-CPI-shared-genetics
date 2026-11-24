import json
from pathlib import Path
import pandas as pd
R=Path(__file__).resolve().parents[1];D=R/'08_smr_heidi_coloc/shared_available4_traits_zero';q=json.loads((D/'MODULE_STATUS.json').read_text());assert q['status']=='SUCCESS_AVAILABLE_4_TRAITS_32_PANELS'
z=pd.read_csv(D/'SMR_ALL_TRAITS_ALL_PANELS.tsv.gz',sep='\t');p=pd.read_csv(D/'SMR_PROBE_TISSUE_PAIR_AUDIT_ALL.tsv',sep='\t');rows=[]
for t,g in z.groupby('trait'):
 rows.append({'trait':t,'n_tests':len(g),'n_global_BH05':int((g.FDR_global_all_traits_panels<.05).sum()),'n_strict_BH_HEIDI05':int(g.primary_pass.sum()),'n_legacy_BH_HEIDI01':int(g.legacy_HEIDI01_pass.sum()),'n_within_trait_BH_HEIDI05_sensitivity':int(((g.FDR_within_trait_all_panels<.05)&(g.p_HEIDI>.05)&(~g.MHC_flag)).sum()),'multiplicity_family':'14591 available4-trait x8-panel hypotheses; missing resources explicitly excluded, not full registered family'})
pd.DataFrame(rows).to_csv(D/'SMR_TRAIT_LEVEL_COUNTS.tsv',sep='\t',index=False)
r=[]
for t,g in p.groupby('trait_pain'):
 aligned=g.allele_alignment.isin(['SAME','FLIPPED','COMPLEMENT','COMPLEMENT_FLIPPED'])&g.eqtl_effect_consistent
 legacy=g.legacy_HEIDI01_pass_LDH&g.legacy_HEIDI01_pass_pain&aligned
 r.append({'pain_trait':t,'n_same_probe_panel_pairs':len(g),'n_strict_shared':int(g.shared_primary_pass.sum()),'n_legacy_HEIDI01_shared_sensitivity':int(legacy.sum()),'scope':'PRIMARY' if t!='MCP_SECONDARY' else 'SECONDARY_NOT_REPLACEMENT'})
pd.DataFrame(r).to_csv(D/'SMR_SHARED_COUNTS_BY_AXIS.tsv',sep='\t',index=False)
z[z.primary_pass].to_csv(D/'SMR_SINGLE_TRAIT_PASS_NOT_SHARED.tsv',sep='\t',index=False)
print(pd.DataFrame(rows).to_string(index=False));print(pd.DataFrame(r).to_string(index=False))
