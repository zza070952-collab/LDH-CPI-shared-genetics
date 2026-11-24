import json
from pathlib import Path
import numpy as np,pandas as pd
P=Path('/analysis/projects/JTM_recent_highIF_direct_clone/P2_LDH_shared_architecture')
N=P/'11_NatCommun_strict_clone_Pro/19_computational_finalize_JTM_20261103'
D=N/'05_DRG'
d=pd.read_csv(D/'DRG_CANDIDATE_DONOR_EXPRESSION.tsv',sep='\t')
assert len(d)==5243 and not d[['dataset','donor','cell_class','gene']].duplicated().any()
meas=d[d.status.str.startswith('MEASURED')]
assert (meas.n_expressing>=0).all() and (meas.n_expressing<=meas.n_cells).all()
assert ((meas.percent_expressing>=0)&(meas.percent_expressing<=100)).all()
assert np.allclose(meas.percent_expressing,100*meas.n_expressing/meas.n_cells)
missing=d[d.status.eq('GENE_NOT_IN_DEPOSITED_MATRIX')]
assert missing.percent_expressing.isna().all() and missing.mean_log1p_CP10K.isna().all()
old=pd.read_csv(P/'06_results_tables/P2_HUMAN_DRG_DONOR_CELLTYPE_SUMMARY.tsv',sep='\t')
old=old.rename(columns={'cell_type':'cell_class','mean_expression':'mean_log1p_CP10K_old','n_expressing':'n_expressing_old','n_cells':'n_cells_old'})
m=d.merge(old,on=['dataset','donor','cell_class','gene'],how='inner')
assert len(m)>0 and (m.n_expressing==m.n_expressing_old).all() and (m.n_cells==m.n_cells_old).all()
assert np.allclose(m.mean_log1p_CP10K,m.mean_log1p_CP10K_old,rtol=0,atol=1e-10)
result={'status':'PASS','unique_expression_rows':len(d),'all_count_proportions_consistent':True,'missing_values_remain_NA':True,'old_successful_target_rows_crosschecked':len(m),'raw_count_detection_exact_agreement':True,'mean_log1p_CP10K_agrees_atol':1e-10}
(D/'DRG_VALIDATION.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
