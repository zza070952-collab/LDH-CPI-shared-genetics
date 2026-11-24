from pathlib import Path
import pandas as pd,json
P=Path('/analysis/projects/JTM_recent_highIF_direct_clone/P2_LDH_shared_architecture')
p=P/'05_results_intermediate/mda/drg/GSE168243_RAW/GSM5134537_prep1.csv.gz'
d=pd.read_csv(p,usecols=[0]).iloc[:,0].astype(str)
s=d.str.replace(r'\..*$','',regex=True)
print(json.dumps({'n_raw_labels':len(d),'n_dotted_labels':int(d.str.contains('.',regex=False).sum()),'n_labels_changed':int(d.ne(s).sum()),'n_collisions_created':len(s)-s.nunique(),'examples_changed':d[d.ne(s)].head(25).to_list()},indent=2))
