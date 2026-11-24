"""Independent multiplicity, alignment, arithmetic and evidence-boundary checks."""
from pathlib import Path
import pandas as pd,numpy as np,hashlib,json,math
N=Path(__file__).resolve().parents[1];D=N/'01_twas';read=lambda p:pd.read_csv(p,sep='\t',low_memory=False)
for p in D.glob('*.tests.tsv.gz'):
    assert hashlib.file_digest(p.open('rb'),'sha256').hexdigest()==p.with_suffix('.sha256').read_text().strip()
r=read(D/'TWAS_ALL_TESTS.tsv.gz');j=read(D/'TWAS_SAME_TISSUE_CONJUNCTION.tsv');s=read(D/'TWAS_CURRENT101_AND_HISTORICAL6.tsv');e=read(N/'02_external/EXTERNAL_FIXED18_ALL_STATUS.tsv')
families=[]
def verify(x,p,q):
    x=x[x[q].notna()].sort_values(p,ascending=False);best=1.;n=len(x)
    for rank,(a,b) in enumerate(zip(x[p],x[q])):
        best=min(best,a*n/(n-rank));assert np.isclose(best,b,rtol=2e-9,atol=2e-13),(q,a,b,best)
    return n
for (fam,tr),frame in r[r.predictor.eq('trait')].groupby(['family','trait']):
    ok=frame[frame.status.eq('ELIGIBLE')];assert len(ok)==frame.q_trait_global4tissue.notna().sum()
    families.append([fam,tr,'trait BH',verify(frame,'P','q_trait_global4tissue')])
    if fam=='en' and tr=='CPI':
        f=ok[ok.q_phi_global4tissue.notna()];assert (f.n_allele_aligned==f.n_weight_snps).all()
        families.append([fam,tr,'matched phi BH',verify(f,'CPI_phi_P','q_phi_global4tissue')])
for fam,f in j.groupby('family'):
    ok=f[f.q_conjunction.notna()]
    assert (ok.subset_hash_CPI==ok.subset_hash_LDH).all() and np.allclose(ok.P_conjunction,np.maximum(ok.P_LDH,ok.P_CPI))
    families.append([fam,'both','common predictor BH',verify(f,'P_conjunction','q_conjunction')])
checks=pd.concat([read(p).assign(source=p.name) for p in D.glob('*.checks.tsv')])
assert set(['DCC','BSN','APEH'])==set(checks.symbol) and np.allclose(checks.official_Z.astype(float),checks.checked_Z.astype(float))
assert len(e)==18 and e.MVP_P.notna().sum()==5 and e.MVP_P.isna().sum()==13 and e.lead_snp.nunique()==18
f=e[e.MVP_P.notna()];assert np.allclose(f.MVP_adjusted_P,f.MVP_P*18) and (f.MVP_adjusted_P<.05).all() and e.exact_direction.eq('UNASSESSABLE').all()
assert len(s)==856 and s.groupby('scope').gene.nunique().to_dict()=={'CURRENT101':101,'HISTORICAL6':6}
for fam in ['mashr','en']:assert (N/'logs'/f'{fam}.exitcode').read_text().strip()=='0'
assert len(list(D.glob('*.tests.tsv.sha256')))==8
pd.DataFrame(families,columns=['family','trait','check','independently_verified_denominator']).to_csv(N/'00_admin/INDEPENDENT_MULTIPLICITY_QC.tsv',sep='\t',index=False)
print('PASS:',len(families),'BH families,',len(checks),'manual Z comparisons; 18 fixed leads, five unsigned exact tests, 13 unavailable; 856 candidate/model/tissue rows; eight model checkpoints.')

