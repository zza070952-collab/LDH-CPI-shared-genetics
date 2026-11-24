"""Summarize completed chr18 models without refitting or modifying old results."""
import hashlib
import json
from pathlib import Path
import pandas as pd

R=Path('/analysis/projects/JTM_recent_highIF_direct_clone/P2_LDH_shared_architecture/11_NatCommun_strict_clone_Pro')
NEW=R/'19_computational_finalize_JTM_20261103'
O=NEW/'03_signal_resolution'

def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as hnd:
        for block in iter(lambda:hnd.read(1024**2),b''):h.update(block)
    return h.hexdigest()

prep=json.loads((O/'CHR18_PREPARED.json').read_text())
for p,h in prep['outputs'].items():assert sha(p)==h
fits=pd.read_csv(O/'CHR18_SUSIE_FIT_QC.tsv',sep='\t')
qc=pd.read_csv(O/'CHR18_RSS_LD_N_DIAGNOSTIC.tsv',sep='\t')
pip=pd.read_csv(O/'CHR18_SUSIE_PIP_ALL.tsv',sep='\t')
cs=pd.read_csv(O/'CHR18_SUSIE_CREDIBLE_SETS_ALL.tsv',sep='\t')
assert len(fits)==4 and fits.converged.all() and qc['passes_fixed_s_0.1'].all()
co=pd.read_csv(O/'CHR18_COLOC_SUSIE_ALL_N_PRIORS.tsv',sep='\t')
assert len(co)==9 and co.status.eq('NO_CREDIBLE_SET_NOT_PP4_ZERO').all()
for name in ['CHR18_COLOC_SUSIE_ALL_N_PRIORS.tsv','CHR18_COLOC_SUSIE_PRIMARY.tsv']:
    x=pd.read_csv(O/name,sep='\t')
    for col in ['nsnps','hit1','hit2','PP.H0.abf','PP.H1.abf','PP.H2.abf','PP.H3.abf','PP.H4.abf']:
        if col not in x:x[col]=pd.NA
    x['scientific_interpretation']='Not estimable because CPI produced no 95% credible set; no evidence for H3 or PP4=0 is inferred.'
    x.to_csv(O/name,sep='\t',index=False,na_rep='NA')
summary={'execution_status':'COMPLETED','scientific_resolution_status':'UNRESOLVED_NO_CPI_CREDIBLE_SET','fixed_locus_id':'LDH_CPI_p05_chr18_001','build':'GRCh37','window':[18,50358109,51055069],'n_harmonized_source':1357,'n_susie_input':1345,'n_reference_missing':12,'n_reference_samples':503,'reference_source':'Existing official MAGMA 1000 Genomes Phase3 EUR503 archive, exact historical SHA256 verified','reference_sha256':prep['signature']['reference_archive_sha256'],'n_fits':4,'n_converged':int(fits.converged.sum()),'LDH_n_CS':1,'CPI_n_CS_by_N':{str(int(r.N)):int(r.n_CS) for r in fits.itertuples() if r.trait=='CPI'},'LDH_CS_size':len(cs),'LDH_CS_coverage':float(cs.coverage.iloc[0]),'LDH_CS_min_abs_corr':float(cs.purity_min_abs_corr.iloc[0]),'RSS_s':{str(r.trait)+'_'+str(int(r.N)):float(r.RSS_s) for r in qc.itertuples()},'n_N_prior_configurations':9,'n_estimable_signal_pairs':0,'coloc_signal_PP4':None,'PP4_unavailable_reason':'No CPI credible set under unchanged coverage .95 and purity min_abs_corr .5, at any source-derived N value','max_PIP_by_fit':{id:{'SNP':g.loc[g.PIP.idxmax(),'SNP'],'BP':int(g.loc[g.PIP.idxmax(),'BP']),'PIP':float(g.PIP.max())} for id,g in pip.groupby('fit_id')},'new_QTL_panels':0,'repeated_model_fits':0,'source_beta_alignment_note':'BETA_CPI was already flipped in historical revision_locus_coloc.R; no second flip used','historical_results_preserved':True}
(O/'CHR18_ANALYSIS_SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n')
text='''# Targeted chromosome 18 signal resolution

Execution completed with exit code 0; scientific signal-pair resolution remains unresolved.

The unchanged GRCh37 interval 18:50,358,109–51,055,069 was examined using all 1,357 already harmonized LDH–CPI variants. An existing, checksum-verified official MAGMA 1000 Genomes Phase 3 European archive provided 503 reference individuals. This reference had previously supported other P2 SuSiE audits; it is distinct from the UK Biobank reference used for LAVA. The original extracted full-reference folder no longer existed, so only the required genotypes were decoded directly from the existing archive. No reference or GWAS data were downloaded. Twelve SNP IDs were absent from the existing reference; all 1,345 remaining variants were retained, without P-value selection, monomorphic-variant exclusion or missing-genotype imputation. The reference allele dosages were matched by SNP, chromosome, position and allele pair, then signed to the frozen LDH effect allele. The CPI beta field in the historical harmonized table had already been signed to LDH by revision_locus_coloc.R; it was not flipped a second time.

Signed LD was symmetric with unit diagonal to numerical precision and positive semidefinite within the existing tolerance. The official susieR RSS consistency diagnostic was applied to each fixed trait/sample-size model. The LDH s estimate was 0.0111114864. CPI s ranged from 0.0120920867 to 0.0120929873 across source-derived N values 109,571, 134,396 and 134,627. All models passed the project's pre-existing s≤0.1 operational rule. This operational rule is not presented as a universal cutoff specified by susieR.

Using coloc 5.2.2 and susieR 0.12.35 under R 4.3.3, SuSiE-RSS was run with L=10, maximum 200 iterations, 95% credible-set coverage, minimum absolute within-set correlation 0.5 and residual-variance estimation disabled for external LD. The fixed LDH N was 388,514. CPI used the minimum, median and maximum per-SNP N of the retained region as explicit scalar-N sensitivity models; these are not exact variable-N likelihoods, and no best-N model was selected. The median-N model was primary. All four fits converged. LDH formed one credible set; CPI formed no credible set at any of its three N values. The LDH set contained {CS_SIZE} SNPs, with attained coverage {CS_COVERAGE} and minimum absolute correlation {CS_PURITY}.

All nine combinations of three CPI N values and three frozen p12 priors (10^-6, 10^-5 and 5×10^-5; p1=p2=10^-4) were retained. No LDH–CPI signal-pair posterior could be estimated because CPI had no credible set. The unavailable PP4 values are NA. They do not indicate PP4=0, a distinct-signal result, or a failed software run. Coverage, purity, N, reference and window were not adjusted in response to this finding.

The existing regional ABF PP4≈0.9835 and LAVA local rg≈0.469 (BH q≈0.0350) remain complementary locus-level evidence from the same two GWAS. The targeted supplement did not resolve whether a particular LDH and CPI signal colocalizes at multi-signal resolution. It adds no gene-specific attribution to DCC. No further chr18 model search is required for the bounded computational manuscript.

Official method guidance consulted on 2026-11-13: https://stephenslab.github.io/susieR/articles/susierss_diagnostic.html and https://chr1swallace.github.io/coloc/articles/a06_SuSiE.html . Their current online package versions were not substituted for the installed historical versions.
'''.replace('{CS_SIZE}',str(len(cs))).replace('{CS_COVERAGE}',f'{float(cs.coverage.iloc[0]):.9f}').replace('{CS_PURITY}',f'{float(cs.purity_min_abs_corr.iloc[0]):.9f}')
(O/'CHR18_INTERPRETATION.md').write_text(text)
manifest=[]
for key in ['source','reference_archive']:
    manifest.append({'role':'historical_input','path':prep['signature'][key],'sha256':prep['signature'][key+'_sha256'],'bytes':Path(prep['signature'][key]).stat().st_size})
for p in sorted(O.glob('CHR18_*')):
    if p.is_file() and p.name!='CHR18_INPUT_OUTPUT_MANIFEST.tsv':manifest.append({'role':'new_output','path':str(p),'sha256':sha(p),'bytes':p.stat().st_size})
for p in sorted((NEW/'scripts').glob('chr18*')):
    if p.is_file():manifest.append({'role':'new_script','path':str(p),'sha256':sha(p),'bytes':p.stat().st_size})
pd.DataFrame(manifest).to_csv(O/'CHR18_INPUT_OUTPUT_MANIFEST.tsv',sep='\t',index=False)
print(json.dumps(summary,indent=2))
