#!/usr/bin/env python3
"""Bounded official Open Targets L2G retrieval; no nearest-gene substitution."""
from pathlib import Path
from collections import defaultdict
import csv,gzip,hashlib,json,os,re,resource,time,urllib.request,urllib.error
R=Path('/analysis/projects/JTM_recent_highIF_direct_clone/P2_LDH_shared_architecture/11_NatCommun_strict_clone_Pro');M=R/'08_smr_heidi_coloc';O=R/'06_fuma_magma_otg/otg_targeted_public';O.mkdir(parents=True,exist_ok=True);C=O/'api_responses';C.mkdir(exist_ok=True)
resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
if hasattr(os,'sched_getaffinity'):os.sched_setaffinity(0,{min(os.sched_getaffinity(0))})
URL='https://api.platform.opentargets.org/api/v4/graphql';calls=[]
def api(query,label,variables=None):
    if (R/'00_admin/PAUSE_REQUESTED').exists():raise SystemExit('PAUSE_REQUESTED')
    payload=dict(query=query,variables=variables or {});key=hashlib.sha256(json.dumps(payload,sort_keys=True).encode()).hexdigest();path=C/(label+'_'+key[:12]+'.json')
    if path.exists():rec=json.loads(path.read_text());assert rec['request']==payload;return rec['response']['data']
    assert sum(f.stat().st_size for f in O.rglob('*') if f.is_file())<200*1024**2,'200MB bound'
    req=urllib.request.Request(URL,data=json.dumps(payload).encode(),headers={'Content-Type':'application/json','User-Agent':'JTM-P2-reproducibility-audit/1.0'})
    for attempt in range(3):
        try:
            with urllib.request.urlopen(req,timeout=75) as h:
                raw=h.read(8*1024**2+1);assert len(raw)<=8*1024**2,'Single-response bound';headers=dict(h.headers);status=h.status
            response=json.loads(raw);assert not response.get('errors'),response.get('errors');break
        except (urllib.error.URLError,TimeoutError) as e:
            if attempt==2:raise
            time.sleep(3*(attempt+1))
    rec=dict(endpoint=URL,retrieved_UTC=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),request=payload,response=response,response_sha256=hashlib.sha256(raw).hexdigest(),http_status=status,headers=headers)
    path.write_text(json.dumps(rec,indent=2));calls.append(str(path));time.sleep(1);return response['data']
def write_tsv(name,rows):
    if not rows:return
    fields=list(dict.fromkeys(k for r in rows for k in r))
    with (O/name).open('w') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rows)
meta=api('{meta{name dataPrefix dataVersion{year month iteration} apiVersion{x y z}}}','release')['meta'];release=meta['dataPrefix'];assert release=='platform2606','Audit a new release explicitly before combining';(O/'RELEASE.json').write_text(json.dumps(meta,indent=2))
GENES=['TGFA','GFPT1','SPOCK2','CHST3','WNT2B','ZNF502','ZNF501'];ids={}
with gzip.open(M/'shared_LDH_CPI_provisional/SMR_ALL_TRAITS_ALL_PANELS.tsv.gz','rt') as f:
    for x in csv.DictReader(f,delimiter='\t'):
        if x['Gene'] in GENES:ids[x['Gene']]='ENSG'+x['probeID'].split('_ENSG')[1]
for r in json.loads((M/'gtex_catalogue_full_cis/GTEX_FULL_CIS_RESOURCE_AUDIT.json').read_text()):ids[r['gene']]=r['ensembl']
assert all(g in ids for g in GENES)
q='query($ids:[String!]!){targets(ensemblIds:$ids){id approvedSymbol genomicLocation{chromosome start end strand}}}'
targets=api(q,'targets',{'ids':[ids[g] for g in GENES]})['targets'];bygene={t['approvedSymbol']:t for t in targets};assert set(bygene)==set(GENES)
write_tsv('TARGET_IDENTIFIERS_GRCH38.tsv',[dict(gene=g,ensembl=ids[g],build='GRCh38',**bygene[g]['genomicLocation']) for g in GENES])
terms=['lumbar disc','disc herniation','intervertebral','sciatica','pain intensity','chronic pain','multisite chronic pain','analgesic','chronic back pain']
studies=set(['FINNGEN_R12_M13_LUMBAR_PROLAPSE','FINNGEN_R12_M13_SCIATICA','GCST008512','GCST90705000']);searchaudit=[]
for term in terms:
    seen=[];idx=0
    while True:
        dat=api('query($q:String!,$i:Int!){search(queryString:$q,entityNames:["study"],page:{index:$i,size:100}){total hits{id name}}}',re.sub('[^a-zA-Z0-9]+','_',term)+'_'+str(idx),{'q':term,'i':idx})['search'];seen.extend(x['id'] for x in dat['hits'])
        if len(seen)>=dat['total']:break
        assert dat['hits'] and idx<10,'Bounded complete pagination';idx+=1
    assert len(set(seen))==len(seen);studies.update(seen);searchaudit.append(dict(query=term,total=dat['total'],retrieved=len(seen),pages=idx+1,status='COMPLETE_SEARCH_INDEX_QUERY'))
write_tsv('STUDY_SEARCH_COVERAGE.tsv',searchaudit)
fields='id traitFromSource nSamples nCases nControls cohorts publicationFirstAuthor pubmedId publicationTitle publicationDate projectId hasSumstats qualityControls analysisFlags'
csfields='studyLocusId studyId chromosome position locusStart locusEnd pValueMantissa pValueExponent finemappingMethod confidence qualityControls variant{id rsIds}'
studyrecords=[];uncovered=[];near={};all_cs=[]
for ix,sid in enumerate(sorted(studies)):
    dat=api('query($s:String!){study(studyId:$s){'+fields+' credibleSets(page:{index:0,size:500}){count rows{'+csfields+'}}}}','study_'+sid,{'s':sid})['study']
    if dat is None:uncovered.append(dict(study_id=sid,status='NOT_COVERED_BY_EXACT_STUDY_ID_IN_RELEASE',release=release));continue
    st={k:v for k,v in dat.items() if k!='credibleSets'};studyrecords.append(st);sets=dat['credibleSets'];rows=sets['rows'];page=0
    while len(rows)<sets['count']:
        page+=1;more=api('query($s:String!,$i:Int!){study(studyId:$s){credibleSets(page:{index:$i,size:500}){count rows{'+csfields+'}}}}','study_'+sid+'_'+str(page),{'s':sid,'i':page})['study']['credibleSets'];assert more['rows'] and page<10;rows.extend(more['rows'])
    assert len(rows)==sets['count'];assert len({x['studyLocusId'] for x in rows})==len(rows)
    for cs in rows:
        all_cs.append(dict(study_id=sid,study_locus_id=cs['studyLocusId'],trait=st['traitFromSource'],chr=cs['chromosome'],bp=cs['position'],build='GRCh38'))
        close=[]
        for g,t in bygene.items():
            loc=t['genomicLocation']
            # Retrieval screen only, never a gene-mapping or evidence criterion.
            if cs['chromosome']==loc['chromosome'] and loc['start']-1000000<=cs['position']<=loc['end']+1000000:close.append(g)
        if close:
            if cs['studyLocusId'] not in near:near[cs['studyLocusId']]=dict(cs=cs,study=st,query_candidate_genes=close)
    print('STUDY',ix+1,len(studies),sid,sets['count'],'near_total',len(near),flush=True)
write_tsv('EXACT_STUDY_NOT_COVERED.tsv',uncovered);(O/'STUDY_METADATA_ALL.json').write_text(json.dumps(studyrecords,indent=2));write_tsv('STUDY_CREDIBLE_SET_INDEX.tsv',all_cs)
def traitclass(trait):
    s=trait.lower()
    if ('disc' in s or 'disk' in s) and any(k in s for k in ['herniat','prolaps']):return 'LDH_RELATED_VERIFY_COHORT_AND_LOCUS'
    if 'sciatica' in s:return 'RADICULAR_AXIS_VERIFY_COHORT_AND_LOCUS'
    if 'pain' in s and 'intensity' in s:return 'PAIN_SEVERITY_VERIFY_COHORT_AND_LOCUS'
    if 'multisite' in s and 'pain' in s:return 'SECONDARY_MULTISITE_PAIN'
    if 'analgesic' in s:return 'TREATMENT_RELATED_90DAY_DEFINITION_NOT_ASSUMED'
    if ('back' in s and 'pain' in s):return 'BACK_PAIN_VERIFY_CHRONICITY_AND_COHORT'
    if 'intervertebral' in s:return 'RELATED_DISC_TRAIT_NOT_AUTOMATICALLY_LDH'
    if 'chronic' in s and 'pain' in s:return 'CHRONIC_PAIN_NOT_AUTOMATICALLY_INTENSITY'
    return 'OTHER_SEARCH_HIT_NOT_PRIMARY_EVIDENCE'
l2grows=[];locusraw=[];locuscoverage=[]
for ix,(lid,obj) in enumerate(sorted(near.items())):
    q='query($l:String!){credibleSet(studyLocusId:$l){studyLocusId l2GPredictions(page:{index:0,size:100}){count rows{studyLocusId score shapBaseValue target{id approvedSymbol} features{name value shapValue}}}}}'
    dat=api(q,'l2g_'+lid,{'l':lid})['credibleSet'];pred=dat['l2GPredictions'];assert len(pred['rows'])==pred['count'],'Need complete L2G pagination';locusraw.append(dict(**obj,l2g=pred));st=obj['study'];cs=obj['cs'];covered=set()
    for pr in pred['rows']:
        target=pr['target'] or {};g=target.get('approvedSymbol')
        if g not in GENES:continue
        assert target['id']==ids[g];covered.add(g)
        positive=[x for x in pr.get('features',[]) if x['shapValue']>0]
        l2grows.append(dict(gene=g,ensembl=ids[g],release=release,api_version='26.6.3',study_id=cs['studyId'],study_locus_id=lid,trait=st['traitFromSource'],trait_class=traitclass(st['traitFromSource']),cohorts=';'.join(st.get('cohorts') or []),study_N=st.get('nSamples'),study_cases=st.get('nCases'),study_author=st.get('publicationFirstAuthor'),pubmed_id=st.get('pubmedId'),build='GRCh38',chr=cs['chromosome'],position=cs['position'],lead_variant=cs['variant']['id'] if cs['variant'] else '',lead_rsids=';'.join(cs['variant']['rsIds'] or []) if cs['variant'] else '',L2G_score=pr['score'],fine_mapping=cs['finemappingMethod'],CS_confidence=cs['confidence'],CS_QC=';'.join(cs.get('qualityControls') or []),L2G_positive_feature_names=';'.join(x['name'] for x in positive),L2G_source='OFFICIAL_PREDICTION_NOT_NEAREST_GENE',formal_status='ANNOTATION_AVAILABLE_REQUIRES_NEW_PROJECT_LOCUS_MATCH',dependence_group=cs['studyId']+'|'+lid,evidence_count_policy='L2G and its source fine-mapped CS are one dependent support group; never count twice',url='https://platform.opentargets.org/credible-set/'+lid))
    for g in obj['query_candidate_genes']:locuscoverage.append(dict(gene=g,study_id=cs['studyId'],study_locus_id=lid,trait=st['traitFromSource'],query_retrieval_window='gene +/-1Mb retrieval only',status='OFFICIAL_L2G_PREDICTION_PRESENT' if g in covered else 'NO_PUBLISHED_L2G_PREDICTION_ABOVE_PLATFORM_DISPLAY_FILTER',not_inferred='No exact zero score, nearest-gene mapping, or biological absence inferred',release=release))
    print('L2G',ix+1,len(near),'candidate_rows',len(l2grows),flush=True)
(O/'L2G_CREDIBLE_SET_RAW_ALL.json').write_text(json.dumps(locusraw,indent=2));write_tsv('OTG_L2G_CANDIDATE_ALL.tsv',l2grows);write_tsv('OTG_L2G_LOCUS_COVERAGE.tsv',locuscoverage)
write_tsv('CURRENT_CANDIDATE_ELIGIBILITY.tsv',[dict(gene=g,official_L2G_prediction_rows=sum(x['gene']==g for x in l2grows),shared_LDH_CPI_SMR='NO_R_VALIDATED',new_conjFDR='PENDING',project_coding_finemapping='PENDING_EXACT_VARIANT_GENE_LINK_AUDIT',n_independent_new_qualifying_layers='NOT_YET_ESTABLISHED',coloc_eligibility='NOT_ELIGIBLE',reason='External L2G requires matching new project locus; single-side SMR is not shared and nearest gene never counts; no duplicate CS/L2G layers') for g in GENES])
metaend=api('{meta{dataPrefix dataVersion{year month iteration}}}','release_end')['meta'];assert metaend['dataPrefix']==release
(O/'OTG_SMALL_RETRIEVAL.done').write_text('SUCCESS_PUBLIC_ANNOTATION_AUDIT_ONLY_NOT_COLOC_ELIGIBILITY\n')
print(json.dumps(dict(release=release,studies=len(studyrecords),near_candidate_credible_sets=len(near),L2G_candidate_rows=len(l2grows),retrieved_bytes=sum(f.stat().st_size for f in O.rglob('*') if f.is_file()),formal_coloc_run=False)))
