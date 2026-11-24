#!/usr/bin/env Rscript
# Formal locus-complete pairwise colocalization after new evidence eligibility.
suppressPackageStartupMessages({library(data.table);library(coloc);library(susieR)})
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=2)stop('Usage: Rscript replan_regulation_coloc.R jobs.tsv output_directory')
jobs<-fread(args[1],na.strings=c('NA',''));outdir<-args[2];dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
stopifnot(as.character(packageVersion('coloc'))=='5.2.2',as.character(packageVersion('susieR'))=='0.12.35')
required<-c('job_id','gene','qtl_panel','axis','comparison','input1','input2','type1','type2','N1','N2','s1','s2','sdY1','sdY2','build','locus_chr','locus_start','locus_end','eligibility_pass','evidence_layer_count','pairing_approved','input1_complete_window','input2_complete_window')
stopifnot(all(required%in%names(jobs)))
 result<-list();sensitivity<-list();align_qc<-list();susie_qc<-list();conditional_sets<-list();ii<-0L
read_input<-function(path,j){
 d<-fread(path);need<-c('SNP','chr','bp','A1','A2','beta','SE','MAF');stopifnot(all(need%in%names(d)))
 d<-d[chr==j$locus_chr & bp>=j$locus_start & bp<=j$locus_end]
 d<-d[is.finite(beta)&is.finite(SE)&SE>0&is.finite(MAF)&MAF>.01&MAF<=.5]
 d[,`:=`(A1=toupper(A1),A2=toupper(A2))]
 d<-d[A1%in%c('A','C','G','T')&A2%in%c('A','C','G','T')&A1!=A2]
 # Drop all duplicates and ambiguous palindromes: QTL metadata has MAF, not EAF.
 duplicate<-duplicated(d$SNP)|duplicated(d$SNP,fromLast=TRUE);d<-d[!duplicate]
 d<-d[!paste0(A1,A2)%in%c('AT','TA','CG','GC')]
 d
}
make_dataset<-function(z,side,j){
 suffix<-paste0('_',side);grab<-function(x)z[[paste0(x,suffix)]]
 nn<-if(paste0('N',suffix)%in%names(z))as.numeric(grab('N')) else as.numeric(j[[paste0('N',side)]])
 if(any(!is.finite(nn)|nn<=0))stop('Invalid source sample size; never replace per-variant N with maximum cohort N')
 ds<-list(snp=z$SNP,position=z$bp_1,beta=grab('beta'),varbeta=grab('SE')^2,MAF=grab('MAF'),type=j[[paste0('type',side)]],N=nn)
 if(ds$type=='cc'){ss<-as.numeric(j[[paste0('s',side)]]);if(!is.finite(ss)||ss<=0||ss>=1)stop('Audited case fraction required');ds$s<-ss}
 if(ds$type=='quant'){
  yy<-as.numeric(j[[paste0('sdY',side)]])
  if(is.finite(yy)&&yy>0)ds$sdY<-yy
  else if(any(!is.finite(ds$N)|ds$N<=0))stop('N/MAF needed to estimate sdY; do not invent trait variance')
 }
 check_dataset(ds);ds
}
for(k in seq_len(nrow(jobs))){
 j<-as.list(jobs[k]);pref<-file.path(outdir,j$job_id)
 if(j$eligibility_pass!='YES'||j$evidence_layer_count<2||j$pairing_approved!='YES')stop(paste('Not yet eligible:',j$job_id))
 if(j$build!='GRCh37'||j$input1_complete_window!='YES'||j$input2_complete_window!='YES')stop('Require matched build and both full tested cis windows')
 if(file.exists(paste0(pref,'.done')))stop('Existing completed coloc job: runner must checksum-audit and skip exact job externally')
 fwrite(data.table(file=c(args[1],j$input1,j$input2),md5=unname(tools::md5sum(c(args[1],j$input1,j$input2)))),paste0(pref,'.input_checksums.tsv'),sep='\t')
 zz1<-read_input(j$input1,j);zz2<-read_input(j$input2,j)
 z<-merge(zz1,zz2,by='SNP',suffixes=c('_1','_2'));z<-z[chr_1==chr_2 & bp_1==bp_2]
 complement<-function(x)chartr('ACGT','TGCA',x)
 same<-z$A1_1==z$A1_2&z$A2_1==z$A2_2
 flip<-z$A1_1==z$A2_2&z$A2_1==z$A1_2
 c_same<-z$A1_1==complement(z$A1_2)&z$A2_1==complement(z$A2_2)
 c_flip<-z$A1_1==complement(z$A2_2)&z$A2_1==complement(z$A1_2)
 z[,alignment:=fifelse(same,'SAME',fifelse(flip,'FLIP',fifelse(c_same,'COMPLEMENT',fifelse(c_flip,'COMPLEMENT_FLIP','MISMATCH'))))]
 z[alignment%in%c('FLIP','COMPLEMENT_FLIP'),beta_2:=-beta_2]
 z<-z[alignment!='MISMATCH'];setorder(z,bp_1,SNP)
 align_qc[[k]]<-data.table(job_id=j$job_id,n1_after_QC=nrow(zz1),n2_after_QC=nrow(zz2),n_overlap=nrow(z),n_flipped=sum(z$alignment%in%c('FLIP','COMPLEMENT_FLIP')),palindrome_policy='ALL_REMOVED_MAF_NOT_EAF')
 fwrite(z,paste0(pref,'.harmonized.tsv.gz'),sep='\t')
 if(nrow(z)<100)stop('Insufficient overlapping full-window SNPs; no PP4 claim')
 d1<-make_dataset(z,1,j);d2<-make_dataset(z,2,j)
 for(p12 in c(1e-6,1e-5,5e-5)){
  rr<-coloc.abf(d1,d2,p1=1e-4,p2=1e-4,p12=p12)
  sm<-as.list(rr$summary)
  one<-data.table(job_id=j$job_id,gene=j$gene,qtl_panel=j$qtl_panel,axis=j$axis,comparison=j$comparison,n_SNP=nrow(z),p1=1e-4,p2=1e-4,p12=p12,PP0=sm[['PP.H0.abf']],PP1=sm[['PP.H1.abf']],PP2=sm[['PP.H2.abf']],PP3=sm[['PP.H3.abf']],PP4=sm[['PP.H4.abf']],lead1=z$SNP[which.max(abs(z$beta_1/z$SE_1))],lead2=z$SNP[which.max(abs(z$beta_2/z$SE_2))],method='coloc.abf',status='SUCCESS')
  sensitivity[[length(sensitivity)+1L]]<-one
  if(p12==1e-5){
   one[,strength:=fifelse(PP4>=.8,'STRONG',fifelse(PP4>=.5,'MODERATE','WEAK'))];result[[length(result)+1L]]<-one
   vr<-as.data.table(rr$results);setorder(vr,-SNP.PP.H4);vr[,H4_conditional_cumulative:=cumsum(SNP.PP.H4)]
   n95<-which(vr$H4_conditional_cumulative>=.95)[1];if(is.na(n95))n95<-nrow(vr)
   conditional_sets[[j$job_id]]<-vr$snp[seq_len(n95)]
   fwrite(vr,paste0(pref,'.variant_posteriors.tsv.gz'),sep='\t')
  }
 }
 # SuSiE is performed only on explicitly aligned and QC-approved LD matrices.
 if(!is.null(j$ld1)&&!is.null(j$ld2)&&!is.na(j$ld1)&&!is.na(j$ld2)&&j$ld_qc_pass=='YES'){
  if(is.null(j$ld1_alleles)||is.null(j$ld2_alleles))stop('LD allele manifests required, not just matching SNP names')
  for(am in c(j$ld1_alleles,j$ld2_alleles)){
   aa<-fread(am);stopifnot(all(c('SNP','A1','A2')%in%names(aa)),!anyDuplicated(aa$SNP))
   aa<-aa[match(z$SNP,SNP)]
   if(anyNA(aa$SNP)||any(aa$A1!=z$A1_1)|any(aa$A2!=z$A2_1))stop('LD must be signed to the common harmonized effect alleles')
  }
  ld1<-readRDS(j$ld1);ld2<-readRDS(j$ld2)
  stopifnot(all(z$SNP%in%rownames(ld1)),all(z$SNP%in%rownames(ld2)))
  d1$LD<-ld1[z$SNP,z$SNP];d2$LD<-ld2[z$SNP,z$SNP]
  if(max(abs(d1$LD-t(d1$LD)))>1e-6||max(abs(d2$LD-t(d2$LD)))>1e-6)stop('Asymmetric LD')
  if(length(unique(d1$N))!=1L||length(unique(d2$N))!=1L)stop('SuSiE requires separately audited effective sample-size model; do not silently collapse variable N')
  d1$N<-unique(d1$N);d2$N<-unique(d2$N)
  ss1<-runsusie(d1,maxit=200,coverage=.95);ss2<-runsusie(d2,maxit=200,coverage=.95)
  saveRDS(ss1,paste0(pref,'.susie1.rds'));saveRDS(ss2,paste0(pref,'.susie2.rds'))
  converged<-isTRUE(ss1$converged)&&isTRUE(ss2$converged)
  if(converged){sc<-coloc.susie(ss1,ss2);if(!is.null(sc$summary))fwrite(as.data.table(sc$summary),paste0(pref,'.coloc_susie.tsv'),sep='\t')}
  susie_qc[[k]]<-data.table(job_id=j$job_id,status=if(converged)'SUCCESS' else 'FAILED_CONVERGENCE',cs1=length(ss1$sets$cs),cs2=length(ss2$sets$cs))
 }else susie_qc[[k]]<-data.table(job_id=j$job_id,status='NOT_STARTED_REQUIRES_ALIGNED_LD',cs1=NA_integer_,cs2=NA_integer_)
 writeLines('0',paste0(pref,'.exitcode'));writeLines('SUCCESS_ABF',paste0(pref,'.done'))
}
fwrite(rbindlist(result,fill=TRUE),file.path(outdir,'TARGETED_COLOC_ALL.tsv'),sep='\t')
fwrite(rbindlist(sensitivity,fill=TRUE),file.path(outdir,'TARGETED_COLOC_PRIOR_SENSITIVITY.tsv'),sep='\t')
fwrite(rbindlist(align_qc,fill=TRUE),file.path(outdir,'TARGETED_COLOC_ALIGNMENT_QC.tsv'),sep='\t')
fwrite(rbindlist(susie_qc,fill=TRUE),file.path(outdir,'TARGETED_COLOC_SUSIE_QC.tsv'),sep='\t')
all_result<-rbindlist(result,fill=TRUE);triples<-list()
for(g in unique(all_result$gene))for(ax in unique(all_result[gene==g,axis]))for(qp in unique(all_result[gene==g & axis==ax,qtl_panel])){
 rr<-all_result[gene==g & axis==ax & qtl_panel==qp];need<-c('LDH_QTL','PAIN_QTL','LDH_PAIN')
 if(all(need%in%rr$comparison)&&!anyDuplicated(rr$comparison)){
  rr<-rr[match(need,comparison)];ss<-Reduce(intersect,conditional_sets[rr$job_id])
  triples[[length(triples)+1L]]<-data.table(gene=g,qtl_panel=qp,axis=ax,PP4_LDH_QTL=rr$PP4[1],PP4_PAIN_QTL=rr$PP4[2],PP4_LDH_PAIN=rr$PP4[3],n_shared_H4_conditional95_SNPs=length(ss),shared_SNPs=paste(ss,collapse=';'),TRIPLE_CONVERGENCE=if(all(rr$PP4>=.5)&&length(ss)>0)'PAIRWISE_SUPPORTED_REQUIRES_FINEMAP_REVIEW' else 'NOT_SUPPORTED',set_definition='ABF_H4_CONDITIONAL_95P_NOT_SUSIE_CREDIBLE_SET')
 }
}
if(length(triples))fwrite(rbindlist(triples),file.path(outdir,'TRIPLE_CONVERGENCE_SUMMARY.tsv'),sep='\t')
writeLines(capture.output(sessionInfo()),file.path(outdir,'software_versions.tsv'))
