args <- commandArgs(TRUE)
root <- args[1]
o <- file.path(root, '03_signal_resolution')
Sys.setenv(OMP_NUM_THREADS='4', OPENBLAS_NUM_THREADS='4', MKL_NUM_THREADS='4')
.libPaths(c(file.path(o,'R_library'), .libPaths()))
suppressPackageStartupMessages({library(data.table);library(coloc);library(susieR)})
stopifnot(as.character(packageVersion('coloc'))=='5.2.2',as.character(packageVersion('susieR'))=='0.12.35')
set.seed(20261103)
writeLines(capture.output(sessionInfo()),file.path(o,'CHR18_R_SESSION_INFO.txt'))
versions <- data.table(component=c('R','coloc','susieR','data.table'),version=c(as.character(getRversion()),as.character(packageVersion('coloc')),as.character(packageVersion('susieR')),as.character(packageVersion('data.table'))))
fwrite(versions,file.path(o,'CHR18_SOFTWARE_VERSIONS.tsv'),sep='\t')
if(file.exists(file.path(o,'CHR18_COMPLETE.done'))) {
  cat('REUSE_EXISTING_COMPLETED_CHR18_ANALYSIS\n')
  quit(status=0)
}
z <- fread(file.path(o,'CHR18_TWO_TRAIT_COMMON_FIXED_WINDOW.tsv'))
x <- as.matrix(read.delim(gzfile(file.path(o,'CHR18_effect_allele_genotypes.tsv.gz')),check.names=FALSE))
stopifnot(nrow(x)==503,identical(colnames(x),z$SNP),nrow(z)>100)
ld <- cor(x);dimnames(ld)<-list(z$SNP,z$SNP)
stopifnot(identical(rownames(ld),colnames(ld)),max(abs(ld-t(ld)))<1e-10,max(abs(diag(ld)-1))<1e-10,all(is.finite(ld)))
ee<-eigen(ld,symmetric=TRUE)
stopifnot(min(ee$values)>-1e-7)
attr(ld,'eigen')<-ee
fwrite(data.table(n_SNP=nrow(z),n_reference_samples=nrow(x),min_eigenvalue=min(ee$values),symmetry_max_error=max(abs(ld-t(ld))),diagonal_max_error=max(abs(diag(ld)-1)),LD_matrix='signed r, aligned to unchanged LDH effect alleles',memory_dense_matrix_bytes=as.numeric(object.size(ld))),file.path(o,'CHR18_LD_MATRIX_QC.tsv'),sep='\t')

# The complete source-derived min/median/max N family is fixed before diagnostics.
nvalues<-list(LDH=unique(c(min(z$N_LDH),median(z$N_LDH),max(z$N_LDH))),CPI=unique(c(min(z$N_CPI),median(z$N_CPI),max(z$N_CPI))))
diagnostics<-list()
for(t in c('LDH','CPI')) for(n in nvalues[[t]]) {
 zz<-z[[paste0('beta_',t)]]/z[[paste0('SE_',t)]]
 s<-estimate_s_rss(zz,ld,n=n)
 diagnostics[[length(diagnostics)+1]]<-data.table(trait=t,N=n,n_SNP=nrow(z),RSS_s=s,passes_fixed_s_0.1=is.finite(s)&&s<=.1,N_role=if(t=='LDH')'CONSTANT_SOURCE_N' else if(n==median(z$N_CPI))'PRIMARY_MEDIAN_N' else 'SOURCE_RANGE_SENSITIVITY',threshold_origin='existing project operational threshold, not a universal official cutoff')
}
qc<-rbindlist(diagnostics)
fwrite(qc,file.path(o,'CHR18_RSS_LD_N_DIAGNOSTIC.tsv'),sep='\t')
if(any(!qc$passes_fixed_s_0.1)) {
 fwrite(data.table(module='targeted_chr18',status='UNRESOLVED_LD',reason='At least one fixed trait/N RSS diagnostic failed pre-existing s<=0.1 threshold; no fit or allele/SNP tuning performed',n_input=nrow(z)),file.path(o,'CHR18_MODULE_STATUS.tsv'),sep='\t')
 writeLines('UNRESOLVED_LD_FIXED_RSS_QC',file.path(o,'CHR18_COMPLETE.done'))
 cat('CHR18_STOP_FIXED_RSS_QC\n')
 quit(status=0)
}

fits<-list();fitrows<-list();csrows<-list();piprows<-list()
for(t in c('LDH','CPI')) for(n in nvalues[[t]]) {
 id<-paste(t,n,sep='_');p<-file.path(o,paste0('CHR18_',id,'.rds'))
 spec<-list(input_md5=tools::md5sum(file.path(o,c('CHR18_TWO_TRAIT_COMMON_FIXED_WINDOW.tsv','CHR18_effect_allele_genotypes.tsv.gz','CHR18_RUN_CONFIG.yaml'))),trait=t,N=n,L=10,coverage=.95,min_abs_corr=.5,maxit=200,estimate_residual_variance=FALSE,versions=versions)
 lock<-paste0(p,'.signature.rds')
 if(file.exists(p)) {
  stopifnot(file.exists(lock),identical(readRDS(lock),spec));ss<-readRDS(p);cat('REUSE_FIT',id,'\n')
 } else {
  dd<-list(beta=z[[paste0('beta_',t)]],varbeta=z[[paste0('SE_',t)]]^2,snp=z$SNP,position=z$BP,N=n,LD=ld,MAF=z[[paste0('MAF_',t)]],type=if(t=='LDH')'cc' else 'quant')
  if(t=='LDH')dd$s<-35290/388514
  ss<-runsusie(dd,maxit=200,repeat_until_convergence=FALSE,L=10,coverage=.95,min_abs_corr=.5,estimate_residual_variance=FALSE)
  saveRDS(ss,p);saveRDS(spec,lock)
 }
 fits[[id]]<-ss
 fitrows[[length(fitrows)+1]]<-data.table(trait=t,fit_id=id,N=n,n_SNP=nrow(z),converged=isTRUE(ss$converged),n_CS=length(ss$sets$cs),iterations=ss$niter,primary_model=t=='LDH'||n==median(z$N_CPI),N_scope=if(t=='LDH')'CONSTANT_SOURCE_N' else 'SCALAR_N_SENSITIVITY_NOT_EXACT_PER_SNP_LIKELIHOOD')
 piprows[[length(piprows)+1]]<-data.table(trait=t,fit_id=id,N=n,SNP=z$SNP,BP=z$BP,PIP=as.numeric(ss$pip))
 if(length(ss$sets$cs))for(i in seq_along(ss$sets$cs)) {
  cs<-ss$sets$cs[[i]]; snps<-z$SNP[cs]
  pur<-ss$sets$purity[i,,drop=FALSE]
  csrows[[length(csrows)+1]]<-data.table(trait=t,fit_id=id,N=n,signal=names(ss$sets$cs)[i],SNP=snps,BP=z$BP[cs],PIP=as.numeric(ss$pip[cs]),coverage=ss$sets$coverage[i],purity_min_abs_corr=pur$min.abs.corr,purity_mean_abs_corr=pur$mean.abs.corr,purity_median_abs_corr=pur$median.abs.corr)
 }
}
fqc<-rbindlist(fitrows,fill=TRUE);fwrite(fqc,file.path(o,'CHR18_SUSIE_FIT_QC.tsv'),sep='\t')
fwrite(rbindlist(piprows),file.path(o,'CHR18_SUSIE_PIP_ALL.tsv'),sep='\t')
if(length(csrows))fwrite(rbindlist(csrows,fill=TRUE),file.path(o,'CHR18_SUSIE_CREDIBLE_SETS_ALL.tsv'),sep='\t') else fwrite(data.table(status='NO_CREDIBLE_SET_NOT_PP4_ZERO'),file.path(o,'CHR18_SUSIE_CREDIBLE_SETS_ALL.tsv'),sep='\t')
all_coloc<-list()
for(a in grep('^LDH_',names(fits),value=TRUE))for(b in grep('^CPI_',names(fits),value=TRUE))for(p12 in c(1e-6,1e-5,5e-5)) {
 if(!isTRUE(fits[[a]]$converged)||!isTRUE(fits[[b]]$converged)) {
  all_coloc[[length(all_coloc)+1]]<-data.table(fit1=a,fit2=b,p1=1e-4,p2=1e-4,p12=p12,status='NOT_ESTIMABLE_NONCONVERGED_MODEL')
  next
 }
 out<-coloc.susie(fits[[a]],fits[[b]],p1=1e-4,p2=1e-4,p12=p12)
 if(!is.null(out$summary)) {
  q<-as.data.table(out$summary);q[,`:=`(fit1=a,fit2=b,p1=1e-4,p2=1e-4,p12=p12,status='SUCCESS')]
  all_coloc[[length(all_coloc)+1]]<-q
 } else all_coloc[[length(all_coloc)+1]]<-data.table(fit1=a,fit2=b,p1=1e-4,p2=1e-4,p12=p12,status='NO_CREDIBLE_SET_NOT_PP4_ZERO')
}
co<-rbindlist(all_coloc,fill=TRUE)
fwrite(co,file.path(o,'CHR18_COLOC_SUSIE_ALL_N_PRIORS.tsv'),sep='\t')
primary<-co[fit2==paste0('CPI_',median(z$N_CPI)) & p12==1e-5]
fwrite(primary,file.path(o,'CHR18_COLOC_SUSIE_PRIMARY.tsv'),sep='\t')
state<-if(all(fqc$converged))'COMPLETED' else 'UNRESOLVED_NONCONVERGENCE'
reason<-if(all(fqc$converged))'Fixed-window diagnostics passed and all four source-derived trait/N fits completed; all available signal pairs and priors retained' else 'At least one fixed model did not converge; no parameter tuning performed'
fwrite(data.table(module='targeted_chr18',status=state,reason=reason,n_input=nrow(z),n_fits=nrow(fqc),n_primary_signal_pairs=sum(primary$status=='SUCCESS')),file.path(o,'CHR18_MODULE_STATUS.tsv'),sep='\t')
writeLines(state,file.path(o,'CHR18_COMPLETE.done'))
cat('CHR18_TARGETED_ANALYSIS_COMPLETE',state,'\n')
