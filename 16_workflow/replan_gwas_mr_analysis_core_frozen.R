root <- '/analysis/projects/JTM_recent_highIF_direct_clone/P2_LDH_shared_architecture/11_NatCommun_strict_clone_Pro'
.libPaths(c(file.path(root,'16_workflow/envs/mr_r43_lib'),.libPaths()))
suppressPackageStartupMessages(library(TwoSampleMR))
stopifnot(as.character(packageVersion('TwoSampleMR'))=='0.7.9')
args <- commandArgs(trailingOnly=TRUE);stage <- if(length(args))args[1] else 'core'
out <- file.path(root,'08_mr/LDH_CPI');lock <- jsonlite::read_json(file.path(root,'16_workflow/replan_gwas_mr_lock.json'),simplifyVector=TRUE)
if(file.exists(file.path(root,'00_admin/PAUSE_REQUESTED')))stop('PAUSE_REQUESTED')
put <- function(x,name) write.table(x,file.path(out,name),sep='\t',row.names=FALSE,quote=FALSE,na='NA')
ids <- c(LDH='FINNGEN_R12_LDH',CPI='UKB_CPI_WANG2025')
harm <- function(direction) {
  x <- read.delim(file.path(out,paste0(direction,'.instrument_pairs.tsv')),check.names=FALSE)
  stopifnot(!anyDuplicated(x$SNP),all(x$status=='READY_FOR_TWOSAMPLEMR_HARMONISE'),all(x$sample_overlap=='LOW'))
  traits <- strsplit(direction,'_to_')[[1]]
  d <- lapply(c('exposure','outcome'),function(t){
    a <- data.frame(SNP=x$SNP,beta=x[[paste0('BETA_',t)]],se=x[[paste0('SE_',t)]],p=x[[paste0('P_',t)]],eaf=x[[paste0('EAF_',t)]],a1=x[[paste0('A1_',t)]],a2=x[[paste0('A2_',t)]],N=x[[paste0('N_',t)]],phenotype=traits[if(t=='exposure')1 else 2])
    format_data(a,type=t,phenotype_col='phenotype',snp_col='SNP',beta_col='beta',se_col='se',pval_col='p',eaf_col='eaf',effect_allele_col='a1',other_allele_col='a2',samplesize_col='N')
  })
  d[[1]]$id.exposure <- ids[[traits[1]]];d[[2]]$id.outcome <- ids[[traits[2]]]
  h <- harmonise_data(d[[1]],d[[2]],action=lock$harmonise_action)
  h$F_stat <- (h$beta.exposure/h$se.exposure)^2
  h$sample_overlap <- 'LOW';h$direction <- direction
  put(h,paste0(direction,'.harmonised_all.tsv'))
  h
}
directions <- c('LDH_to_CPI','CPI_to_LDH')
if(stage=='core') {
  all <- list();hets <- list();pleio <- list();qc <- list()
  for(direction in directions) {
    set.seed(lock$seed);h <- harm(direction);d <- h[h$mr_keep,,drop=FALSE]
    stopifnot(nrow(d)>0,all(is.finite(d$beta.exposure)),all(is.finite(d$beta.outcome)),all(d$se.exposure>0),all(d$se.outcome>0))
    parameters <- default_parameters();parameters$nboot <- lock$weighted_median_bootstrap
    methods <- if(nrow(d)==1)'mr_wald_ratio' else c('mr_ivw_mre','mr_weighted_median','mr_egger_regression')
    res <- mr(d,parameters=parameters,method_list=methods)
    res$direction <- direction;res$primary <- res$method %in% c('Inverse variance weighted (multiplicative random effects)','Wald ratio')
    res$P_Bonferroni5 <- pmin(1,res$pval*5);res$P_Bonferroni10_sensitivity <- pmin(1,res$pval*10)
    res$primary_alpha <- .05/5;res$passes_primary_threshold <- res$primary & res$pval<.05/5
    res$OR <- if(direction=='CPI_to_LDH')exp(res$b)else NA_real_
    res$CI_lower <- res$b-1.96*res$se;res$CI_upper <- res$b+1.96*res$se
    res$OR_CI_lower <- if(direction=='CPI_to_LDH')exp(res$CI_lower)else NA_real_
    res$OR_CI_upper <- if(direction=='CPI_to_LDH')exp(res$CI_upper)else NA_real_
    res$effect_scale <- if(direction=='LDH_to_CPI')lock$LDH_to_CPI_effect else lock$CPI_to_LDH_effect
    res$sample_overlap <- 'LOW';res$Steiger_status <- lock$Steiger
    put(res,paste0(direction,'.MR_ALL.tsv'));all[[direction]] <- res
    if(nrow(d)>=2) {
      q <- TwoSampleMR::mr_ivw_mre(d$beta.exposure,d$beta.outcome,d$se.exposure,d$se.outcome,parameters)
      hets[[direction]] <- data.frame(direction=direction,method='mr_ivw_mre',Q=q$Q,Q_df=q$Q_df,Q_P=q$Q_pval)
      loo <- do.call(rbind,lapply(c(d$SNP,'ALL_INSTRUMENTS'),function(snp){
        z <- if(snp=='ALL_INSTRUMENTS')d else d[d$SNP!=snp,,drop=FALSE]
        r <- if(nrow(z)==1)TwoSampleMR::mr_wald_ratio(z$beta.exposure,z$beta.outcome,z$se.exposure,z$se.outcome,parameters)else TwoSampleMR::mr_ivw_mre(z$beta.exposure,z$beta.outcome,z$se.exposure,z$se.outcome,parameters)
        data.frame(direction=direction,omitted_SNP=snp,nsnp=nrow(z),beta=r$b,SE=r$se,P=r$pval,method='Official TwoSampleMR mr_ivw_mre called after each omission')
      }))
      put(loo,paste0(direction,'.LEAVE_ONE_OUT.tsv'))
    }
    if(nrow(d)>=3){p <- mr_pleiotropy_test(d);p$direction <- direction;pleio[[direction]] <- p}
    qc[[direction]] <- data.frame(direction=direction,n_input=nrow(h),n_harmonised_keep=nrow(d),n_removed=sum(!h$mr_keep),F_min=min(d$F_stat),F_mean=mean(d$F_stat),N_exposure_min=min(d$samplesize.exposure),N_exposure_max=max(d$samplesize.exposure),N_outcome_min=min(d$samplesize.outcome),N_outcome_max=max(d$samplesize.outcome),Steiger_status=lock$Steiger)
    saveRDS(d,file.path(out,paste0(direction,'.harmonised_keep.rds')))
  }
  put(do.call(rbind,all),'MR_CORE_ALL.tsv');put(do.call(rbind,hets),'MR_HETEROGENEITY.tsv');put(do.call(rbind,pleio),'MR_EGGER_INTERCEPT.tsv');put(do.call(rbind,qc),'MR_ANALYSIS_QC.tsv')
  capture.output(sessionInfo(),file=file.path(out,'MR_CORE_SESSION_INFO.txt'))
  jsonlite::write_json(list(status='SUCCESS_CORE_ESTIMATION',Steiger=lock$Steiger,full_module_status='PARTIAL_UNTIL_MR_PRESSO_AND_OTHER_AXES'),file.path(out,'MR_CORE_STATUS.json'),auto_unbox=TRUE,pretty=TRUE)
} else if(stage %in% directions) {
  direction <- stage;d <- readRDS(file.path(out,paste0(direction,'.harmonised_keep.rds')))
  row.names(d) <- d$SNP
  if(nrow(d)<4)stop('MR-PRESSO requires at least4 IVs; record NOT_APPLICABLE instead of silently substituting')
  p <- MRPRESSO::mr_presso(BetaOutcome='beta.outcome',BetaExposure='beta.exposure',SdOutcome='se.outcome',SdExposure='se.exposure',data=as.data.frame(d),OUTLIERtest=TRUE,DISTORTIONtest=TRUE,SignifThreshold=.05,NbDistribution=10000,seed=lock$seed)
  saveRDS(p,file.path(out,paste0(direction,'.MR_PRESSO.raw.rds')))
  capture.output(str(p),print(p),file=file.path(out,paste0(direction,'.MR_PRESSO.raw.txt')))
  main <- p[['Main MR results']];main$direction <- direction;main$use <- 'SENSITIVITY_ONLY_NOT_PRIMARY_IVW';put(main,paste0(direction,'.MR_PRESSO_MAIN.tsv'))
  gt <- p[['MR-PRESSO results']][['Global Test']]
  jsonlite::write_json(list(direction=direction,status='SUCCESS',n_IV=nrow(d),NbDistribution=10000,seed=lock$seed,global_test=gt,outlier_SNP_mapping=d$SNP,note='Raw canonical MR-PRESSO finite-permutation P retained including boundary strings; corrected MR is sensitivity only'),file.path(out,paste0(direction,'.MR_PRESSO_STATUS.json')),auto_unbox=TRUE,pretty=TRUE)
} else stop('Unknown stage')
