library(data.table);library(coloc);library(susieR)
stopifnot(as.character(packageVersion('coloc'))=='5.2.2',as.character(packageVersion('susieR'))=='0.12.35')
root<-commandArgs(TRUE)[1];base<-file.path(root,'08_smr_heidi_coloc/multisignal_ld503_audit_zero')
for(gene in c('MST1','COL7A1')){
 d<-file.path(base,gene);o<-file.path(d,'susie_N_sensitivity');dir.create(o,showWarnings=FALSE)
 z<-fread(file.path(d,'THREE_TRAIT_COMMON_FULL_WINDOW.tsv'));qc<-fread(file.path(d,'RSS_LD_N_DIAGNOSTIC.tsv'));stopifnot(all(qc$s<=.1))
 x<-as.matrix(read.delim(gzfile(file.path(d,'effect_allele_genotypes.tsv.gz')),check.names=FALSE));stopifnot(identical(colnames(x),z$SNP),nrow(x)==503);ld<-cor(x);dimnames(ld)<-list(z$SNP,z$SNP)
 spec<-list(input_md5=tools::md5sum(file.path(d,c('THREE_TRAIT_COMMON_FULL_WINDOW.tsv','effect_allele_genotypes.tsv.gz'))),versions=c(coloc=as.character(packageVersion('coloc')),susieR=as.character(packageVersion('susieR'))),N=lapply(c('LDH','CPI','QTL'),function(t)unique(c(min(z[[paste0('N_',t)]]),median(z[[paste0('N_',t)]]),max(z[[paste0('N_',t)]])))),L=10,coverage=.95,maxit=200,estimate_residual_variance=FALSE,min_abs_corr=.5,scope='Explicit scalar-N sensitivity, not an exact per-SNP-N multi-signal likelihood; no best-N selection; primary ABF retains original vector N')
 lock<-file.path(o,'BEFORE_FIT_LOCK.rds');if(file.exists(lock))stopifnot(identical(readRDS(lock),spec)) else saveRDS(spec,lock)
 # runsusie uses z, LD and N; its generic quant dataset validator additionally
 # requires source MAF when sdY is unknown. Do not invent a phenotype SD.
 panel<-if(gene=='MST1')'Nerve_Tibial' else 'Brain_Cortex'
 qid<-paste0(gene,'_CPI_',panel,'_LDH_QTL');cid<-paste0(gene,'_CPI_COMMON_WINDOW_LDH_PAIN')
 qp<-file.path(root,'08_smr_heidi_coloc/new_candidates_coloc_zero',qid,paste0(qid,'.harmonized.tsv.gz'));cp<-file.path(root,'08_smr_heidi_coloc/new_candidates_coloc_zero',cid,paste0(cid,'.harmonized.tsv.gz'))
 qq<-fread(qp);cc<-fread(cp);qq<-qq[match(z$SNP,SNP)];cc<-cc[match(z$SNP,SNP)];stopifnot(identical(qq$SNP,z$SNP),identical(cc$SNP,z$SNP))
 maf<-list(LDH=qq$MAF_1,CPI=cc$MAF_2,QTL=qq$MAF_2);stopifnot(all(vapply(maf,function(v)all(is.finite(v)&v>0&v<=.5),logical(1))))
 fwrite(data.table(source=c(qp,cp),md5=unname(tools::md5sum(c(qp,cp))),scope='Original source MAF for generic validator; no sdY invention; RSS z/R/N inputs and completed fits unchanged'),file.path(o,'MAF_VALIDATOR_PROVENANCE.tsv'),sep='\t')
 fits<-list();fitrows<-list()
 for(t in c('LDH','CPI','QTL')){
  nn<-z[[paste0('N_',t)]];values<-unique(c(min(nn),median(nn),max(nn)))
  for(n in values){
   id<-paste(t,n,sep='_');p<-file.path(o,paste0(id,'.rds'));sig<-list(spec=spec,trait=t,N=n);mp<-paste0(p,'.signature.rds')
   if(file.exists(p)){stopifnot(file.exists(mp),identical(readRDS(mp),sig));ss<-readRDS(p);cat('REUSE_FIT',gene,id,'\n')}
   else{
    dd<-list(beta=z[[paste0('beta_',t)]],varbeta=z[[paste0('SE_',t)]]^2,snp=z$SNP,position=z$BP,N=n,LD=ld,MAF=maf[[t]],type=if(t=='LDH')'cc' else 'quant');if(t=='LDH')dd$s<-35290/388514
    ss<-runsusie(dd,maxit=200,repeat_until_convergence=FALSE,L=10,coverage=.95,min_abs_corr=.5,estimate_residual_variance=FALSE)
    stopifnot(isTRUE(ss$converged));saveRDS(ss,p);saveRDS(sig,mp)
   }
   fits[[id]]<-ss;fwrite(data.table(SNP=names(ss$pip),PIP=ss$pip),file.path(o,paste0(id,'.PIP.tsv')),sep='\t')
   cs<-rbindlist(lapply(seq_along(ss$sets$cs),function(i)data.table(signal=names(ss$sets$cs)[i],SNP=names(ss$sets$cs[[i]]))),fill=TRUE);if(nrow(cs))fwrite(cs,file.path(o,paste0(id,'.CS95.tsv')),sep='\t')
   fitrows[[length(fitrows)+1]]<-data.table(gene=gene,trait=t,N=n,n_SNP=nrow(z),converged=ss$converged,n_CS=length(ss$sets$cs),N_scope=if(length(unique(nn))==1)'EXACT_CONSTANT_N' else 'SCALAR_N_RANGE_SENSITIVITY')
  }
 }
 ans<-list()
 for(pair in list(c('LDH','QTL'),c('CPI','QTL'),c('LDH','CPI'))){
  aa<-grep(paste0('^',pair[1],'_'),names(fits),value=TRUE);bb<-grep(paste0('^',pair[2],'_'),names(fits),value=TRUE)
  for(a in aa)for(b in bb)for(p12 in c(1e-6,1e-5,5e-5)){
   res<-coloc.susie(fits[[a]],fits[[b]],p1=1e-4,p2=1e-4,p12=p12)
   if(!is.null(res$summary)){rr<-as.data.table(res$summary);rr[,`:=`(gene=gene,fit1=a,fit2=b,p12=p12,status='SUCCESS_CONDITIONAL_N_SENSITIVITY')];ans[[length(ans)+1]]<-rr}
   else ans[[length(ans)+1]]<-data.table(gene=gene,fit1=a,fit2=b,p12=p12,status='NO_CREDIBLE_SET_NOT_PP4_ZERO')
  }
 }
 fwrite(rbindlist(fitrows,fill=TRUE),file.path(o,'SUSIE_FIT_QC.tsv'),sep='\t');fwrite(rbindlist(ans,fill=TRUE),file.path(o,'COLOC_SUSIE_ALL_N_PRIORS.tsv'),sep='\t');writeLines('SUCCESS_N_SENSITIVITY_ONLY_NOT_EXACT_PER_SNP_N',file.path(o,'COMPLETE.done'));cat(gene,'ALL_N_SENSITIVITY_COMPLETE\n');rm(x,ld,fits);gc()
}
