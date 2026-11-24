# Mother-study SCAVENGE algorithm with explicit data gates and reproducible age block.
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=5)stop("Usage: Rscript replan_spinal_scavenge.R SE.rds PIP_hg38.bed output_dir seed isolated_R_library")
.libPaths(c(args[5],.libPaths()));seed<-as.integer(args[4]);set.seed(seed)
suppressPackageStartupMessages({library(SCAVENGE);library(chromVAR);library(gchromVAR);library(SummarizedExperiment);library(Matrix);library(data.table);library(BSgenome.Hsapiens.UCSC.hg38)})
BiocParallel::register(BiocParallel::SerialParam())
stopifnot(as.character(packageVersion("SCAVENGE"))=="1.0.2")
od<-args[3];dir.create(od,recursive=TRUE,showWarnings=FALSE)
se<-readRDS(args[1]);stopifnot(inherits(se,"RangedSummarizedExperiment"),all(genome(rowRanges(se))%in%c("hg38",NA)),"cell_cluster"%in%colnames(colData(se)),"donor"%in%colnames(colData(se)))
pip<-fread(args[2],header=FALSE);if(ncol(pip)<5||!nrow(pip))stop("True fine-mapped PIP BED required, >=5 columns")
if(any(!is.finite(pip[[5]])|pip[[5]]<=0|pip[[5]]>1))stop("PIP values not in (0,1]")
qc<-list(seed=seed,input_cells=ncol(se),input_peaks=nrow(se),input_donors=length(unique(se$donor)))
se<-se[rowSums(assay(se,"counts"))>=10,colSums(assay(se,"counts"))>1000]
se<-se[,!se$cell_cluster%in%c("GLIALPROG","OUT")]
cells<-unlist(lapply(split(colnames(se),se$cell_cluster),function(v)if(length(v)<=700)v else sample(v,700)),use.names=FALSE)
# Exactly the mother's Adult block; insufficient peaks are a reported deviation.
np<-min(104791,nrow(se));peaks<-sample(rownames(se),np);se<-se[peaks,cells]
qc$peak_subsample_target<-104791;qc$peak_subsample_actual<-np;qc$selected_cells<-ncol(se);qc$deviation<-if(np<104791)"ALL_AVAILABLE_PEAKS_LT_TARGET" else "NONE"
se<-addGCBias(se,genome=BSgenome.Hsapiens.UCSC.hg38)
se<-se[rowSums(assay(se,"counts"))>0,]
# Installed official importer uses one-based makeGRangesFromDataFrame internally.
# Convert standard BED zero-based start to the importer's expected one-based point.
pip_onebased<-copy(pip);pip_onebased[[2]]<-pip_onebased[[2]]+1L
stopifnot(all(pip_onebased[[2]]==pip_onebased[[3]]))
point_file<-file.path(od,"PIP_ONEBASED_FOR_GCHROMVAR.tsv")
fwrite(pip_onebased,point_file,sep="\t",col.names=FALSE)
weights<-importBedScore(rowRanges(se),point_file,colidx=5)
point_gr<-GenomicRanges::GRanges(pip_onebased[[1]],IRanges::IRanges(pip_onebased[[2]],pip_onebased[[3]]))
hits<-GenomicRanges::findOverlaps(rowRanges(se),point_gr);expected_weights<-numeric(nrow(se))
point_sums<-tapply(pip_onebased[[5]][S4Vectors::subjectHits(hits)],S4Vectors::queryHits(hits),sum)
expected_weights[as.integer(names(point_sums))]<-point_sums
stopifnot(max(abs(as.numeric(assay(weights,"weights"))-expected_weights))<1e-12)
qc$BED_zero_to_importer_one_based<-"VERIFIED_EXACT_POINT_WEIGHTS"
stopifnot(inherits(weights,"SummarizedExperiment"),"weights"%in%assayNames(weights))
weight_values<-assay(weights,"weights")
stopifnot(nrow(weight_values)==nrow(se),ncol(weight_values)==1L,all(is.finite(weight_values)))
qc$positive_weight_peaks<-sum(weight_values>0)
if(!any(weight_values>0))stop("No fine-mapped variants overlap selected snATAC peaks")
bg<-getBackgroundPeaks(se,niterations=50)
dev<-computeWeightedDeviations(object=se,weights=weights,background_peaks=bg)
z<-as.numeric(t(assays(dev)[["z"]]));if(length(z)!=ncol(se)||any(!is.finite(z)))stop("Nonfinite or mismatched chromVAR Z")
sidx<-seedindex(z_score=z,percent_cut=.05);sf<-cal_scalefactor(z_score=z,percent_cut=.01)
fwrite(data.table(cell_id=colnames(se),z=z,seed_cell=sidx),file.path(od,"CHROMVAR_Z_AND_SEEDS.tsv.gz"),sep="\t",compress="gzip")
if(!any(sidx)||!is.finite(sf)||sf<=0)stop("NO_POSITIVE_SEED_OR_SCALE: preserve real chromVAR result; no forced seeds")
tf<-tfidf(bmat=assay(se,"counts"),mat_binary=TRUE,TF=TRUE,log_TF=TRUE)
lsi<-do_lsi(mat=tf,dims=30);knn<-getmutualknn(lsimat=lsi,num_k=30)
np_score<-randomWalk_sparse(intM=knn,queryCells=rownames(knn)[sidx],gamma=.05)
stopifnot(length(np_score)==ncol(se),all(is.finite(np_score)));connected<-np_score!=0
if(!any(connected)||diff(range(np_score[connected]))==0)stop("UNIDENTIFIABLE_NETWORK_RANGE: no forced TRS")
trs<-rep(NA_real_,length(np_score));trs[connected]<-max_min_scale(capOutlierQuantile(x=np_score[connected],q_ceiling=.95))*sf
stopifnot(all(is.finite(trs[connected])))
cells<-as.data.frame(colData(se));cells$cell_id<-colnames(se);cells$z<-z;cells$seed_cell<-as.logical(sidx);cells$network_score<-np_score;cells$connected<-connected;cells$TRS<-trs;cells$enriched<-connected&trs>2;cells$random_seed<-seed
fwrite(cells,file.path(od,"SCAVENGE_CELL_LEVEL.tsv.gz"),sep="\t",compress="gzip")
# Cell-level hypergeometric enrichment mirrors SCAVENGE; donor numbers are disclosed.
cd<-as.data.table(cells[connected,]);N<-nrow(cd);K<-sum(cd$enriched)
sm<-cd[,.(n_cells=.N,n_enriched=sum(enriched),n_donors=uniqueN(donor),median_TRS=median(TRS),mean_TRS=mean(TRS)),by=cell_cluster]
allct<-as.data.table(cells)[,.(n_reference_donors=uniqueN(donor),n_input_cells=.N),by=cell_cluster];sm<-merge(allct,sm,by="cell_cluster",all.x=TRUE)
sm[is.na(n_cells),`:=`(n_cells=0L,n_enriched=0L,n_donors=0L)];sm[,test_status:=ifelse(n_cells>0,"TESTED_CONNECTED_CELLS","NOT_TESTED_NO_CONNECTED_CELLS")]
sm[,P:=ifelse(n_cells>0,phyper(n_enriched-1,K,N-K,n_cells,lower.tail=FALSE),NA_real_)];sm[,FDR:=p.adjust(P,"BH",n=nrow(sm))];sm[,random_seed:=seed]
fwrite(sm,file.path(od,"SCAVENGE_CELLTYPE_SUMMARY.tsv"),sep="\t")
ds<-cd[,.(n_cells=.N,n_enriched=sum(enriched),median_TRS=median(TRS)),by=.(donor,cell_cluster)];fwrite(ds,file.path(od,"SCAVENGE_DONOR_DESCRIPTIVE.tsv"),sep="\t")
qc$connected_cells<-N;qc$unconnected_cells<-sum(!connected);qc$enriched_cells<-K
writeLines(unlist(Map(function(k,v)paste(k,v,sep="\t"),names(qc),qc)),file.path(od,"SCAVENGE_QC.tsv"))
writeLines(capture.output(sessionInfo()),file.path(od,"software_versions.txt"))
writeLines("0",file.path(od,"scavenge.exitcode"));writeLines(format(Sys.time(),tz="UTC"),file.path(od,"scavenge.done"))
