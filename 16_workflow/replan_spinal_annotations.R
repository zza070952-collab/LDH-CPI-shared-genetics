args<-commandArgs(trailingOnly=TRUE);if(length(args)!=3)stop("Usage: Rscript annotations.R SE_Adult.rds hg38ToHg19.chain outputdir")
suppressPackageStartupMessages({library(SummarizedExperiment);library(GenomicRanges);library(rtracklayer);library(Matrix);library(data.table)})
se<-readRDS(args[1]);od<-args[3];dir.create(od,recursive=TRUE,showWarnings=FALSE)
stopifnot(all(se$analysis=="AUTHOR_ADULT_PRIMARY"),all(se$annotation_source=="UCSC author Lineage"))
keep<-colSums(assay(se,"counts"))>1000;cd<-as.data.table(as.data.frame(colData(se)));cd$eligible_library_counts<-keep;cd$library_counts<-colSums(assay(se,"counts"));fwrite(cd,file.path(od,"SLDSC_CELL_QC.tsv.gz"),sep="\t",compress="gzip")
se<-se[,keep];mat<-assay(se,"counts");gr<-rowRanges(se);auto<-as.character(seqnames(gr))%in%paste0("chr",1:22)
ch<-import.chain(args[2]);mapped<-liftOver(gr,ch);nm<-elementNROWS(mapped);unique<-nm==1;g19<-unlist(mapped[unique],use.names=FALSE);origin<-which(unique);ok<-as.character(seqnames(g19))%in%paste0("chr",1:22)&width(g19)>=.95*width(gr)[origin]&width(g19)<=2*width(gr)[origin]&width(g19)>0
valid<-rep(FALSE,length(gr));valid[origin[ok]]<-TRUE
mapped_chr<-rep(NA_character_,length(gr));mapped_start<-mapped_end<-rep(NA_integer_,length(gr));mapped_chr[origin]<-as.character(seqnames(g19));mapped_start[origin]<-start(g19)-1L;mapped_end[origin]<-end(g19)
mp<-data.table(peak=names(gr),chr38=as.character(seqnames(gr)),start38=start(gr)-1L,end38=end(gr),n_mapped_intervals=nm,unique_valid=valid,chr19=mapped_chr,start19=mapped_start,end19=mapped_end)
fwrite(mp,file.path(od,"SLDSC_PEAK_LIFTOVER_ALL.tsv.gz"),sep="\t",compress="gzip")
res<-list();donors<-list()
for(ct in sort(unique(as.character(se$cell_cluster)))){
  ix<-which(se$cell_cluster==ct);tot<-rowSums(mat[,ix,drop=FALSE]);present<-tot>=10;safe<-gsub("[^A-Za-z0-9_]","_",ct)
  b38<-data.table(chr=as.character(seqnames(gr))[present&auto],start=start(gr)[present&auto]-1L,end=end(gr)[present&auto],peak=names(gr)[present&auto]);setorder(b38,chr,start,end)
  fwrite(b38,file.path(od,paste0(safe,".hg38.bed")),sep="\t",col.names=FALSE)
  b19<-mp[present&auto&valid,.(chr=chr19,start=start19,end=end19,peak)];dup<-duplicated(b19,by=c("chr","start","end"));ndup<-sum(dup);b19<-b19[!dup];setorder(b19,chr,start,end);fwrite(b19,file.path(od,paste0(safe,".hg19.bed")),sep="\t",col.names=FALSE)
  dns<-unique(as.character(se$donor[ix]));ndetect<-integer(nrow(mat))
  for(dn in dns){ii<-which(se$cell_cluster==ct&se$donor==dn);dc<-rowSums(mat[,ii,drop=FALSE]);ndetect<-ndetect+(dc>=10);donors[[length(donors)+1]]<-data.table(cell_type=ct,donor=dn,n_cells=length(ii),matrix_library_counts=sum(dc),peaks_ge10_counts=sum(dc>=10),autosomal_peaks_ge10_counts=sum(dc>=10&auto))}
  fwrite(data.table(peak=names(gr),cell_type=ct,total_matrix_counts=tot,n_donors_peak_ge10=ndetect,n_donors=length(dns)),file.path(od,paste0(safe,".peak_donor_detection.tsv.gz")),sep="\t",compress="gzip")
  res[[length(res)+1]]<-data.table(annotation_id=safe,cell_type=ct,source="PMID37824647",scope="STRICT_CLONE_BENCHMARK_CORTEX",segment="cortex",n_donors=length(dns),n_cells=length(ix),peak_matrix_count_threshold=10,n_peaks_all=sum(present),n_peaks_autosomal_hg38=nrow(b38),n_unmapped=sum(present&nm==0),n_multimap=sum(present&nm>1),n_single_invalid=sum(present&unique&!valid),n_duplicate_mapped=ndup,n_peaks_hg19=nrow(b19),mapping_rate=nrow(b19)/max(1,sum(present&auto)),eligible_S_LDSC=ct!="OUT",bed_hg38=file.path(od,paste0(safe,".hg38.bed")),bed_hg19=file.path(od,paste0(safe,".hg19.bed")))
}
fwrite(rbindlist(res),file.path(od,"SLDSC_ANNOTATION_MANIFEST.tsv"),sep="\t");fwrite(rbindlist(donors),file.path(od,"SLDSC_ANNOTATION_DONOR_QC.tsv"),sep="\t")
writeLines(capture.output(sessionInfo()),file.path(od,"software_versions.txt"));writeLines(format(Sys.time(),tz="UTC"),file.path(od,"annotations.done"))
