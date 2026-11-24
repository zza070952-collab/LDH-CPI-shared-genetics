#!/usr/bin/env Rscript

root <- "/analysis/projects/JTM_recent_highIF_direct_clone/P2_LDH_shared_architecture/11_NatCommun_strict_clone_Pro/18_targeted_LAVA_20260913"
.libPaths(c(file.path(root, "02_reference/R_library"), .libPaths()))
suppressPackageStartupMessages(library(LAVA))

deps <- c("LAVA", "data.table", "matrixsampling", "cpp11", "R.utils")
dep.versions <- data.frame(
  component=c("R", deps),
  version=c(as.character(getRversion()), vapply(deps, function(x) as.character(packageVersion(x)), character(1))),
  stringsAsFactors=FALSE
)
write.table(dep.versions, file.path(root, "00_admin/R_DEPENDENCY_VERSIONS.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

info.file <- file.path(root, "01_inputs/input.info.tsv")
overlap.primary <- file.path(root, "01_inputs/sample.overlap.primary.txt")
overlap.zero <- file.path(root, "01_inputs/sample.overlap.zero_sensitivity.txt")
ref.prefix <- file.path(root, "02_reference/analysis_reference/lava-ukb-v1.1")
locus.file <- file.path(root, "01_inputs/selected_blocks.loci")

capture_eval <- function(expr) {
  warns <- character()
  value <- tryCatch(
    withCallingHandlers(expr, warning=function(w) {
      warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning")
    }),
    error=function(e) structure(list(message=conditionMessage(e)), class="captured_error")
  )
  list(value=value, warnings=unique(warns))
}

collapse_warnings <- function(x) if (length(x)) paste(x, collapse=" | ") else ""

message("Processing primary sample-overlap input")
primary.cap <- capture_eval(process.input(info.file, overlap.primary, ref.prefix, c("LDH", "CPI")))
if (inherits(primary.cap$value, "captured_error")) stop(primary.cap$value$message)
input.primary <- primary.cap$value
loci <- read.loci(locus.file)

univ.rows <- list()
bivar.rows <- list()
locus.objects <- list()

for (i in seq_len(nrow(loci))) {
  block <- as.character(loci$LOC[i])
  message("Primary block ", block, " (", i, "/", nrow(loci), ")")
  pc <- capture_eval(process.locus(loci[i, ], input.primary))
  locus <- pc$value
  if (inherits(locus, "captured_error") || is.null(locus)) {
    reason <- if (inherits(locus, "captured_error")) locus$message else "process.locus returned NULL"
    for (phen in c("LDH", "CPI")) {
      univ.rows[[length(univ.rows)+1]] <- data.frame(
        lava_block_id=as.integer(block), phenotype=phen, chr=loci$CHR[i], start=loci$START[i], end=loci$STOP[i],
        n_snps=NA_integer_, n_pcs=NA_integer_, local_variance=NA_real_, local_h2=NA_real_,
        local_h2_latent=NA_real_, SE=NA_real_, P=NA_real_, status="TECHNICAL_NA",
        warning=paste(c(reason, pc$warnings), collapse=" | "), stringsAsFactors=FALSE
      )
    }
    bivar.rows[[length(bivar.rows)+1]] <- data.frame(
      lava_block_id=as.integer(block), chr=loci$CHR[i], start=loci$START[i], end=loci$STOP[i], n_snps=NA_integer_, n_pcs=NA_integer_,
      rg_local=NA_real_, rho_lower=NA_real_, rho_upper=NA_real_, SE_approx_from_95CI=NA_real_, P=NA_real_,
      status="TECHNICAL_NA", warning=paste(c(reason, pc$warnings), collapse=" | "), stringsAsFactors=FALSE
    )
    next
  }
  locus.objects[[block]] <- locus
  available.phenos <- intersect(c("LDH", "CPI"), locus$phenos)
  missing.phenos <- setdiff(c("LDH", "CPI"), available.phenos)
  if (length(available.phenos) == 0) {
    for (phen in c("LDH", "CPI")) {
      univ.rows[[length(univ.rows)+1]] <- data.frame(
        lava_block_id=as.integer(block), phenotype=phen, chr=locus$chr, start=locus$start, end=locus$stop,
        n_snps=locus$n.snps, n_pcs=locus$K, local_variance=NA_real_, local_h2=NA_real_,
        local_h2_latent=NA_real_, SE=NA_real_, P=NA_real_, status="NOT_ESTIMABLE_NEGATIVE_OR_INVALID_VARIANCE",
        warning="process.locus removed both phenotypes before univariate analysis", stringsAsFactors=FALSE
      )
    }
    bivar.rows[[length(bivar.rows)+1]] <- data.frame(
      lava_block_id=as.integer(block), chr=locus$chr, start=locus$start, end=locus$stop, n_snps=locus$n.snps, n_pcs=locus$K,
      rg_local=NA_real_, rho_lower=NA_real_, rho_upper=NA_real_, SE_approx_from_95CI=NA_real_, P=NA_real_,
      status="NOT_TESTED_NOT_BOTH_ELIGIBLE", warning="both phenotypes removed by process.locus", stringsAsFactors=FALSE
    )
    next
  }
  uc <- capture_eval(run.univ(locus, phenos=available.phenos, var=TRUE, cap.estimates=FALSE))
  if (inherits(uc$value, "captured_error") || is.null(uc$value)) {
    reason <- if (inherits(uc$value, "captured_error")) uc$value$message else "run.univ returned NULL"
    for (phen in c("LDH", "CPI")) {
      univ.rows[[length(univ.rows)+1]] <- data.frame(
        lava_block_id=as.integer(block), phenotype=phen, chr=locus$chr, start=locus$start, end=locus$stop,
        n_snps=locus$n.snps, n_pcs=locus$K, local_variance=NA_real_, local_h2=NA_real_,
        local_h2_latent=NA_real_, SE=NA_real_, P=NA_real_, status="TECHNICAL_NA",
        warning=paste(c(reason, pc$warnings, uc$warnings), collapse=" | "), stringsAsFactors=FALSE
      )
    }
    bivar.rows[[length(bivar.rows)+1]] <- data.frame(
      lava_block_id=as.integer(block), chr=locus$chr, start=locus$start, end=locus$stop, n_snps=locus$n.snps, n_pcs=locus$K,
      rg_local=NA_real_, rho_lower=NA_real_, rho_upper=NA_real_, SE_approx_from_95CI=NA_real_, P=NA_real_,
      status="TECHNICAL_NA", warning=paste(c(reason, pc$warnings, uc$warnings), collapse=" | "), stringsAsFactors=FALSE
    )
    next
  }
  u <- uc$value
  for (j in seq_len(nrow(u))) {
    raw.h2 <- as.numeric(u$h2.obs[j])
    latent <- if ("h2.latent" %in% names(u)) as.numeric(u$h2.latent[j]) else NA_real_
    status <- if (is.finite(raw.h2) && is.finite(u$p[j])) "SUCCESS" else "TECHNICAL_NA"
    warning <- paste(c(pc$warnings, uc$warnings, if (is.finite(raw.h2) && raw.h2 < 0) "negative raw local h2; eligibility uses returned P and h2 is retained uncapped" else NULL), collapse=" | ")
    univ.rows[[length(univ.rows)+1]] <- data.frame(
      lava_block_id=as.integer(block), phenotype=as.character(u$phen[j]), chr=locus$chr, start=locus$start, end=locus$stop,
      n_snps=locus$n.snps, n_pcs=locus$K, local_variance=as.numeric(u$var[j]), local_h2=raw.h2,
      local_h2_latent=latent, SE=NA_real_, P=as.numeric(u$p[j]), status=status,
      warning=warning, stringsAsFactors=FALSE
    )
  }
  if (length(missing.phenos)) {
    for (phen in missing.phenos) {
      univ.rows[[length(univ.rows)+1]] <- data.frame(
        lava_block_id=as.integer(block), phenotype=phen, chr=locus$chr, start=locus$start, end=locus$stop,
        n_snps=locus$n.snps, n_pcs=locus$K, local_variance=NA_real_, local_h2=NA_real_,
        local_h2_latent=NA_real_, SE=NA_real_, P=NA_real_, status="NOT_ESTIMABLE_NEGATIVE_OR_INVALID_VARIANCE",
        warning="phenotype removed by process.locus before univariate analysis", stringsAsFactors=FALSE
      )
    }
  }
  uu <- do.call(rbind, univ.rows)
  cur <- subset(uu, lava_block_id == as.integer(block))
  eligible <- nrow(cur) == 2 && all(cur$status == "SUCCESS") && all(cur$P < 0.05)
  if (!eligible) {
    bivar.rows[[length(bivar.rows)+1]] <- data.frame(
      lava_block_id=as.integer(block), chr=locus$chr, start=locus$start, end=locus$stop, n_snps=locus$n.snps, n_pcs=locus$K,
      rg_local=NA_real_, rho_lower=NA_real_, rho_upper=NA_real_, SE_approx_from_95CI=NA_real_, P=NA_real_,
      status="NOT_TESTED_NOT_BOTH_ELIGIBLE", warning=collapse_warnings(c(pc$warnings, uc$warnings)), stringsAsFactors=FALSE
    )
    next
  }
  set.seed(2026091300 + as.integer(block))
  bc <- capture_eval(run.bivar(locus, phenos=c("LDH", "CPI"), adap.thresh=c(1e-4, 1e-6), p.values=TRUE, CIs=TRUE, param.lim=1.25, cap.estimates=TRUE))
  if (inherits(bc$value, "captured_error") || is.null(bc$value)) {
    reason <- if (inherits(bc$value, "captured_error")) bc$value$message else "run.bivar returned NULL"
    bivar.rows[[length(bivar.rows)+1]] <- data.frame(
      lava_block_id=as.integer(block), chr=locus$chr, start=locus$start, end=locus$stop, n_snps=locus$n.snps, n_pcs=locus$K,
      rg_local=NA_real_, rho_lower=NA_real_, rho_upper=NA_real_, SE_approx_from_95CI=NA_real_, P=NA_real_,
      status="TECHNICAL_NA", warning=paste(c(reason, pc$warnings, uc$warnings, bc$warnings), collapse=" | "), stringsAsFactors=FALSE
    )
  } else {
    b <- bc$value[1, ]
    se.approx <- if (all(is.finite(c(b$rho.lower, b$rho.upper)))) (b$rho.upper - b$rho.lower)/(2*1.96) else NA_real_
    bivar.rows[[length(bivar.rows)+1]] <- data.frame(
      lava_block_id=as.integer(block), chr=locus$chr, start=locus$start, end=locus$stop, n_snps=locus$n.snps, n_pcs=locus$K,
      rg_local=as.numeric(b$rho), rho_lower=as.numeric(b$rho.lower), rho_upper=as.numeric(b$rho.upper),
      SE_approx_from_95CI=se.approx, P=as.numeric(b$p), status=ifelse(is.finite(b$rho) && is.finite(b$p), "SUCCESS", "NUMERICAL_BOUNDARY_NA"),
      warning=collapse_warnings(c(pc$warnings, uc$warnings, bc$warnings)), stringsAsFactors=FALSE
    )
  }
}

univ <- do.call(rbind, univ.rows)
bivar <- do.call(rbind, bivar.rows)
bivar$FDR_BH <- NA_real_
tested <- which(bivar$status == "SUCCESS" & is.finite(bivar$P))
if (length(tested)) bivar$FDR_BH[tested] <- p.adjust(bivar$P[tested], method="BH")

write.table(univ, file.path(root, "04_univariate/LAVA_UNIVARIATE_ALL.tsv"), sep="\t", quote=FALSE, row.names=FALSE, na="NA")
write.table(bivar, file.path(root, "05_bivariate/LAVA_BIVARIATE_ALL.tsv"), sep="\t", quote=FALSE, row.names=FALSE, na="NA")

# Zero-overlap sensitivity: rerun only the frozen primary-eligible blocks.
message("Processing zero-overlap sensitivity input")
zero.cap <- capture_eval(process.input(info.file, overlap.zero, ref.prefix, c("LDH", "CPI")))
if (inherits(zero.cap$value, "captured_error")) stop(zero.cap$value$message)
input.zero <- zero.cap$value
sens.rows <- list()
for (block in as.character(bivar$lava_block_id[bivar$status == "SUCCESS"])) {
  i <- which(as.character(loci$LOC) == block)[1]
  pc <- capture_eval(process.locus(loci[i, ], input.zero))
  if (inherits(pc$value, "captured_error") || is.null(pc$value)) {
    sens.rows[[length(sens.rows)+1]] <- data.frame(lava_block_id=as.integer(block), rg_local_zero_overlap=NA_real_, rho_lower_zero_overlap=NA_real_, rho_upper_zero_overlap=NA_real_, P_zero_overlap=NA_real_, status="TECHNICAL_NA", warning=collapse_warnings(pc$warnings))
    next
  }
  set.seed(2026091300 + as.integer(block))
  bc <- capture_eval(run.bivar(pc$value, phenos=c("LDH", "CPI"), adap.thresh=c(1e-4, 1e-6), p.values=TRUE, CIs=TRUE, param.lim=1.25, cap.estimates=TRUE))
  if (inherits(bc$value, "captured_error") || is.null(bc$value)) {
    sens.rows[[length(sens.rows)+1]] <- data.frame(lava_block_id=as.integer(block), rg_local_zero_overlap=NA_real_, rho_lower_zero_overlap=NA_real_, rho_upper_zero_overlap=NA_real_, P_zero_overlap=NA_real_, status="TECHNICAL_NA", warning=paste(c(if(inherits(bc$value,"captured_error")) bc$value$message else "NULL",pc$warnings,bc$warnings),collapse=" | "))
  } else {
    b <- bc$value[1,]
    sens.rows[[length(sens.rows)+1]] <- data.frame(lava_block_id=as.integer(block), rg_local_zero_overlap=as.numeric(b$rho), rho_lower_zero_overlap=as.numeric(b$rho.lower), rho_upper_zero_overlap=as.numeric(b$rho.upper), P_zero_overlap=as.numeric(b$p), status="SUCCESS", warning=collapse_warnings(c(pc$warnings,bc$warnings)))
  }
}
sens <- if(length(sens.rows)) do.call(rbind,sens.rows) else data.frame(lava_block_id=integer(),rg_local_zero_overlap=double(),rho_lower_zero_overlap=double(),rho_upper_zero_overlap=double(),P_zero_overlap=double(),status=character(),warning=character())
write.table(sens, file.path(root, "05_bivariate/LAVA_SAMPLE_OVERLAP_SENSITIVITY.tsv"), sep="\t", quote=FALSE, row.names=FALSE, na="NA")

qc <- data.frame(
  scope=c("GLOBAL_COMMON_ALIGNED","GLOBAL_UNALIGNABLE"),
  phenotype=c("LDH_AND_CPI","LDH_AND_CPI"),
  lava_block_id=NA_integer_,
  n_input=NA_integer_, n_reference_matched=NA_integer_, n_allele_aligned=c(length(input.primary$analysis.snps),NA_integer_),
  n_ambiguous_removed=c(NA_integer_,length(input.primary$unalignable.snps)), n_duplicate_removed=0L,
  n_final=c(length(input.primary$analysis.snps),NA_integer_),
  status=c("SUCCESS","AUDIT_COUNT"), warning=c(collapse_warnings(primary.cap$warnings),"")
)
block.qc <- unique(univ[,c("lava_block_id","n_snps","status","warning")])
block.qc$scope <- "LAVA_BLOCK"
block.qc$phenotype <- "LDH_AND_CPI"
block.qc$n_input <- NA_integer_; block.qc$n_reference_matched <- block.qc$n_snps
block.qc$n_allele_aligned <- block.qc$n_snps; block.qc$n_ambiguous_removed <- NA_integer_
block.qc$n_duplicate_removed <- 0L; block.qc$n_final <- block.qc$n_snps
block.qc$n_snps <- NULL
block.qc <- block.qc[,names(qc)]
qc <- rbind(qc,block.qc)
write.table(qc, file.path(root, "03_harmonized/LAVA_HARMONIZATION_QC.tsv"), sep="\t", quote=FALSE, row.names=FALSE, na="NA")

writeLines("0", file.path(root,"00_admin/lava_analysis.exitcode"))
writeLines(paste0("SUCCESS_",nrow(loci),"_UNIQUE_BLOCKS_",length(tested),"_BIVARIATE_TESTS"), file.path(root,"00_admin/lava_analysis.done"))
cat("LAVA_PRIMARY_COMPLETE blocks=",nrow(loci)," tested=",length(tested),"\n",sep="")
