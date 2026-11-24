# Lumbar disc herniation and chronic pain intensity

Analysis code for the study of shared genetic architecture between lumbar disc herniation (LDH) and chronic pain intensity (CPI). The primary pair uses FinnGen R12 M13_LUMBAR_PROLAPSE and the UK Biobank full-sample chronic pain intensity GWAS (Zenodo 17219195). Long-term analgesic use and sciatica are distinct phenotypes.

## Analysis modules

| Directory | Analysis |
| --- | --- |
| 16_workflow | Global LDSC, MiXeR, conditional and conjunctional FDR, LD-based locus definition, MAGMA, Open Targets annotations, bidirectional MR, shared SMR/HEIDI, colocalization, cortical S-LDSC and SCAVENGE |
| 18_targeted_LAVA_20260913/scripts | Local heritability and eligible local genetic correlations in the 18 fixed LDH–CPI regions using LAVA v0.1.5 and the UK Biobank European LD reference v1.1 |
| signal_resolution/scripts | Chromosome 18 variant alignment, signed LD checks, SuSiE fits and sensitivity summaries |
| drg/scripts | Donor-level human DRG reference expression, identifier harmonization and prior genomic evidence |
| 21_two_finite_enhancements_20261122/07_scripts | PredictDB expression models, S-PrediXcan, external pain/LDH comparisons and multiplicity checks |
| figures/scripts | Rendering and assembly of computational Figures 1–5 |

## Using the scripts

These are module drivers from the study. They expect the original relative directory layout, input tables and module-specific software. The source snapshot alone does not execute the entire workflow on a clean machine. Machine-specific home prefixes have been replaced with /analysis. Adapt input, interpreter, reference and output paths in a disposable analysis tree before running a driver. Some drivers write or replace intermediates. Read the module input checks before execution.

Use the source accessions and allele definitions in manuscript Additional file 2. Obtain any restricted resources directly from their providers. Place European LD matrices, identifiers and variant ordering in the locations expected by each module. External binaries and vendor libraries are installed separately. Important dependencies include LDSC, MiXeR protocol 1.3, pleioFDR with Octave, MAGMA 1.10, SMR 1.03, R with LAVA, coloc and susieR, Python with NumPy, pandas and SciPy, and MetaXcan at fc735bf01a1cf43f25cef466439dcfa04712dd50. Figure drivers consume source tables under figures/source_data. Run them in a working copy with those tables present.

The code archive contains study-authored analysis scripts and this file. SHA256SUMS.tsv identifies each script. The manuscript reports inferential limits for boundary mixture estimates and signal comparisons lacking a credible set.
