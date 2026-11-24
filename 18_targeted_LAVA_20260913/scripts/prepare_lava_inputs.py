#!/usr/bin/env python3
from __future__ import annotations

import gzip
import hashlib
import json
import math
from pathlib import Path

import pandas as pd

REPLAN = Path("/analysis/projects/JTM_recent_highIF_direct_clone/P2_LDH_shared_architecture/11_NatCommun_strict_clone_Pro")
REVISION = REPLAN / "17_evidence_tier_revision_20260912"
ROOT = REPLAN / "18_targeted_LAVA_20260913"
LAVA_SRC = ROOT / "02_reference/LAVA_v0.1.5"
LDH = REPLAN / "02_qc/gwas/FINNGEN_R12_LDH.GRCh37.QC.tsv.gz"
CPI = REPLAN / "02_qc/gwas/UKB_CPI_WANG2025.GRCh37.QC.tsv.gz"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


fixed = pd.read_csv(REVISION / "01_locus_evidence/LOCUS_EVIDENCE_18.tsv", sep="\t")
assert len(fixed) == 18 and fixed["locus_id"].is_unique
frozen = pd.DataFrame(
    {
        "locus_id": fixed["locus_id"],
        "chr": fixed["chr"].astype(int),
        "start": fixed["start"].astype(int),
        "end": fixed["end"].astype(int),
        "lead_snp": fixed["lead_SNP"],
        "lead_conjFDR": fixed["conjFDR"],
        "conjFDR_tier": fixed["conjFDR"].apply(lambda x: "LT_0.01" if x < 0.01 else "LT_0.05"),
        "existing_coloc_tier": fixed["evidence_class"],
        "regional_PP4": fixed["PP4"],
        "signal_level_PP4": fixed["locus_id"].map({"LDH_CPI_p05_chr3_001": 0.93967409258028}),
        "technical_status": fixed["trait_trait_coloc_available"],
    }
)
frozen.to_csv(ROOT / "01_inputs/FIXED_18_LOCI_FOR_LAVA.tsv", sep="\t", index=False, na_rep="NA")

blocks = pd.read_csv(
    LAVA_SRC / "support_data/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile",
    sep=r"\s+",
)
blocks.columns = ["lava_block_id", "chr", "block_start", "block_end"]
maps = []
for row in frozen.itertuples(index=False):
    hit = blocks[
        (blocks["chr"] == row.chr)
        & (blocks["block_end"] >= row.start)
        & (blocks["block_start"] <= row.end)
    ].copy()
    assert len(hit) >= 1, row.locus_id
    for b in hit.itertuples(index=False):
        maps.append(
            {
                "fixed_locus_id": row.locus_id,
                "chr": row.chr,
                "fixed_start": row.start,
                "fixed_end": row.end,
                "lead_snp": row.lead_snp,
                "lava_block_id": int(b.lava_block_id),
                "block_start": int(b.block_start),
                "block_end": int(b.block_end),
                "overlap_bp": max(0, min(row.end, b.block_end) - max(row.start, b.block_start) + 1),
            }
        )
mapping = pd.DataFrame(maps)
mapping["n_blocks_for_fixed_locus"] = mapping.groupby("fixed_locus_id")["lava_block_id"].transform("nunique")
mapping["n_fixed_loci_in_lava_block"] = mapping.groupby("lava_block_id")["fixed_locus_id"].transform("nunique")
mapping["fixed_locus_crosses_blocks"] = mapping["n_blocks_for_fixed_locus"] > 1
mapping["shared_lava_block"] = mapping["n_fixed_loci_in_lava_block"] > 1
mapping["mapping_rule"] = "All official LAVA GRCh37 blocks intersecting the unchanged fixed conjFDR interval"
mapping.to_csv(ROOT / "01_inputs/FIXED_LOCUS_TO_LAVA_BLOCK_MAP.tsv", sep="\t", index=False)

selected = blocks[blocks["lava_block_id"].isin(mapping["lava_block_id"].unique())].copy()
selected = selected.sort_values(["chr", "block_start"])
selected[["lava_block_id", "chr", "block_start", "block_end"]].to_csv(
    ROOT / "01_inputs/SELECTED_UNIQUE_LAVA_BLOCKS.tsv", sep="\t", index=False
)
selected.rename(
    columns={"lava_block_id": "LOC", "chr": "CHR", "block_start": "START", "block_end": "STOP"}
)[["LOC", "CHR", "START", "STOP"]].to_csv(
    ROOT / "01_inputs/selected_blocks.loci", sep=" ", index=False
)

intervals: dict[int, list[tuple[int, int]]] = {}
for b in selected.itertuples(index=False):
    intervals.setdefault(int(b.chr), []).append((int(b.block_start), int(b.block_end)))


def extract_sumstats(source: Path, destination: Path) -> dict:
    selected_rows = []
    input_rows = 0
    for chunk in pd.read_csv(source, sep="\t", compression="gzip", chunksize=400_000, low_memory=False):
        input_rows += len(chunk)
        mask = pd.Series(False, index=chunk.index)
        for chrom, ranges in intervals.items():
            cmask = chunk["CHR"].eq(chrom)
            if not cmask.any():
                continue
            rmask = pd.Series(False, index=chunk.index)
            for start, stop in ranges:
                rmask |= cmask & chunk["BP"].between(start, stop)
            mask |= rmask
        if mask.any():
            out = chunk.loc[mask, ["SNP", "A1", "A2", "N", "BETA", "P", "CHR", "BP"]].copy()
            selected_rows.append(out)
    df = pd.concat(selected_rows, ignore_index=True)
    duplicate_rows = int(df.duplicated("SNP", keep=False).sum())
    df = df.sort_values(["CHR", "BP", "SNP"]).drop_duplicates("SNP", keep="first")
    nonfinite = (~pd.to_numeric(df["N"], errors="coerce").map(math.isfinite)) | (~pd.to_numeric(df["BETA"], errors="coerce").map(math.isfinite)) | (~pd.to_numeric(df["P"], errors="coerce").map(math.isfinite))
    invalid = int(nonfinite.sum())
    df = df.loc[~nonfinite]
    df.to_csv(destination, sep="\t", index=False, compression={"method": "gzip", "compresslevel": 6})
    return {
        "source": str(source),
        "source_sha256": sha256(source),
        "n_input_qc_table": input_rows,
        "n_selected_block_rows_before_dedup": int(sum(len(x) for x in selected_rows)),
        "n_duplicate_rows_flagged": duplicate_rows,
        "n_nonfinite_removed": invalid,
        "n_written": len(df),
        "output": str(destination),
        "output_sha256": sha256(destination),
    }


ldh_out = ROOT / "03_harmonized/LDH.selected_blocks.lava.tsv.gz"
cpi_out = ROOT / "03_harmonized/CPI.selected_blocks.lava.tsv.gz"
ldh_qc = extract_sumstats(LDH, ldh_out)
cpi_qc = extract_sumstats(CPI, cpi_out)

pd.DataFrame(
    [
        {"phenotype": "LDH", "cases": 35290, "controls": 353224, "filename": str(ldh_out)},
        {"phenotype": "CPI", "cases": 1, "controls": 0, "filename": str(cpi_out)},
    ]
).to_csv(ROOT / "01_inputs/input.info.tsv", sep="\t", index=False)

# LDSC covariance intercept matrix -> sampling correlation, as documented by LAVA.
diag_ldh = 1.1547
diag_cpi = 1.0223
cross = 0.0048
sampling_r = cross / math.sqrt(diag_ldh * diag_cpi)
overlap_matrix = pd.DataFrame([[1.0, sampling_r], [sampling_r, 1.0]], index=["LDH", "CPI"], columns=["LDH", "CPI"])
overlap_matrix.to_csv(ROOT / "01_inputs/sample.overlap.primary.txt", sep="\t")
pd.DataFrame([[1.0, 0.0], [0.0, 1.0]], index=["LDH", "CPI"], columns=["LDH", "CPI"]).to_csv(
    ROOT / "01_inputs/sample.overlap.zero_sensitivity.txt", sep="\t"
)
pd.DataFrame(
    [
        {
            "trait1": "LDH",
            "trait2": "CPI",
            "LDH_LDSC_intercept": diag_ldh,
            "CPI_LDSC_intercept": diag_cpi,
            "cross_trait_gcov_intercept": cross,
            "cross_trait_gcov_intercept_SE": 0.0067,
            "sampling_correlation_primary": sampling_r,
            "zero_overlap_sensitivity": 0.0,
            "overlap_risk": "LOW_NOT_PROVEN_ZERO",
        }
    ]
).to_csv(ROOT / "01_inputs/LAVA_SAMPLE_OVERLAP_INPUT.tsv", sep="\t", index=False)

(ROOT / "01_inputs/LAVA_SAMPLE_OVERLAP_METHOD.md").write_text(
    f"""# LAVA sample-overlap method\n\nThe primary residual sampling correlation is {sampling_r:.8f}. It was obtained using the official LAVA tutorial transformation of the frozen LDSC intercept covariance matrix: 0.0048 / sqrt(1.1547 x 1.0223). The cross-trait intercept SE is 0.0067, so this near-zero estimate is uncertain; the pair is described as low-overlap, not proven zero-overlap. A zero residual-correlation matrix is run only as a sensitivity analysis.\n""",
    encoding="utf-8",
)

pre_qc = {"LDH": ldh_qc, "CPI": cpi_qc, "selected_unique_blocks": len(selected), "fixed_to_block_rows": len(mapping)}
(ROOT / "03_harmonized/PREPROCESS_QC.json").write_text(json.dumps(pre_qc, indent=2), encoding="utf-8")
print(json.dumps({"fixed_loci": len(frozen), "mapping_rows": len(mapping), "unique_blocks": len(selected), "sampling_r": sampling_r, "LDH_rows": ldh_qc["n_written"], "CPI_rows": cpi_qc["n_written"]}, indent=2))

