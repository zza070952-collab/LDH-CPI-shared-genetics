#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
import re
import subprocess
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np
import pandas as pd

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["font.size"] = 6.5
plt.rcParams["axes.linewidth"] = 0.7
plt.rcParams["axes.spines.top"] = False
plt.rcParams["axes.spines.right"] = False

REPLAN = Path("/analysis/projects/JTM_recent_highIF_direct_clone/P2_LDH_shared_architecture/11_NatCommun_strict_clone_Pro")
REV = REPLAN / "17_evidence_tier_revision_20260912"
ROOT = REPLAN / "18_targeted_LAVA_20260913"


def read(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", low_memory=False)


def finite(x) -> bool:
    try:
        return math.isfinite(float(x))
    except Exception:
        return False


def fmt(x, digits=4) -> str:
    return f"{float(x):.{digits}g}" if finite(x) else "NA"


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


fixed = read(ROOT / "01_inputs/FIXED_18_LOCI_FOR_LAVA.tsv")
mapping = read(ROOT / "01_inputs/FIXED_LOCUS_TO_LAVA_BLOCK_MAP.tsv")
univ = read(ROOT / "04_univariate/LAVA_UNIVARIATE_ALL.tsv")
bivar = read(ROOT / "05_bivariate/LAVA_BIVARIATE_ALL.tsv")
sens = read(ROOT / "05_bivariate/LAVA_SAMPLE_OVERLAP_SENSITIVITY.tsv")
assert len(fixed) == 18 and fixed.locus_id.is_unique
assert len(univ) == 36 and univ.lava_block_id.nunique() == 18
assert len(bivar) == 18 and bivar.lava_block_id.is_unique


def add_fixed_locus_reporting_fields(frame: pd.DataFrame) -> pd.DataFrame:
    """Annotate unique block-level outputs without duplicating statistical tests."""
    report_fields = ["fixed_locus_id", "fixed_locus_count", "shared_lava_block"]
    out = frame.drop(columns=report_fields, errors="ignore").copy()
    block_loci = mapping.groupby("lava_block_id").fixed_locus_id.agg(lambda x: ";".join(sorted(set(x)))).to_dict()
    out.insert(1, "fixed_locus_id", out.lava_block_id.map(block_loci))
    out.insert(2, "fixed_locus_count", out.fixed_locus_id.str.count(";").add(1))
    out.insert(3, "shared_lava_block", out.fixed_locus_count.gt(1))
    return out

up = univ.pivot(index="lava_block_id", columns="phenotype", values=["local_h2", "local_h2_latent", "P", "status", "n_snps", "n_pcs"])
elig = []
for block in sorted(univ.lava_block_id.unique()):
    u = univ[univ.lava_block_id.eq(block)].set_index("phenotype")
    technical = len(u) != 2 or u.status.eq("TECHNICAL_NA").all()
    ldh = "LDH" in u.index and u.loc["LDH", "status"] == "SUCCESS" and finite(u.loc["LDH", "P"]) and float(u.loc["LDH", "P"]) < 0.05
    cpi = "CPI" in u.index and u.loc["CPI", "status"] == "SUCCESS" and finite(u.loc["CPI", "P"]) and float(u.loc["CPI", "P"]) < 0.05
    if technical:
        category = "TECHNICAL_NA"
    elif ldh and cpi:
        category = "BOTH_ELIGIBLE"
    elif ldh:
        category = "LDH_ONLY"
    elif cpi:
        category = "CPI_ONLY"
    else:
        category = "NEITHER"
    elig.append({
        "lava_block_id": int(block), "LDH_local_h2": u.loc["LDH", "local_h2"] if "LDH" in u.index else np.nan,
        "LDH_local_h2_latent": u.loc["LDH", "local_h2_latent"] if "LDH" in u.index else np.nan,
        "LDH_P": u.loc["LDH", "P"] if "LDH" in u.index else np.nan,
        "CPI_local_h2": u.loc["CPI", "local_h2"] if "CPI" in u.index else np.nan,
        "CPI_P": u.loc["CPI", "P"] if "CPI" in u.index else np.nan,
        "eligibility": category,
        "n_snps": u.n_snps.dropna().iloc[0] if u.n_snps.notna().any() else np.nan,
        "n_pcs": u.n_pcs.dropna().iloc[0] if u.n_pcs.notna().any() else np.nan,
    })
elig = pd.DataFrame(elig)
add_fixed_locus_reporting_fields(elig).to_csv(ROOT / "04_univariate/LAVA_BIVARIATE_ELIGIBILITY.tsv", sep="\t", index=False, na_rep="NA")
bivar = bivar.drop(columns=["fixed_locus_id", "fixed_locus_count", "shared_lava_block"], errors="ignore")
add_fixed_locus_reporting_fields(bivar).to_csv(ROOT / "05_bivariate/LAVA_BIVARIATE_ALL.tsv", sep="\t", index=False, na_rep="NA")

block = elig.merge(bivar, on="lava_block_id", how="left", suffixes=("_univ", "_bivar")).merge(sens, on="lava_block_id", how="left")
bm = mapping.merge(block, on="lava_block_id", how="left")

def aggregate_locus(frow):
    m = bm[bm.fixed_locus_id.eq(frow.locus_id)].copy()
    cats = set(m.eligibility.dropna())
    if "BOTH_ELIGIBLE" in cats:
        category = "BOTH_ELIGIBLE"
    elif cats == {"TECHNICAL_NA"} or not cats:
        category = "TECHNICAL_NA"
    elif "LDH_ONLY" in cats and "CPI_ONLY" not in cats:
        category = "LDH_ONLY"
    elif "CPI_ONLY" in cats and "LDH_ONLY" not in cats:
        category = "CPI_ONLY"
    else:
        category = "NEITHER"
    tested = m[(m.status_x == "SUCCESS") & pd.to_numeric(m.FDR_BH, errors="coerce").notna()].copy()
    if len(tested):
        best = tested.sort_values(["FDR_BH", "P", "lava_block_id"]).iloc[0]
        rg, lo, hi, p, q = best.rg_local, best.rho_lower, best.rho_upper, best.P, best.FDR_BH
        selected_block = int(best.lava_block_id)
        lava_status = "SUCCESS"
        p0 = best.P_zero_overlap if "P_zero_overlap" in best else np.nan
        rg0 = best.rg_local_zero_overlap if "rg_local_zero_overlap" in best else np.nan
    else:
        rg = lo = hi = p = q = p0 = rg0 = np.nan
        selected_block = np.nan
        lava_status = "NOT_TESTED_NOT_BOTH_ELIGIBLE" if category != "TECHNICAL_NA" else "TECHNICAL_NA"
    strong_coloc = finite(frow.regional_PP4) and float(frow.regional_PP4) >= 0.80
    if finite(frow.signal_level_PP4) and float(frow.signal_level_PP4) >= 0.80:
        strong_coloc = True
    supportive_coloc = finite(frow.regional_PP4) and 0.50 <= float(frow.regional_PP4) < 0.80
    lava_sig = finite(q) and float(q) < 0.05 and finite(rg)
    lava_concordant = lava_sig and float(rg) > 0
    original = str(frow.existing_coloc_tier)
    if strong_coloc and lava_concordant:
        tier = "TIER_S"
    elif lava_sig and supportive_coloc:
        tier = "TIER_A"
    elif original.startswith("A_"):
        tier = "TIER_A"
    elif original.startswith("B_"):
        tier = "TIER_B"
    elif lava_sig:
        tier = "TIER_B"
    elif original.startswith("E_"):
        tier = "TIER_E"
    else:
        tier = "TIER_C"
    interp = (
        f"Original coloc={original}; LAVA blocks={','.join(m.lava_block_id.astype(int).astype(str))}; "
        f"eligibility={category}; selected local rg={fmt(rg)}, P={fmt(p)}, BH q={fmt(q)}. "
        "LAVA is locus-level evidence and is not gene-specific."
    )
    return {
        "locus_id": frow.locus_id, "lead_snp": frow.lead_snp, "conjFDR": frow.lead_conjFDR,
        "strict_conjFDR_flag": frow.conjFDR_tier == "LT_0.01", "regional_PP4": frow.regional_PP4,
        "signal_level_PP4": frow.signal_level_PP4, "coloc_tier_before_LAVA": original,
        "lava_block_ids": ";".join(m.lava_block_id.astype(int).astype(str)),
        "shared_lava_block": bool(m.shared_lava_block.any()), "fixed_locus_crosses_blocks": bool(m.fixed_locus_crosses_blocks.any()),
        "LDH_local_h2": ";".join(fmt(x,5) for x in m.LDH_local_h2), "LDH_local_h2_P": ";".join(fmt(x,5) for x in m.LDH_P),
        "CPI_local_h2": ";".join(fmt(x,5) for x in m.CPI_local_h2), "CPI_local_h2_P": ";".join(fmt(x,5) for x in m.CPI_P),
        "LAVA_eligible": category, "selected_lava_block": selected_block,
        "local_rg": rg, "local_rg_lower": lo, "local_rg_upper": hi, "local_rg_P": p, "local_rg_FDR": q,
        "local_rg_zero_overlap": rg0, "local_rg_P_zero_overlap": p0,
        "local_rg_direction_vs_positive_global": ("CONCORDANT_POSITIVE" if finite(rg) and float(rg) > 0 else "OPPOSING_NEGATIVE" if finite(rg) and float(rg) < 0 else "NOT_ESTIMATED"),
        "LAVA_status": lava_status, "evidence_tier_original": original,
        "LAVA_ENHANCED_TIER": tier, "interpretation": interp,
    }

integrated = pd.DataFrame([aggregate_locus(r) for r in fixed.itertuples(index=False)])
assert len(integrated) == 18 and integrated.locus_id.is_unique
integrated.to_csv(ROOT / "06_integration/LAVA_COLOC_INTEGRATED_18.tsv", sep="\t", index=False, na_rep="NA")
integrated.to_csv(ROOT / "08_tables/Supplementary_LAVA_18_loci.tsv", sep="\t", index=False, na_rep="NA")

# Candidate ranking: LAVA support is inherited from the locus and never treated as gene-specific.
cmat = read(REV / "04_gene_evidence/CANDIDATE_CLAIM_EVIDENCE_MATRIX.tsv")
focus = cmat[cmat.gene.isin(["DCC", "COL7A1", "APEH"])].copy()
assert len(focus) == 3
focus = focus.merge(integrated[["locus_id", "LAVA_ENHANCED_TIER", "LAVA_eligible", "local_rg", "local_rg_P", "local_rg_FDR"]], on="locus_id", how="left")
tier_order = {"TIER_S": 5, "TIER_A": 4, "TIER_B": 3, "TIER_C": 2, "TIER_E": 1}
candidate_tie = {"DCC": 3, "COL7A1": 2, "APEH": 1}
focus["tier_order_internal"] = focus.LAVA_ENHANCED_TIER.map(tier_order)
focus["candidate_tiebreak_internal"] = focus.gene.map(candidate_tie)
focus = focus.sort_values(["tier_order_internal", "candidate_tiebreak_internal"], ascending=False).reset_index(drop=True)
focus.insert(0, "rank_after_LAVA", np.arange(1, len(focus)+1))
focus["rank_before_LAVA"] = focus.gene.map({"DCC": 1, "COL7A1": 2, "APEH": 3})
focus["rank_changed"] = focus.rank_after_LAVA.ne(focus.rank_before_LAVA)
focus["LAVA_LOCUS_SUPPORT"] = focus.apply(lambda r: f"{r.LAVA_ENHANCED_TIER}; eligibility={r.LAVA_eligible}; local rg={fmt(r.local_rg)}, BH q={fmt(r.local_rg_FDR)}; locus-level only", axis=1)
focus["ranking_boundary"] = "LAVA changes locus confidence only; ranking also preserves L2G/QTL/candidate-coloc/human-reference/testability evidence and does not infer gene causality."
focus = focus.drop(columns=["tier_order_internal", "candidate_tiebreak_internal"])
focus.to_csv(ROOT / "06_integration/FINAL_CANDIDATE_RANKING_AFTER_LAVA.tsv", sep="\t", index=False, na_rep="NA")

# Manifests.
git_sha = subprocess.check_output(["/analysis/miniconda3/bin/git", "-C", str(ROOT / "02_reference/LAVA_v0.1.5"), "rev-parse", "HEAD"], text=True).strip()
deps = read(ROOT / "00_admin/R_DEPENDENCY_VERSIONS.tsv")
version_rows = [
    {"component": "LAVA", "version": "0.1.5", "tag_or_commit": git_sha, "source": "https://github.com/josefin-werme/LAVA/releases/tag/v0.1.5", "sha256": "GIT_COMMIT_SHA_IS_CANONICAL"},
]
for r in deps.itertuples(index=False):
    version_rows.append({"component": r.component, "version": r.version, "tag_or_commit": "NA", "source": "R runtime/package library", "sha256": "NA"})
pd.DataFrame(version_rows).to_csv(ROOT / "00_admin/LAVA_VERSION_MANIFEST.tsv", sep="\t", index=False)

url_map = {
    "lava-ukb-v1.1_chr1-2.zip": "https://vu.data.surfsara.nl/index.php/s/7NBVIvtPRdu7Qhz/download",
    "lava-ukb-v1.1_chr3-4.zip": "https://vu.data.surfsara.nl/index.php/s/fy6ITboMojHrQXr/download",
    "lava-ukb-v1.1_chr5-6.zip": "https://vu.data.surfsara.nl/index.php/s/mRz31q0lq7KMcuI/download",
    "lava-ukb-v1.1_chr7-9.zip": "https://vu.data.surfsara.nl/index.php/s/89RXxPN2BlOqxLs/download",
    "lava-ukb-v1.1_chr10-12.zip": "https://vu.data.surfsara.nl/index.php/s/3dU1L2Hap43xuCs/download",
    "lava-ukb-v1.1_chr17-23.zip": "https://vu.data.surfsara.nl/index.php/s/U8eg5XfTr8qPzPp/download",
}
ref_rows = []
for p in sorted((ROOT / "02_reference/archives").glob("*.zip")):
    ref_rows.append({"record_type": "OFFICIAL_ARCHIVE", "file": str(p.relative_to(ROOT)), "chromosome": "BUNDLE", "bytes": p.stat().st_size, "sha256": sha(p), "version": "UKB_EUR_v1.1", "N_reference": 100000, "source": url_map[p.name], "used": True})
for chrom in sorted(mapping.chr.unique()):
    for ext in ("info", "bcor"):
        p = ROOT / f"02_reference/analysis_reference/lava-ukb-v1.1_chr{int(chrom)}.{ext}"
        assert p.exists() and p.stat().st_size > 0
        ref_rows.append({"record_type": "ANALYSIS_FILE", "file": str(p.relative_to(ROOT)), "chromosome": int(chrom), "bytes": p.stat().st_size, "sha256": sha(p), "version": "UKB_EUR_v1.1", "N_reference": 100000, "source": "extracted from official archive", "used": True})
pd.DataFrame(ref_rows).to_csv(ROOT / "00_admin/LAVA_REFERENCE_MANIFEST.tsv", sep="\t", index=False)

# Figure source data.
figdata = integrated.copy()
figdata.to_csv(ROOT / "07_figures/Figure3_LAVA_updated_source_data.tsv", sep="\t", index=False, na_rep="NA")

labels = [x.replace("LDH_CPI_p05_", "") for x in figdata.locus_id]
n = len(figdata)
fig = plt.figure(figsize=(7.2, 7.7), constrained_layout=False)
gs = fig.add_gridspec(1, 2, width_ratios=[3.7, 1.55], wspace=0.42, left=0.18, right=0.98, top=0.94, bottom=0.08)
ax = fig.add_subplot(gs[0, 0])

# Panel F: a normalized, fully annotated evidence ledger.
columns = ["conjFDR", "strict", "regional\nPP4", "signal\nPP4", "h2\neligible", "local rg", "BH q", "tier"]
vals = np.zeros((n, len(columns)))
ann = np.empty((n, len(columns)), dtype=object)
tier_num = {"TIER_S": 1.0, "TIER_A": 0.8, "TIER_B": 0.6, "TIER_C": 0.35, "TIER_E": 0.1}
elig_num = {"BOTH_ELIGIBLE": 1.0, "LDH_ONLY": 0.7, "CPI_ONLY": 0.55, "NEITHER": 0.25, "TECHNICAL_NA": 0.0}
for i, r in figdata.iterrows():
    vals[i, 0] = min(-math.log10(float(r.conjFDR))/6, 1)
    ann[i, 0] = f"{float(r.conjFDR):.1e}"
    vals[i, 1] = 1.0 if bool(r.strict_conjFDR_flag) else 0.15; ann[i, 1] = "yes" if bool(r.strict_conjFDR_flag) else "no"
    vals[i, 2] = float(r.regional_PP4) if finite(r.regional_PP4) else 0; ann[i, 2] = fmt(r.regional_PP4, 3)
    vals[i, 3] = float(r.signal_level_PP4) if finite(r.signal_level_PP4) else 0; ann[i, 3] = fmt(r.signal_level_PP4, 3)
    vals[i, 4] = elig_num.get(r.LAVA_eligible, 0); ann[i, 4] = {"BOTH_ELIGIBLE":"both","LDH_ONLY":"LDH","CPI_ONLY":"CPI","NEITHER":"none","TECHNICAL_NA":"NA"}.get(r.LAVA_eligible,"NA")
    vals[i, 5] = (float(r.local_rg)+1)/2 if finite(r.local_rg) else 0; ann[i, 5] = fmt(r.local_rg, 3)
    vals[i, 6] = min(-math.log10(float(r.local_rg_FDR))/3, 1) if finite(r.local_rg_FDR) and float(r.local_rg_FDR)>0 else 0; ann[i, 6] = fmt(r.local_rg_FDR, 3)
    vals[i, 7] = tier_num[r.LAVA_ENHANCED_TIER]; ann[i, 7] = r.LAVA_ENHANCED_TIER.replace("TIER_", "")
im = ax.imshow(vals, cmap="Blues", vmin=0, vmax=1, aspect="auto")
ax.set_xticks(range(len(columns))); ax.set_xticklabels(columns, rotation=35, ha="right")
ax.set_yticks(range(n)); ax.set_yticklabels(labels)
ax.tick_params(length=0)
for i in range(n):
    for j in range(len(columns)):
        color = "white" if vals[i, j] > 0.62 else "#252525"
        ax.text(j, i, ann[i, j], ha="center", va="center", fontsize=5.1, color=color)
for s in ax.spines.values(): s.set_visible(False)
ax.set_title("f  Fixed-locus evidence ledger", loc="left", fontweight="bold", fontsize=8)

# Panel G: all official blocks intersecting chr3/chr18, preserving NA.
axg = fig.add_subplot(gs[0, 1])
focus_ids = ["LDH_CPI_p05_chr3_001", "LDH_CPI_p05_chr18_001"]
submap = bm[bm.fixed_locus_id.isin(focus_ids)].copy().sort_values(["fixed_locus_id", "lava_block_id"])
submap["plot_label"] = submap.apply(lambda r: f"{r.fixed_locus_id.replace('LDH_CPI_p05_','')}\nblock {int(r.lava_block_id)}", axis=1)
y = np.arange(len(submap))[::-1]
for yi, (_, r) in zip(y, submap.iterrows()):
    if r.status_x == "SUCCESS" and finite(r.rg_local):
        axg.plot([r.rho_lower, r.rho_upper], [yi, yi], color="#0F4D92", lw=1.5)
        axg.plot(r.rg_local, yi, "o", color="#B64342" if finite(r.FDR_BH) and r.FDR_BH < 0.05 else "#0F4D92", ms=4)
        axg.text(1.03, yi, f"q={fmt(r.FDR_BH,2)}", va="center", fontsize=5)
    else:
        axg.plot(0, yi, marker="x", color="#8F8F8F", ms=4)
        axg.text(0.05, yi, r.eligibility.replace("_", " ").lower(), va="center", fontsize=5, color="#606060")
axg.axvline(0, color="#8F8F8F", lw=0.8, ls="--")
axg.set_xlim(-1.05, 1.30)
axg.set_yticks(y); axg.set_yticklabels(submap.plot_label)
axg.set_xlabel("Local genetic correlation (95% CI)")
axg.set_title("g  chr3 and chr18", loc="left", fontweight="bold", fontsize=8)
fig.suptitle("Local genetic-correlation evidence across 18 frozen LDH–CPI regions", fontsize=9, fontweight="bold", y=0.985)
fig.text(0.18, 0.015, "Eligibility: LDH and CPI local-h² P<0.05 in the same official LAVA block. BH q covers all tested blocks; NA is not zero.", fontsize=5.5)
for ext, kwargs in [("pdf", {}), ("png", {"dpi": 600}), ("svg", {})]:
    fig.savefig(ROOT / f"07_figures/Figure3_LAVA_updated.{ext}", bbox_inches="tight", **kwargs)
plt.close(fig)

# Manuscript documents.
counts = integrated.LAVA_eligible.value_counts().to_dict()
tier_counts = integrated.LAVA_ENHANCED_TIER.value_counts().to_dict()
nominal = int(((bivar.status == "SUCCESS") & (pd.to_numeric(bivar.P, errors="coerce") < 0.05)).sum())
fdr = int(((bivar.status == "SUCCESS") & (pd.to_numeric(bivar.FDR_BH, errors="coerce") < 0.05)).sum())

def locus_text(lid):
    r = integrated[integrated.locus_id.eq(lid)].iloc[0]
    return f"blocks {r.lava_block_ids}; LDH local h2={r.LDH_local_h2} (P={r.LDH_local_h2_P}); CPI local h2={r.CPI_local_h2} (P={r.CPI_local_h2_P}); eligibility={r.LAVA_eligible}; selected local rg={fmt(r.local_rg)} (95% CI {fmt(r.local_rg_lower)} to {fmt(r.local_rg_upper)}), P={fmt(r.local_rg_P)}, BH q={fmt(r.local_rg_FDR)}; {r.LAVA_ENHANCED_TIER}."

methods = f"""# LAVA methods draft\n\nWe conducted a bounded local genetic-correlation analysis of the 18 previously fixed non-MHC LDH-CPI conjunctional-FDR regions. Region definitions, lead variants and conjunctional-FDR thresholds were not changed. Fixed regions were mapped to every intersecting block in the official LAVA GRCh37 partition; unique official blocks were analysed once and then mapped back to fixed locus IDs. Two chr2 loci shared one LAVA block, while the chr18 locus crossed two blocks.\n\nLAVA v0.1.5 (commit {git_sha}) was used with the official UK Biobank European binary LD reference v1.1 based on 100,000 unrelated participants. Summary statistics were the checksum-preserved FinnGen R12 LDH GWAS (35,290 cases and 353,224 controls) and Wang 2025 UK Biobank CPI GWAS. LDH was treated as binary and CPI as continuous. No population prevalence was supplied; the v0.1.5 no-ascertainment transformation was retained and latent-scale values are not interpreted as population-liability estimates.\n\nThe LDSC intercept covariance matrix yielded a primary sampling correlation of 0.004416; zero correlation was evaluated only as sensitivity. For every selected block, run.univ() was applied to both traits. Bivariate run.bivar() analysis was restricted to blocks in which both univariate P values were below 0.05. All actually tested block-level P values formed one Benjamini-Hochberg family. LAVA v0.1.5 does not return a native standard error for univariate h2; that field is therefore NA. Bivariate uncertainty is reported using native 95% confidence limits, with an explicitly labelled approximate SE derived from CI width only for tabular compatibility.\n"""
(ROOT / "09_manuscript/LAVA_METHODS_DRAFT.md").write_text(methods, encoding="utf-8")

results = f"""# LAVA results draft\n\nThe 18 frozen regions mapped to 18 unique official LAVA blocks through 19 fixed-locus-to-block links. At fixed-locus level, {counts.get('BOTH_ELIGIBLE',0)} were BOTH_ELIGIBLE, {counts.get('LDH_ONLY',0)} were LDH_ONLY, {counts.get('CPI_ONLY',0)} were CPI_ONLY, {counts.get('NEITHER',0)} were NEITHER and {counts.get('TECHNICAL_NA',0)} were TECHNICAL_NA. Among the {int((bivar.status == 'SUCCESS').sum())} tested unique blocks, {nominal} had nominal local-rg P<0.05 and {fdr} passed BH q<0.05.\n\nCHR3: {locus_text('LDH_CPI_p05_chr3_001')} This is interpreted alongside regional PP4=0.8409 and the existing LDH-CPI SuSiE signal-pair PP4=0.939674; LAVA does not identify COL7A1, APEH or another gene as causal.\n\nCHR18: {locus_text('LDH_CPI_p05_chr18_001')} This is interpreted alongside regional PP4=0.9835. The previous multi-signal uncertainty is retained, and any LAVA support applies to the locus rather than DCC.\n\nFinal LAVA-enhanced tiers were S={tier_counts.get('TIER_S',0)}, A={tier_counts.get('TIER_A',0)}, B={tier_counts.get('TIER_B',0)}, C={tier_counts.get('TIER_C',0)} and E={tier_counts.get('TIER_E',0)}. Shared loci therefore showed heterogeneous support across complementary local genetic methods; gene prioritization remained separately dependent on mapping, QTL/coloc scope, human reference tissue evidence and experimental testability.\n"""
(ROOT / "09_manuscript/LAVA_RESULTS_DRAFT.md").write_text(results, encoding="utf-8")

limitations = """# LAVA limitations\n\n- The primary CPI GWAS is a population chronic-pain burden score, not pain severity conditional on MRI-confirmed LDH.\n- Fixed conjFDR regions were mapped to official LD blocks. A block can extend beyond a fixed region; chr18 crosses two blocks and two chr2 loci share one block. These are follow-up tests, not independent replication.\n- Bivariate testing is conditional on both local-h2 P values being below 0.05, so non-eligible regions are not local-rg=0.\n- The LDSC-derived sampling correlation is near zero but uncertain; the pairing is low-overlap rather than proven zero-overlap.\n- LAVA v0.1.5 local-h2 output has no native SE. Bivariate confidence intervals arise from the package procedure, and CI-derived approximate SE is descriptive only.\n- A significant local rg is locus-level covariance evidence. It does not establish a causal variant, effector gene, shared expression mediator or disc-to-DRG mechanism.\n- Non-significant or technically unresolved LAVA results do not negate the existing conjunctional-FDR or colocalization evidence.\n- MVP replication is NOT_AVAILABLE / NOT_PURSUED and is not a blocker for computational freezing.\n"""
(ROOT / "09_manuscript/LAVA_LIMITATIONS.md").write_text(limitations, encoding="utf-8")

rank_order = " > ".join(focus.sort_values("rank_after_LAVA").gene)
preferred = "chr3" if tier_order[integrated.loc[integrated.locus_id.eq('LDH_CPI_p05_chr3_001'),'LAVA_ENHANCED_TIER'].iloc[0]] >= tier_order[integrated.loc[integrated.locus_id.eq('LDH_CPI_p05_chr18_001'),'LAVA_ENHANCED_TIER'].iloc[0]] else "chr18"
decision = f"""# Final scope and decision after LAVA\n\nCOMPUTATION_FROZEN = YES\n\nLAVA_STATUS = SUCCESS\n\nLAVA_VERSION = 0.1.5\n\nLD_REFERENCE = LAVA UK Biobank European binary LD v1.1, N=100,000\n\nMVP_REPLICATION = NOT_AVAILABLE / NOT_PURSUED\n\nMVP_BLOCKS_FREEZE = NO\n\nHISTORICAL_P2_HARD_GO = FAIL\n\nREADY_FOR_WETLAB = NO\n\nSTRICT_SHARED_SMR_HEIDI = 0\n\nESTABLISHED_REGULATORY_TRIPLES = 0\n\nLAVA_ENHANCED_TIERS = S:{tier_counts.get('TIER_S',0)}; A:{tier_counts.get('TIER_A',0)}; B:{tier_counts.get('TIER_B',0)}; C:{tier_counts.get('TIER_C',0)}; E:{tier_counts.get('TIER_E',0)}\n\nCANDIDATE_ORDER_AFTER_LAVA = {rank_order}\n\nPREFERRED_LOCUS_FOR_EXPERIMENTAL_DISCUSSION = {preferred}\n\nThe final computation is frozen. No new GWAS, QTL panel, single-cell atlas, machine-learning model, network analysis, MR, colocalization or conjunctional-FDR run is authorized. LAVA changes locus confidence only and does not convert DCC, COL7A1 or APEH into a causal gene. Historical decisions and all negative or unresolved results remain preserved.\n"""
(ROOT / "09_manuscript/FINAL_SCOPE_AND_DECISION_AFTER_LAVA.md").write_text(decision, encoding="utf-8")

outputs = [
    ROOT / "00_admin/LAVA_VERSION_MANIFEST.tsv", ROOT / "00_admin/LAVA_REFERENCE_MANIFEST.tsv",
    ROOT / "01_inputs/FIXED_18_LOCI_FOR_LAVA.tsv", ROOT / "01_inputs/FIXED_LOCUS_TO_LAVA_BLOCK_MAP.tsv",
    ROOT / "03_harmonized/LAVA_HARMONIZATION_QC.tsv", ROOT / "01_inputs/LAVA_SAMPLE_OVERLAP_INPUT.tsv",
    ROOT / "04_univariate/LAVA_UNIVARIATE_ALL.tsv", ROOT / "04_univariate/LAVA_BIVARIATE_ELIGIBILITY.tsv",
    ROOT / "05_bivariate/LAVA_BIVARIATE_ALL.tsv", ROOT / "06_integration/LAVA_COLOC_INTEGRATED_18.tsv",
    ROOT / "06_integration/FINAL_CANDIDATE_RANKING_AFTER_LAVA.tsv", ROOT / "09_manuscript/LAVA_METHODS_DRAFT.md",
    ROOT / "09_manuscript/LAVA_RESULTS_DRAFT.md", ROOT / "09_manuscript/LAVA_LIMITATIONS.md",
    ROOT / "09_manuscript/FINAL_SCOPE_AND_DECISION_AFTER_LAVA.md", ROOT / "07_figures/Figure3_LAVA_updated.pdf",
    ROOT / "07_figures/Figure3_LAVA_updated.png", ROOT / "07_figures/Figure3_LAVA_updated.svg",
    ROOT / "07_figures/Figure3_LAVA_updated_source_data.tsv", ROOT / "08_tables/Supplementary_LAVA_18_loci.tsv",
]
manifest = pd.DataFrame([{"file": str(p.relative_to(REPLAN)), "bytes": p.stat().st_size, "sha256": sha(p)} for p in outputs])
manifest.to_csv(ROOT / "00_admin/FINAL_OUTPUT_MANIFEST.tsv", sep="\t", index=False)
summary = {
    "computation_frozen": True, "lava_version": "0.1.5", "reference": "UKB_EUR_v1.1_N100000",
    "fixed_loci": 18, "unique_blocks": 18, "mapping_rows": 19,
    "eligibility_counts": counts, "bivariate_tested_blocks": int((bivar.status == "SUCCESS").sum()),
    "nominal_p_lt_0.05": nominal, "BH_q_lt_0.05": fdr, "tier_counts": tier_counts,
    "candidate_order": list(focus.sort_values("rank_after_LAVA").gene), "preferred_locus": preferred,
    "MVP_replication": "NOT_AVAILABLE_NOT_PURSUED", "ready_for_wetlab": False,
}
(ROOT / "00_admin/COMPUTATION_FROZEN.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
(ROOT / "00_admin/lava_finalize.done").write_text("COMPUTATION_FROZEN=YES\n", encoding="utf-8")
print(json.dumps(summary, indent=2))
