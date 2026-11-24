#!/usr/bin/env python3
"""Add frozen fixed-locus identifiers to completed block-level LAVA tables.

This is a reporting-only repair. It does not invoke LAVA or alter any statistic.
"""
import csv
import os
import tempfile

ROOT = "/analysis/projects/JTM_recent_highIF_direct_clone/P2_LDH_shared_architecture/11_NatCommun_strict_clone_Pro/18_targeted_LAVA_20260913"
MAP = os.path.join(ROOT, "01_inputs", "FIXED_LOCUS_TO_LAVA_BLOCK_MAP.tsv")
TARGETS = [
    os.path.join(ROOT, "04_univariate", "LAVA_UNIVARIATE_ALL.tsv"),
    os.path.join(ROOT, "04_univariate", "LAVA_BIVARIATE_ELIGIBILITY.tsv"),
    os.path.join(ROOT, "05_bivariate", "LAVA_BIVARIATE_ALL.tsv"),
]

def read_tsv(path):
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return reader.fieldnames, list(reader)

_, mapping = read_tsv(MAP)
by_block = {}
for row in mapping:
    by_block.setdefault(str(row["lava_block_id"]), []).append(row["fixed_locus_id"])

for path in TARGETS:
    fields, data = read_tsv(path)
    report_fields = ["fixed_locus_id", "fixed_locus_count", "shared_lava_block"]
    fields = [f for f in fields if f not in report_fields]
    insert_at = fields.index("lava_block_id") + 1
    fields[insert_at:insert_at] = report_fields
    for row in data:
        ids = sorted(set(by_block[str(row["lava_block_id"])]))
        row["fixed_locus_id"] = ";".join(ids)
        row["fixed_locus_count"] = str(len(ids))
        row["shared_lava_block"] = "TRUE" if len(ids) > 1 else "FALSE"
    fd, tmp = tempfile.mkstemp(prefix=".enrich.", dir=os.path.dirname(path), text=True)
    try:
        with os.fdopen(fd, "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n", extrasaction="ignore")
            writer.writeheader()
            writer.writerows(data)
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)

print("REPORTING_REPAIR_COMPLETE: fixed_locus_id added; no statistics changed")
