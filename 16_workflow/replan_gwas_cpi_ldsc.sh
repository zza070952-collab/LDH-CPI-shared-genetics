#!/usr/bin/env bash
set -Eeuo pipefail
R=/analysis/projects/JTM_recent_highIF_direct_clone/P2_LDH_shared_architecture/11_NatCommun_strict_clone_Pro
cd "$R"
exec > "$R/logs/gwas/cpi_ldsc_dispatch.log" 2>&1
if test -e "$R/00_admin/PAUSE_REQUESTED"; then printf 'PAUSE_REQUESTED; no new stage\n'; exit 75; fi
test -s "$R/02_qc/gwas/UKB_CPI_WANG2025.full_qc.done"
test -s "$R/02_qc/gwas/UKB_CPI_WANG2025.done"
/analysis/miniconda3/bin/python "$R/16_workflow/replan_gwas_register_cpi_ldsc.py"
set +e
/analysis/miniconda3/bin/python "$R/16_workflow/replan_run_task.py" ldsc_cpi
rc=$?
set -e
/analysis/miniconda3/bin/python "$R/16_workflow/replan_gwas_ldsc_publish.py"
printf '%s\n' "$rc" > "$R/logs/gwas/cpi_ldsc_dispatch.exitcode"
exit "$rc"
