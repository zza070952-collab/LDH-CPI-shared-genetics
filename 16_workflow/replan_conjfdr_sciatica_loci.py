"""Exact same preregistered LD locus algorithm applied to the new LOW sciatica core."""
import difflib
from pathlib import Path
from replan_run_task import digest
R=Path(__file__).resolve().parents[1];base=R/'16_workflow/replan_conjfdr_loci.py'
assert digest(base)=='5ac907d856be6e2a613285439d530e7f6fda052b1564f31f787909cdee1b9c8d'
src=base.read_text();new=src.replace('UKB_CPI_WANG2025','PANUKB_SCIATICA_1476_EUR').replace('LDH_CPI_low_overlap','LDH_SCIATICA_PANUKB_low_overlap').replace('conjfdr_cpi','conjfdr_pan_sciatica').replace('CPI','SCIATICA')
o=R/'05_conjfdr/LDH_SCIATICA_PANUKB_low_overlap';o.mkdir(exist_ok=True)
(o/'LOCUS_SCOPE_ONLY.diff').write_text(''.join(difflib.unified_diff(src.splitlines(True),new.splitlines(True),fromfile=str(base),tofile='SCIATICA_loci')))
exec(compile(new,str(base),'exec'),{'__file__':__file__,'__name__':'__main__'})
