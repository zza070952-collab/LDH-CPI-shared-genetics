"""Read the frozen MATLAB output on Xiyou. Write display-bin representatives to stdout."""
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
from scipy.io import loadmat

pair = Path(sys.argv[1])
core = pair / "results/result_core.mat"
with core.open("rb") as source:
    h = hashlib.sha256()
    for chunk in iter(lambda: source.read(2**20), b""):
        h.update(chunk)
qc = json.loads((pair / "CORE_QC.json").read_text())
assert h.hexdigest() == qc["result_sha256"]
m = loadmat(core, variable_names=["fdrmat", "chrnumvec", "posvec"])
f, chrom, bp = [m[k].ravel() for k in ("fdrmat","chrnumvec","posvec")]
finite = np.isfinite(f) & (f > 0) & (f <= 1)
keep = finite & (chrom >= 1) & (chrom <= 22) & ~((chrom == 6) & (bp >= 25000000) & (bp <= 34000000))
assert (keep & (f < .05)).sum() == 762
i = np.flatnonzero(keep & (f >= .05))
# One true variant coordinate per 2-Mb by 0.05-log10 unit display bin.
bins = np.stack((chrom[i].astype(int), bp[i].astype(int)//2000000, (-np.log10(f[i])/.05).astype(int)), axis=1)
_, first, counts = np.unique(bins, axis=0, return_index=True, return_counts=True)
print("CHR\tBP\tconjFDR\tvariants_in_bin\tcore_sha256")
for j, count in zip(i[first], counts):
    print(f"{int(chrom[j])}\t{int(bp[j])}\t{f[j]:.14g}\t{count}\t{h.hexdigest()}")
