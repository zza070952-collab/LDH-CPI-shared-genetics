"""Regenerate the five main figures from the bundled frozen plotting tables."""
from pathlib import Path
import hashlib, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.colors import LinearSegmentedColormap

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/"source_data"
BLUE="#356A99"; COPPER="#B75E32"; TEAL="#087F78"; RED="#A14352"
INK="#24313B"; GRAY="#687581"; LIGHT="#E4E9EB"; PALE="#F2F5F5"
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":7.2,
 "axes.labelsize":7.2,"axes.edgecolor":GRAY,"axes.labelcolor":INK,"text.color":INK,
 "xtick.color":INK,"ytick.color":INK,"axes.spines.right":False,"axes.spines.top":False,
 "axes.linewidth":.55,"xtick.major.width":.55,"ytick.major.width":.55,
 "xtick.major.size":2.5,"ytick.major.size":2.5,"xtick.labelsize":6.7,"ytick.labelsize":6.7,
 "legend.frameon":False,"legend.fontsize":6.5,"pdf.fonttype":42,"svg.fonttype":"none",
 "svg.hashsalt":"ldh-cpi-figures-20261122","figure.facecolor":"white","savefig.facecolor":"white"})
def read(name):return pd.read_csv(DATA/name,sep="\t")
def fig(height):return plt.figure(figsize=(183/25.4,height/25.4))
def axis(f,x,y,w,h):return f.add_axes([x,y,w,h])
def title(f,x,y,letter,name):
    f.text(x,y,letter,fontweight="bold",fontsize=10.2,va="baseline",ha="left",color="#161C21")
    f.text(x+.032,y,name,fontweight="bold",fontsize=8.25,va="baseline",ha="left")
def text(f,x,y,s,**kw):return f.text(x,y,s,**kw)
def ruler(f,x0,x1,y):
    f.add_artist(Line2D([x0,x1],[y,y],transform=f.transFigure,color=LIGHT,lw=.65))
NAMES = ["workflow", "genomewide", "shared_regions", "local_evidence", "candidate_genes"]
def save(f, i):
    for ext in ("pdf", "svg", "png"):
        f.savefig(ROOT / f"Figure_{i}_{NAMES[i-1]}_strictTemplate.{ext}", dpi=600)
    plt.close(f)
def ticks_off(ax):ax.tick_params(length=0)
def point_legend(color,shape,label,open_point=False):
    return Line2D([],[],ls="",marker=shape,markersize=4.3,markeredgewidth=.9,
                  markerfacecolor="white" if open_point else color,markeredgecolor=color,label=label)
leads=read("Figure_2_fixed_leads_and_MVP.tsv")
ledger=read("Figure_3_all_18_ledger.tsv")
lava=read("Figure_3_LAVA_all_tested.tsv")
regions=pd.read_csv(DATA/"hg19.chrom.sizes",sep="\t",header=None,names=["chrom","bp"]).set_index("chrom").bp
lengths=np.array([regions[f"chr{i}"] for i in range(1,23)])/1e6
offsets=np.r_[0,np.cumsum(lengths+15)[:-1]]
def xgen(c,p):return offsets[np.array(c,dtype=int)-1]+np.array(p)/1e6
assert len(leads)==18 and list(leads.lead_snp)==list(ledger.lead_snp)
