"""Figure 3. Genome-wide landscape and aligned per-lead evidence."""
from plot_core import *

def build():
    f=fig(203)
    title(f,.025,.973,"a","Shared association landscape")
    ax=axis(f,.094,.671,.875,.269)
    bg=read("GENOME_DISPLAY_BACKGROUND.tsv")
    for c,p in bg.groupby("CHR"):
        ax.scatter(xgen(p.CHR,p.BP),-np.log10(p.conjFDR),s=.85,c="#A6B1BA" if c%2 else "#CAD1D5",lw=0,zorder=1)
    sig=read("SIGNIFICANT_SNP_INPUT.tsv")
    ax.scatter(xgen(sig.CHR,sig.BP),-np.log10(sig.conjFDR),s=2.8,color="#6C92A4",lw=0,zorder=2)
    ax.axhline(-np.log10(.05),ls="--",color=GRAY,lw=.65,zorder=0)
    for r in leads.itertuples():
        xp=xgen(r.chromosome,r.lead_position);yp=-np.log10(r.conjFDR)
        priority=r.lead_snp in ("rs11921590","rs62100775")
        ax.scatter(xp,yp,marker="D" if priority else "o",s=19 if priority else 11,facecolor="white",edgecolor=TEAL if priority else BLUE,lw=.8,zorder=4)
        if priority:
            ax.annotate(r.lead_snp,(xp,yp),xytext=(9,7),textcoords="offset points",fontsize=6.8,color=TEAL,arrowprops=dict(arrowstyle="-",lw=.55,color=TEAL))
    ax.set(xlim=(0,sum(lengths+15)),ylim=(0,6.1),yticks=range(0,7),xticks=offsets+lengths/2,xticklabels=[str(i) if i<17 or i in (18,20,22) else "" for i in range(1,23)],ylabel="−log₁₀(conjFDR)",xlabel="Chromosome (GRCh37)")
    ax.legend(handles=[point_legend(BLUE,"o","Regional lead",True),point_legend(TEAL,"D","Positive local rg and regional PP4 ≥ 0.80",True)],loc="upper left",ncol=1,handletextpad=.4)
    title(f,.025,.592,"b","Discovery leads")
    title(f,.517,.592,"c","MVP exact SNP")
    title(f,.770,.592,"d","Prior reports")
    colx=[.037,.102,.247,.348,.450]
    for x,h in zip(colx,["Chr","Lead SNP","Mb","conjFDR","Sign"]):
        text(f,x,.555,h,fontweight="bold",fontsize=6.7)
    ax=axis(f,.552,.118,.173,.422)
    ax.axvline(-np.log10(.05/18),ls=":",lw=.75,color=GRAY)
    ax.set(ylim=(17.5,-.5),xlim=(0,15.5),xticks=[0,5,10,15],yticks=[],xlabel="−log₁₀(MVP P)")
    ax.spines["left"].set_visible(False)
    labels={"PUBLISHED_MATCH_BOTH_PRIMARY_STUDIES":("Both primary",TEAL),"PUBLISHED_PRIMARY_LDH_LEAD":("LDH only",BLUE),"PUBLISHED_PRIMARY_CPI_REGION":("CPI only",COPPER),"PUBLISHED_RELATED_TRAIT_REGION":("Related trait","#64788B"),"NO_MATCH_IN_BOUNDED_PUBLISHED_TABLES":("No match",GRAY)}
    for i,r in enumerate(leads.itertuples()):
        y=.540-(i+.5)*.422/18
        if i%2==0:
            f.add_artist(Rectangle((.029,y-.0114),.939,.0228,transform=f.transFigure,fc="#F2F5F6",lw=0,zorder=-2))
        for x,v in zip(colx[:4],[r.chromosome,r.lead_snp,f"{r.lead_position/1e6:.3f}",f"{r.conjFDR:.4g}"]):
            text(f,x,y,str(v),va="center",fontsize=6.75)
        neg=r.effect_direction=="DISCORDANT"
        text(f,.454,y,"−" if neg else "+",va="center",ha="center",color=RED if neg else BLUE,fontsize=8,fontweight="bold")
        if pd.isna(r.MVP_P):
            ax.text(.10,i,"Unavailable",fontsize=6.1,va="center",color=GRAY,bbox=dict(fc="white",pad=.35,ec="none"))
        else:
            val=-np.log10(r.MVP_P);ax.hlines(i,0,val,color="#DAB6A2",lw=1)
            ax.plot(val,i,"D",ms=3.5,color=COPPER)
        desc,col=labels[ledger.iloc[i].novelty_category]
        text(f,.793,y,desc,va="center",fontsize=6.7,color=col)
    text(f,.035,.045,"+  Concordant discovery effects     −  Discordant effects      |      MVP export: 5 exact matches, 13 unavailable.",fontsize=6.8)
    text(f,.035,.022,"Prior categories refer to a bounded regional search. Exact MVP effect estimates were unsigned.",fontsize=6.6,color=GRAY)
    save(f,3)

if __name__=="__main__":
    build()
