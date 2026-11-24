"""Figure 4. Cross-trait regional and block-level evidence."""
from plot_core import *
def build():
    f=fig(190);reg=read("Figure_3_regional_GWAS.tsv")
    for k,c in enumerate([3,18]):
        rec=ledger[(ledger.chromosome==c)&ledger.locus_id.str.endswith("_001")].iloc[0]
        title(f,.025+.50*k,.965,"ab"[k],f"Chromosome {c}  ·  {rec.lead_snp}")
        a=axis(f,.094+k*.495,.685,.375,.230);r=reg[reg.locus_id==rec.locus_id]
        x=r.position_GRCh37/1e6
        for side,tr,col in [(1,"LDH",BLUE),(-1,"CPI",COPPER)]:
            y=-side*np.log10(r[tr+"_P"])
            a.scatter(x,y,s=2.8,alpha=.78,color=col,linewidths=0)
            sel=r.SNP==rec.lead_snp
            a.scatter(x[sel],y[sel],s=14,facecolors="white",edgecolors=col,marker="D",linewidths=.8,zorder=5)
        a.axvline(rec.lead_position/1e6,color=GRAY,ls=":",lw=.8);a.axhline(0,color="#A7B1B9",lw=.6)
        a.set_xlim(rec.start/1e6,rec.end/1e6);a.set_ylim(-13.9,13.9)
        a.set_yticks([-12,-8,-4,0,4,8,12],["12","8","4","0","4","8","12"]);a.set_xlabel("Position (Mb, GRCh37)")
        if k==0:a.set_ylabel("−log₁₀(GWAS P)")
        a.text(.02,.93,"LDH",transform=a.transAxes,color=BLUE,fontweight="bold",fontsize=6.5);a.text(.02,.04,"CPI",transform=a.transAxes,color=COPPER,fontweight="bold",fontsize=6.5)
    title(f,.025,.584,"c","Local genetic correlation")
    title(f,.565,.584,"d","LDH–CPI colocalization")
    a=axis(f,.147,.385,.255,.153)
    a.axvline(0,ls=":",lw=.75,color=GRAY)
    for i,r in enumerate(lava.itertuples()):
        y=4-i;col=RED if r.rg_local<0 else TEAL if r.FDR_BH<.05 else GRAY
        a.errorbar(r.rg_local,y,xerr=[[r.rg_local-r.rho_lower],[r.rho_upper-r.rg_local]],fmt="o",color=col,ms=4,capsize=2,lw=1)
        a.text(1.04,y,f"{r.FDR_BH:.4f}",transform=a.get_yaxis_transform(),va="center",fontsize=6.5)
    a.text(1.04,1.07,"BH q",transform=a.transAxes,fontsize=6.6,fontweight="bold")
    a.set(xlim=(-1.05,1.05),ylim=(-.55,4.55),xticks=[-1,-.5,0,.5,1],yticks=range(4,-1,-1),yticklabels=[f"Chr{r.chr}  /  {r.lava_block_id}" for r in lava.itertuples()],xlabel="Local genetic correlation (95% CI)")
    a.spines["left"].set_visible(False);a.tick_params(axis="y",length=0)
    a=axis(f,.706,.385,.249,.153)
    row3=ledger[ledger.lead_snp=="rs11921590"].iloc[0];row18=ledger[ledger.lead_snp=="rs62100775"].iloc[0]
    for y,v,shape in [(3,row3.PP4,"o"),(2,row3.PP4_signal,"D"),(1,row18.PP4,"o")]:
        a.hlines(y,0,v,color=LIGHT,lw=1);a.plot(v,y,marker=shape,color=TEAL,ms=4.6)
        a.annotate(f"{v:.3f}",(v,y),xytext=(-5,1),ha="right",va="bottom",textcoords="offset points",fontsize=6.4)
    a.text(.03,0,"No CPI credible set",va="center",color=GRAY,fontsize=6.5)
    a.set(ylim=(-.55,3.55),yticks=[3,2,1,0],yticklabels=["Chr3 region","Signal pair¹","Chr18 region","Signal pair"],xlim=(0,1.04),xticks=[0,.25,.5,.75,1],xlabel="Posterior probability (PP4)")
    a.spines["left"].set_visible(False);a.tick_params(axis="y",length=0)

    title(f,.025,.295,"e","Evidence in all 18 fixed discovery regions")
    ax=axis(f,.171,.087,.802,.170)
    names=["−log₁₀(conjFDR)","Regional PP4","Signal-pair PP4","Local rg","Local BH q","Local h²"]
    def fmt(v):
        if pd.isna(v):return "—"
        if 0<v<.01:return "<.01"
        return f"{v:.2f}".replace("0.",".").replace("-","−")
    for i,r in enumerate(ledger.itertuples()):
        test=r.eligibility=="BOTH_ELIGIBLE"
        h2={"LDH_ONLY":"L","CPI_ONLY":"C","BOTH_ELIGIBLE":"B"}[r.eligibility]
        rows=[(f"{-np.log10(r.conjFDR):.1f}",r.conjFDR<.01,False),(fmt(r.PP4),r.PP4>=.8,False),(fmt(r.PP4_signal),r.PP4_signal>=.8,False),(fmt(r.local_rg),r.local_rg_FDR<.05,r.local_rg<0),(fmt(r.local_rg_FDR),r.local_rg_FDR<.05,r.local_rg<0),(h2,test,False)]
        for j,(s,highlight,negative) in enumerate(rows):
            col=(RED if negative else TEAL) if highlight else "#E5ECEE" if s!="—" else "#F5F7F7"
            ax.add_patch(Rectangle((i-.5,j-.5),1,1,fc=col,ec="white",lw=.55))
            ax.text(i,j,s,va="center",ha="center",fontsize=5.85,color="white" if highlight else INK)
    reglabels=[r.locus_id.split("chr")[1].replace("_00",".") for r in ledger.itertuples()]
    ax.set(xlim=(-.5,17.5),ylim=(5.5,-.5),xticks=range(18),xticklabels=reglabels,yticks=range(6),yticklabels=names)
    ax.tick_params(labelsize=6.15,length=0)
    for sp in ax.spines.values():sp.set_visible(False)
    text(f,.027,.041,"Local h²: B, both traits; L, LDH only; C, CPI only.  —  No estimate available. Columns identify chromosome.region.",fontsize=6.6)
    text(f,.027,.018,"¹ The chromosome 3 signal-pair estimate comes from the overlapping COL7A1-labelled window.",fontsize=6.6)
    save(f,4)

if __name__=="__main__":
    build()

