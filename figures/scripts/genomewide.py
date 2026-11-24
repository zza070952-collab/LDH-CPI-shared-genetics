"""Figure 2. Statistics, fitted boundary diagrams and reciprocal conditional enrichment."""
from plot_core import *
from matplotlib.patches import Circle

def build():
    f=fig(178)
    title(f,.029,.967,"a","Observed-scale SNP heritability")
    title(f,.524,.967,"b","Cross-trait genetic correlation")
    g=read("Figure_1_LDSC.tsv").sort_values("display_order")
    primary=g.iloc[0]
    ax=axis(f,.109,.799,.332,.124)
    for x,(v,se,col) in enumerate([(primary.LDH_h2,primary.LDH_h2_SE,BLUE),(primary.pain_h2,primary.pain_h2_SE,COPPER)]):
        ax.vlines(x,0,v,color=col,lw=.8)
        ax.errorbar(x,v,yerr=1.96*se,marker="o",color=col,ms=4,capsize=2,lw=.8)
    ax.set(xlim=(-.5,1.5),ylim=(0,.13),xticks=[0,1],xticklabels=["LDH","CPI"],yticks=[0,.05,.10],ylabel="SNP heritability")
    ax=axis(f,.683,.800,.273,.123)
    for i,r in enumerate(g.itertuples()):
        ax.errorbar(r.rg,2-i,xerr=1.96*r.se,fmt="D" if i==0 else "o",color=TEAL if i==0 else GRAY,ms=4.5,capsize=2,lw=.95)
    ax.set(xlim=(-.04,1),ylim=(-.55,2.55),yticks=range(3),yticklabels=["Sciatica","Multisite pain","CPI"],xticks=[0,.5,1],xlabel="Genetic correlation (95% CI)")
    ax.axvline(0,color=GRAY,ls=":",lw=.65);ax.spines["left"].set_visible(False);ax.tick_params(axis="y",length=0)
    title(f,.029,.724,"c","MiXeR overlap at the LDH complete-sharing boundary")
    m=read("Figure_1_MiXeR.tsv")
    for k,r in enumerate(m.itertuples()):
        ax=axis(f,.035+.319*k,.536,.306,.163);ax.set_aspect("equal");ax.axis("off")
        ax.set(xlim=(0,2.58),ylim=(0,1.53))
        outer=.66*np.sqrt(r.nc_pain_total/m.nc_pain_total.max())
        inner=.66*np.sqrt(r.shared_nc12/m.nc_pain_total.max())
        cx=.70;cy=.73
        ax.add_patch(Circle((cx,cy),outer,fc="#F2DAC9",ec=COPPER,lw=.85))
        ax.add_patch(Circle((cx-outer+inner+.04,cy),inner,fc="#B8D9D9",ec=BLUE,lw=.85))
        ax.text(1.43,1.17,f"Seed {r.seed}",fontsize=6.9,fontweight="bold")
        ax.text(1.43,.87,f"Shared  {r.shared_nc12:,.0f}",fontsize=6.7)
        ax.text(1.43,.59,f"CPI only  {r.pain_unique_nc2:,.0f}",fontsize=6.7)
        ax.text(1.43,.30,f"Dice  {r.Dice:.3f}",fontsize=6.7,color=TEAL,fontweight="bold")
    text(f,.050,.529,"Blue outline, LDH; copper outline, CPI. All three fits place fewer than 0.001 variants in the LDH-only component.",fontsize=6.55)
    title(f,.029,.499,"d","Reciprocal conditional Q–Q enrichment")
    qq=read("LDH_CPI_low_overlap.tsv")
    stamp=json.loads((DATA/"LDH_CPI_low_overlap.done.json").read_text())
    assert hashlib.sha256((DATA/"LDH_CPI_low_overlap.tsv").read_bytes()).hexdigest()==stamp["output_sha256"]
    assert qq.n_preserved_prunes.eq(500).all()
    colors=[GRAY,"#92BEC8","#4879A0",TEAL]
    handles=[]
    for k,(direction,label) in enumerate([("LDH_given_CPI","LDH conditioned on CPI"),("CPI_given_LDH","CPI conditioned on LDH")]):
        a=axis(f,.108+k*.485,.306,.365,.155)
        a.plot([0,7],[0,7],ls=":",color="#9DA9AF",lw=.8)
        for j,(t,d) in enumerate(qq[qq.direction.eq(direction)].groupby("conditioning_logP_threshold")):
            h=d[d.mean_conditional_survival.gt(0)]
            line,=a.plot(-np.log10(h.mean_conditional_survival),h.nominal_neglogP,color=colors[j],lw=1.1,label=["All","≤ 0.1","≤ 0.01","≤ 0.001"][j])
            if k==0:handles.append(line)
        a.set(xlim=(0,6.5),ylim=(0,7.2),xticks=[0,2,4,6],yticks=[0,2,4,6],xlabel="Empirical −log₁₀(quantile)")
        a.set_title(label,loc="left",color=[BLUE,COPPER][k],fontsize=7.3,pad=2)
        if k==0:a.set_ylabel("Nominal −log₁₀(P)")
    f.legend(handles=handles,loc="center right",bbox_to_anchor=(.977,.508),ncol=4,handlelength=1.35,columnspacing=.85,title="Conditioning P",title_fontsize=6.5)

    title(f,.029,.220,"e","Discovery variants within the 18 fixed regions")
    ax=axis(f,.11,.073,.847,.108)
    snps=read("SIGNIFICANT_SNP_INPUT.tsv")
    for c,part in snps.groupby("CHR"):
        ax.scatter(xgen(part.CHR,part.BP),-np.log10(part.conjFDR),s=4.1,c=[GRAY if c%2 else "#90A2AC"],lw=0)
    ax.scatter(xgen(leads.chromosome,leads.lead_position),-np.log10(leads.conjFDR),s=15,facecolors="white",edgecolors=TEAL,lw=.8)
    ax.axhline(-np.log10(.05),color=GRAY,ls="--",lw=.7)
    ax.set(xlim=(0,sum(lengths+15)),ylim=(1.13,6),yticks=[2,4,6],xticks=offsets+lengths/2,xticklabels=[str(c) if c<17 or c in [18,20,22] else "" for c in range(1,23)],xlabel="Chromosome (GRCh37)",ylabel="−log₁₀(conjFDR)")
    ax.legend(handles=[point_legend(GRAY,"o","Passing SNP"),point_legend(TEAL,"o","Fixed lead",True)],loc="upper left",ncol=2,handletextpad=.3,borderaxespad=.15)
    save(f,2)

if __name__=="__main__":
    build()

