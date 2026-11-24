"""Figure 5. Matched tissue predictors, donor-level detection and follow-up rationale."""
from plot_core import *
def build():
    f=fig(165);s=read("Figure_4_expression_and_status.tsv")
    genes=s.sort_values("display_order").gene.drop_duplicates().tolist()
    ts=["Nerve_Tibial","Brain_Cortex","Brain_Spinal_cord_cervical_c-1","Whole_Blood"]
    title(f,.025,.968,"a","Candidate evidence across scales")
    for x,head,detail in [(.057,"Shared region","18 fixed intervals"),(.291,"External pain","5 exact MVP leads"),(.518,"Predicted expression","4 GTEx tissues"),(.779,"Tissue and function","2 adult DRG datasets")]:
        text(f,x,.923,head,fontsize=7.6,fontweight="bold",color=TEAL);text(f,x,.893,detail,fontsize=6.8)
    ruler(f,.029,.976,.870)
    title(f,.025,.834,"b","Predicted expression in the same tissue")
    for k,(fam,label) in enumerate([("mashr","MASHR  ·  primary"),("en","Elastic net  ·  sensitivity")]):
        ax=axis(f,.145+k*.455,.485,.332,.277)
        ax.set_title(label,loc="left",fontsize=7.5,pad=6)
        for ix in range(len(ts)):
            ax.axvline(ix,color=LIGHT,lw=.5,zorder=0)
        for iy,g in enumerate(genes):
            ax.axhline(iy,color=LIGHT,lw=.45,zorder=0)
            for ix,t in enumerate(ts):
                d=s[(s.gene==g)&(s.family==fam)&(s.tissue==t)].iloc[0]
                if d.status_LDH=="NO_MODEL":ax.text(ix,iy,"—",fontsize=7,color="#99A3A8",va="center",ha="center");continue
                if pd.isna(d.q_conjunction):ax.plot(ix,iy,"x",ms=3.0,c=GRAY,mew=.65);continue
                joint=d.q_conjunction<.05
                separate=max(d.trait_q_trait_global4tissue_LDH,d.trait_q_trait_global4tissue_CPI)<.05
                # Circle area reflects the unadjusted conjunction P on a log scale.
                area=12+8*(-np.log10(d.P_conjunction))
                col=TEAL if joint and d.same_tissue_Z_direction=="SAME_SIGN" else RED if joint else BLUE if separate else "#A9B7C0"
                ax.scatter(ix,iy,s=area,marker="D" if joint else "o",facecolors="white" if separate and not joint else col,edgecolors=col,linewidths=1.05 if separate else .55,zorder=3)
        ax.set(xlim=(-.45,3.45),ylim=(len(genes)-.48,-.5),xticks=range(4),xticklabels=["Tibial nerve","Cortex","Spinal cord","Blood"],yticks=range(len(genes)),yticklabels=genes if k==0 else [])
        if k==0:
            for tt in ax.get_yticklabels():tt.set_fontstyle("italic")
        ticks_off(ax)
        for sp in ax.spines.values():sp.set_visible(False)
    for x,col,shape,label,open_ in [(.063,TEAL,"D","Conjunction q < 0.05",False),(.324,BLUE,"o","Both trait q < 0.05 only",True),(.615,"#A9B7C0","o","Other eligible",False),(.804,GRAY,"x","Ineligible",False)]:
        f.add_artist(Line2D([x],[.422],transform=f.transFigure,ls="",marker=shape,ms=5,mec=col,mew=.9,mfc="white" if open_ else col))
        text(f,x+.013,.422,label,va="center",fontsize=6.7)
    text(f,.053,.387,"—  No model",fontsize=6.7,color=GRAY);text(f,.222,.387,"Marker area increases with −log₁₀(common-predictor P).",fontsize=6.7)
    ruler(f,.029,.976,.363)
    title(f,.025,.326,"c","Adult DRG transcript detection")
    title(f,.555,.326,"d","DCC functional follow-up")
    a=axis(f,.127,.089,.37,.193)
    drg=read("DRG_CANDIDATE_CLASS_SUMMARY.tsv")
    series=["DRG Neuron","Satellite glia","Myelinating Schwann cell","NeuN-enriched neuronal nuclei"]
    mat=drg[drg.cell_class.isin(series)&drg.gene.isin(genes)].pivot(index="gene",columns="cell_class",values="median_percent_expressing").loc[genes,series]
    cm=LinearSegmentedColormap.from_list("detection",["#F1F5F5","#98C6BD","#327D77","#174F4F"])
    vals=mat.to_numpy();upper=max(40,np.ceil(np.nanmax(vals)/10)*10)
    a.imshow(vals,aspect="auto",cmap=cm,vmin=0,vmax=upper,interpolation="nearest")
    for i in range(len(genes)):
        for j in range(4):
            v=vals[i,j];a.text(j,i,"—" if pd.isna(v) else f"{v:.1f}",ha="center",va="center",color="white" if v>upper*.54 else INK,fontsize=6.1)
    a.set(yticks=range(len(genes)),yticklabels=genes,xticks=range(4),xticklabels=["Neuron","Satellite\nglia","Schwann","NeuN+¹"])
    a.set_xticks(np.arange(-.5,4,1),minor=True);a.set_yticks(np.arange(-.5,7,1),minor=True)
    a.grid(which="minor",color="white",lw=.65);a.tick_params(which="minor",length=0)
    a.axvline(2.5,color="white",lw=2.4);ticks_off(a)
    for tt in a.get_yticklabels():tt.set_fontstyle("italic")
    for sp in a.spines.values():sp.set_visible(False)
    dc=ledger[ledger.lead_snp=="rs62100775"].iloc[0];exact=leads[leads.lead_snp==dc.lead_snp].iloc[0]
    items=[("Chr18 region",f"PP4 {dc.PP4:.3f}  ·  local q {dc.local_rg_FDR:.3f}"),
           ("External pain",f"{dc.lead_snp}  ·  MVP P {exact.MVP_P:.2g}"),
           ("Receptor context","Adult neuronal and glial expression"),
           ("Functional rationale","Netrin receptor implicated in pain biology")]
    for i,(h,v) in enumerate(items):
        y=.278-i*.053;text(f,.592,y,h,fontsize=6.7,fontweight="bold",color=TEAL)
        text(f,.592,y-.023,v,fontsize=6.95)
    text(f,.592,.055,"Shared expression mediation remains unresolved.",fontsize=6.8,fontstyle="italic")
    text(f,.038,.016,"Median donor-level detection (%).  ¹ GSE168243 NeuN-enriched nuclei. Other columns: GSE201586.",fontsize=6.4)
    save(f,5)


if __name__=="__main__":
    build()

