"""Figure 1. Four-tier study diagram modeled on the JTM workflow."""
from plot_core import *
from matplotlib.patches import FancyArrowPatch

def build():
    f = fig(165)
    def box(x,y,w,h,heading,detail="",fill="#E8F0F3",edge="#637984",fs=7.5):
        f.add_artist(Rectangle((x,y),w,h,transform=f.transFigure,facecolor=fill,edgecolor=edge,lw=.6,zorder=.5))
        text(f,x+w/2,y+h*(.68 if detail else .5),heading,ha="center",va="center",fontsize=fs,fontweight="bold",zorder=3)
        if detail:
            text(f,x+w/2,y+h*.31,detail,ha="center",va="center",fontsize=6.45,linespacing=1.22,zorder=3)
    def arrow(x1,y1,x2,y2):
        f.add_artist(FancyArrowPatch((x1,y1),(x2,y2),transform=f.transFigure,arrowstyle="-|>",mutation_scale=7,lw=.72,color="#617580",zorder=1))
    def connector(points):
        x,y=zip(*points)
        f.add_artist(Line2D(x,y,transform=f.transFigure,lw=.65,color="#617580"))
    def tier(y,letter,heading):
        title(f,.028,y,letter,heading)
    # Shared frame and uninterrupted horizontal rules follow the four-stage template.
    f.add_artist(Rectangle((.015,.018),.970,.967,transform=f.transFigure,fill=False,edgecolor=GRAY,lw=.75))
    for y in (.674,.415,.173):
        f.add_artist(Line2D([.015,.985],[y,y],transform=f.transFigure,color=GRAY,lw=.65))
    tier(.952,"a","Discovery cohorts and genome-wide sharing")
    box(.037,.862,.433,.065,"Lumbar disc herniation  ·  FinnGen R12","35,290 cases and 353,224 controls",fill="#DDEAF4",edge=BLUE)
    box(.531,.862,.431,.065,"Chronic pain intensity  ·  UK Biobank","134,627 participants",fill="#F7E8DD",edge=COPPER)
    connector([(.25,.862),(.25,.847),(.749,.847),(.749,.862)])
    arrow(.50,.847,.50,.829)
    for x,head in [(.04,"LDSC"),(.351,"MiXeR"),(.663,"Conditional / conjunctional FDR")]:
        box(x,.799,.298,.031,head,fill="#C5DCE0",fs=7.2)
    for x in (.189,.812):
        connector([(.50,.840),(x,.840)]);arrow(x,.840,x,.830)
    # Miniatures in the first tier use the actual estimates and passing variants.
    a=axis(f,.066,.716,.241,.050);r=read("Figure_1_LDSC.tsv").iloc[0]
    a.errorbar(r.rg,0,xerr=1.96*r.se,fmt="D",color=TEAL,capsize=2,ms=4)
    a.axvline(0,lw=.5,ls=":",color=GRAY);a.set(xlim=(-.05,.8),ylim=(-1,1),yticks=[],xticks=[0,.4,.8]);a.tick_params(labelsize=5.5)
    a.spines["left"].set_visible(False)
    a=axis(f,.374,.708,.247,.073);mx=read("Figure_1_MiXeR.tsv")
    for i,r in enumerate(mx.itertuples()):
        a.vlines(i,0,r.Dice,color="#8DB6B5",lw=1.3);a.plot(i,r.Dice,"o",color=TEAL,ms=3.5)
    a.set(xlim=(-.55,2.55),ylim=(0,.45),yticks=[],xticks=[]);a.spines["left"].set_visible(False)
    a=axis(f,.686,.708,.254,.069);s=read("SIGNIFICANT_SNP_INPUT.tsv")
    a.scatter(xgen(s.CHR,s.BP),-np.log10(s.conjFDR),s=1,c=TEAL,lw=0)
    a.set(xlim=(0,sum(lengths+15)),ylim=(1,6),xticks=[],yticks=[]);a.spines["left"].set_visible(False)
    for x,label in [(.186,"Genome-wide correlation"),(.498,"Overlap across three seeds"),(.812,"Regions shared by both traits")]:
        text(f,x,.682,label,ha="center",fontsize=6.65)
    tier(.641,"b","Follow-up of fixed association regions")
    box(.036,.491,.217,.084,"18 fixed regions","LD-based grouping\nof discovery variants")
    arrow(.253,.534,.294,.534)
    connector([(.294,.470),(.294,.587)])
    for y,heading,sub in [(.587,"LAVA","Local heritability and correlation"),(.528,"Colocalization","Regional and signal-pair posteriors"),(.470,"MVP pain intensity","Exact-SNP external associations")]:
        arrow(.294,y,.328,y);box(.328,y-.022,.225,.044,heading,fill="#C5DCE0",fs=7)
        arrow(.553,y,.602,y)
        text(f,.619,y,sub,va="center",fontsize=7.1)
    text(f,.50,.426,"Prior LDH and pain reports provide context for the discovery regions.",ha="center",fontsize=6.8)
    tier(.382,"c","Gene-level associations and adult sensory tissue")
    box(.035,.257,.234,.087,"Gene mapping","Positional candidates\nand MAGMA context")
    arrow(.269,.301,.305,.301);connector([(.305,.244),(.305,.334)])
    for y,head,sub in [(.334,"S-PrediXcan","MASHR and elastic-net models"),(.246,"Adult human DRG","Donor-level transcript detection")]:
        arrow(.305,y,.342,y);box(.342,y-.034,.277,.068,head,sub,fill="#DDEBED")
        arrow(.619,y,.656,y)
    connector([(.656,.246),(.656,.334)]);arrow(.656,.290,.697,.290)
    box(.697,.247,.266,.086,"Candidate follow-up","Region, expression and\nreceptor function",fill="#E4F1EE",edge=TEAL)
    text(f,.50,.186,"Gene attribution and regional association are evaluated separately.",ha="center",fontsize=6.8)
    tier(.139,"d","Interpretation and functional follow-up")
    box(.034,.047,.240,.067,"Shared architecture","LDH and population CPI")
    arrow(.274,.080,.310,.080)
    box(.310,.047,.264,.067,"Prioritized regions","Chromosomes 3 and 18",fill="#E4F1EE",edge=TEAL)
    arrow(.574,.080,.609,.080)
    box(.609,.037,.352,.088,"DCC receptor follow-up","Shared expression mediation\nremains unresolved",fill="#E4F1EE",edge=TEAL)
    save(f,1)

if __name__ == "__main__":
    build()
