#!/bin/bash
#SBATCH --nodes=1
#SBATCH --job-name="MAFFT_iqtree"
#SBATCH --account="aggstrategic1"
#SBATCH --partition="batch"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=128GB
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=james.o\'dwyer@agriculture.vic.gov.au
#SBATCH --time=12:0:00

# Build a taxon-matched whole-genome baseline tree and quantitatively compare it
# with the existing DP>=10 partitioned core-CDS tree.


module load MAFFT/7.526-GCC-13.3.0-with-extensions

module load IQ-TREE/2.3.6-gompi-2024a

#module load Python/3.12.3-GCCcore-13.3.0

#module load Python-bundle-PyPI


eval "$(conda shell.bash hook)"
conda activate /home/vidh72t/.conda/envs/oshv1_annot


MASTER_DIR="/group/sequencing/assembly/James/OsHV_annot/Align_raws/dp10_complete_CDS_sets/mafft_iqtree_partitioned_dp10_added_SRA_nodups"
CDS_TREE="${MASTER_DIR}/figures/OsHV_dp10_partitioned_midpoint_rooted.treefile"
GENOMES="/group/sequencing/assembly/James/OsHV_annot/Align_raws/dp10_complete_CDS_sets/all_genomes_compined_forbasetree/All_genomes_fastas.fasta"
OUTDIR="${MASTER_DIR}/whole_genome_vs_dp10_cds_comparison"

THREADS="${THREADS:-24}"
MAFFT_EXE="${MAFFT_EXE:-mafft}"
IQTREE_EXE="${IQTREE_EXE:-iqtree2}"
PYTHON_EXE="${PYTHON_EXE:-python3}"

mkdir -p "${OUTDIR}"/{alignment,tree,figures,metrics,logs}

TAXA_FILE="${OUTDIR}/logs/dp10_cds_tree_taxa.txt"
FILTERED_FASTA="${OUTDIR}/alignment/whole_genomes_taxon_matched_unaligned.fasta"
ALIGNED_FASTA="${OUTDIR}/alignment/whole_genomes_taxon_matched_aligned.fasta"
BASE_PREFIX="${OUTDIR}/tree/whole_genome_baseline"
BASE_TREE="${BASE_PREFIX}.treefile"
BASE_MID="${OUTDIR}/tree/whole_genome_baseline_midpoint_rooted.treefile"
METRICS_TSV="${OUTDIR}/metrics/whole_genome_vs_dp10_cds_metrics.tsv"
PAIRWISE_TSV="${OUTDIR}/metrics/pairwise_patristic_distances.tsv"
SPLITS_TSV="${OUTDIR}/metrics/support_filtered_split_metrics.tsv"
SUMMARY_PNG="${OUTDIR}/figures/whole_genome_vs_dp10_cds_multi_metric_summary.png"
SCATTER_PNG="${OUTDIR}/figures/patristic_distance_comparison.png"
TREE_PNG="${OUTDIR}/figures/taxon_matched_tree_comparison.png"

for command in "${PYTHON_EXE}" "${MAFFT_EXE}" "${IQTREE_EXE}"; do
  command -v "${command}" >/dev/null 2>&1 || { echo "ERROR: executable not found: ${command}" >&2; exit 1; }
done
[[ -s "${CDS_TREE}" ]] || { echo "ERROR: CDS tree not found: ${CDS_TREE}" >&2; exit 1; }
[[ -s "${GENOMES}" ]] || { echo "ERROR: genome FASTA not found: ${GENOMES}" >&2; exit 1; }

# Extract exact taxa from the existing DP>=10 CDS tree.
"${PYTHON_EXE}" - "${CDS_TREE}" "${TAXA_FILE}" <<'PY'
import sys
from Bio import Phylo

tree = Phylo.read(sys.argv[1], "newick")
taxa = sorted(t.name for t in tree.get_terminals())
if len(taxa) != len(set(taxa)):
    raise SystemExit("Duplicate tip labels found in DP>=10 CDS tree")
with open(sys.argv[2], "w") as out:
    out.write("\n".join(taxa) + "\n")
print(f"Extracted {len(taxa)} taxa from DP>=10 CDS tree")
PY

# Match the unaligned genome FASTA to exactly those taxa without requiring seqkit.
"${PYTHON_EXE}" - "${GENOMES}" "${TAXA_FILE}" "${FILTERED_FASTA}" "${OUTDIR}/logs/taxon_matching_report.tsv" <<'PY'
import sys
from Bio import SeqIO

genomes, taxa_file, output, report = sys.argv[1:]
wanted = {x.strip() for x in open(taxa_file) if x.strip()}
records = list(SeqIO.parse(genomes, "fasta"))
by_id = {}
for rec in records:
    if rec.id in by_id:
        raise SystemExit(f"Duplicate FASTA identifier: {rec.id}")
    by_id[rec.id] = rec
missing = sorted(wanted - set(by_id))
extra = sorted(set(by_id) - wanted)
if missing:
    raise SystemExit("Taxa in CDS tree missing from whole-genome FASTA: " + ", ".join(missing))
selected = [by_id[name] for name in sorted(wanted)]
SeqIO.write(selected, output, "fasta")
with open(report, "w") as out:
    out.write("category\ttaxon\n")
    for name in sorted(wanted): out.write(f"retained\t{name}\n")
    for name in extra: out.write(f"excluded_not_in_cds_tree\t{name}\n")
print(f"Retained {len(selected)} whole genomes; excluded {len(extra)} records not present in CDS tree")
PY

# Align the taxon-matched whole genomes.
"${MAFFT_EXE}" --auto --thread "${THREADS}" "${FILTERED_FASTA}" \
  > "${ALIGNED_FASTA}" 2> "${OUTDIR}/logs/mafft_whole_genome.log"

# Sanity-check aligned taxa and equal alignment lengths.
"${PYTHON_EXE}" - "${ALIGNED_FASTA}" "${TAXA_FILE}" <<'PY'
import sys
from Bio import SeqIO
records = list(SeqIO.parse(sys.argv[1], "fasta"))
wanted = {x.strip() for x in open(sys.argv[2]) if x.strip()}
seen = {r.id for r in records}
lengths = {len(r.seq) for r in records}
assert seen == wanted, f"Aligned FASTA taxa differ: missing={wanted-seen}, extra={seen-wanted}"
assert len(lengths) == 1, "Aligned sequences have unequal lengths"
print(f"Alignment verified: {len(records)} taxa, {next(iter(lengths))} columns")
PY

# Infer a whole-genome baseline with equivalent support statistics.
"${IQTREE_EXE}" \
  -s "${ALIGNED_FASTA}" \
  -st DNA \
  -m MFP \
  -B 1000 \
  --alrt 1000 \
  -T "${THREADS}" \
  --prefix "${BASE_PREFIX}" \
  --redo 2>&1 | tee "${OUTDIR}/logs/iqtree_whole_genome.log"

[[ -s "${BASE_TREE}" ]] || { echo "ERROR: IQ-TREE did not produce ${BASE_TREE}" >&2; exit 1; }

# Compute taxon-matched RF, support-filtered splits, patristic metrics,
# nearest-neighbour retention, quartet concordance, and figures.
"${PYTHON_EXE}" - \
  "${BASE_TREE}" "${CDS_TREE}" "${BASE_MID}" "${METRICS_TSV}" "${PAIRWISE_TSV}" \
  "${SPLITS_TSV}" "${SUMMARY_PNG}" "${SCATTER_PNG}" "${TREE_PNG}" <<'PY'
import sys, copy, csv, itertools, math, re
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio import Phylo
from scipy.stats import pearsonr, spearmanr

(base_path, cds_path, base_mid_path, metrics_path, pairwise_path,
 splits_path, summary_png, scatter_png, tree_png) = sys.argv[1:]

def tips(t): return {x.name for x in t.get_terminals()}
def prune(t, keep):
    t=copy.deepcopy(t)
    for x in list(t.get_terminals()):
        if x.name not in keep: t.prune(x)
    return t

def support(clade):
    raw = clade.name if clade.name not in (None, "") else clade.confidence
    if raw is None: return None
    values=[float(x) for x in re.findall(r"\d+(?:\.\d+)?", str(raw))]
    if not values: return None
    values=[x/100 if x>1 else x for x in values]
    return min(values[:2]) if len(values)>=2 else values[0]

def canonical(desc, all_taxa):
    other=all_taxa-desc
    if len(desc)<2 or len(other)<2: return None
    a,b=frozenset(desc),frozenset(other)
    if len(a)<len(b): return a
    if len(b)<len(a): return b
    return min((a,b), key=lambda x: tuple(sorted(x)))

def splitmap(tree):
    all_taxa=tips(tree); out={}
    for clade in tree.find_clades(order="postorder"):
        if clade is tree.root or clade.is_terminal(): continue
        sp=canonical({x.name for x in clade.get_terminals()}, all_taxa)
        if sp is not None: out[sp]=support(clade)
    return out

def split_metrics(a,b):
    shared=len(a&b); union=len(a|b); rf=len(a-b)+len(b-a); den=len(a)+len(b)
    return dict(rf=rf, normalized_rf=rf/den if den else 0.0,
                shared=shared, baseline=len(a), cds=len(b),
                jaccard=shared/union if union else 1.0,
                precision=shared/len(b) if b else (1.0 if not a else 0.0),
                recall=shared/len(a) if a else (1.0 if not b else 0.0))

def distances(tree, labels):
    pairs=[]; vals=[]
    for i,a in enumerate(labels):
        for b in labels[i+1:]: pairs.append((a,b)); vals.append(tree.distance(a,b))
    return pairs,np.asarray(vals,float)

def nearest(tree, labels, k):
    out={}
    for a in labels:
        ds=sorted((tree.distance(a,b),b) for b in labels if b!=a)
        cutoff=ds[min(k,len(ds))-1][0]; tol=max(1e-12,abs(cutoff)*1e-9)
        out[a]={b for d,b in ds if d<=cutoff+tol}
    return out

def quartet_state(tree,q):
    a,b,c,d=q
    sums={"ab|cd":tree.distance(a,b)+tree.distance(c,d),
          "ac|bd":tree.distance(a,c)+tree.distance(b,d),
          "ad|bc":tree.distance(a,d)+tree.distance(b,c)}
    z=sorted(sums.items(),key=lambda x:x[1]); scale=max(1.0,max(abs(v) for v in sums.values()))
    return "unresolved" if abs(z[1][1]-z[0][1])<=1e-12*scale else z[0][0]

base=Phylo.read(base_path,"newick"); cds=Phylo.read(cds_path,"newick")
shared=tips(base)&tips(cds)
base=prune(base,shared); cds=prune(cds,shared)
assert tips(base)==tips(cds)
labels=sorted(shared)
base_mid=copy.deepcopy(base); base_mid.root_at_midpoint(); Phylo.write(base_mid,base_mid_path,"newick")

bm,cm=splitmap(base),splitmap(cds)
all_split=split_metrics(set(bm),set(cm))
threshold_results={}
for th in (0.5,0.8,0.9):
    threshold_results[th]=split_metrics({s for s,v in bm.items() if v is not None and v>=th},
                                        {s for s,v in cm.items() if v is not None and v>=th})

pairs,db=distances(base,labels); pairs2,dc=distances(cds,labels); assert pairs==pairs2
scale=float(np.dot(dc,db)/np.dot(dc,dc)); dcs=dc*scale
spr=spearmanr(db,dc); prs=pearsonr(db,dcs)
rmse=float(np.sqrt(np.mean((db-dcs)**2))); nrmse=rmse/float(np.mean(db))
mae=float(np.mean(np.abs(db-dcs)))
sst=float(np.sum((db-np.mean(db))**2)); r2=1-float(np.sum((db-dcs)**2))/sst

b1,c1=nearest(base,labels,1),nearest(cds,labels,1); b3,c3=nearest(base,labels,3),nearest(cds,labels,3)
top1=float(np.mean([bool(b1[x]&c1[x]) for x in labels]))
top3j=float(np.mean([len(b3[x]&c3[x])/len(b3[x]|c3[x]) for x in labels]))
base_nn_top3=float(np.mean([bool(b1[x]&c3[x]) for x in labels]))

con=dis=unres=0
for q in itertools.combinations(labels,4):
    x,y=quartet_state(base,q),quartet_state(cds,q)
    if "unresolved" in (x,y): unres+=1
    elif x==y: con+=1
    else: dis+=1
comparable=con+dis; total=comparable+unres
qcon=con/comparable if comparable else float("nan")

metrics={
 "shared_taxa":len(labels), "pairwise_distances":len(db),
 "normalized_rf_all":all_split["normalized_rf"],
 "patristic_scale_factor":scale, "patristic_spearman_r":float(spr.statistic),
 "patristic_spearman_p":float(spr.pvalue), "patristic_scaled_pearson_r":float(prs.statistic),
 "patristic_scaled_pearson_p":float(prs.pvalue), "patristic_scaled_r2":r2,
 "patristic_scaled_rmse":rmse, "patristic_scaled_nrmse_mean":nrmse,
 "patristic_scaled_mae":mae, "nearest_top1_agreement":top1,
 "nearest_top3_jaccard_mean":top3j, "baseline_nearest_in_cds_top3":base_nn_top3,
 "quartets_total":total, "quartets_comparable":comparable,
 "quartet_concordance":qcon, "quartet_discordance":dis/comparable if comparable else float("nan"),
 "quartet_unresolved_fraction":unres/total if total else float("nan")}
for th,res in threshold_results.items():
    for k,v in res.items(): metrics[f"support_{th}_{k}"]=v
with open(metrics_path,"w",newline="") as f:
    w=csv.writer(f,delimiter="\t"); w.writerow(["metric","value"]); w.writerows(metrics.items())
with open(pairwise_path,"w",newline="") as f:
    w=csv.writer(f,delimiter="\t"); w.writerow(["taxon_1","taxon_2","baseline_distance","cds_distance","scaled_cds_distance"])
    for (a,b),x,y,z in zip(pairs,db,dc,dcs): w.writerow([a,b,x,y,z])
with open(splits_path,"w",newline="") as f:
    fields=["threshold","rf","normalized_rf","shared","baseline","cds","jaccard","precision","recall"]
    w=csv.DictWriter(f,fieldnames=fields,delimiter="\t"); w.writeheader()
    w.writerow({"threshold":"all",**all_split})
    for th,res in threshold_results.items(): w.writerow({"threshold":th,**res})

# Distance scatter
fig,ax=plt.subplots(1,2,figsize=(12,5))
ax[0].scatter(db,dc,s=14,alpha=.55); ax[0].set(xlabel="Whole-genome baseline patristic distance",ylabel="DP>=10 CDS patristic distance",title=f"Raw distances\nSpearman r={spr.statistic:.3f}")
lim=max(db.max(),dcs.max()); ax[1].scatter(db,dcs,s=14,alpha=.55); ax[1].plot([0,lim],[0,lim],color="black",lw=1)
ax[1].set(xlabel="Whole-genome baseline patristic distance",ylabel="Scaled DP>=10 CDS patristic distance",title=f"Scaled distances\nPearson r={prs.statistic:.3f}; NRMSE={nrmse:.3f}")
fig.tight_layout(); fig.savefig(scatter_png,dpi=300,bbox_inches="tight"); plt.close(fig)

# Multi-metric summary
fig,ax=plt.subplots(2,2,figsize=(11,8))
items=[("Patristic distance ordering",float(spr.statistic),"Spearman r",(0,1.05)),
       ("Scaled branch-length agreement",1-nrmse,"1 - scaled NRMSE",(0,1.05)),
       ("Quartet topology agreement",qcon,"Quartet concordance",(0,1.05)),
       (">=80% supported split recovery",threshold_results[0.8]["precision"],"Precision",(0,1.05))]
for a,(title,val,label,ylim) in zip(ax.flat,items):
    a.bar([0],[val],width=.55); a.set_xticks([0],[label]); a.set_ylim(*ylim); a.set_title(title); a.text(0,val+0.025,f"{val:.3f}",ha="center")
fig.suptitle("Whole-genome baseline versus DP>=10 CDS phylogeny",fontweight="bold"); fig.tight_layout(); fig.savefig(summary_png,dpi=300,bbox_inches="tight"); plt.close(fig)

# Side-by-side tree plots, midpoint rooted only for display
bplot=copy.deepcopy(base); cplot=copy.deepcopy(cds); bplot.root_at_midpoint(); cplot.root_at_midpoint(); bplot.ladderize(); cplot.ladderize()
fig,axes=plt.subplots(1,2,figsize=(20,max(8,.3*len(labels))))
Phylo.draw(bplot,axes=axes[0],do_show=False,show_confidence=False,label_func=lambda c:c.name if c.is_terminal() else None); axes[0].set_title("Whole-genome baseline")
Phylo.draw(cplot,axes=axes[1],do_show=False,show_confidence=False,label_func=lambda c:c.name if c.is_terminal() else None); axes[1].set_title("DP>=10 complete-CDS tree")
fig.tight_layout(); fig.savefig(tree_png,dpi=300,bbox_inches="tight"); plt.close(fig)
print(f"Shared taxa: {len(labels)}")
print(f"Spearman r: {spr.statistic:.4f}")
print(f"Quartet concordance: {qcon:.4f}")
print(f">=80% supported-split precision/recall: {threshold_results[0.8]['precision']:.4f}/{threshold_results[0.8]['recall']:.4f}")
PY

echo "Completed whole-genome versus DP>=10 CDS comparison."
echo "Baseline tree: ${BASE_TREE}"
echo "Metrics: ${METRICS_TSV}"
echo "Support-filtered splits: ${SPLITS_TSV}"
echo "Figures: ${OUTDIR}/figures"
