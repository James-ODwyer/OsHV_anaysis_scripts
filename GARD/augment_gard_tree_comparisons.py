#!/usr/bin/env python3
"""Augment a GARD partition_tree_summary.tsv with complementary tree metrics.

For each successfully compared partition this script:
  * uses the taxon-matched baseline and partition Newick trees;
  * re-checks and, if necessary, prunes both trees to their exact shared taxa;
  * compares all pairwise patristic distances using Spearman correlation;
  * scales partition distances to baseline distances by least squares through 0;
  * reports scaled Pearson r, R2, RMSE and normalized RMSE;
  * reports nearest-neighbour top-1 and top-3 overlap;
  * reports quartet concordance, discordance and unresolved proportions;
  * retains existing RF / support-filtered split statistics from the input TSV;
  * creates per-partition distance scatterplots and a multi-metric summary plot.

Quartet comparison is topology-only and evaluates every quartet by default.
For much larger trees, use --max-quartets to reproducibly subsample quartets.
"""
from __future__ import annotations
import argparse
import copy
import csv
import itertools
import math
import random
from pathlib import Path
from typing import Dict, List, Sequence, Set, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from Bio import Phylo
from scipy.stats import pearsonr, spearmanr


def read_tree(path: str):
    return Phylo.read(path, "newick")


def taxa(tree) -> Set[str]:
    return {tip.name for tip in tree.get_terminals()}


def prune_exact(tree, keep: Set[str]):
    out = copy.deepcopy(tree)
    for tip in list(out.get_terminals()):
        if tip.name not in keep:
            out.prune(tip)
    return out


def matched_trees(base_path: str, part_path: str):
    base, part = read_tree(base_path), read_tree(part_path)
    shared = taxa(base) & taxa(part)
    if len(shared) < 4:
        raise ValueError(f"Only {len(shared)} shared taxa")
    base, part = prune_exact(base, shared), prune_exact(part, shared)
    if taxa(base) != taxa(part):
        raise RuntimeError("Taxon matching failed")
    return base, part, sorted(shared)


def patristic_vector(tree, labels: Sequence[str]):
    pairs, values = [], []
    for i, a in enumerate(labels):
        for b in labels[i + 1:]:
            pairs.append((a, b))
            values.append(tree.distance(a, b))
    return pairs, np.asarray(values, dtype=float)


def safe_correlation(func, x, y):
    if len(x) < 3 or np.allclose(x, x[0]) or np.allclose(y, y[0]):
        return float("nan"), float("nan")
    result = func(x, y)
    return float(result.statistic), float(result.pvalue)


def pairwise_metrics(base, part, labels: Sequence[str]):
    pairs, db = patristic_vector(base, labels)
    pairs2, dp = patristic_vector(part, labels)
    if pairs != pairs2:
        raise RuntimeError("Pairwise taxon order differs")
    denom = float(np.dot(dp, dp))
    scale = float(np.dot(dp, db) / denom) if denom > 0 else float("nan")
    scaled = dp * scale if math.isfinite(scale) else np.full_like(dp, np.nan)
    spearman_r, spearman_p = safe_correlation(spearmanr, db, dp)
    pearson_r, pearson_p = safe_correlation(pearsonr, db, scaled)
    residual = db - scaled
    rmse = float(np.sqrt(np.mean(residual ** 2)))
    mean_base = float(np.mean(db))
    nrmse_mean = rmse / mean_base if mean_base > 0 else float("nan")
    sst = float(np.sum((db - mean_base) ** 2))
    r2 = 1.0 - float(np.sum(residual ** 2)) / sst if sst > 0 else float("nan")
    mae = float(np.mean(np.abs(residual)))
    return {
        "patristic_n_pairs": len(db),
        "patristic_scale_factor": scale,
        "patristic_spearman_r": spearman_r,
        "patristic_spearman_p": spearman_p,
        "patristic_scaled_pearson_r": pearson_r,
        "patristic_scaled_pearson_p": pearson_p,
        "patristic_scaled_r2": r2,
        "patristic_scaled_rmse": rmse,
        "patristic_scaled_nrmse_mean": nrmse_mean,
        "patristic_scaled_mae": mae,
    }, db, dp, scaled


def nearest_sets(tree, labels: Sequence[str], k: int) -> Dict[str, Set[str]]:
    output = {}
    for a in labels:
        distances = sorted((tree.distance(a, b), b) for b in labels if b != a)
        if not distances:
            output[a] = set()
            continue
        # Include ties at the kth distance so arbitrary ordering does not penalize polytomies.
        cutoff = distances[min(k, len(distances)) - 1][0]
        tolerance = max(1e-12, abs(cutoff) * 1e-9)
        output[a] = {b for d, b in distances if d <= cutoff + tolerance}
    return output


def nearest_metrics(base, part, labels: Sequence[str]):
    b1, p1 = nearest_sets(base, labels, 1), nearest_sets(part, labels, 1)
    b3, p3 = nearest_sets(base, labels, 3), nearest_sets(part, labels, 3)
    exact1 = np.mean([bool(b1[x] & p1[x]) for x in labels])
    overlap3 = np.mean([
        len(b3[x] & p3[x]) / len(b3[x] | p3[x]) if (b3[x] | p3[x]) else 1.0
        for x in labels
    ])
    baseline_nn_in_part_top3 = np.mean([bool(b1[x] & p3[x]) for x in labels])
    return {
        "nearest_top1_agreement": float(exact1),
        "nearest_top3_jaccard_mean": float(overlap3),
        "baseline_nearest_in_partition_top3": float(baseline_nn_in_part_top3),
    }


def quartet_state(tree, quartet: Tuple[str, str, str, str], tol=1e-12):
    a, b, c, d = quartet
    sums = {
        "ab|cd": tree.distance(a, b) + tree.distance(c, d),
        "ac|bd": tree.distance(a, c) + tree.distance(b, d),
        "ad|bc": tree.distance(a, d) + tree.distance(b, c),
    }
    ordered = sorted(sums.items(), key=lambda kv: kv[1])
    scale = max(1.0, max(abs(v) for v in sums.values()))
    if abs(ordered[1][1] - ordered[0][1]) <= tol * scale:
        return "unresolved"
    return ordered[0][0]


def quartet_metrics(base, part, labels: Sequence[str], max_quartets: int, seed: int):
    total_possible = math.comb(len(labels), 4)
    if max_quartets > 0 and total_possible > max_quartets:
        rng = random.Random(seed)
        sampled = set()
        while len(sampled) < max_quartets:
            sampled.add(tuple(sorted(rng.sample(list(labels), 4))))
        quartets = sorted(sampled)
        sampled_flag = 1
    else:
        quartets = itertools.combinations(labels, 4)
        sampled_flag = 0
    concordant = discordant = unresolved_either = comparable = 0
    n = 0
    for quartet in quartets:
        n += 1
        left, right = quartet_state(base, quartet), quartet_state(part, quartet)
        if left == "unresolved" or right == "unresolved":
            unresolved_either += 1
            continue
        comparable += 1
        if left == right:
            concordant += 1
        else:
            discordant += 1
    return {
        "quartets_possible": total_possible,
        "quartets_evaluated": n,
        "quartets_sampled": sampled_flag,
        "quartets_comparable": comparable,
        "quartet_concordance": concordant / comparable if comparable else float("nan"),
        "quartet_discordance": discordant / comparable if comparable else float("nan"),
        "quartet_unresolved_either_fraction": unresolved_either / n if n else float("nan"),
    }


def scatter_plot(db, dp, scaled, metrics, outpath: Path, title: str):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    axes[0].scatter(db, dp, s=13, alpha=0.55)
    axes[0].set_xlabel("Baseline patristic distance")
    axes[0].set_ylabel("Partition patristic distance")
    axes[0].set_title(f"Raw distances\nSpearman r={metrics['patristic_spearman_r']:.3f}")
    axes[1].scatter(db, scaled, s=13, alpha=0.55)
    limit = max(float(np.max(db)), float(np.max(scaled))) if len(db) else 1
    axes[1].plot([0, limit], [0, limit], color="black", linewidth=1)
    axes[1].set_xlabel("Baseline patristic distance")
    axes[1].set_ylabel("Scaled partition patristic distance")
    axes[1].set_title(
        f"Least-squares scaled\nPearson r={metrics['patristic_scaled_pearson_r']:.3f}; "
        f"NRMSE={metrics['patristic_scaled_nrmse_mean']:.3f}"
    )
    fig.suptitle(title, fontweight="bold")
    fig.tight_layout()
    fig.savefig(outpath, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def summary_plot(rows: List[dict], outpath: Path):
    usable = [r for r in rows if r.get("extended_metric_status") == "calculated"]
    if not usable:
        return
    x = np.arange(len(usable))
    labels = [f"P{r['partition']}" for r in usable]
    fig, axes = plt.subplots(2, 2, figsize=(13, 9), sharex=True)
    axes[0,0].plot(x, [float(r["patristic_spearman_r"]) for r in usable], marker="o")
    axes[0,0].set_ylabel("Spearman r")
    axes[0,0].set_ylim(-0.05, 1.05)
    axes[0,0].set_title("Pairwise patristic-distance ordering")
    axes[0,1].plot(x, [float(r["patristic_scaled_nrmse_mean"]) for r in usable], marker="o")
    axes[0,1].set_ylabel("Scaled NRMSE / mean baseline distance")
    axes[0,1].set_ylim(bottom=0)
    axes[0,1].set_title("Scaled branch-length discrepancy")
    axes[1,0].plot(x, [float(r["quartet_concordance"]) for r in usable], marker="o")
    axes[1,0].set_ylabel("Quartet concordance")
    axes[1,0].set_ylim(-0.05, 1.05)
    axes[1,0].set_title("Topology agreement across quartets")
    axes[1,1].plot(x, [float(r["support_0_8_precision"]) if r.get("support_0_8_precision") else np.nan for r in usable], marker="o", label="Supported-split precision >=0.8")
    axes[1,1].plot(x, [float(r["support_0_8_recall"]) if r.get("support_0_8_recall") else np.nan for r in usable], marker="o", label="Supported-split recall >=0.8")
    axes[1,1].set_ylabel("Proportion")
    axes[1,1].set_ylim(-0.05, 1.05)
    axes[1,1].set_title("Supported split recovery")
    axes[1,1].legend(frameon=False, fontsize=8)
    for ax in axes[1,:]:
        ax.set_xticks(x, labels, rotation=45)
        ax.set_xlabel("GARD partition")
    fig.tight_layout()
    fig.savefig(outpath, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary", required=True, type=Path, help="partition_tree_summary.tsv")
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--max-quartets", type=int, default=0, help="0 evaluates all quartets")
    ap.add_argument("--seed", type=int, default=20260702)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    figdir = args.outdir / "pairwise_distance_figures"
    figdir.mkdir(exist_ok=True)

    with args.summary.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        original_fields = list(reader.fieldnames or [])

    new_fields = [
        "extended_metric_status", "extended_metric_details", "metric_shared_taxa",
        "patristic_n_pairs", "patristic_scale_factor", "patristic_spearman_r",
        "patristic_spearman_p", "patristic_scaled_pearson_r", "patristic_scaled_pearson_p",
        "patristic_scaled_r2", "patristic_scaled_rmse", "patristic_scaled_nrmse_mean",
        "patristic_scaled_mae", "nearest_top1_agreement", "nearest_top3_jaccard_mean",
        "baseline_nearest_in_partition_top3", "quartets_possible", "quartets_evaluated",
        "quartets_sampled", "quartets_comparable", "quartet_concordance",
        "quartet_discordance", "quartet_unresolved_either_fraction", "patristic_figure"
    ]

    for row in rows:
        if not row.get("baseline_taxon_matched_tree") or not row.get("partition_taxon_matched_tree"):
            row["extended_metric_status"] = "skipped"
            row["extended_metric_details"] = "Matched tree paths unavailable"
            continue
        try:
            base, part, labels = matched_trees(row["baseline_taxon_matched_tree"], row["partition_taxon_matched_tree"])
            pm, db, dp, scaled = pairwise_metrics(base, part, labels)
            nm = nearest_metrics(base, part, labels)
            qm = quartet_metrics(base, part, labels, args.max_quartets, args.seed + int(row["partition"]))
            row.update({k: str(v) for k, v in {**pm, **nm, **qm}.items()})
            row["metric_shared_taxa"] = str(len(labels))
            row["extended_metric_status"] = "calculated"
            row["extended_metric_details"] = ""
            fig = figdir / f"partition_{int(row['partition']):02d}_patristic_comparison.png"
            scatter_plot(db, dp, scaled, pm, fig, f"GARD partition {row['partition']}: matched-tree distances")
            row["patristic_figure"] = str(fig)
        except Exception as error:
            row["extended_metric_status"] = "failed"
            row["extended_metric_details"] = str(error)

    out_tsv = args.outdir / "partition_tree_summary_extended.tsv"
    fields = original_fields + [f for f in new_fields if f not in original_fields]
    with out_tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    summary_plot(rows, args.outdir / "multi_metric_partition_summary.png")
    print(f"Wrote {out_tsv}")
    print(f"Wrote figures under {args.outdir}")


if __name__ == "__main__":
    main()
