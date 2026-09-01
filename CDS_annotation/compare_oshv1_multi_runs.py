#!/usr/bin/env python3

import os
import re
import glob
import csv
import argparse
from itertools import combinations

ORF_REGEX = re.compile(r'\bORF[_\-\s]?[A-Za-z0-9.]+\b', re.IGNORECASE)

def normalize_gene_name(name):
    if name is None:
        return None
    name = str(name).strip()
    if not name or name.upper() == "NA":
        return None

    m = ORF_REGEX.search(name)
    if m:
        x = m.group(0).upper()
        x = x.replace("ORF_", "ORF").replace("ORF-", "ORF").replace(" ", "")
        return x

    # fallback for non-ORF names
    name = re.sub(r'\s+', '_', name)
    return name

def parse_attrs(attr_str):
    attrs = {}
    for item in attr_str.strip().split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            attrs[k] = v
    return attrs

def extract_gene_from_attrs(attrs):
    for key in ["gene", "Name", "locus_tag", "product", "note", "ID"]:
        if key in attrs:
            val = normalize_gene_name(attrs[key])
            if val:
                return val
    return None

def load_gff_gene_set(gff_file):
    genes = set()
    with open(gff_file) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue
            seqid, source, feature_type, start, end, score, strand, phase, attrs_str = parts
            if feature_type not in {"gene", "mRNA", "CDS"}:
                continue
            attrs = parse_attrs(attrs_str)
            g = extract_gene_from_attrs(attrs)
            if g:
                genes.add(g)
    return genes

def find_final_gffs(run_dir):
    pattern = os.path.join(run_dir, "05_final", "*", "*.final.gff3")
    gffs = sorted(glob.glob(pattern))
    sample_to_file = {}
    for gff in gffs:
        sample = os.path.basename(gff).replace(".final.gff3", "")
        sample_to_file[sample] = gff
    return sample_to_file

def ensure_dir(path):
    os.makedirs(path, exist_ok=True)

def resolve_primary_reference_gff(run_dir, label):
    """
    Try to resolve the primary reference GFF from:
    run_dir/01_reference/ref_gffs/

    Handles:
      label.gff3
      label*.gff3
    """
    ref_dir = os.path.join(run_dir, "01_reference", "ref_gffs")
    if not os.path.isdir(ref_dir):
        raise FileNotFoundError(f"Missing reference GFF directory: {ref_dir}")

    exact = os.path.join(ref_dir, f"{label}.gff3")
    if os.path.isfile(exact):
        return exact

    matches = sorted(glob.glob(os.path.join(ref_dir, f"{label}*.gff3")))
    if len(matches) == 1:
        return matches[0]
    elif len(matches) > 1:
        raise RuntimeError(
            f"Multiple reference GFF files match label '{label}' in {ref_dir}: {matches}\n"
            f"Please use a more specific label or rename files."
        )

    raise FileNotFoundError(
        f"Could not find primary reference GFF for label '{label}' in {ref_dir}.\n"
        f"Looked for {label}.gff3 and {label}*.gff3"
    )

def write_per_sample_missing(out_tsv, samples, reference_genes, sample_gene_sets):
    with open(out_tsv, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["sample", "missing_count", "missing_genes"])
        for sample in sorted(samples):
            present = sample_gene_sets.get(sample, set())
            missing = sorted(reference_genes - present)
            writer.writerow([sample, len(missing), ",".join(missing)])

def write_missing_frequency(out_tsv, samples, reference_genes, sample_gene_sets):
    with open(out_tsv, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["gene", "missing_in_n_samples", "samples_missing"])
        for gene in sorted(reference_genes):
            missing_samples = [
                sample for sample in sorted(samples)
                if gene not in sample_gene_sets.get(sample, set())
            ]
            writer.writerow([gene, len(missing_samples), ",".join(missing_samples)])

def write_reference_presence_matrix(out_tsv, samples, reference_genes, sample_gene_sets):
    genes_sorted = sorted(reference_genes)
    with open(out_tsv, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["sample"] + genes_sorted)
        for sample in sorted(samples):
            present = sample_gene_sets.get(sample, set())
            row = [sample] + [1 if gene in present else 0 for gene in genes_sorted]
            writer.writerow(row)

def write_pairwise_difference_lists(out_tsv, samples, a_gene_sets, b_gene_sets, label_a, label_b):
    with open(out_tsv, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow([
            "sample",
            f"genes_only_in_{label_a}_count",
            f"genes_only_in_{label_a}",
            f"genes_only_in_{label_b}_count",
            f"genes_only_in_{label_b}"
        ])
        for sample in sorted(samples):
            genes_a = a_gene_sets.get(sample, set())
            genes_b = b_gene_sets.get(sample, set())
            only_a = sorted(genes_a - genes_b)
            only_b = sorted(genes_b - genes_a)
            writer.writerow([
                sample,
                len(only_a),
                ",".join(only_a),
                len(only_b),
                ",".join(only_b)
            ])

def write_pairwise_difference_frequency(out_tsv, samples, a_gene_sets, b_gene_sets, label_a, label_b):
    all_genes = set()
    for sample in samples:
        all_genes.update(a_gene_sets.get(sample, set()))
        all_genes.update(b_gene_sets.get(sample, set()))

    with open(out_tsv, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow([
            "gene",
            f"present_only_in_{label_a}_n_samples",
            f"samples_only_in_{label_a}",
            f"present_only_in_{label_b}_n_samples",
            f"samples_only_in_{label_b}",
        ])

        for gene in sorted(all_genes):
            only_a_samples = []
            only_b_samples = []
            for sample in sorted(samples):
                in_a = gene in a_gene_sets.get(sample, set())
                in_b = gene in b_gene_sets.get(sample, set())
                if in_a and not in_b:
                    only_a_samples.append(sample)
                elif in_b and not in_a:
                    only_b_samples.append(sample)

            writer.writerow([
                gene,
                len(only_a_samples),
                ",".join(only_a_samples),
                len(only_b_samples),
                ",".join(only_b_samples),
            ])

def write_multirun_gene_support_matrix(out_tsv, samples, run_labels, run_sample_gene_sets):
    all_genes = set()
    for label in run_labels:
        for sample in samples:
            all_genes.update(run_sample_gene_sets[label].get(sample, set()))

    with open(out_tsv, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["sample", "gene"] + run_labels + ["pattern", "n_runs_present"])

        for sample in sorted(samples):
            for gene in sorted(all_genes):
                statuses = []
                for label in run_labels:
                    statuses.append(1 if gene in run_sample_gene_sets[label].get(sample, set()) else 0)
                pattern = "".join(str(x) for x in statuses)
                writer.writerow([sample, gene] + statuses + [pattern, sum(statuses)])

def write_multirun_gene_pattern_frequency(out_tsv, samples, run_labels, run_sample_gene_sets):
    all_genes = set()
    for label in run_labels:
        for sample in samples:
            all_genes.update(run_sample_gene_sets[label].get(sample, set()))

    with open(out_tsv, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["gene"] + run_labels + ["pattern", "n_samples_with_pattern", "samples"])

        for gene in sorted(all_genes):
            pattern_to_samples = {}
            for sample in sorted(samples):
                statuses = []
                for label in run_labels:
                    statuses.append(1 if gene in run_sample_gene_sets[label].get(sample, set()) else 0)
                pattern = "".join(str(x) for x in statuses)
                pattern_to_samples.setdefault(pattern, []).append(sample)

            for pattern, pl_samples in sorted(pattern_to_samples.items()):
                bits = [int(x) for x in pattern]
                writer.writerow([gene] + bits + [pattern, len(pl_samples), ",".join(pl_samples)])

def write_summary_counts(out_tsv, samples, run_labels, run_sample_gene_sets, run_reference_genes):
    header = ["sample"]
    for label in run_labels:
        header.append(f"{label}_gene_count")
    for label in run_labels:
        header.append(f"missing_from_{label}_reference_count")

    with open(out_tsv, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(header)

        for sample in sorted(samples):
            row = [sample]
            for label in run_labels:
                row.append(len(run_sample_gene_sets[label].get(sample, set())))
            for label in run_labels:
                missing = run_reference_genes[label] - run_sample_gene_sets[label].get(sample, set())
                row.append(len(missing))
            writer.writerow(row)

def parse_run_arg(run_arg):
    """
    Expect LABEL=PATH
    """
    if "=" not in run_arg:
        raise ValueError(f"Run argument must be LABEL=PATH, got: {run_arg}")
    label, path = run_arg.split("=", 1)
    label = label.strip()
    path = path.strip()
    if not label or not path:
        raise ValueError(f"Invalid run argument: {run_arg}")
    return label, path

def main():
    parser = argparse.ArgumentParser(
        description="Compare OsHV-1 annotations across multiple primary-reference runs."
    )
    parser.add_argument(
        "--run",
        action="append",
        required=True,
        help="Run specification in the form LABEL=PATH. Use multiple times."
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory"
    )

    args = parser.parse_args()

    runs = {}
    for run_arg in args.run:
        label, path = parse_run_arg(run_arg)
        if label in runs:
            raise ValueError(f"Duplicate label provided: {label}")
        runs[label] = path

    if len(runs) < 2:
        raise ValueError("Please provide at least two --run arguments.")
    ensure_dir(args.outdir)

    run_labels = list(runs.keys())

    # Load primary reference genes and sample gene sets
    run_reference_genes = {}
    run_sample_gene_sets = {}
    run_final_gffs = {}

    all_samples = set()

    for label, run_dir in runs.items():
        if not os.path.isdir(run_dir):
            raise FileNotFoundError(f"Run directory not found: {run_dir}")

        primary_ref_gff = resolve_primary_reference_gff(run_dir, label)
        reference_genes = load_gff_gene_set(primary_ref_gff)
        run_reference_genes[label] = reference_genes

        final_gffs = find_final_gffs(run_dir)
        run_final_gffs[label] = final_gffs

        sample_gene_sets = {}
        for sample, gff_file in final_gffs.items():
            sample_gene_sets[sample] = load_gff_gene_set(gff_file)

        run_sample_gene_sets[label] = sample_gene_sets
        all_samples.update(final_gffs.keys())

    all_samples = sorted(all_samples)

    # ------------------------------------------------------------------
    # 1) For each run: genes present in that primary reference but absent
    # ------------------------------------------------------------------
    per_reference_dir = os.path.join(args.outdir, "per_reference_missing")
    ensure_dir(per_reference_dir)

    for label in run_labels:
        reference_genes = run_reference_genes[label]
        sample_gene_sets = run_sample_gene_sets[label]

        write_per_sample_missing(
            os.path.join(per_reference_dir, f"missing_from_{label}_reference_by_sample.tsv"),
            all_samples, reference_genes, sample_gene_sets
        )

        write_missing_frequency(
            os.path.join(per_reference_dir, f"missing_from_{label}_reference_gene_frequency.tsv"),
            all_samples, reference_genes, sample_gene_sets
        )

        write_reference_presence_matrix(
            os.path.join(per_reference_dir, f"presence_matrix_vs_{label}_reference.tsv"),
            all_samples, reference_genes, sample_gene_sets
        )

    # ------------------------------------------------------------------
    # 2) Pairwise comparison of runs
    # ------------------------------------------------------------------
    pairwise_dir = os.path.join(args.outdir, "pairwise_run_differences")
    ensure_dir(pairwise_dir)

    for label_a, label_b in combinations(run_labels, 2):
        gene_sets_a = run_sample_gene_sets[label_a]
        gene_sets_b = run_sample_gene_sets[label_b]

        write_pairwise_difference_lists(
            os.path.join(pairwise_dir, f"{label_a}_vs_{label_b}_per_sample_gene_differences.tsv"),
            all_samples, gene_sets_a, gene_sets_b, label_a, label_b
        )

        write_pairwise_difference_frequency(
            os.path.join(pairwise_dir, f"{label_a}_vs_{label_b}_gene_difference_frequency.tsv"),
            all_samples, gene_sets_a, gene_sets_b, label_a, label_b
        )

    # ------------------------------------------------------------------
    # 3) Multi-run support matrices
    # ------------------------------------------------------------------
    multirun_dir = os.path.join(args.outdir, "multirun_support")
    ensure_dir(multirun_dir)

    write_multirun_gene_support_matrix(
        os.path.join(multirun_dir, "sample_gene_multirun_support_matrix.tsv"),
        all_samples, run_labels, run_sample_gene_sets
    )

    write_multirun_gene_pattern_frequency(
        os.path.join(multirun_dir, "gene_multirun_pattern_frequency.tsv"),
        all_samples, run_labels, run_sample_gene_sets
    )

    # ------------------------------------------------------------------
    # 4) Summary counts
    # ------------------------------------------------------------------
    write_summary_counts(
        os.path.join(args.outdir, "summary_counts.tsv"),
        all_samples, run_labels, run_sample_gene_sets, run_reference_genes
    )

    print("Done.")
    print(f"Output directory: {args.outdir}")
    print("")
    print("Created:")
    print(f"  {per_reference_dir}/")
    print(f"  {pairwise_dir}/")
    print(f"  {multirun_dir}/")
    print(f"  {os.path.join(args.outdir, 'summary_counts.tsv')}")

if __name__ == "__main__":
    main()