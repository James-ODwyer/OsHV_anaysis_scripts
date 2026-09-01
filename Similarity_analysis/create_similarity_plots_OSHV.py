#!/usr/bin/env python3
"""Create sliding-window nucleotide-identity profiles from an aligned FASTA.

The script compares one or more focal genomes with selected reference genomes.
For each reference and alignment window, it calculates the mean pairwise identity
across all focal genomes. Pairwise identity is calculated only from positions at
which both sequences contain resolved nucleotides (A, C, G, or T). Alignment
positions containing gaps, Ns, or other ambiguous symbols in either sequence are
excluded from the denominator.

Outputs
-------
1. A long-format TSV containing per-focal and mean window statistics.
2. A PNG, PDF, and SVG similarity-profile figure.
3. A summary TSV with whole-alignment identity statistics.

Coordinates are 1-based inclusive alignment coordinates. If the input alignment
contains gaps, alignment coordinates will not necessarily equal coordinates in
any individual unaligned genome.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from collections import OrderedDict
from pathlib import Path
from typing import Iterable, Sequence

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

RESOLVED = frozenset("ACGT")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate sliding-window nucleotide-identity profiles comparing "
            "focal OsHV-1 genomes with selected reference genomes."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--alignment", required=True, type=Path,
        help="Aligned multi-FASTA containing focal and reference genomes.",
    )
    parser.add_argument(
        "--focal-samples", required=True, type=Path,
        help="Text file containing one focal sequence identifier per line.",
    )
    parser.add_argument(
        "--references", required=True, type=Path,
        help=(
            "Reference specification file. Use one identifier per line, or "
            "two tab-separated columns: sequence_identifier and display_label."
        ),
    )
    parser.add_argument(
        "--outdir", required=True, type=Path,
        help="Directory for tables and figures.",
    )
    parser.add_argument(
        "--window-size", type=int, default=1000,
        help="Sliding-window width in alignment columns.",
    )
    parser.add_argument(
        "--step-size", type=int, default=100,
        help="Distance between consecutive window starts in alignment columns.",
    )
    parser.add_argument(
        "--min-comparable-sites", type=int, default=100,
        help=(
            "Minimum resolved pairwise sites required for a focal-reference "
            "comparison within a window. Comparisons below this are recorded as NA."
        ),
    )
    parser.add_argument(
        "--min-focal-comparisons", type=int, default=1,
        help="Minimum valid focal-reference comparisons required for a window mean.",
    )
    parser.add_argument(
        "--identity-scale", choices=("percent", "fraction"), default="percent",
        help="Report identity as 0-100 percent or 0-1 fraction.",
    )
    parser.add_argument(
        "--smooth-windows", type=int, default=1,
        help=(
            "Centred rolling mean width applied only to plotted mean profiles. "
            "A value of 1 disables smoothing."
        ),
    )
    parser.add_argument(
        "--low-identity-threshold", type=float, default=None,
        help="Optional horizontal reference line on the selected identity scale.",
    )
    parser.add_argument(
        "--title", default="OsHV-1 whole-genome similarity profiles",
        help="Figure title.",
    )
    parser.add_argument(
        "--dpi", type=int, default=600,
        help="Resolution of the PNG output.",
    )
    parser.add_argument(
        "--figure-width", type=float, default=12.0,
        help="Figure width in inches.",
    )
    parser.add_argument(
        "--panel-height", type=float, default=2.4,
        help="Height in inches allocated to each reference panel.",
    )
    parser.add_argument(
        "--show-individual", action="store_true",
        help="Draw faint profiles for individual focal genomes behind the mean.",
    )
    parser.add_argument(
        "--overwrite", action="store_true",
        help="Permit replacement of existing output files.",
    )
    return parser.parse_args()


def read_fasta(path: Path) -> OrderedDict[str, str]:
    records: OrderedDict[str, str] = OrderedDict()
    name: str | None = None
    chunks: list[str] = []

    def store_record() -> None:
        nonlocal name, chunks
        if name is None:
            return
        if name in records:
            raise ValueError(f"Duplicate FASTA identifier: {name}")
        records[name] = "".join(chunks).upper()

    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                store_record()
                name = line[1:].split()[0]
                if not name:
                    raise ValueError("Encountered an empty FASTA identifier")
                chunks = []
            else:
                if name is None:
                    raise ValueError("Sequence data occurred before the first FASTA header")
                chunks.append(re.sub(r"\s+", "", line))
    store_record()

    if not records:
        raise ValueError(f"No FASTA records found in {path}")
    lengths = {len(sequence) for sequence in records.values()}
    if len(lengths) != 1:
        detail = ", ".join(f"{name}:{len(seq)}" for name, seq in records.items())
        raise ValueError(f"Input is not an alignment; sequence lengths differ: {detail}")
    return records


def read_name_list(path: Path) -> list[str]:
    names: list[str] = []
    with path.open(encoding="utf-8-sig") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            names.append(line.split("\t", 1)[0].strip())
    if not names:
        raise ValueError(f"No identifiers found in {path}")
    if len(names) != len(set(names)):
        raise ValueError(f"Duplicate identifiers found in {path}")
    return names


def read_references(path: Path) -> list[tuple[str, str]]:
    references: list[tuple[str, str]] = []
    with path.open(encoding="utf-8-sig") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            fields = [field.strip() for field in line.split("\t")]
            identifier = fields[0]
            label = fields[1] if len(fields) > 1 and fields[1] else identifier
            references.append((identifier, label))
    if not references:
        raise ValueError(f"No reference identifiers found in {path}")
    identifiers = [identifier for identifier, _ in references]
    if len(identifiers) != len(set(identifiers)):
        raise ValueError(f"Duplicate reference identifiers found in {path}")
    return references


def validate_names(
    records: OrderedDict[str, str],
    focal_names: Sequence[str],
    references: Sequence[tuple[str, str]],
) -> None:
    available = set(records)
    missing_focal = sorted(set(focal_names) - available)
    missing_reference = sorted({identifier for identifier, _ in references} - available)
    if missing_focal or missing_reference:
        messages: list[str] = []
        if missing_focal:
            messages.append("missing focal identifiers: " + ", ".join(missing_focal))
        if missing_reference:
            messages.append("missing reference identifiers: " + ", ".join(missing_reference))
        raise ValueError("; ".join(messages))


def identity_stats(left: str, right: str) -> tuple[float, int, int]:
    comparable = 0
    matches = 0
    for a, b in zip(left, right):
        if a in RESOLVED and b in RESOLVED:
            comparable += 1
            if a == b:
                matches += 1
    identity = matches / comparable if comparable else math.nan
    return identity, comparable, matches


def rolling_nanmean(values: np.ndarray, width: int) -> np.ndarray:
    if width <= 1:
        return values.copy()
    result = np.full(values.shape, np.nan, dtype=float)
    half = width // 2
    for index in range(len(values)):
        start = max(0, index - half)
        end = min(len(values), index + half + 1)
        section = values[start:end]
        if np.isfinite(section).any():
            result[index] = np.nanmean(section)
    return result


def scale_identity(value: float, identity_scale: str) -> float:
    if not math.isfinite(value):
        return math.nan
    return value * 100.0 if identity_scale == "percent" else value


def make_windows(length: int, width: int, step: int) -> Iterable[tuple[int, int]]:
    if width > length:
        yield 1, length
        return
    starts = list(range(1, length - width + 2, step))
    final_start = length - width + 1
    if starts[-1] != final_start:
        starts.append(final_start)
    for start in starts:
        yield start, min(length, start + width - 1)


def ensure_outputs_available(paths: Sequence[Path], overwrite: bool) -> None:
    existing = [str(path) for path in paths if path.exists()]
    if existing and not overwrite:
        raise FileExistsError(
            "Output files already exist. Use --overwrite to replace them: "
            + ", ".join(existing)
        )


def write_tables_and_collect_profiles(
    records: OrderedDict[str, str],
    focal_names: Sequence[str],
    references: Sequence[tuple[str, str]],
    args: argparse.Namespace,
    window_table: Path,
    summary_table: Path,
) -> dict[str, dict[str, np.ndarray]]:
    alignment_length = len(next(iter(records.values())))
    windows = list(make_windows(alignment_length, args.window_size, args.step_size))
    profiles: dict[str, dict[str, np.ndarray]] = {}

    with window_table.open("w", newline="", encoding="utf-8") as output:
        fieldnames = [
            "reference_id", "reference_label", "window_start", "window_end",
            "window_midpoint", "window_length", "focal_id", "identity",
            "comparable_sites", "matching_sites", "valid_comparison",
            "mean_identity", "median_identity", "minimum_identity",
            "maximum_identity", "n_valid_focal_comparisons",
        ]
        writer = csv.DictWriter(output, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()

        for reference_id, reference_label in references:
            mean_values: list[float] = []
            midpoint_values: list[float] = []
            individual: dict[str, list[float]] = {name: [] for name in focal_names}
            reference_sequence = records[reference_id]

            for start, end in windows:
                midpoint = (start + end) / 2.0
                midpoint_values.append(midpoint)
                values: list[float] = []
                per_focal_rows: list[dict[str, object]] = []
                ref_window = reference_sequence[start - 1:end]

                for focal_id in focal_names:
                    focal_window = records[focal_id][start - 1:end]
                    identity, comparable, matches = identity_stats(focal_window, ref_window)
                    valid = comparable >= args.min_comparable_sites
                    reported_identity = scale_identity(identity, args.identity_scale) if valid else math.nan
                    individual[focal_id].append(reported_identity)
                    if valid:
                        values.append(reported_identity)
                    per_focal_rows.append({
                        "reference_id": reference_id,
                        "reference_label": reference_label,
                        "window_start": start,
                        "window_end": end,
                        "window_midpoint": midpoint,
                        "window_length": end - start + 1,
                        "focal_id": focal_id,
                        "identity": "" if not valid else f"{reported_identity:.8f}",
                        "comparable_sites": comparable,
                        "matching_sites": matches,
                        "valid_comparison": int(valid),
                    })

                enough = len(values) >= args.min_focal_comparisons
                mean_identity = float(np.mean(values)) if enough else math.nan
                median_identity = float(np.median(values)) if enough else math.nan
                minimum_identity = float(np.min(values)) if enough else math.nan
                maximum_identity = float(np.max(values)) if enough else math.nan
                mean_values.append(mean_identity)

                for row in per_focal_rows:
                    row.update({
                        "mean_identity": "" if not enough else f"{mean_identity:.8f}",
                        "median_identity": "" if not enough else f"{median_identity:.8f}",
                        "minimum_identity": "" if not enough else f"{minimum_identity:.8f}",
                        "maximum_identity": "" if not enough else f"{maximum_identity:.8f}",
                        "n_valid_focal_comparisons": len(values),
                    })
                    writer.writerow(row)

            profiles[reference_id] = {
                "label": np.asarray([reference_label], dtype=object),
                "midpoints": np.asarray(midpoint_values, dtype=float),
                "mean": np.asarray(mean_values, dtype=float),
                "individual_names": np.asarray(list(focal_names), dtype=object),
                "individual": np.asarray([individual[name] for name in focal_names], dtype=float),
            }

    with summary_table.open("w", newline="", encoding="utf-8") as output:
        fieldnames = [
            "reference_id", "reference_label", "focal_id", "identity",
            "comparable_sites", "matching_sites", "alignment_length",
        ]
        writer = csv.DictWriter(output, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for reference_id, reference_label in references:
            aggregate_values: list[float] = []
            for focal_id in focal_names:
                identity, comparable, matches = identity_stats(records[focal_id], records[reference_id])
                value = scale_identity(identity, args.identity_scale)
                aggregate_values.append(value)
                writer.writerow({
                    "reference_id": reference_id,
                    "reference_label": reference_label,
                    "focal_id": focal_id,
                    "identity": f"{value:.8f}",
                    "comparable_sites": comparable,
                    "matching_sites": matches,
                    "alignment_length": alignment_length,
                })
            writer.writerow({
                "reference_id": reference_id,
                "reference_label": reference_label,
                "focal_id": "MEAN_ACROSS_FOCAL_GENOMES",
                "identity": f"{np.mean(aggregate_values):.8f}",
                "comparable_sites": "",
                "matching_sites": "",
                "alignment_length": alignment_length,
            })
    return profiles


def plot_profiles(
    profiles: dict[str, dict[str, np.ndarray]],
    references: Sequence[tuple[str, str]],
    args: argparse.Namespace,
    output_stem: Path,
) -> None:
    n_panels = len(references)
    figure_height = max(3.5, args.panel_height * n_panels)
    fig, axes = plt.subplots(
        n_panels, 1,
        figsize=(args.figure_width, figure_height),
        sharex=True,
        squeeze=False,
    )
    axes_flat = axes[:, 0]
    cmap = plt.get_cmap("tab20")

    for panel_index, ((reference_id, reference_label), axis) in enumerate(zip(references, axes_flat)):
        profile = profiles[reference_id]
        x = profile["midpoints"]
        mean = profile["mean"]
        plotted_mean = rolling_nanmean(mean, args.smooth_windows)

        if args.show_individual:
            for focal_index, values in enumerate(profile["individual"]):
                axis.plot(
                    x, values,
                    color=cmap(focal_index % 20),
                    linewidth=0.45,
                    alpha=0.22,
                    zorder=1,
                )

        axis.plot(x, plotted_mean, color="#1769AA", linewidth=1.35, zorder=3)
        axis.fill_between(x, plotted_mean, alpha=0.12, color="#1769AA", zorder=2)
        if args.low_identity_threshold is not None:
            axis.axhline(
                args.low_identity_threshold,
                color="#C62828",
                linestyle="--",
                linewidth=0.8,
                zorder=0,
            )
        axis.set_title(f"{reference_label} ({reference_id})", loc="left", fontsize=10, fontweight="bold")
        axis.set_ylabel("Identity (%)" if args.identity_scale == "percent" else "Identity")
        axis.grid(axis="y", linewidth=0.35, alpha=0.35)
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
        if args.identity_scale == "fraction":
            axis.set_ylim(0, 1.01)
        else:
            finite = plotted_mean[np.isfinite(plotted_mean)]
            if finite.size:
                lower = max(0.0, math.floor((float(np.min(finite)) - 2.0) / 5.0) * 5.0)
                axis.set_ylim(lower, 100.5)

    axes_flat[-1].set_xlabel("Alignment position (bp)")
    subtitle = (
        f"Window = {args.window_size:,} bp; step = {args.step_size:,} bp; "
        f"minimum comparable sites = {args.min_comparable_sites:,}"
    )
    fig.suptitle(args.title, fontsize=14, fontweight="bold", y=0.995)
    fig.text(0.5, 0.975, subtitle, ha="center", va="top", fontsize=9)
    fig.tight_layout(rect=(0.04, 0.035, 0.99, 0.945))
    fig.savefig(output_stem.with_suffix(".png"), dpi=args.dpi, bbox_inches="tight", facecolor="white")
    fig.savefig(output_stem.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    fig.savefig(output_stem.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main() -> int:
    args = parse_args()
    if args.window_size < 1 or args.step_size < 1:
        raise ValueError("--window-size and --step-size must be positive integers")
    if args.min_comparable_sites < 1:
        raise ValueError("--min-comparable-sites must be a positive integer")
    if args.min_comparable_sites > args.window_size:
        raise ValueError("--min-comparable-sites cannot exceed --window-size")
    if args.min_focal_comparisons < 1:
        raise ValueError("--min-focal-comparisons must be a positive integer")
    if args.smooth_windows < 1:
        raise ValueError("--smooth-windows must be at least 1")

    args.outdir.mkdir(parents=True, exist_ok=True)
    window_table = args.outdir / "OsHV_similarity_windows.tsv"
    summary_table = args.outdir / "OsHV_whole_alignment_identity.tsv"
    figure_stem = args.outdir / "OsHV_similarity_profiles"
    outputs = [
        window_table,
        summary_table,
        figure_stem.with_suffix(".png"),
        figure_stem.with_suffix(".pdf"),
        figure_stem.with_suffix(".svg"),
    ]
    ensure_outputs_available(outputs, args.overwrite)

    records = read_fasta(args.alignment)
    focal_names = read_name_list(args.focal_samples)
    references = read_references(args.references)
    validate_names(records, focal_names, references)
    if args.min_focal_comparisons > len(focal_names):
        raise ValueError("--min-focal-comparisons exceeds the number of focal genomes")

    profiles = write_tables_and_collect_profiles(
        records,
        focal_names,
        references,
        args,
        window_table,
        summary_table,
    )
    plot_profiles(profiles, references, args, figure_stem)

    print(f"Alignment records: {len(records)}")
    print(f"Alignment length: {len(next(iter(records.values()))):,} columns")
    print(f"Focal genomes: {len(focal_names)}")
    print(f"Reference genomes: {len(references)}")
    print(f"Window table: {window_table}")
    print(f"Whole-alignment summary: {summary_table}")
    print(f"Figures: {figure_stem}.png/.pdf/.svg")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValueError, FileNotFoundError, FileExistsError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(2)
