#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import math
import re

INPUT_DIR = Path(
    "/group/sequencing/assembly/James/OsHV_annot/Align_raws/bbmapaligned_to_final_genomes_for_SRA_upload_1to1_2"
)

OUTDIR = INPUT_DIR / "coverage_plots"
OUTDIR.mkdir(exist_ok=True)

files = sorted(INPUT_DIR.glob("*_perbase_sta*"))

if len(files) == 0:
    raise RuntimeError("No perbase_stats files found")

coverage_data = []

# ---------------------------------------------------------------------
# Individual coverage plots
# ---------------------------------------------------------------------

for file in files:

    sample = re.sub(
        r"_perbase_stats.*$",
        "",
        file.name
    )

    df = pd.read_csv(
        file,
        sep="\t",
        comment="#",
        header=None,
        names=["RefName", "Pos", "Coverage"]
    )

    coverage_data.append((sample, df))

    plt.figure(figsize=(10, 3))

    plt.plot(
        df["Pos"],
        df["Coverage"],
        linewidth=0.35,
        color="steelblue"
    )

    # 5x assembly threshold
    plt.axhline(
        y=5,
        color="red",
        linestyle="--",
        linewidth=1,
        label="5x coverage threshold"
    )

    plt.title(sample)
    plt.xlabel("Genome position (bp)")
    plt.ylabel("Read depth")

    plt.legend(
        loc="upper right",
        fontsize=8,
        frameon=False
    )

    plt.tight_layout()

    plt.savefig(
        OUTDIR / f"{sample}_coverage.png",
        dpi=300
    )

    plt.close()

# ---------------------------------------------------------------------
# Combined multi-panel figure
# ---------------------------------------------------------------------

n = len(coverage_data)

ncols = 3
nrows = math.ceil(n / ncols)

fig, axes = plt.subplots(
    nrows,
    ncols,
    figsize=(18, 4 * nrows),
    sharex=False,
    sharey=False
)

axes = axes.flatten()

for ax, (sample, df) in zip(axes, coverage_data):

    ax.plot(
        df["Pos"],
        df["Coverage"],
        linewidth=0.30,
        color="darkgreen"
    )

    # 5x assembly threshold
    ax.axhline(
        y=5,
        color="red",
        linestyle="--",
        linewidth=0.8
    )

    ax.set_title(
        sample,
        fontsize=10,
        fontweight="bold"
    )

    ax.tick_params(
        axis="both",
        labelsize=7
    )

    ax.set_xlabel("Position")
    ax.set_ylabel("Depth")

# Remove unused panels
for ax in axes[len(coverage_data):]:
    fig.delaxes(ax)

# Figure-level legend
from matplotlib.lines import Line2D

legend_elements = [
    Line2D(
        [0],
        [0],
        color="red",
        linestyle="--",
        linewidth=1,
        label="5x coverage threshold"
    )
]

fig.legend(
    handles=legend_elements,
    loc="upper center",
    ncol=1,
    frameon=False
)

plt.tight_layout(rect=[0, 0, 1, 0.97])

combined_file = OUTDIR / "OsHV_all_samples_coverage.png"

plt.savefig(
    combined_file,
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print(f"Written individual plots to: {OUTDIR}")
print(f"Combined plot: {combined_file}")