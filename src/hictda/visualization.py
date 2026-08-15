"""Publication-oriented visualization functions."""

from __future__ import annotations

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from persim import plot_diagrams


def _save(fig, path):
    if path is not None:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=300, bbox_inches="tight")
    return fig


def plot_contact_map(
    matrix,
    chromosome,
    start,
    end,
    output=None,
):
    matrix = np.asarray(matrix)
    fig, ax = plt.subplots(figsize=(7, 6))

    log_matrix = np.log1p(matrix)
    cmap = LinearSegmentedColormap.from_list(
        "hictda_hic_red", ["white", "red", "darkred"]
    )
    im = ax.imshow(log_matrix, cmap=cmap, interpolation="nearest")

    ax.set_title(
        f"Hi-C Contact Map (Chr{chromosome}: "
        f"{start/1e6:.1f}–{end/1e6:.1f} Mb)"
    )

    ticks = np.linspace(0, matrix.shape[0] - 1, 5)
    labels = np.linspace(start / 1e6, end / 1e6, 5)
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{x:.1f} Mb" for x in labels])
    ax.set_yticks(ticks)
    ax.set_yticklabels([f"{x:.1f} Mb" for x in labels])
    ax.set_xlabel("Genomic position")
    ax.set_ylabel("Genomic position")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="log(1 + contact)")

    fig.tight_layout()
    return _save(fig, output)


def plot_persistence_diagram(diagrams, output=None):
    fig, ax = plt.subplots(figsize=(6, 6))
    plot_diagrams(diagrams, show=False, ax=ax, title="Persistence Diagram")
    ax.set_xlabel("Birth (distance)")
    ax.set_ylabel("Death (distance)")
    fig.tight_layout()
    return _save(fig, output)


def plot_barcode(diagrams, output=None):
    colors = ["#1f77b4", "#ff7f0e"]
    labels = ["$H_0$ (Connected Components)", "$H_1$ (Loops)"]

    max_finite = 0.0
    for dgm in diagrams:
        if len(dgm):
            finite_deaths = dgm[np.isfinite(dgm[:, 1]), 1]
            if finite_deaths.size:
                max_finite = max(max_finite, float(np.max(finite_deaths)))

    inf_cap = max_finite * 1.1 if max_finite > 0 else 1.0

    fig, ax = plt.subplots(figsize=(9, 6))
    y_pos = 0

    for dim, dgm in enumerate(diagrams):
        if len(dgm) == 0:
            continue

        lifespans = []
        for birth, death in dgm:
            plotted_death = float(death) if np.isfinite(death) else inf_cap
            lifespans.append((float(birth), plotted_death))

        lifespans.sort(key=lambda x: x[1] - x[0], reverse=True)

        added_label = False
        for birth, death in lifespans:
            ax.hlines(
                y=y_pos,
                xmin=birth,
                xmax=death,
                color=colors[dim] if dim < len(colors) else None,
                linewidth=1.5,
                label=labels[dim] if not added_label and dim < len(labels) else None,
            )
            added_label = True
            y_pos += 1

    ax.set_xlabel("Filtration Value (Distance)")
    ax.set_ylabel("Topological Features (ranked by lifespan)")
    ax.set_title("Persistence Barcode")
    ax.set_yticks([])
    if y_pos:
        ax.legend(loc="lower right")

    fig.tight_layout()
    return _save(fig, output)


def plot_h1_distribution(persistence, chromosome, start, end, output=None):
    persistence = np.asarray(persistence)
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(
        persistence,
        bins=30,
        edgecolor="black",
        alpha=0.75,
    )
    ax.set_title(
        f"$H_1$ Persistence Distribution "
        f"(Chr{chromosome}: {start/1e6:.1f}–{end/1e6:.1f} Mb)"
    )
    ax.set_xlabel("Persistence lifespan (death − birth)")
    ax.set_ylabel("Number of $H_1$ features")
    ax.grid(axis="y", linestyle="--", alpha=0.7)
    fig.tight_layout()
    return _save(fig, output)


def plot_h1_boxplot(persistence, chromosome, start, end, output=None):
    persistence = np.asarray(persistence)
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.boxplot(
        persistence,
        patch_artist=True,
        boxprops=dict(facecolor="#ff7f0e", color="black"),
        medianprops=dict(color="black", linewidth=2),
        whiskerprops=dict(color="black"),
        capprops=dict(color="black"),
        flierprops=dict(
            marker="o",
            markerfacecolor="red",
            markersize=5,
            linestyle="none",
        ),
    )
    ax.set_title(
        f"$H_1$ Persistence Box Plot "
        f"(Chr{chromosome}: {start/1e6:.1f}–{end/1e6:.1f} Mb)"
    )
    ax.set_ylabel("Persistence value (death − birth)")
    ax.set_xticks([1])
    ax.set_xticklabels([r"$H_1$ features"])
    ax.grid(axis="y", linestyle="--", alpha=0.7)
    fig.tight_layout()
    return _save(fig, output)


def plot_summary(
    matrix,
    diagrams,
    persistence,
    metadata,
    stats,
    output=None,
):
    """Generate a compact six-panel summary figure."""
    fig = plt.figure(figsize=(20, 12))
    grid = fig.add_gridspec(2, 3)

    ax1 = fig.add_subplot(grid[0, 0])
    log_matrix = np.log1p(matrix)
    cmap = LinearSegmentedColormap.from_list(
        "hictda_summary", ["white", "red", "darkred"]
    )
    im = ax1.imshow(log_matrix, cmap=cmap, interpolation="nearest")
    ax1.set_title("Hi-C Contact Map")
    fig.colorbar(im, ax=ax1, fraction=0.046, pad=0.04)

    ax2 = fig.add_subplot(grid[0, 1])
    plot_diagrams(diagrams, show=False, ax=ax2, title="Persistence Diagram")
    ax2.set_xlabel("Birth")
    ax2.set_ylabel("Death")

    ax3 = fig.add_subplot(grid[0, 2])
    max_finite = 0.0
    for dgm in diagrams:
        if len(dgm):
            finite = dgm[np.isfinite(dgm[:, 1]), 1]
            if finite.size:
                max_finite = max(max_finite, float(np.max(finite)))
    cap = max_finite * 1.1 if max_finite else 1.0
    y = 0
    for dim, dgm in enumerate(diagrams[:2]):
        for birth, death in dgm:
            d = float(death) if np.isfinite(death) else cap
            ax3.hlines(y, birth, d, linewidth=1.2)
            y += 1
    ax3.set_title("Persistence Barcode")
    ax3.set_xlabel("Distance")
    ax3.set_yticks([])

    ax4 = fig.add_subplot(grid[1, 0])
    ax4.hist(persistence, bins=30, edgecolor="black", alpha=0.75)
    ax4.set_title("$H_1$ Persistence Distribution")
    ax4.set_xlabel("Persistence")
    ax4.set_ylabel("Number of features")

    ax5 = fig.add_subplot(grid[1, 1])
    ax5.boxplot(persistence, patch_artist=True)
    ax5.set_title("$H_1$ Persistence Box Plot")
    ax5.set_ylabel("Persistence")

    ax6 = fig.add_subplot(grid[1, 2])
    ax6.axis("off")
    summary = (
        "$H_1$ PERSISTENCE SUMMARY\n\n"
        f"Chromosome: {metadata['chromosome']}\n"
        f"Region: {metadata['start']/1e6:.2f}–{metadata['end']/1e6:.2f} Mb\n"
        f"Resolution: {metadata['resolution']:,} bp\n"
        f"Normalization: {metadata['normalization']}\n\n"
        f"Number of $H_1$ features: {stats['n_h1']}\n"
        f"Mean persistence: {stats['mean_persistence']:.6g}\n"
        f"Median persistence: {stats['median_persistence']:.6g}\n"
        f"Standard deviation: {stats['std_persistence']:.6g}"
    )
    ax6.text(
        0.05, 0.95, summary,
        va="top",
        ha="left",
        fontsize=12,
        family="monospace",
    )

    fig.suptitle("HiCTDA Persistent-Homology Analysis", fontsize=16)
    fig.tight_layout()
    return _save(fig, output)
