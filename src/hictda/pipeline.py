"""High-level HiCTDA analysis pipeline."""

from __future__ import annotations

from dataclasses import dataclass, asdict
from pathlib import Path
import json
import numpy as np

from .io import load_hic, extract_matrix
from .topology import contact_to_distance, compute_persistence
from .statistics import extract_h1_persistence, h1_statistics
from .visualization import (
    plot_contact_map,
    plot_persistence_diagram,
    plot_barcode,
    plot_h1_distribution,
    plot_h1_boxplot,
    plot_summary,
)


@dataclass
class AnalysisResult:
    """Container returned by analyze_hic."""

    matrix: np.ndarray
    distance_matrix: np.ndarray
    diagrams: list[np.ndarray]
    h1_persistence: np.ndarray
    statistics: dict
    metadata: dict
    figures: dict

    def save_arrays(self, output_dir: str | Path) -> None:
        """Save numerical arrays to disk."""
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)

        np.save(out / "contact_matrix.npy", self.matrix)
        np.save(out / "distance_matrix.npy", self.distance_matrix)
        np.save(out / "h1_persistence.npy", self.h1_persistence)

        for dim, diagram in enumerate(self.diagrams):
            np.save(out / f"persistence_diagram_h{dim}.npy", diagram)


def _json_safe(obj):
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, Path):
        return str(obj)
    raise TypeError(f"Object of type {type(obj)} is not JSON serializable")


def _save_csv(persistence, path):
    with open(path, "w", encoding="utf-8") as f:
        f.write("h1_feature,persistence\n")
        for i, value in enumerate(persistence, start=1):
            f.write(f"{i},{value:.17g}\n")


def analyze_hic(
    hic_file,
    chromosome: str,
    start: int | None = None,
    end: int | None = None,
    resolution: int = 50_000,
    normalization: str = "KR",
    maxdim: int = 1,
    distance_method: str = "inverse",
    epsilon: float = 1e-5,
    output_dir: str | Path | None = None,
    save_arrays: bool = True,
) -> AnalysisResult:
    """
    Run the complete HiCTDA pipeline.

    Normalization is mandatory. If start/end are omitted, the entire
    chromosome is analyzed.
    """
    if not normalization or str(normalization).upper() == "NONE":
        raise ValueError(
            "Normalization is mandatory. Specify a supported normalized "
            "Hi-C matrix type, e.g. normalization='KR'."
        )

    hic = load_hic(hic_file)

    matrix, metadata = extract_matrix(
        hic=hic,
        chromosome=chromosome,
        start=start,
        end=end,
        resolution=resolution,
        normalization=normalization,
    )

    distance_matrix = contact_to_distance(
        matrix,
        method=distance_method,
        epsilon=epsilon,
    )

    ph = compute_persistence(
        distance_matrix,
        maxdim=maxdim,
    )
    diagrams = ph["dgms"]

    persistence = extract_h1_persistence(diagrams)
    stats = h1_statistics(diagrams)

    metadata.update({
        "hic_file": str(hic_file),
        "distance_method": distance_method,
        "epsilon": epsilon,
        "maxdim": maxdim,
    })

    figures = {}

    if output_dir is not None:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)

        figures["contact_map"] = plot_contact_map(
            matrix,
            metadata["chromosome"],
            metadata["start"],
            metadata["end"],
            out / "contact_map.png",
        )
        figures["persistence_diagram"] = plot_persistence_diagram(
            diagrams,
            out / "persistence_diagram.png",
        )
        figures["barcode"] = plot_barcode(
            diagrams,
            out / "persistence_barcode.png",
        )
        figures["h1_distribution"] = plot_h1_distribution(
            persistence,
            metadata["chromosome"],
            metadata["start"],
            metadata["end"],
            out / "h1_distribution.png",
        )
        figures["h1_boxplot"] = plot_h1_boxplot(
            persistence,
            metadata["chromosome"],
            metadata["start"],
            metadata["end"],
            out / "h1_boxplot.png",
        )
        figures["summary"] = plot_summary(
            matrix,
            diagrams,
            persistence,
            metadata,
            stats,
            out / "hictda_summary.png",
        )

        _save_csv(persistence, out / "h1_persistence.csv")

        summary = {
            "metadata": metadata,
            "statistics": stats,
        }
        with open(out / "summary.json", "w", encoding="utf-8") as f:
            json.dump(summary, f, indent=2, default=_json_safe)

        if save_arrays:
            result_for_save = AnalysisResult(
                matrix=matrix,
                distance_matrix=distance_matrix,
                diagrams=diagrams,
                h1_persistence=persistence,
                statistics=stats,
                metadata=metadata,
                figures=figures,
            )
            result_for_save.save_arrays(out)

    return AnalysisResult(
        matrix=matrix,
        distance_matrix=distance_matrix,
        diagrams=diagrams,
        h1_persistence=persistence,
        statistics=stats,
        metadata=metadata,
        figures=figures,
    )
