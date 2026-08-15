"""Hi-C file access, metadata validation, and matrix extraction."""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Union

import hicstraw
import numpy as np


PathLike = Union[str, Path]


def load_hic(hic_file: PathLike):
    """Open a local path or URL understood by hic-straw."""
    return hicstraw.HiCFile(str(hic_file))


def list_chromosomes(hic) -> list[dict]:
    """Return chromosomes available in the .hic file."""
    return [
        {"name": str(c.name), "length": int(c.length)}
        for c in hic.getChromosomes()
        if str(c.name).lower() != "all"
    ]


def _find_chromosome(hic, chromosome: str):
    """Find a chromosome by name, accepting both '3' and 'chr3'."""
    target = str(chromosome).lower().removeprefix("chr")

    for c in hic.getChromosomes():
        name = str(c.name).lower().removeprefix("chr")

        if name == target:
            return c

    available = ", ".join(
        c["name"] for c in list_chromosomes(hic)
    )

    raise ValueError(
        f"Chromosome '{chromosome}' is not available in the supplied .hic file. "
        f"Available chromosomes: {available}"
    )


def available_resolutions(hic, chromosome: Optional[str] = None) -> list[int]:
    """Return resolutions exposed by hic-straw for the file."""
    # hic-straw exposes base-pair resolutions through the file's attributes.
    # The chromosome argument is accepted for a stable public API; resolution
    # availability is generally file-wide in .hic files.
    resolutions = []
    for attr in ("getResolutions", "getBpResolutions"):
        if hasattr(hic, attr):
            try:
                vals = getattr(hic, attr)()
                resolutions = [int(x) for x in vals]
                break
            except Exception:
                pass
    return sorted(set(resolutions))


def validate_region(
    hic,
    chromosome: str,
    start: Optional[int],
    end: Optional[int],
    resolution: int,
) -> tuple[str, int, int, int]:
    """Validate and canonicalize chromosome and genomic coordinates."""
    chrom_obj = _find_chromosome(hic, chromosome)
    chrom_name = str(chrom_obj.name)
    chrom_length = int(chrom_obj.length)

    if not isinstance(resolution, int) or resolution <= 0:
        raise ValueError("resolution must be a positive integer in base pairs.")

    if start is None:
        start = 0
    if end is None:
        end = chrom_length

    if not isinstance(start, int) or not isinstance(end, int):
        raise TypeError("start and end must be integers in base pairs.")

    if start < 0 or end < 0:
        raise ValueError("start and end must be non-negative.")
    if start >= end:
        raise ValueError("The genomic interval must satisfy start < end.")
    if end > chrom_length:
        raise ValueError(
            f"end={end:,} exceeds chromosome {chrom_name} length "
            f"{chrom_length:,} bp."
        )

    return chrom_name, start, end, resolution


def extract_matrix(
    hic,
    chromosome: str,
    start: Optional[int],
    end: Optional[int],
    resolution: int,
    normalization: str = "KR",
) -> tuple[np.ndarray, dict]:
    """Extract a normalized observed intra-chromosomal Hi-C matrix."""
    chrom_name, start, end, resolution = validate_region(
        hic, chromosome, start, end, resolution
    )

    if not normalization or str(normalization).upper() == "NONE":
        raise ValueError(
            "Normalization is mandatory for the high-level HiCTDA analysis."
        )
    normalization = str(normalization).upper()

    try:
        mzd = hic.getMatrixZoomData(
            chrom_name,
            chrom_name,
            "observed",
            normalization,
            "BP",
            resolution,
        )
    except Exception as exc:
        raise ValueError(
            f"Could not obtain a normalized Hi-C matrix for chromosome "
            f"{chrom_name}, resolution {resolution}, normalization "
            f"{normalization}. Check that the requested normalization and "
            f"resolution are present in the .hic file."
        ) from exc

    matrix = np.asarray(
        mzd.getRecordsAsMatrix(start, end, start, end),
        dtype=float,
    )

    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError(
            f"Expected a square intra-chromosomal matrix; received shape "
            f"{matrix.shape}."
        )

    matrix = np.nan_to_num(matrix, nan=0.0, posinf=0.0, neginf=0.0)
    matrix[matrix < 0] = 0.0

    metadata = {
        "chromosome": chrom_name,
        "start": start,
        "end": end,
        "resolution": resolution,
        "normalization": normalization,
        "n_bins": int(matrix.shape[0]),
        "matrix_shape": [int(x) for x in matrix.shape],
    }
    return matrix, metadata
