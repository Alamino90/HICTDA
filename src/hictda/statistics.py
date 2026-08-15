"""Persistence statistics."""

from __future__ import annotations

import numpy as np


def extract_h1_persistence(diagrams: list[np.ndarray]) -> np.ndarray:
    """Return finite H1 persistence lifetimes."""
    if len(diagrams) < 2:
        return np.empty(0, dtype=float)

    h1 = np.asarray(diagrams[1], dtype=float)
    if h1.size == 0:
        return np.empty(0, dtype=float)

    finite = h1[np.isfinite(h1[:, 1])]
    if finite.size == 0:
        return np.empty(0, dtype=float)

    persistence = finite[:, 1] - finite[:, 0]
    persistence = persistence[np.isfinite(persistence)]
    return persistence[persistence >= 0]


def h1_statistics(diagrams: list[np.ndarray]) -> dict:
    """Calculate H1 count, mean, median, and standard deviation."""
    persistence = extract_h1_persistence(diagrams)
    n = int(persistence.size)

    return {
        "n_h1": n,
        "mean_persistence": float(np.mean(persistence)) if n else 0.0,
        "median_persistence": float(np.median(persistence)) if n else 0.0,
        "std_persistence": float(np.std(persistence)) if n else 0.0,
    }
