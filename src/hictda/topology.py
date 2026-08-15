"""Distance construction and persistent homology."""

from __future__ import annotations

import numpy as np
from ripser import ripser


def contact_to_distance(
    matrix: np.ndarray,
    method: str = "inverse",
    epsilon: float = 1e-5,
) -> np.ndarray:
    """Convert a Hi-C contact-frequency matrix to a dissimilarity matrix."""
    matrix = np.asarray(matrix, dtype=float)

    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("matrix must be square.")
    if not np.all(np.isfinite(matrix)):
        raise ValueError("matrix contains non-finite values.")
    if np.any(matrix < 0):
        raise ValueError("contact frequencies must be non-negative.")
    if epsilon <= 0:
        raise ValueError("epsilon must be positive.")

    method = method.lower()

    if method == "inverse":
        distance = 1.0 / (matrix + epsilon)
    elif method == "negative_log":
        distance = -np.log(matrix + epsilon)
        distance -= np.min(distance)
    elif method == "max_minus_contact":
        distance = np.max(matrix) - matrix
    else:
        raise ValueError(
            f"Unknown distance method '{method}'. "
            "Supported methods: inverse, negative_log, max_minus_contact."
        )

    distance = np.asarray(distance, dtype=float)
    distance = (distance + distance.T) / 2.0
    np.fill_diagonal(distance, 0.0)

    if not np.all(np.isfinite(distance)):
        raise ValueError("Distance transformation produced non-finite values.")

    return distance


def compute_persistence(
    distance_matrix: np.ndarray,
    maxdim: int = 1,
) -> dict:
    """Compute Vietoris-Rips persistent homology from a distance matrix."""
    if maxdim < 0:
        raise ValueError("maxdim must be >= 0.")

    distance_matrix = np.asarray(distance_matrix, dtype=float)

    if distance_matrix.ndim != 2:
        raise ValueError("distance_matrix must be two-dimensional.")
    if distance_matrix.shape[0] != distance_matrix.shape[1]:
        raise ValueError("distance_matrix must be square.")
    if not np.all(np.isfinite(distance_matrix)):
        raise ValueError("distance_matrix contains non-finite values.")

    # Numerical symmetrization.
    distance_matrix = (distance_matrix + distance_matrix.T) / 2.0
    np.fill_diagonal(distance_matrix, 0.0)

    return ripser(
        distance_matrix,
        distance_matrix=True,
        maxdim=maxdim,
    )
