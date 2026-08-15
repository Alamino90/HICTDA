"""HiCTDA: persistent-homology analysis of Hi-C contact maps."""

from .pipeline import AnalysisResult, analyze_hic
from .io import load_hic, list_chromosomes, available_resolutions, extract_matrix
from .topology import contact_to_distance, compute_persistence
from .statistics import h1_statistics

__version__ = "0.1.1"

__all__ = [
    "AnalysisResult",
    "analyze_hic",
    "load_hic",
    "list_chromosomes",
    "available_resolutions",
    "extract_matrix",
    "contact_to_distance",
    "compute_persistence",
    "h1_statistics",
]
