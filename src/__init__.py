# Ana Carsi 2024
# Project: Review on Adaptive Canonical Correlation Analysis

from .cca_functions import (
    cca_on_stiefel,
    cca_with_tracking,
    canonical_correlation_heatmap,
    optimization_path_plot,
)

__all__ = [
    "cca_on_stiefel",
    "cca_with_tracking",
    "canonical_correlation_heatmap",
    "optimization_path_plot",
]
