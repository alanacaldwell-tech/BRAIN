"""
brainpy.mie — Multi-Ion Event (MIE) artifact correction for mass spectrometry.

Quick start
-----------
    from brainpy.mie import process_dmt

    results = process_dmt("experiment.dmt")
    rescued = results[results["IsEstimated"]]   # reconstructed missing peaks

Lower-level access
------------------
    from brainpy.mie import read_ions, find_isotopic_clusters, classify, reconstruct
"""
from .pipeline import process_dmt
from .reader import read_ions
from .envelope import bin_ions, find_isotopic_clusters
from .averagine import averagine_composition, theoretical_envelope
from .classifier import classify
from .reconstructor import reconstruct
from .peaks import ObservedPeak, IsotopicCluster, ClassificationResult, ReconstructedPeak

__all__ = [
    "process_dmt",
    "read_ions",
    "bin_ions",
    "find_isotopic_clusters",
    "averagine_composition",
    "theoretical_envelope",
    "classify",
    "reconstruct",
    "ObservedPeak",
    "IsotopicCluster",
    "ClassificationResult",
    "ReconstructedPeak",
]
