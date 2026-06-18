from dataclasses import dataclass, field
from typing import List, Optional


@dataclass
class ObservedPeak:
    mz: float
    intensity: float
    charge: int
    isotope_index: int  # index relative to the first observed peak in the cluster


@dataclass
class IsotopicCluster:
    peaks: List[ObservedPeak]
    charge: int
    gap_indices: List[int]   # cluster-relative indices where isotope peaks are absent
    neutral_mass: float      # estimated from the lowest-m/z observed peak


@dataclass
class ClassificationResult:
    """
    Outcome of single-species-vs-overlap hypothesis test for a gapped cluster.

    hypothesis: 'mie'       — single species, apex removed by multi-ion event
                'overlap'   — two distinct co-eluting species
                'ambiguous' — insufficient evidence to decide
                'no_gap'    — cluster has no missing peaks; classification not needed
    confidence: 0–1 probability-like score for the stated hypothesis
    fit_r_squared: R² of the single-species averagine fit over observed peaks
    gap_at_apex: whether any gap falls at the theoretical intensity maximum
    n_consecutive_gaps: longest run of consecutive missing isotope peaks
    notes: human-readable caveats (e.g. reconstruction reliability warnings)
    """
    hypothesis: str
    confidence: float
    fit_r_squared: float
    gap_at_apex: bool
    n_consecutive_gaps: int
    notes: str = ""


@dataclass
class ReconstructedPeak:
    """A missing peak whose intensity has been estimated from the surviving envelope."""
    mz: float
    charge: int
    estimated_intensity: float
    uncertainty: float    # 1-sigma; inf when the gap is out of the theoretical range
    isotope_index: int    # cluster-relative index matching the gap
    confidence: float     # inherits from ClassificationResult.confidence
