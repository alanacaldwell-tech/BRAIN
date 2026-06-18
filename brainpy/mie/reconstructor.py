"""
Reconstruct missing isotope peak intensities for MIE-affected clusters.

For each gap index in a classified cluster the reconstructor:

  1. Generates (or reuses) the averagine theoretical envelope.
  2. Finds the best alignment of observed peaks to the theoretical series
     (same alignment used by the classifier).
  3. Fits a single scale factor *s* from the surviving observed peaks via
     ordinary least squares.
  4. Predicts the missing peak intensity as  s × theoretical[gap].
  5. Estimates a 1-sigma prediction interval from the fit residuals,
     widened by a penalty that grows with the number of consecutive gaps
     (the more consecutive peaks are missing, the less constrained the fit).

Reconstruction is skipped for 'overlap' clusters (two real species) and
returns an empty list for clusters that have no gaps.
"""
import numpy as np
from typing import List

from .averagine import theoretical_envelope
from .classifier import _fit_alignment
from .peaks import IsotopicCluster, ClassificationResult, ReconstructedPeak

ISOTOPE_SPACING = 1.003355  # Da


def reconstruct(
    cluster: IsotopicCluster,
    classification: ClassificationResult,
    npeaks: int = 25,
) -> List[ReconstructedPeak]:
    """
    Return a list of ReconstructedPeak objects for every gap in *cluster*.

    For 'overlap' clusters an empty list is returned because no single-species
    reconstruction is appropriate.

    For 'ambiguous' clusters the peaks are still estimated but inherit the
    reduced confidence score from *classification*.
    """
    if classification.hypothesis == "overlap":
        return []
    if not cluster.gap_indices:
        return []

    theoretical = theoretical_envelope(cluster.neutral_mass, cluster.charge, npeaks=npeaks)
    if not theoretical:
        return []

    theo_offset, scale, _r2 = _fit_alignment(cluster.peaks, theoretical)

    # --- Residuals of the fit (observed peaks only) -----------------------
    obs = {p.isotope_index: p.intensity for p in cluster.peaks}
    residuals = []
    for obs_idx, obs_int in obs.items():
        theo_idx = obs_idx + theo_offset
        if 0 <= theo_idx < len(theoretical):
            pred = scale * theoretical[theo_idx].intensity
            residuals.append(obs_int - pred)

    residual_std = float(np.std(residuals)) if len(residuals) >= 2 else 0.0
    n_obs = max(len(residuals), 1)

    # Prediction interval base: std × sqrt(1 + 1/n)
    pred_interval_base = residual_std * np.sqrt(1.0 + 1.0 / n_obs)

    # Penalty multiplier grows with consecutive gaps
    n_consec = classification.n_consecutive_gaps
    consec_penalty = 1.0 + 0.5 * max(0, n_consec - 1)

    delta = ISOTOPE_SPACING / cluster.charge
    first_mz = cluster.peaks[0].mz

    reconstructed: List[ReconstructedPeak] = []

    for gap_idx in cluster.gap_indices:
        theo_idx = gap_idx + theo_offset

        # m/z of the missing peak (cluster indices are relative to peaks[0])
        missing_mz = first_mz + gap_idx * delta

        if 0 <= theo_idx < len(theoretical):
            theo_intensity = theoretical[theo_idx].intensity
            estimated = max(0.0, scale * theo_intensity)
            uncertainty = float(pred_interval_base * consec_penalty)
        else:
            # Gap is outside the theoretical range — extrapolation only
            estimated = 0.0
            uncertainty = float("inf")

        reconstructed.append(
            ReconstructedPeak(
                mz=missing_mz,
                charge=cluster.charge,
                estimated_intensity=estimated,
                uncertainty=uncertainty,
                isotope_index=gap_idx,
                confidence=classification.confidence,
            )
        )

    return reconstructed
