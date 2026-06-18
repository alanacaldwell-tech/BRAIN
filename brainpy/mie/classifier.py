"""
Classify isotopic clusters with gaps as MIE artifacts or true overlapping species.

Decision logic
--------------
Three features drive the classification:

1. gap_at_apex  — Is the missing peak where the theoretical intensity maximum
                  would be?  This is the strongest single indicator of MIE:
                  the apex is the most abundant peak and therefore the most
                  likely to trigger a multi-ion coincidence.

2. fit_r_squared — How well does a single-species averagine envelope fit the
                   *surviving* peaks?  A high R² with a gap at the apex is the
                   canonical MIE signature.  A low R² suggests two distinct
                   envelopes are present.

3. n_consecutive_gaps — More consecutive missing peaks mean less evidence and
                        lower reconstruction confidence.

Scoring formula (additive, then thresholded):
  mie_score  = 0.50 × gap_at_apex
             + [0.35 / 0.20 / 0.05] × r² tier
             - 0.05 × n_consecutive_gaps   (capped at 0.15)

  mie_score ≥ 0.65  →  'mie'
  mie_score ≤ 0.30  →  'overlap'
  otherwise         →  'ambiguous'
"""
import numpy as np
from typing import List, Tuple

from .averagine import theoretical_envelope
from .peaks import IsotopicCluster, ClassificationResult


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _fit_alignment(
    observed_peaks,
    theoretical_peaks,
) -> Tuple[int, float, float]:
    """
    Find which theoretical peak best aligns with the first observed peak, and
    fit a single scale factor s such that observed ≈ s × theoretical.

    Returns (theo_offset, scale, r_squared).

    theo_offset is the index into *theoretical_peaks* that corresponds to
    observed isotope_index 0.
    """
    obs = {p.isotope_index: p.intensity for p in observed_peaks}
    n_theo = len(theoretical_peaks)
    if n_theo == 0 or not obs:
        return 0, 1.0, 0.0

    best_r2 = -np.inf
    best_offset = 0
    best_scale = 1.0

    # Try each possible alignment (don't need to check more than n_theo offsets)
    for offset in range(min(n_theo, 15)):
        pairs = []
        for obs_idx, obs_int in obs.items():
            theo_idx = obs_idx + offset
            if 0 <= theo_idx < n_theo:
                pairs.append((obs_int, theoretical_peaks[theo_idx].intensity))

        if len(pairs) < 2:
            continue

        obs_arr = np.array([p[0] for p in pairs])
        theo_arr = np.array([p[1] for p in pairs])

        denom = float(np.sum(theo_arr ** 2))
        if denom == 0:
            continue
        scale = float(np.dot(obs_arr, theo_arr)) / denom

        pred = scale * theo_arr
        ss_res = float(np.sum((obs_arr - pred) ** 2))
        ss_tot = float(np.sum((obs_arr - obs_arr.mean()) ** 2))
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0

        if r2 > best_r2:
            best_r2 = r2
            best_offset = offset
            best_scale = scale

    return best_offset, best_scale, float(best_r2)


def _apex_in_gaps(gap_indices: List[int], theoretical_peaks, theo_offset: int) -> bool:
    """Return True if any gap index maps to the theoretical intensity maximum."""
    if not gap_indices or not theoretical_peaks:
        return False
    apex_theo_idx = int(np.argmax([p.intensity for p in theoretical_peaks]))
    # Convert from theoretical index to cluster index
    apex_cluster_idx = apex_theo_idx - theo_offset
    return apex_cluster_idx in set(gap_indices)


def _max_consecutive_gaps(gap_indices: List[int]) -> int:
    """Return the length of the longest run of consecutive missing isotope indices."""
    if not gap_indices:
        return 0
    gaps = sorted(gap_indices)
    max_run = current_run = 1
    for i in range(1, len(gaps)):
        if gaps[i] == gaps[i - 1] + 1:
            current_run += 1
            max_run = max(max_run, current_run)
        else:
            current_run = 1
    return max_run


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def classify(cluster: IsotopicCluster, npeaks: int = 25) -> ClassificationResult:
    """
    Classify *cluster* as an MIE artifact, a true spectral overlap, or ambiguous.

    Clusters without any gaps are returned immediately as 'no_gap'.
    """
    if not cluster.gap_indices:
        return ClassificationResult(
            hypothesis="no_gap",
            confidence=1.0,
            fit_r_squared=1.0,
            gap_at_apex=False,
            n_consecutive_gaps=0,
        )

    theoretical = theoretical_envelope(cluster.neutral_mass, cluster.charge, npeaks=npeaks)
    if not theoretical:
        return ClassificationResult(
            hypothesis="ambiguous",
            confidence=0.0,
            fit_r_squared=0.0,
            gap_at_apex=False,
            n_consecutive_gaps=_max_consecutive_gaps(cluster.gap_indices),
            notes="Could not generate theoretical envelope for this mass/charge.",
        )

    theo_offset, _scale, r2 = _fit_alignment(cluster.peaks, theoretical)
    apex_gap = _apex_in_gaps(cluster.gap_indices, theoretical, theo_offset)
    n_consec = _max_consecutive_gaps(cluster.gap_indices)

    # --- Scoring -----------------------------------------------------------
    mie_score = 0.0

    if apex_gap:
        mie_score += 0.50

    if r2 >= 0.85:
        mie_score += 0.35
    elif r2 >= 0.70:
        mie_score += 0.20
    elif r2 >= 0.50:
        mie_score += 0.05

    # Consecutive gaps reduce confidence (penalise up to -0.15)
    mie_score -= min(0.15, n_consec * 0.05)

    # Reconstruction reliability decays with consecutive missing peaks
    recon_confidence = max(0.0, 1.0 - (n_consec - 1) * 0.25)

    # --- Decision ----------------------------------------------------------
    notes_parts = []

    if mie_score >= 0.65:
        hypothesis = "mie"
        confidence = min(1.0, mie_score) * recon_confidence
    elif mie_score <= 0.30:
        hypothesis = "overlap"
        confidence = min(1.0, max(0.0, 1.0 - mie_score))
    else:
        hypothesis = "ambiguous"
        confidence = recon_confidence * (0.5 + abs(mie_score - 0.475))
        notes_parts.append("borderline score — treat reconstruction with caution")

    if n_consec >= 2:
        notes_parts.append(
            f"{n_consec} consecutive missing peaks: reconstruction uncertainty is elevated"
        )
    if r2 < 0.50:
        notes_parts.append("poor single-species fit (R²={:.2f})".format(r2))
    if not apex_gap and hypothesis == "mie":
        notes_parts.append("gap is not at the apex — verify MIE assignment")

    return ClassificationResult(
        hypothesis=hypothesis,
        confidence=confidence,
        fit_r_squared=r2,
        gap_at_apex=apex_gap,
        n_consecutive_gaps=n_consec,
        notes="; ".join(notes_parts),
    )
