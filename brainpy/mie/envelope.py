"""
Isotopic cluster detection with gap awareness.

Clusters are built from individual (mz, intensity) pairs by linking peaks
that are separated by integer multiples of the charge-state-specific isotope
spacing (1.003355 / z Da).  Gaps of up to *max_gap* consecutive missing
isotopes are tolerated so that MIE-induced holes do not split one cluster
into two.
"""
from collections import defaultdict
from typing import List, Tuple

from brainpy import PROTON

from .peaks import IsotopicCluster, ObservedPeak

ISOTOPE_SPACING = 1.003355  # Da  (¹³C – ¹²C mass difference)


# ---------------------------------------------------------------------------
# Ion binning
# ---------------------------------------------------------------------------

def bin_ions(mz_values, bin_width: float = 0.005) -> List[Tuple[float, float]]:
    """
    Aggregate individual ion detection events into (mz, count) peaks.

    *bin_width* should be well below the isotope spacing (1/z Da) for the
    charge states of interest — 0.005 Da works for z up to ~200.

    Returns a list of (bin_centre_mz, ion_count) sorted ascending by m/z.
    """
    bins: dict = {}
    for mz in mz_values:
        key = round(float(mz) / bin_width) * bin_width
        bins[key] = bins.get(key, 0.0) + 1.0
    return sorted(bins.items())


# ---------------------------------------------------------------------------
# Cluster detection
# ---------------------------------------------------------------------------

def find_isotopic_clusters(
    peaks: List[Tuple[float, float]],
    charge: int,
    max_gap: int = 3,
    ppm_tolerance: float = 20.0,
) -> List[IsotopicCluster]:
    """
    Group (mz, intensity) pairs into isotopic clusters for *charge*.

    Two peaks are joined in the same cluster when their m/z difference is
    within *ppm_tolerance* ppm of k × (ISOTOPE_SPACING / charge) for
    k ∈ {1 … max_gap}.  This allows clusters to span gaps produced by
    multi-ion event removal.

    Each returned IsotopicCluster carries:
      • peaks          — observed peaks with cluster-relative isotope indices
      • gap_indices    — cluster indices where no peak was detected
      • neutral_mass   — estimated from the lowest-m/z observed peak
    """
    if not peaks or charge <= 0:
        return []

    peaks = sorted(peaks, key=lambda p: p[0])
    n = len(peaks)
    delta = ISOTOPE_SPACING / charge

    # --- Union-Find --------------------------------------------------------
    parent = list(range(n))

    def _find(x: int) -> int:
        root = x
        while parent[root] != root:
            root = parent[root]
        while parent[x] != root:
            parent[x], x = root, parent[x]
        return root

    def _union(x: int, y: int) -> None:
        parent[_find(x)] = _find(y)

    for i in range(n):
        mz_i = peaks[i][0]
        for j in range(i + 1, n):
            mz_j = peaks[j][0]
            diff = mz_j - mz_i
            # Stop scanning forward once we're definitely past max_gap steps
            if diff > (max_gap + 0.5) * delta:
                break
            # Tolerance: use the larger of ppm-derived or 0.01 Da floor
            tol = max(mz_i * ppm_tolerance * 1e-6, 0.01)
            for k in range(1, max_gap + 1):
                if abs(diff - k * delta) < tol:
                    _union(i, j)
                    break

    # --- Group components --------------------------------------------------
    components: dict = defaultdict(list)
    for i in range(n):
        components[_find(i)].append(i)

    clusters: List[IsotopicCluster] = []
    for indices in components.values():
        if len(indices) < 2:
            continue

        comp_peaks = [peaks[i] for i in sorted(indices)]
        mz_start = comp_peaks[0][0]

        # Assign isotope indices relative to the first observed peak
        indexed = []
        for mz, intensity in comp_peaks:
            idx = round((mz - mz_start) / delta)
            indexed.append((mz, float(intensity), idx))

        max_idx = max(p[2] for p in indexed)
        observed_set = {p[2] for p in indexed}
        gap_indices = sorted(set(range(max_idx + 1)) - observed_set)

        # Neutral mass: mz of the lowest-m/z peak back-calculated to mass.
        # NOTE: swap this line for the formula from branch claude/awesome-planck-25y17c
        # once that branch is available.
        neutral_mass = mz_start * charge - charge * PROTON

        observed_peaks = [
            ObservedPeak(
                mz=mz,
                intensity=intensity,
                charge=charge,
                isotope_index=idx,
            )
            for mz, intensity, idx in indexed
        ]

        clusters.append(
            IsotopicCluster(
                peaks=observed_peaks,
                charge=charge,
                gap_indices=gap_indices,
                neutral_mass=neutral_mass,
            )
        )

    return clusters
