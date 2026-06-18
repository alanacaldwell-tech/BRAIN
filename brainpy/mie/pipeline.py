"""
End-to-end pipeline: .dmt file → corrected isotopic peak table.

Typical usage
-------------
    from brainpy.mie import process_dmt

    df = process_dmt("experiment.dmt")
    # All observed peaks plus reconstructed MIE peaks in one DataFrame.

    # Only the reconstructed (rescued) peaks:
    rescued = df[df["IsEstimated"]]

    # Write corrected data back to SQLite:
    import sqlite3
    conn = sqlite3.connect("corrected.dmt")
    df.to_sql("CorrectedIon", conn, if_exists="replace", index=False)
    conn.close()
"""
import pandas as pd
from typing import Tuple

from .reader import read_ions
from .envelope import bin_ions, find_isotopic_clusters
from .classifier import classify
from .reconstructor import reconstruct


def process_dmt(
    filepath: str,
    charge_range: Tuple[int, int] = (1, 10),
    bin_width: float = 0.005,
    max_gap: int = 3,
    ppm_tolerance: float = 20.0,
    min_cluster_peaks: int = 3,
    npeaks: int = 25,
) -> pd.DataFrame:
    """
    Read a .dmt SQLite file and return a corrected isotopic peak DataFrame.

    Parameters
    ----------
    filepath : str
        Path to the .dmt file.
    charge_range : (int, int)
        Inclusive range of charge states to process.
    bin_width : float
        m/z bin width in Da used to aggregate individual ion events into peaks.
        Should be << 1/max_charge.  Default 0.005 Da.
    max_gap : int
        Maximum number of consecutive missing isotope peaks that are still
        treated as belonging to the same cluster.  Default 3.
    ppm_tolerance : float
        m/z tolerance in ppm for assigning peaks to isotope positions.
    min_cluster_peaks : int
        Clusters with fewer observed peaks than this (and no gaps) are skipped.
    npeaks : int
        Number of theoretical isotope peaks to generate per cluster.

    Returns
    -------
    pd.DataFrame with columns:

    Mz                  — observed m/z (NaN for estimated peaks)
    Charge              — charge state
    Intensity           — observed ion count (NaN for estimated peaks)
    NeutralMass         — estimated neutral mass of the cluster
    IsotopicIndex       — 0-based index within the isotopic series
    IsEstimated         — True for reconstructed MIE peaks
    EstimatedIntensity  — predicted ion count (populated for estimated peaks)
    Uncertainty         — 1-sigma prediction interval of EstimatedIntensity
    Confidence          — 0–1 score for the stated hypothesis / reconstruction
    Hypothesis          — 'mie', 'overlap', 'ambiguous', or 'no_gap'
    FitRSquared         — R² of the single-species averagine fit
    GapAtApex           — whether the gap sits at the theoretical envelope apex
    Notes               — human-readable caveats
    """
    ions_df = read_ions(filepath)
    records = []

    for charge in range(charge_range[0], charge_range[1] + 1):
        subset = ions_df[ions_df["Charge"] == charge]
        if subset.empty:
            continue

        if "Intensity" in subset.columns:
            # Pre-aggregated peaks: use Intensity directly
            peaks = list(zip(subset["Mz"].values, subset["Intensity"].values))
        else:
            # Individual ion events: bin by m/z to derive peak intensities
            peaks = bin_ions(subset["Mz"].values, bin_width=bin_width)

        clusters = find_isotopic_clusters(
            peaks,
            charge,
            max_gap=max_gap,
            ppm_tolerance=ppm_tolerance,
        )

        for cluster in clusters:
            n_obs = len(cluster.peaks)
            if n_obs < min_cluster_peaks and not cluster.gap_indices:
                continue

            classification = classify(cluster, npeaks=npeaks)
            recon_peaks = reconstruct(cluster, classification, npeaks=npeaks)

            # --- Observed peaks -------------------------------------------
            for peak in cluster.peaks:
                records.append(
                    {
                        "Mz": peak.mz,
                        "Charge": peak.charge,
                        "Intensity": peak.intensity,
                        "NeutralMass": cluster.neutral_mass,
                        "IsotopicIndex": peak.isotope_index,
                        "IsEstimated": False,
                        "EstimatedIntensity": None,
                        "Uncertainty": None,
                        "Confidence": classification.confidence,
                        "Hypothesis": classification.hypothesis,
                        "FitRSquared": classification.fit_r_squared,
                        "GapAtApex": classification.gap_at_apex,
                        "Notes": classification.notes,
                    }
                )

            # --- Reconstructed (rescued) peaks ----------------------------
            for rec in recon_peaks:
                records.append(
                    {
                        "Mz": None,
                        "Charge": rec.charge,
                        "Intensity": None,
                        "NeutralMass": cluster.neutral_mass,
                        "IsotopicIndex": rec.isotope_index,
                        "IsEstimated": True,
                        "EstimatedIntensity": rec.estimated_intensity,
                        "Uncertainty": rec.uncertainty if rec.uncertainty != float("inf") else None,
                        "Confidence": rec.confidence,
                        "Hypothesis": classification.hypothesis,
                        "FitRSquared": classification.fit_r_squared,
                        "GapAtApex": classification.gap_at_apex,
                        "Notes": classification.notes,
                    }
                )

    return pd.DataFrame(records)
