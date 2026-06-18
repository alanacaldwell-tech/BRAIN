"""
Tests for brainpy.mie — Multi-Ion Event correction pipeline.

Synthetic data is generated directly rather than reading .dmt files so that
the tests run without external fixtures.
"""
import sqlite3
import tempfile
import os
import math

import pytest

from brainpy.mie import (
    bin_ions,
    find_isotopic_clusters,
    averagine_composition,
    theoretical_envelope,
    classify,
    reconstruct,
    process_dmt,
    ObservedPeak,
    IsotopicCluster,
)
from brainpy import PROTON

ISOTOPE_SPACING = 1.003355


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_dmt(ions, tmp_dir):
    """Write a minimal .dmt SQLite file with the given (Mz, Charge) rows."""
    path = os.path.join(tmp_dir, "test.dmt")
    conn = sqlite3.connect(path)
    conn.execute("CREATE TABLE Ion (Mz REAL, Charge INTEGER)")
    conn.executemany("INSERT INTO Ion VALUES (?, ?)", ions)
    conn.commit()
    conn.close()
    return path


def _cluster_from_peaks(mz_list, charge, gap_indices=None):
    """Build an IsotopicCluster directly from a list of m/z values."""
    neutral_mass = mz_list[0] * charge - charge * PROTON
    peaks = [
        ObservedPeak(mz=mz, intensity=1.0, charge=charge, isotope_index=i)
        for i, mz in enumerate(mz_list)
    ]
    return IsotopicCluster(
        peaks=peaks,
        charge=charge,
        gap_indices=gap_indices or [],
        neutral_mass=neutral_mass,
    )


# ---------------------------------------------------------------------------
# bin_ions
# ---------------------------------------------------------------------------

class TestBinIons:
    def test_groups_nearby_values(self):
        # Values chosen so they land cleanly inside two 0.01-Da bins
        # 100.001–100.003 → all round to 100.00; 101.001–101.002 → both round to 101.00
        mz_values = [100.001, 100.002, 100.003, 101.001, 101.002]
        result = bin_ions(mz_values, bin_width=0.01)
        assert len(result) == 2
        counts = [r[1] for r in result]
        assert counts == [3.0, 2.0]

    def test_empty_input(self):
        assert bin_ions([], bin_width=0.005) == []

    def test_single_ion(self):
        result = bin_ions([500.123], bin_width=0.005)
        assert len(result) == 1
        assert result[0][1] == 1.0


# ---------------------------------------------------------------------------
# find_isotopic_clusters
# ---------------------------------------------------------------------------

class TestFindIsotopicClusters:
    def _make_peaks(self, start_mz, charge, n_peaks, missing=None):
        """Generate a clean isotopic series, optionally with gaps."""
        delta = ISOTOPE_SPACING / charge
        all_peaks = [(start_mz + i * delta, float(n_peaks - i)) for i in range(n_peaks)]
        if missing:
            all_peaks = [p for i, p in enumerate(all_peaks) if i not in missing]
        return all_peaks

    def test_clean_cluster_detected(self):
        peaks = self._make_peaks(500.0, charge=2, n_peaks=5)
        clusters = find_isotopic_clusters(peaks, charge=2)
        assert len(clusters) == 1
        assert len(clusters[0].peaks) == 5
        assert clusters[0].gap_indices == []

    def test_gap_detected(self):
        # Remove the middle peak (index 2)
        peaks = self._make_peaks(500.0, charge=2, n_peaks=5, missing={2})
        clusters = find_isotopic_clusters(peaks, charge=2, max_gap=3)
        assert len(clusters) == 1
        assert 2 in clusters[0].gap_indices

    def test_two_separate_clusters(self):
        # Two clusters far apart in m/z
        c1 = self._make_peaks(500.0, charge=1, n_peaks=4)
        c2 = self._make_peaks(800.0, charge=1, n_peaks=4)
        clusters = find_isotopic_clusters(c1 + c2, charge=1)
        assert len(clusters) == 2

    def test_empty_input(self):
        assert find_isotopic_clusters([], charge=1) == []

    def test_invalid_charge_skipped(self):
        peaks = [(500.0, 10.0), (501.003, 8.0)]
        assert find_isotopic_clusters(peaks, charge=0) == []


# ---------------------------------------------------------------------------
# averagine_composition / theoretical_envelope
# ---------------------------------------------------------------------------

class TestAveragine:
    def test_composition_positive_counts(self):
        comp = averagine_composition(1000.0)
        for el, count in comp.items():
            assert count >= 1, f"Element {el} has non-positive count {count}"

    def test_scales_with_mass(self):
        comp_small = averagine_composition(500.0)
        comp_large = averagine_composition(2000.0)
        assert comp_large["C"] > comp_small["C"]

    def test_envelope_normalised(self):
        peaks = theoretical_envelope(1000.0, charge=2)
        assert peaks, "Should return non-empty peak list"
        max_intensity = max(p.intensity for p in peaks)
        assert abs(max_intensity - 1.0) < 1e-9

    def test_envelope_mz_ordered(self):
        peaks = theoretical_envelope(1000.0, charge=2)
        mz_values = [p.mz for p in peaks]
        assert mz_values == sorted(mz_values)


# ---------------------------------------------------------------------------
# classify
# ---------------------------------------------------------------------------

class TestClassify:
    def _make_mie_cluster(self, neutral_mass=5000.0, charge=3):
        """
        Simulate a MIE artifact using BRAIN's own theoretical intensities so
        the single-species fit is perfect and the gap is at the true apex.

        A large molecule (default 5 kDa) is used so the apex sits well inside
        the isotopic envelope (not at index 0), making it detectable as a gap.
        """
        theoretical = theoretical_envelope(neutral_mass, charge, npeaks=20)
        if not theoretical:
            pytest.skip("Could not generate theoretical envelope")

        apex_idx = max(range(len(theoretical)), key=lambda i: theoretical[i].intensity)

        peaks = []
        gap_indices = []
        for i, theo_peak in enumerate(theoretical):
            if i == apex_idx:
                gap_indices.append(i)
            else:
                peaks.append(
                    ObservedPeak(
                        mz=theo_peak.mz,
                        intensity=theo_peak.intensity,
                        charge=charge,
                        isotope_index=i,
                    )
                )

        # neutral_mass from the lowest-m/z observed peak (as the pipeline does)
        first_mz = min(p.mz for p in peaks)
        computed_mass = first_mz * charge - charge * PROTON

        return IsotopicCluster(
            peaks=peaks,
            charge=charge,
            gap_indices=gap_indices,
            neutral_mass=computed_mass,
        )

    def test_no_gap_returns_no_gap(self):
        cluster = _cluster_from_peaks(
            [500.0, 500.0 + ISOTOPE_SPACING, 500.0 + 2 * ISOTOPE_SPACING],
            charge=1,
        )
        result = classify(cluster)
        assert result.hypothesis == "no_gap"
        assert result.confidence == 1.0

    def test_mie_cluster_classified(self):
        cluster = self._make_mie_cluster()
        result = classify(cluster)
        # Should be mie or at worst ambiguous — not overlap
        assert result.hypothesis in ("mie", "ambiguous"), (
            f"Expected mie/ambiguous, got {result.hypothesis} "
            f"(R²={result.fit_r_squared:.2f}, gap_at_apex={result.gap_at_apex})"
        )

    def test_confidence_between_0_and_1(self):
        cluster = self._make_mie_cluster()
        result = classify(cluster)
        assert 0.0 <= result.confidence <= 1.0

    def test_consecutive_gaps_noted(self):
        # Two consecutive gaps → note should appear
        delta = ISOTOPE_SPACING / 2
        start = 700.0
        neutral_mass = start * 2 - 2 * PROTON
        peaks = [
            ObservedPeak(mz=start, intensity=0.3, charge=2, isotope_index=0),
            ObservedPeak(mz=start + 3 * delta, intensity=0.3, charge=2, isotope_index=3),
        ]
        cluster = IsotopicCluster(
            peaks=peaks, charge=2, gap_indices=[1, 2], neutral_mass=neutral_mass
        )
        result = classify(cluster)
        assert result.n_consecutive_gaps == 2


# ---------------------------------------------------------------------------
# reconstruct
# ---------------------------------------------------------------------------

class TestReconstruct:
    def test_no_reconstruction_for_overlap(self):
        from brainpy.mie.peaks import ClassificationResult
        cluster = _cluster_from_peaks(
            [500.0, 500.0 + ISOTOPE_SPACING], charge=1, gap_indices=[1]
        )
        classification = ClassificationResult(
            hypothesis="overlap",
            confidence=0.9,
            fit_r_squared=0.4,
            gap_at_apex=False,
            n_consecutive_gaps=1,
        )
        assert reconstruct(cluster, classification) == []

    def test_estimated_intensity_positive(self):
        from brainpy.mie.peaks import ClassificationResult
        delta = ISOTOPE_SPACING / 2
        start = 600.0
        neutral_mass = start * 2 - 2 * PROTON
        # 5-peak cluster with middle peak (index 2) missing
        peaks = [
            ObservedPeak(mz=start + i * delta, intensity=float(3 - abs(i - 2.5)),
                         charge=2, isotope_index=i)
            for i in [0, 1, 3, 4]
        ]
        cluster = IsotopicCluster(
            peaks=peaks, charge=2, gap_indices=[2], neutral_mass=neutral_mass
        )
        classification = ClassificationResult(
            hypothesis="mie",
            confidence=0.85,
            fit_r_squared=0.90,
            gap_at_apex=True,
            n_consecutive_gaps=1,
        )
        result = reconstruct(cluster, classification)
        assert len(result) == 1
        assert result[0].estimated_intensity > 0
        assert result[0].confidence == 0.85

    def test_no_gaps_returns_empty(self):
        from brainpy.mie.peaks import ClassificationResult
        cluster = _cluster_from_peaks([500.0, 501.003, 502.006], charge=1)
        classification = ClassificationResult(
            hypothesis="no_gap",
            confidence=1.0,
            fit_r_squared=1.0,
            gap_at_apex=False,
            n_consecutive_gaps=0,
        )
        assert reconstruct(cluster, classification) == []


# ---------------------------------------------------------------------------
# process_dmt (integration)
# ---------------------------------------------------------------------------

class TestProcessDmt:
    def test_reads_and_returns_dataframe(self, tmp_path):
        # Build a simple 3-peak z=1 series
        delta = ISOTOPE_SPACING
        start = 500.0
        ions = [(start + i * delta, 1) for i in range(3) for _ in range(10)]
        path = _make_dmt(ions, str(tmp_path))

        df = process_dmt(path, charge_range=(1, 1))
        assert not df.empty
        assert "IsEstimated" in df.columns
        assert "Hypothesis" in df.columns

    def test_mie_peak_rescued(self, tmp_path):
        """
        Simulate a 5 kDa peptide at charge 3 whose apex isotope peak was
        discarded by the instrument due to multi-ion coincidences.

        For molecules of this size the BRAIN averagine apex falls several
        isotope steps in from the monoisotopic peak, so the missing peak is
        internal to the cluster and clearly detectable as a gap.
        """
        from brainpy.mie.averagine import theoretical_envelope as theo_env

        neutral_mass = 5000.0
        charge = 3
        theoretical = theo_env(neutral_mass, charge, npeaks=20)
        if not theoretical:
            pytest.skip("No theoretical peaks generated")

        apex_idx = max(range(len(theoretical)), key=lambda i: theoretical[i].intensity)

        # Emit 100 ions at every isotope position except the apex
        ions = []
        for i, peak in enumerate(theoretical):
            if i == apex_idx:
                continue
            for _ in range(100):
                ions.append((peak.mz, charge))

        path = _make_dmt(ions, str(tmp_path))
        df = process_dmt(path, charge_range=(charge, charge), max_gap=3, min_cluster_peaks=2)
        estimated = df[df["IsEstimated"] == True]
        assert len(estimated) >= 1, (
            f"Expected at least one rescued MIE peak "
            f"(apex_idx={apex_idx}, hypotheses seen: {df['Hypothesis'].unique()})"
        )

    def test_empty_table_returns_empty_df(self, tmp_path):
        path = _make_dmt([], str(tmp_path))
        df = process_dmt(path)
        assert df.empty

    def test_required_columns_present(self, tmp_path):
        delta = ISOTOPE_SPACING
        ions = [(500.0 + i * delta, 1) for i in range(4) for _ in range(20)]
        path = _make_dmt(ions, str(tmp_path))
        df = process_dmt(path, charge_range=(1, 1))
        expected_cols = {
            "Mz", "Charge", "Intensity", "NeutralMass", "IsotopicIndex",
            "IsEstimated", "EstimatedIntensity", "Uncertainty", "Confidence",
            "Hypothesis", "FitRSquared", "GapAtApex", "Notes",
        }
        assert expected_cols.issubset(set(df.columns))
