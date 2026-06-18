namespace BrainMie.Core;

using BrainMie.Core.Data;

/// <summary>
/// Classifies gapped isotopic clusters as MIE artifact, true spectral overlap, or ambiguous.
///
/// Scoring logic
/// -------------
/// Three features contribute to a continuous mie_score in [−0.15, 0.85]:
///
///   +0.50  gap_at_apex        Gap sits at the theoretical intensity maximum (strongest signal).
///   +0.05–0.35  fit_r²        Single-species averagine fit quality over surviving peaks.
///   −0.05 × n_consec          Penalty for each consecutive missing peak (cap −0.15).
///
///   mie_score ≥ 0.65  →  "mie"
///   mie_score ≤ 0.30  →  "overlap"
///   otherwise         →  "ambiguous"
/// </summary>
public static class Classifier
{
    /// <summary>
    /// Classify <paramref name="cluster"/> and return a <see cref="ClassificationResult"/>.
    /// </summary>
    public static ClassificationResult Classify(IsotopicCluster cluster, int nPeaks = 25)
    {
        var theoretical = IsotopeDistribution.ComputeEnvelope(
            cluster.NeutralMass, cluster.Charge, nPeaks);

        if (theoretical.Count == 0)
            return new ClassificationResult(
                "ambiguous", 0.0, 0.0, false,
                MaxConsecutiveGaps(cluster.GapIndices),
                "Could not generate theoretical envelope for this mass/charge.");

        var (offset, _, r2) = FitAlignment(cluster.Peaks, theoretical);

        if (cluster.GapIndices.Count == 0)
            return new ClassificationResult("no_gap", 1.0, r2, false, 0);
        bool apexInGap      = ApexInGaps(cluster.GapIndices, theoretical, offset);
        int  nConsec        = MaxConsecutiveGaps(cluster.GapIndices);

        // ---- Scoring ---------------------------------------------------------
        double mieScore = 0.0;
        if (apexInGap)   mieScore += 0.50;
        if (r2 >= 0.85)  mieScore += 0.35;
        else if (r2 >= 0.70) mieScore += 0.20;
        else if (r2 >= 0.50) mieScore += 0.05;
        mieScore -= Math.Min(0.15, nConsec * 0.05);

        double reconConfidence = Math.Max(0.0, 1.0 - (nConsec - 1) * 0.25);

        // ---- Decision --------------------------------------------------------
        var notes = new List<string>();
        string hypothesis;
        double confidence;

        if (mieScore >= 0.65)
        {
            hypothesis = "mie";
            confidence = Math.Min(1.0, mieScore) * reconConfidence;
        }
        else if (mieScore <= 0.30)
        {
            hypothesis = "overlap";
            confidence = Math.Min(1.0, Math.Max(0.0, 1.0 - mieScore));
        }
        else
        {
            hypothesis = "ambiguous";
            confidence = reconConfidence * (0.5 + Math.Abs(mieScore - 0.475));
            notes.Add("borderline score — treat reconstruction with caution");
        }

        if (nConsec >= 2)
            notes.Add($"{nConsec} consecutive missing peaks: reconstruction uncertainty is elevated");
        if (r2 < 0.50)
            notes.Add($"poor single-species fit (R²={r2:F2})");
        if (!apexInGap && hypothesis == "mie")
            notes.Add("gap is not at the apex — verify MIE assignment");

        return new ClassificationResult(
            hypothesis, confidence, r2, apexInGap, nConsec,
            string.Join("; ", notes));
    }

    // ---------------------------------------------------------------------------
    // Suppression detection
    // ---------------------------------------------------------------------------

    /// <summary>
    /// Identify observed peaks whose intensity is below
    /// <paramref name="suppressionThreshold"/> × the fitted theoretical expectation.
    /// Uses a two-pass robust fit: first pass includes all peaks; suppressed peaks are
    /// excluded from the second pass so they don't drag the scale factor down.
    /// </summary>
    /// <returns>
    /// List of (IsotopeIndex, ObservedIntensity, CorrectedIntensity) for every peak
    /// judged to be MIE-suppressed. Empty when fewer than 2 non-suppressed peaks
    /// remain for the second-pass fit.
    /// </returns>
    public static List<(int IsotopeIndex, double ObservedIntensity, double CorrectedIntensity)>
        DetectSuppressed(
            IsotopicCluster cluster,
            double suppressionThreshold = 0.5,
            int nPeaks = 25)
    {
        var theoretical = IsotopeDistribution.ComputeEnvelope(
            cluster.NeutralMass, cluster.Charge, nPeaks);
        if (theoretical.Count == 0) return [];

        // First pass: fit with all peaks
        var (offset1, scale1, _) = FitAlignment(cluster.Peaks, theoretical);
        if (scale1 <= 0) return [];

        var suppressedIdx = cluster.Peaks
            .Where(p =>
            {
                int ti = p.IsotopeIndex + offset1;
                if (ti < 0 || ti >= theoretical.Count) return false;
                double expected = scale1 * theoretical[ti].Intensity;
                return expected > 0 && p.Intensity / expected < suppressionThreshold;
            })
            .Select(p => p.IsotopeIndex)
            .ToHashSet();

        if (suppressedIdx.Count == 0) return [];

        // Second pass: refit without suppressed peaks for an unbiased scale estimate
        var unsuppressed = cluster.Peaks
            .Where(p => !suppressedIdx.Contains(p.IsotopeIndex))
            .ToList();
        if (unsuppressed.Count < 2) return [];

        var (offset2, scale2, _) = FitAlignment(unsuppressed, theoretical);

        return cluster.Peaks
            .Where(p => suppressedIdx.Contains(p.IsotopeIndex))
            .Select(p =>
            {
                int ti = p.IsotopeIndex + offset2;
                double corrected = (ti >= 0 && ti < theoretical.Count)
                    ? Math.Max(0.0, scale2 * theoretical[ti].Intensity)
                    : p.Intensity;
                return (p.IsotopeIndex, p.Intensity, corrected);
            })
            .ToList();
    }

    // ---------------------------------------------------------------------------
    // Internal helpers (internal so Reconstructor can reuse FitAlignment)
    // ---------------------------------------------------------------------------

    /// <summary>
    /// Find the best alignment of observed peaks to the theoretical series and fit
    /// a single scale factor s such that observed ≈ s × theoretical.
    /// </summary>
    /// <returns>(theoOffset, scale, R²) where theoOffset is the theoretical index
    /// that corresponds to cluster isotope index 0.</returns>
    internal static (int Offset, double Scale, double R2) FitAlignment(
        List<ObservedPeak> observed,
        List<(double Mz, double Intensity)> theoretical)
    {
        var obs    = observed.ToDictionary(p => p.IsotopeIndex, p => p.Intensity);
        int nTheo  = theoretical.Count;

        double bestR2    = double.NegativeInfinity;
        int    bestOff   = 0;
        double bestScale = 1.0;

        for (int offset = 0; offset < Math.Min(nTheo, 15); offset++)
        {
            var pairs = obs
                .Where(kvp => kvp.Key + offset < nTheo)
                .Select(kvp => (Obs: kvp.Value, Theo: theoretical[kvp.Key + offset].Intensity))
                .ToList();

            if (pairs.Count < 2) continue;

            double denom = pairs.Sum(p => p.Theo * p.Theo);
            if (denom == 0) continue;
            double scale = pairs.Sum(p => p.Obs * p.Theo) / denom;

            double meanObs = pairs.Average(p => p.Obs);
            double ssRes   = pairs.Sum(p => Math.Pow(p.Obs - scale * p.Theo, 2));
            double ssTot   = pairs.Sum(p => Math.Pow(p.Obs - meanObs, 2));
            double r2      = ssTot > 0 ? 1.0 - ssRes / ssTot : 1.0;

            if (r2 > bestR2) { bestR2 = r2; bestOff = offset; bestScale = scale; }
        }

        return (bestOff, bestScale, bestR2 == double.NegativeInfinity ? 0.0 : bestR2);
    }

    private static bool ApexInGaps(
        List<int> gapIndices,
        List<(double Mz, double Intensity)> theoretical,
        int offset)
    {
        if (gapIndices.Count == 0 || theoretical.Count == 0) return false;

        int apexTheoIdx = theoretical
            .Select((p, i) => (p.Intensity, Index: i))
            .MaxBy(t => t.Intensity)
            .Index;

        return gapIndices.Contains(apexTheoIdx - offset);
    }

    private static int MaxConsecutiveGaps(List<int> gaps)
    {
        if (gaps.Count == 0) return 0;
        var sorted     = gaps.Order().ToList();
        int maxRun = 1, run = 1;
        for (int i = 1; i < sorted.Count; i++)
        {
            run = sorted[i] == sorted[i - 1] + 1 ? run + 1 : 1;
            if (run > maxRun) maxRun = run;
        }
        return maxRun;
    }
}
