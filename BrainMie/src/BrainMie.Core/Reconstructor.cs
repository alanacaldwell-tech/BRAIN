namespace BrainMie.Core;

using BrainMie.Core.Data;

/// <summary>
/// Reconstructs the expected intensity of MIE-discarded isotope peaks.
///
/// For each gap the reconstructor:
///   1. Fits a scale factor s over surviving observed peaks  (s × theoretical ≈ observed).
///   2. Predicts the missing intensity as  s × theoretical[gap].
///   3. Estimates a 1-sigma prediction interval from the fit residuals,
///      widened by a penalty that grows with the number of consecutive gaps.
///
/// "overlap" clusters are returned as-is (empty list); "ambiguous" clusters
/// are still reconstructed but carry a reduced confidence score.
/// </summary>
public static class Reconstructor
{
    /// <summary>
    /// Estimate intensities for every gap in <paramref name="cluster"/>.
    /// Returns an empty list when the cluster is classified as overlap or has no gaps.
    /// </summary>
    public static List<ReconstructedPeak> Reconstruct(
        IsotopicCluster     cluster,
        ClassificationResult classification,
        int                  nPeaks = 25)
    {
        if (classification.Hypothesis == "overlap" || cluster.GapIndices.Count == 0)
            return [];

        var theoretical = IsotopeDistribution.ComputeEnvelope(
            cluster.NeutralMass, cluster.Charge, nPeaks);
        if (theoretical.Count == 0) return [];

        var (offset, scale, _) = Classifier.FitAlignment(cluster.Peaks, theoretical);

        // Residuals over the observed (non-gap) peaks
        var obs = cluster.Peaks.ToDictionary(p => p.IsotopeIndex, p => p.Intensity);
        var residuals = obs
            .Where(kvp => kvp.Key + offset < theoretical.Count)
            .Select(kvp => kvp.Value - scale * theoretical[kvp.Key + offset].Intensity)
            .ToList();

        double residualStd      = residuals.Count >= 2 ? StdDev(residuals) : 0.0;
        int    nObs             = Math.Max(residuals.Count, 1);
        double predIntervalBase = residualStd * Math.Sqrt(1.0 + 1.0 / nObs);
        double consecPenalty    = 1.0 + 0.5 * Math.Max(0, classification.NConsecutiveGaps - 1);

        double delta   = Constants.IsotopeSpacing / cluster.Charge;
        double firstMz = cluster.Peaks.Min(p => p.Mz);

        var result = new List<ReconstructedPeak>();
        foreach (int gapIdx in cluster.GapIndices)
        {
            double missingMz = firstMz + gapIdx * delta;
            int    theoIdx   = gapIdx + offset;

            if (theoIdx >= 0 && theoIdx < theoretical.Count)
            {
                double estimated   = Math.Max(0.0, scale * theoretical[theoIdx].Intensity);
                double uncertainty = predIntervalBase * consecPenalty;
                result.Add(new ReconstructedPeak(
                    missingMz, cluster.Charge, estimated, uncertainty,
                    gapIdx, classification.Confidence));
            }
            else
            {
                // Gap falls outside theoretical range — no meaningful estimate
                result.Add(new ReconstructedPeak(
                    missingMz, cluster.Charge, 0.0, double.PositiveInfinity,
                    gapIdx, 0.0));
            }
        }
        return result;
    }

    private static double StdDev(List<double> values)
    {
        double mean = values.Average();
        return Math.Sqrt(values.Sum(v => Math.Pow(v - mean, 2)) / values.Count);
    }
}
