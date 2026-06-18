namespace BrainMie.Core;

using BrainMie.Core.Data;

/// <summary>
/// End-to-end pipeline: .dmt file → corrected isotopic peak list.
/// </summary>
public static class Pipeline
{
    /// <summary>
    /// Read a .dmt SQLite file, detect MIE artifacts across all charge states in
    /// <paramref name="chargeRange"/>, and return a flat list of processing rows
    /// containing both observed and reconstructed (rescued) peaks.
    /// </summary>
    /// <param name="filePath">Path to the .dmt file.</param>
    /// <param name="chargeRange">Inclusive (min, max) charge state range. Default (1, 10).</param>
    /// <param name="binWidth">m/z bin width in Da for ion-event aggregation. Default 0.02 Da.</param>
    /// <param name="maxGap">Max consecutive missing isotope peaks per cluster. Default 3.</param>
    /// <param name="ppmTolerance">m/z matching tolerance in ppm. Default 10 ppm.</param>
    /// <param name="minClusterPeaks">Minimum observed peaks to retain any cluster. Default 8.</param>
    /// <param name="nPeaks">Number of theoretical isotope peaks to generate. Default 25.</param>
    /// <param name="suppressionThreshold">
    /// Fraction of fitted theoretical intensity below which an observed peak is
    /// considered MIE-suppressed (present but anomalously low). Default 0.5 (50%).
    /// </param>
    /// <param name="minFitR2">
    /// Minimum R² of the single-species averagine fit required to keep a cluster.
    /// Clusters below this threshold are unlikely to be single-species proteoforms.
    /// Default 0.5.
    /// </param>
    public static List<ProcessingRow> ProcessDmt(
        string             filePath,
        (int Min, int Max) chargeRange          = default,
        double             binWidth             = 0.02,
        int                maxGap               = 3,
        double             ppmTolerance         = 10.0,
        int                minClusterPeaks      = 8,
        int                nPeaks               = 25,
        double             suppressionThreshold = 0.5,
        double             minFitR2             = 0.5)
    {
        if (chargeRange == default) chargeRange = (1, 10);

        var ions    = DmtReader.ReadIons(filePath);
        var results = new List<ProcessingRow>();

        for (int charge = chargeRange.Min; charge <= chargeRange.Max; charge++)
        {
            var chargeIons = ions.Where(i => i.Charge == charge).ToList();
            if (chargeIons.Count == 0) continue;

            List<(double Mz, double Intensity)> peaks;
            if (chargeIons[0].Intensity.HasValue)
                peaks = chargeIons.Select(i => (i.Mz, i.Intensity!.Value)).ToList();
            else
            {
                peaks = EnvelopeDetector.BinIons(chargeIons.Select(i => i.Mz), binWidth);
                peaks = EnvelopeDetector.PickLocalMaxima(peaks);
            }

            var clusters = EnvelopeDetector.FindClusters(
                peaks, charge, maxGap, ppmTolerance);

            foreach (var cluster in clusters)
            {
                // Require minimum isotope count for all clusters, gapped or not.
                if (cluster.Peaks.Count < minClusterPeaks) continue;

                // Require complete isotope coverage within ±2 of the apex.
                // A gap this close to the most intense peak makes correction unreliable.
                int        apexIdx     = cluster.Peaks.MaxBy(p => p.Intensity)!.IsotopeIndex;
                var        obsSet      = cluster.Peaks.Select(p => p.IsotopeIndex).ToHashSet();
                bool       gapNearApex = Enumerable.Range(apexIdx - 2, 5)
                                             .Where(i => i >= 0)
                                             .Any(i => !obsSet.Contains(i));
                if (gapNearApex) continue;

                var classification = Classifier.Classify(cluster, nPeaks);

                // Require a minimum averagine fit quality; poor fits indicate noise or
                // multi-species overlap rather than a single proteoform.
                if (classification.FitRSquared < minFitR2) continue;
                var suppressed     = Classifier.DetectSuppressed(
                                         cluster, suppressionThreshold, nPeaks);

                var suppressedMap  = suppressed.ToDictionary(
                                         s => s.IsotopeIndex,
                                         s => s.CorrectedIntensity);

                // CorrectedIonCount: replace suppressed-peak observed intensities with
                // their fitted theoretical values; completely absent peaks are not rescued.
                double clusterIonCount   = cluster.Peaks.Sum(p => p.Intensity);
                double suppressionDelta  = suppressed.Sum(
                    s => s.CorrectedIntensity - s.ObservedIntensity);
                double correctedIonCount = clusterIonCount + suppressionDelta;

                // Observed peaks
                foreach (var peak in cluster.Peaks)
                {
                    bool   isSuppressed      = suppressedMap.TryGetValue(peak.IsotopeIndex, out double ci);
                    double? correctedIntensity = isSuppressed ? ci : null;

                    results.Add(new ProcessingRow(
                        Mz:                 peak.Mz,
                        Charge:             peak.Charge,
                        Intensity:          peak.Intensity,
                        NeutralMass:        cluster.NeutralMass,
                        IsotopicIndex:      peak.IsotopeIndex,
                        IsEstimated:        false,
                        IsSuppressed:       isSuppressed,
                        EstimatedIntensity: null,
                        CorrectedIntensity: correctedIntensity,
                        Uncertainty:        null,
                        Confidence:         classification.Confidence,
                        Hypothesis:         classification.Hypothesis,
                        FitRSquared:        classification.FitRSquared,
                        GapAtApex:          classification.GapAtApex,
                        Notes:              classification.Notes,
                        ClusterIonCount:    clusterIonCount,
                        CorrectedIonCount:  correctedIonCount));
                }
            }
        }

        return results;
    }
}
