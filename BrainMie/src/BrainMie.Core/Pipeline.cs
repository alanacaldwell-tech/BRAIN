namespace BrainMie.Core;

using BrainMie.Core.Data;

/// <summary>
/// End-to-end pipeline: .dmt file → corrected isotopic peak list.
/// </summary>
public static class Pipeline
{
    /// <summary>
    /// Read a .dmt SQLite file, convert all ions to neutral mass, detect MIE artifacts
    /// across all charge states in one unified mass spectrum, and return a flat list of
    /// processing rows containing observed and suppression-corrected peaks.
    /// </summary>
    /// <param name="filePath">Path to the .dmt file.</param>
    /// <param name="binWidth">Mass bin width in Da for ion-event aggregation. Default 0.02 Da.</param>
    /// <param name="maxGap">Max consecutive missing isotope peaks per cluster. Default 3.</param>
    /// <param name="massMatchTol">
    /// Absolute mass tolerance in Da for matching peaks to isotope positions (replaces
    /// the former per-charge ppm parameter). Default 0.2 Da.
    /// </param>
    /// <param name="minClusterPeaks">Minimum observed peaks to retain any cluster. Default 8.</param>
    /// <param name="nPeaks">Number of theoretical isotope peaks to generate. Default 25.</param>
    /// <param name="suppressionThreshold">
    /// Fraction of fitted theoretical intensity below which an observed peak is
    /// considered MIE-suppressed. Default 0.5 (50%).
    /// </param>
    /// <param name="minFitR2">
    /// Minimum R² of the single-species averagine fit required to keep a cluster.
    /// Default 0.5.
    /// </param>
    /// <returns>
    /// A tuple of the corrected row list and the original raw ion list (needed by
    /// <see cref="ProteoformAggregator"/> to derive charge-state annotations).
    /// </returns>
    public static (List<ProcessingRow> Rows, List<(double Mz, int Charge, double? Intensity)> Ions)
        ProcessDmt(
            string filePath,
            double binWidth             = 0.02,
            int    maxGap               = 3,
            double massMatchTol         = 0.2,
            int    minClusterPeaks      = 8,
            int    nPeaks               = 25,
            double suppressionThreshold = 0.5,
            double minFitR2             = 0.5)
    {
        var ions = DmtReader.ReadIons(filePath);

        // Convert every ion to neutral mass immediately.
        // mass = mz × charge − charge × HMass
        List<(double Mass, double Intensity)> peaks;
        if (ions.Count > 0 && ions[0].Intensity.HasValue)
        {
            // File already stores per-peak intensities — use them directly.
            peaks = ions
                .Select(i => (Mass: i.Mz * i.Charge - i.Charge * Constants.HMass, Intensity: i.Intensity!.Value))
                .OrderBy(p => p.Mass)
                .ToList();
            peaks = EnvelopeDetector.PickLocalMaxima(peaks);
        }
        else
        {
            var massValues = ions.Select(i => i.Mz * i.Charge - i.Charge * Constants.HMass);
            peaks = EnvelopeDetector.BinIons(massValues, binWidth);
            peaks = EnvelopeDetector.PickLocalMaxima(peaks);
        }

        var clusters = EnvelopeDetector.FindClusters(peaks, maxGap, massMatchTol);

        var results = new List<ProcessingRow>();

        foreach (var cluster in clusters)
        {
            // Minimum isotope count
            if (cluster.Peaks.Count < minClusterPeaks) continue;

            // Classify first — needed for both R² filter and theory-based apex index.
            var classification = Classifier.Classify(cluster, nPeaks);
            if (classification.FitRSquared < minFitR2) continue;

            // Use the theoretical averagine apex (adjusted by fit offset) for the
            // completeness check. This is robust to MIE suppression of the centroid
            // peak, which would shift the observed intensity maximum to a flanking peak.
            var theoretical  = IsotopeDistribution.ComputeMassEnvelope(cluster.NeutralMass, nPeaks);
            if (theoretical.Count == 0) continue;
            int theoApexTi   = theoretical.Select((t, i) => (t.Intensity, i)).MaxBy(x => x.Intensity)!.i;
            int apexIdx      = theoApexTi - classification.FitOffset;
            var obsSet       = cluster.Peaks.Select(p => p.IsotopeIndex).ToHashSet();
            bool gapNearApex = Enumerable.Range(apexIdx - 2, 5)
                                   .Where(i => i >= 0)
                                   .Any(i => !obsSet.Contains(i));
            if (gapNearApex) continue;

            var suppressed    = Classifier.DetectSuppressed(cluster, suppressionThreshold, nPeaks);
            var suppressedMap = suppressed.ToDictionary(s => s.IsotopeIndex, s => s.CorrectedIntensity);

            double clusterIonCount   = cluster.Peaks.Sum(p => p.Intensity);
            double suppressionDelta  = suppressed.Sum(s => s.CorrectedIntensity - s.ObservedIntensity);
            double correctedIonCount = clusterIonCount + suppressionDelta;

            foreach (var peak in cluster.Peaks)
            {
                bool    isSuppressed       = suppressedMap.TryGetValue(peak.IsotopeIndex, out double ci);
                double? correctedIntensity = isSuppressed ? ci : null;

                results.Add(new ProcessingRow(
                    IsotopeMass:        peak.Mass,
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

        return (results, ions);
    }
}
