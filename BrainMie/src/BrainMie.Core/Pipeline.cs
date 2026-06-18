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
    /// <param name="binWidth">m/z bin width in Da for ion-event aggregation. Default 0.005 Da.</param>
    /// <param name="maxGap">Max consecutive missing isotope peaks per cluster. Default 3.</param>
    /// <param name="ppmTolerance">m/z matching tolerance in ppm. Default 20 ppm.</param>
    /// <param name="minClusterPeaks">Minimum observed peaks required to process a gap-free cluster. Default 3.</param>
    /// <param name="nPeaks">Number of theoretical isotope peaks to generate. Default 25.</param>
    public static List<ProcessingRow> ProcessDmt(
        string       filePath,
        (int Min, int Max) chargeRange   = default,
        double       binWidth            = 0.005,
        int          maxGap              = 3,
        double       ppmTolerance        = 20.0,
        int          minClusterPeaks     = 3,
        int          nPeaks              = 25)
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
                peaks = EnvelopeDetector.BinIons(chargeIons.Select(i => i.Mz), binWidth);

            var clusters = EnvelopeDetector.FindClusters(
                peaks, charge, maxGap, ppmTolerance);

            foreach (var cluster in clusters)
            {
                if (cluster.Peaks.Count < minClusterPeaks && cluster.GapIndices.Count == 0)
                    continue;

                var classification = Classifier.Classify(cluster, nPeaks);
                var reconPeaks     = Reconstructor.Reconstruct(cluster, classification, nPeaks);

                double clusterIonCount   = cluster.Peaks.Sum(p => p.Intensity);
                double correctedIonCount = clusterIonCount + reconPeaks.Sum(r => r.EstimatedIntensity);

                // Observed peaks
                foreach (var peak in cluster.Peaks)
                {
                    results.Add(new ProcessingRow(
                        peak.Mz, peak.Charge, peak.Intensity,
                        cluster.NeutralMass, peak.IsotopeIndex,
                        IsEstimated:        false,
                        EstimatedIntensity: null,
                        Uncertainty:        null,
                        Confidence:         classification.Confidence,
                        Hypothesis:         classification.Hypothesis,
                        FitRSquared:        classification.FitRSquared,
                        GapAtApex:          classification.GapAtApex,
                        Notes:              classification.Notes,
                        ClusterIonCount:    clusterIonCount,
                        CorrectedIonCount:  correctedIonCount));
                }

                // Reconstructed (rescued) peaks
                foreach (var rec in reconPeaks)
                {
                    results.Add(new ProcessingRow(
                        Mz:                 null,
                        Charge:             rec.Charge,
                        Intensity:          null,
                        NeutralMass:        cluster.NeutralMass,
                        IsotopicIndex:      rec.IsotopeIndex,
                        IsEstimated:        true,
                        EstimatedIntensity: rec.EstimatedIntensity,
                        Uncertainty:        double.IsInfinity(rec.Uncertainty) ? null : rec.Uncertainty,
                        Confidence:         rec.Confidence,
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
