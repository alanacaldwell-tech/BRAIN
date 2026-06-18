namespace BrainMie.Core.Data;

/// <summary>
/// One row in the pipeline output table. Observed peaks populate Mz and Intensity;
/// reconstructed peaks populate EstimatedIntensity and Uncertainty instead.
/// ClusterIonCount is the sum of raw observed ion events across the cluster.
/// CorrectedIonCount adds the reconstructed intensities for MIE gap peaks.
/// Both values are repeated on every row belonging to the same cluster.
/// </summary>
public record ProcessingRow(
    double?  Mz,
    int      Charge,
    double?  Intensity,
    double   NeutralMass,
    int      IsotopicIndex,
    bool     IsEstimated,
    double?  EstimatedIntensity,
    double?  Uncertainty,
    double   Confidence,
    string   Hypothesis,
    double   FitRSquared,
    bool     GapAtApex,
    string   Notes,
    double   ClusterIonCount,
    double   CorrectedIonCount);
