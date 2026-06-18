namespace BrainMie.Core.Data;

/// <summary>
/// One row in the pipeline output table.
///
/// Observed peaks:      Mz and Intensity are populated.
/// Gap-reconstructed:   IsEstimated=true; EstimatedIntensity and Uncertainty populated.
/// Suppressed peaks:    Observed but below 50% of theoretical expectation (MIE partial
///                      suppression). IsSuppressed=true; CorrectedIntensity holds the
///                      fitted expected value; Intensity holds the raw observed value.
///
/// ClusterIonCount  = sum of raw observed ions for the cluster.
/// CorrectedIonCount = observed (non-suppressed) + corrected (suppressed) + estimated (gaps).
/// Both are repeated on every row of the same cluster.
/// </summary>
public record ProcessingRow(
    double?  Mz,
    int      Charge,
    double?  Intensity,
    double   NeutralMass,
    int      IsotopicIndex,
    bool     IsEstimated,
    bool     IsSuppressed,
    double?  EstimatedIntensity,
    double?  CorrectedIntensity,
    double?  Uncertainty,
    double   Confidence,
    string   Hypothesis,
    double   FitRSquared,
    bool     GapAtApex,
    string   Notes,
    double   ClusterIonCount,
    double   CorrectedIonCount);
