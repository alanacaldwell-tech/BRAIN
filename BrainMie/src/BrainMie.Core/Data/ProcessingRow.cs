namespace BrainMie.Core.Data;

/// <summary>
/// One row in the pipeline output table. Observed peaks populate Mz and Intensity;
/// reconstructed peaks populate EstimatedIntensity and Uncertainty instead.
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
    string   Notes);
