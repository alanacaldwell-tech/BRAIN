namespace BrainMie.Core.Data;

/// <summary>
/// Outcome of the MIE-vs-overlap hypothesis test for a gapped isotopic cluster.
/// </summary>
/// <param name="Hypothesis">
///   "mie"       – single species whose apex was removed by a multi-ion event.<br/>
///   "overlap"   – two distinct co-eluting species.<br/>
///   "ambiguous" – insufficient evidence to decide.<br/>
///   "no_gap"    – cluster has no missing peaks; classification not needed.
/// </param>
/// <param name="Confidence">0–1 probability-like score for the stated hypothesis.</param>
/// <param name="FitRSquared">R² of the single-species averagine fit over observed peaks.</param>
/// <param name="GapAtApex">Whether any gap falls at the theoretical intensity maximum.</param>
/// <param name="NConsecutiveGaps">Longest run of consecutive missing isotope peaks.</param>
/// <param name="Notes">Human-readable caveats (e.g. reconstruction reliability warnings).</param>
public record ClassificationResult(
    string Hypothesis,
    double Confidence,
    double FitRSquared,
    bool GapAtApex,
    int NConsecutiveGaps,
    string Notes = "");
