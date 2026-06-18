namespace BrainMie.Core.Data;

/// <summary>
/// One row in the proteoform-level summary: all isotopic clusters whose neutral
/// masses fall within <c>massTolerance</c> Da of each other are merged into a
/// single entry, summing ion counts across every contributing charge state.
/// </summary>
public record ProteoformSummary(
    double       NeutralMass,
    string       ChargeStates,
    int          ClusterCount,
    double       TotalObservedIonCount,
    double       TotalCorrectedIonCount,
    string       Hypothesis,
    double       MaxConfidence);
