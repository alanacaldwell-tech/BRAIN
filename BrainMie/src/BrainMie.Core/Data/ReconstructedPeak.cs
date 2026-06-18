namespace BrainMie.Core.Data;

/// <summary>A missing isotope peak whose intensity has been estimated from the surviving envelope.</summary>
/// <param name="Mz">Computed m/z of the missing peak.</param>
/// <param name="Charge">Charge state.</param>
/// <param name="EstimatedIntensity">Predicted ion count derived from the averagine scale factor.</param>
/// <param name="Uncertainty">
///   1-sigma prediction interval. <see cref="double.PositiveInfinity"/> when the gap falls
///   outside the theoretical range and no meaningful estimate is possible.
/// </param>
/// <param name="IsotopeIndex">Cluster-relative isotope index matching the originating gap.</param>
/// <param name="Confidence">Confidence inherited from the parent <see cref="ClassificationResult"/>.</param>
public record ReconstructedPeak(
    double Mz,
    int Charge,
    double EstimatedIntensity,
    double Uncertainty,
    int IsotopeIndex,
    double Confidence);
