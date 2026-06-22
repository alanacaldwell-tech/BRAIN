namespace BrainMie.Core.Data;

/// <summary>A single ion peak observed in the neutral-mass spectrum.</summary>
/// <param name="Mass">Neutral mass of this specific isotope peak.</param>
/// <param name="Intensity">Ion count or signal intensity.</param>
/// <param name="IsotopeIndex">0-based index within the parent isotopic cluster (0 = first observed peak).</param>
public record ObservedPeak(double Mass, double Intensity, int IsotopeIndex);
