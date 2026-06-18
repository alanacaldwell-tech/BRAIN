namespace BrainMie.Core.Data;

/// <summary>A single ion peak observed in the spectrum.</summary>
/// <param name="Mz">Measured mass-to-charge ratio.</param>
/// <param name="Intensity">Ion count or signal intensity.</param>
/// <param name="Charge">Charge state.</param>
/// <param name="IsotopeIndex">0-based index within the parent isotopic cluster (0 = first observed peak).</param>
public record ObservedPeak(double Mz, double Intensity, int Charge, int IsotopeIndex);
