namespace BrainMie.Core.Data;

/// <summary>
/// A group of peaks that belong to the same isotopic envelope, potentially with gaps
/// where multi-ion events caused peaks to be discarded by the instrument.
/// </summary>
/// <param name="Peaks">Observed peaks in the cluster, sorted by IsotopeIndex.</param>
/// <param name="GapIndices">Cluster-relative isotope indices where peaks are absent.</param>
/// <param name="NeutralMass">Intensity-weighted centroid neutral mass of the cluster.</param>
public record IsotopicCluster(
    List<ObservedPeak> Peaks,
    List<int> GapIndices,
    double NeutralMass);
