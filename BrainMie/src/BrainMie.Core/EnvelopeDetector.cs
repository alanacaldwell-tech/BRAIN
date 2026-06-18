namespace BrainMie.Core;

using BrainMie.Core.Data;

/// <summary>
/// Bins individual ion detection events into peaks, then groups peaks into
/// isotopic clusters — tolerating gaps of up to <c>maxGap</c> consecutive
/// missing isotopes so that MIE-induced holes do not split one cluster into two.
/// </summary>
public static class EnvelopeDetector
{
    /// <summary>
    /// Aggregate individual ion m/z values into (centre_mz, count) peaks by binning.
    /// </summary>
    /// <param name="mzValues">Raw m/z values from the Ion table (one per detection event).</param>
    /// <param name="binWidth">Bin width in Da. Should be much smaller than the isotope spacing (1/charge Da).</param>
    public static List<(double Mz, double Intensity)> BinIons(
        IEnumerable<double> mzValues, double binWidth = 0.005)
    {
        var bins = new Dictionary<long, double>();
        foreach (double mz in mzValues)
        {
            // Math.Round uses banker's rounding (round-half-to-even), matching Python's round().
            long key = (long)Math.Round(mz / binWidth);
            bins[key] = bins.GetValueOrDefault(key) + 1.0;
        }
        return bins
            .OrderBy(kvp => kvp.Key)
            .Select(kvp => (kvp.Key * binWidth, kvp.Value))
            .ToList();
    }

    /// <summary>
    /// Reduce a binned peak list to local maxima, rejecting noise bins.
    ///
    /// A bin qualifies as a peak if:
    ///   1. Its intensity exceeds <paramref name="minProminenceFraction"/> × the
    ///      spectrum maximum (with a hard floor of 2 ions), and
    ///   2. It is strictly greater than both of its immediate list-neighbours
    ///      (list-index adjacency on the sorted sparse bin list).
    ///
    /// Because the bin list is sparse, list-index adjacency is intentional: two
    /// real isotope peaks that are far apart in m/z are not adjacent in the list
    /// and do not interfere with each other's local-max test, even at high charge
    /// states where isotope spacing approaches the bin width.
    /// </summary>
    /// <param name="bins">Sorted (by m/z) output of <see cref="BinIons"/>.</param>
    /// <param name="minProminenceFraction">Fraction of the spectrum maximum below
    /// which a bin is unconditionally rejected. Default 0.1 %.</param>
    public static List<(double Mz, double Intensity)> PickLocalMaxima(
        List<(double Mz, double Intensity)> bins,
        double minProminenceFraction = 0.001)
    {
        if (bins.Count == 0) return [];

        double maxIntensity = bins.Max(b => b.Intensity);
        double threshold    = Math.Max(maxIntensity * minProminenceFraction, 2.0);

        var peaks = new List<(double Mz, double Intensity)>();
        for (int i = 0; i < bins.Count; i++)
        {
            if (bins[i].Intensity < threshold) continue;

            double prev = i > 0              ? bins[i - 1].Intensity : 0.0;
            double next = i < bins.Count - 1 ? bins[i + 1].Intensity : 0.0;

            if (bins[i].Intensity >= prev && bins[i].Intensity > next)
                peaks.Add(bins[i]);
        }
        return peaks;
    }

    /// <summary>
    /// Group (mz, intensity) peaks into isotopic clusters for a given charge state.
    /// </summary>
    /// <param name="peaks">Observed peaks sorted (or unsorted) by m/z.</param>
    /// <param name="charge">Charge state to use for expected isotope spacing.</param>
    /// <param name="maxGap">Maximum number of consecutive missing isotopes still treated as the same cluster.</param>
    /// <param name="ppmTolerance">m/z tolerance in parts-per-million for matching peaks to isotope positions.</param>
    public static List<IsotopicCluster> FindClusters(
        List<(double Mz, double Intensity)> peaks,
        int    charge,
        int    maxGap        = 3,
        double ppmTolerance  = 20.0)
    {
        if (peaks.Count == 0 || charge <= 0) return [];

        var sorted = peaks.OrderBy(p => p.Mz).ToList();
        int    n     = sorted.Count;
        double delta = Constants.IsotopeSpacing / charge;

        // ---- Union-Find with path compression --------------------------------
        int[] parent = Enumerable.Range(0, n).ToArray();

        int Find(int x)
        {
            while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; }
            return x;
        }

        void Union(int x, int y) => parent[Find(x)] = Find(y);

        for (int i = 0; i < n; i++)
        {
            for (int j = i + 1; j < n; j++)
            {
                double diff = sorted[j].Mz - sorted[i].Mz;
                if (diff > (maxGap + 0.5) * delta) break;

                double tol = Math.Max(sorted[i].Mz * ppmTolerance * 1e-6, 0.01);
                for (int k = 1; k <= maxGap; k++)
                {
                    if (Math.Abs(diff - k * delta) < tol)
                    {
                        Union(i, j);
                        break;
                    }
                }
            }
        }

        // ---- Group components ------------------------------------------------
        var components = new Dictionary<int, List<int>>();
        for (int i = 0; i < n; i++)
        {
            int root = Find(i);
            if (!components.ContainsKey(root)) components[root] = [];
            components[root].Add(i);
        }

        var clusters = new List<IsotopicCluster>();
        foreach (var indices in components.Values)
        {
            if (indices.Count < 2) continue;

            var compPeaks = indices.OrderBy(i => sorted[i].Mz)
                                   .Select(i => sorted[i])
                                   .ToList();
            double mzStart = compPeaks[0].Mz;

            // Assign integer isotope index relative to the first observed peak.
            // Merge any two binned peaks that round to the same index (sum intensities,
            // intensity-weighted average m/z) so downstream dictionaries never see duplicates.
            var indexed = compPeaks
                .Select(p => (p.Mz, p.Intensity, IsotopeIndex: (int)Math.Round((p.Mz - mzStart) / delta)))
                .GroupBy(p => p.IsotopeIndex)
                .Select(g =>
                {
                    double totalIntensity = g.Sum(p => p.Intensity);
                    double weightedMz     = g.Sum(p => p.Mz * p.Intensity) / totalIntensity;
                    return (Mz: weightedMz, Intensity: totalIntensity, IsotopeIndex: g.Key);
                })
                .OrderBy(p => p.IsotopeIndex)
                .ToList();

            int            maxIdx      = indexed.Max(p => p.IsotopeIndex);
            HashSet<int>   observedSet = indexed.Select(p => p.IsotopeIndex).ToHashSet();
            List<int>      gapIndices  = Enumerable.Range(0, maxIdx + 1)
                                                   .Where(i => !observedSet.Contains(i))
                                                   .ToList();

            double neutralMass = mzStart * charge - charge * Constants.HMass;

            var observedPeaks = indexed
                .Select(p => new ObservedPeak(p.Mz, p.Intensity, charge, p.IsotopeIndex))
                .ToList();

            clusters.Add(new IsotopicCluster(observedPeaks, charge, gapIndices, neutralMass));
        }

        return clusters;
    }
}
