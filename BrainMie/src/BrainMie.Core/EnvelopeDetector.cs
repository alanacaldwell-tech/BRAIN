namespace BrainMie.Core;

using BrainMie.Core.Data;

/// <summary>
/// Bins individual ion detection events into peaks in neutral-mass space, then groups
/// peaks into isotopic clusters — tolerating gaps of up to <c>maxGap</c> consecutive
/// missing isotopes so that MIE-induced holes do not split one cluster into two.
///
/// Working in mass space (rather than m/z) keeps the isotope spacing constant at
/// 1.003355 Da regardless of charge state, giving ~50 mass bins per isotope peak at
/// 0.02 Da bin width. This makes each isotope peak robustly identifiable as a
/// local maximum in the sparse bin list, even when the centroid peak is MIE-suppressed.
/// </summary>
public static class EnvelopeDetector
{
    /// <summary>
    /// Aggregate individual neutral-mass values into (centre_mass, count) peaks by binning.
    /// </summary>
    /// <param name="massValues">Neutral-mass values (one per detection event).</param>
    /// <param name="binWidth">Bin width in Da. Should be much smaller than the isotope spacing (1.003355 Da).</param>
    public static List<(double Mass, double Intensity)> BinIons(
        IEnumerable<double> massValues, double binWidth = 0.02)
    {
        var bins = new Dictionary<long, double>();
        foreach (double mass in massValues)
        {
            long key = (long)Math.Round(mass / binWidth);
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
    /// real isotope peaks are not adjacent in the list and do not interfere with
    /// each other's local-max test. In mass space with 0.02 Da bins and 1.003355 Da
    /// isotope spacing there are ~50 bins per isotope peak, so suppressed peaks are
    /// still reliably detected as local maxima within their own bin cluster.
    /// </summary>
    /// <param name="bins">Sorted (by mass) output of <see cref="BinIons"/>.</param>
    /// <param name="minProminenceFraction">Fraction of the spectrum maximum below
    /// which a bin is unconditionally rejected. Default 0.1 %.</param>
    public static List<(double Mass, double Intensity)> PickLocalMaxima(
        List<(double Mass, double Intensity)> bins,
        double minProminenceFraction = 0.001)
    {
        if (bins.Count == 0) return [];

        double maxIntensity = bins.Max(b => b.Intensity);
        double threshold    = Math.Max(maxIntensity * minProminenceFraction, 2.0);

        var peaks = new List<(double Mass, double Intensity)>();
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
    /// Group (mass, intensity) peaks into isotopic clusters using the constant
    /// neutral-mass isotope spacing (1.003355 Da).
    /// </summary>
    /// <param name="peaks">Observed peaks in neutral-mass space (any order).</param>
    /// <param name="maxGap">Maximum number of consecutive missing isotopes still treated as the same cluster.</param>
    /// <param name="massMatchTol">
    /// Absolute mass tolerance in Da for matching peaks to isotope positions.
    /// Default 0.2 Da (~20% of the isotope spacing), suitable for Orbitrap data across
    /// the full 10–200 kDa range. This replaces the former per-charge ppm tolerance.
    /// </param>
    public static List<IsotopicCluster> FindClusters(
        List<(double Mass, double Intensity)> peaks,
        int    maxGap       = 3,
        double massMatchTol = 0.2)
    {
        if (peaks.Count == 0) return [];

        var sorted = peaks.OrderBy(p => p.Mass).ToList();
        int    n     = sorted.Count;
        double delta = Constants.IsotopeSpacing;

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
                double diff = sorted[j].Mass - sorted[i].Mass;
                if (diff > (maxGap + 0.5) * delta) break;

                for (int k = 1; k <= maxGap; k++)
                {
                    if (Math.Abs(diff - k * delta) < massMatchTol)
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

            var compPeaks = indices.OrderBy(i => sorted[i].Mass)
                                   .Select(i => sorted[i])
                                   .ToList();
            double massStart = compPeaks[0].Mass;

            // Assign integer isotope index relative to the first observed peak.
            // Merge any two binned peaks that round to the same index (sum intensities,
            // intensity-weighted average mass) so downstream dictionaries never see duplicates.
            var indexed = compPeaks
                .Select(p => (p.Mass, p.Intensity, IsotopeIndex: (int)Math.Round((p.Mass - massStart) / delta)))
                .GroupBy(p => p.IsotopeIndex)
                .Select(g =>
                {
                    double totalIntensity = g.Sum(p => p.Intensity);
                    double weightedMass   = g.Sum(p => p.Mass * p.Intensity) / totalIntensity;
                    return (Mass: weightedMass, Intensity: totalIntensity, IsotopeIndex: g.Key);
                })
                .OrderBy(p => p.IsotopeIndex)
                .ToList();

            int          maxIdx      = indexed.Max(p => p.IsotopeIndex);
            HashSet<int> observedSet = indexed.Select(p => p.IsotopeIndex).ToHashSet();
            List<int>    gapIndices  = Enumerable.Range(0, maxIdx + 1)
                                                 .Where(i => !observedSet.Contains(i))
                                                 .ToList();

            // Centroid neutral mass: intensity-weighted average across all isotope peaks.
            double totalInt    = indexed.Sum(p => p.Intensity);
            double neutralMass = indexed.Sum(p => p.Mass * p.Intensity) / totalInt;

            var observedPeaks = indexed
                .Select(p => new ObservedPeak(p.Mass, p.Intensity, p.IsotopeIndex))
                .ToList();

            clusters.Add(new IsotopicCluster(observedPeaks, gapIndices, neutralMass));
        }

        return clusters;
    }
}
