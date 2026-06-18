namespace BrainMie.Core;

/// <summary>
/// Computes theoretical aggregated isotopic distributions using polynomial convolution,
/// analogous to the BRAIN algorithm (Dittwald et al., 2013).
///
/// Each element's isotope pattern is represented as a probability polynomial over
/// integer neutron counts. For n atoms of an element the polynomial is raised to
/// the nth power via binary exponentiation. All element polynomials are then
/// convolved to give the full distribution.
///
/// Intensities are normalised so the most abundant peak equals 1.0.
/// </summary>
public static class IsotopeDistribution
{
    // Natural isotope abundances by element.
    // Each entry: (extra_neutrons_relative_to_lightest_isotope, fractional_abundance)
    private static readonly Dictionary<string, (int Shift, double Abundance)[]> IsotopeTable = new()
    {
        ["C"] = [(0, 0.9893),   (1, 0.0107)],
        ["H"] = [(0, 0.999885), (1, 0.000115)],
        ["N"] = [(0, 0.99632),  (1, 0.00368)],
        ["O"] = [(0, 0.99757),  (1, 0.00038), (2, 0.00205)],
        ["S"] = [(0, 0.9502),   (1, 0.0075),  (2, 0.0421),  (4, 0.0002)],
    };

    /// <summary>
    /// Compute the theoretical isotopic envelope for a molecule of <paramref name="neutralMass"/> Da
    /// observed at charge state <paramref name="charge"/>.
    /// </summary>
    /// <param name="neutralMass">Neutral monoisotopic mass in Da.</param>
    /// <param name="charge">Charge state (must be > 0).</param>
    /// <param name="nPeaks">Maximum number of isotope peaks to return.</param>
    /// <returns>
    /// List of (Mz, NormalisedIntensity) tuples ordered by m/z. The most abundant
    /// peak has intensity 1.0. Returns an empty list if the mass is invalid.
    /// </returns>
    public static List<(double Mz, double Intensity)> ComputeEnvelope(
        double neutralMass, int charge, int nPeaks = 25)
    {
        if (neutralMass <= 0 || charge <= 0) return [];

        var composition = Averagine.GetComposition(neutralMass);
        double[] dist = ComputeDistribution(composition, nPeaks);

        double maxInt = dist.Length > 0 ? dist.Max() : 0;
        if (maxInt <= 0) return [];

        double mz0   = (neutralMass + charge * Constants.HMass) / charge;
        double step  = Constants.IsotopeSpacing / charge;

        var peaks = new List<(double Mz, double Intensity)>(dist.Length);
        for (int k = 0; k < dist.Length; k++)
            peaks.Add((mz0 + k * step, dist[k] / maxInt));

        return peaks;
    }

    // ---------------------------------------------------------------------------
    // Internal polynomial math
    // ---------------------------------------------------------------------------

    private static double[] ComputeDistribution(Dictionary<string, int> composition, int maxPeaks)
    {
        double[] result = [1.0];

        foreach (var (element, count) in composition)
        {
            if (!IsotopeTable.TryGetValue(element, out var isotopes)) continue;
            double[] elemPoly = BuildElementPoly(isotopes);
            double[] elemDist = PolyPower(elemPoly, count, maxPeaks);
            result = Truncate(Convolve(result, elemDist), maxPeaks);
        }

        return result;
    }

    private static double[] BuildElementPoly((int Shift, double Abundance)[] isotopes)
    {
        int maxShift = isotopes.Max(i => i.Shift);
        double[] poly = new double[maxShift + 1];
        foreach (var (shift, abundance) in isotopes)
            poly[shift] = abundance;
        return poly;
    }

    /// <summary>Discrete convolution of two probability arrays.</summary>
    private static double[] Convolve(double[] a, double[] b)
    {
        double[] result = new double[a.Length + b.Length - 1];
        for (int i = 0; i < a.Length; i++)
            for (int j = 0; j < b.Length; j++)
                result[i + j] += a[i] * b[j];
        return result;
    }

    /// <summary>
    /// Raise a probability polynomial to the <paramref name="n"/>th power using
    /// binary (repeated squaring) exponentiation, truncating at <paramref name="maxLen"/>
    /// at each step to keep intermediate arrays small.
    /// </summary>
    private static double[] PolyPower(double[] poly, int n, int maxLen)
    {
        if (n == 0) return [1.0];
        if (n == 1) return Truncate(poly, maxLen);

        double[] half   = PolyPower(poly, n / 2, maxLen);
        double[] result = Truncate(Convolve(half, half), maxLen);
        if (n % 2 != 0)
            result = Truncate(Convolve(result, poly), maxLen);
        return result;
    }

    private static double[] Truncate(double[] a, int maxLen) =>
        a.Length <= maxLen ? a : a[..maxLen];
}
