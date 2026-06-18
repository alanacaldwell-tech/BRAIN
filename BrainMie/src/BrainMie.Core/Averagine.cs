namespace BrainMie.Core;

/// <summary>
/// Scales the standard peptide averagine template (Senko et al., 1995) to a target neutral mass
/// to produce an approximate elemental composition when the true composition is unknown.
/// </summary>
public static class Averagine
{
    private static readonly Dictionary<string, double> Template = new()
    {
        ["C"] = 4.9384,
        ["H"] = 7.7583,
        ["N"] = 1.3577,
        ["O"] = 1.4773,
        ["S"] = 0.0417,
    };

    private const double AverageMass = 111.1254; // Da per residue unit

    /// <summary>
    /// Returns an integer elemental composition scaled to <paramref name="neutralMass"/>.
    /// Each count is at least 1.
    /// </summary>
    public static Dictionary<string, int> GetComposition(double neutralMass)
    {
        double n = neutralMass / AverageMass;
        return Template.ToDictionary(
            kvp => kvp.Key,
            kvp => Math.Max(1, (int)Math.Round(kvp.Value * n)));
    }
}
