namespace BrainMie.Core;

using BrainMie.Core.Data;

/// <summary>
/// Groups per-cluster rows into proteoform-level summaries by merging clusters
/// whose neutral masses fall within <paramref name="massTolerance"/> Da of each
/// other.  Ion counts are summed across all contributing charge states.
///
/// Grouping uses a greedy single-linkage scan (sort by mass, start a new group
/// whenever the gap from the current group anchor exceeds the tolerance).  This
/// matches the behaviour of the ProteoformAnalyzer reference pipeline.
/// </summary>
public static class ProteoformAggregator
{
    public static List<ProteoformSummary> Aggregate(
        List<ProcessingRow> rows,
        double massTolerance = 1.0)
    {
        // One record per (NeutralMass, Charge) cluster — skip estimated rows to
        // avoid counting the same cluster twice (they share identical header fields).
        var clusters = rows
            .Where(r => !r.IsEstimated)
            .GroupBy(r => (r.NeutralMass, r.Charge))
            .Select(g =>
            {
                var first = g.First();
                return (
                    NeutralMass:        first.NeutralMass,
                    Charge:             first.Charge,
                    ObservedIons:       first.ClusterIonCount,
                    CorrectedIons:      first.CorrectedIonCount,
                    Hypothesis:         first.Hypothesis,
                    Confidence:         first.Confidence
                );
            })
            .OrderBy(c => c.NeutralMass)
            .ToList();

        if (clusters.Count == 0) return [];

        var summaries = new List<ProteoformSummary>();

        // Greedy grouping: anchor = first cluster in each group.
        int start = 0;
        while (start < clusters.Count)
        {
            double anchor = clusters[start].NeutralMass;
            int end = start;
            while (end + 1 < clusters.Count &&
                   clusters[end + 1].NeutralMass - anchor <= massTolerance)
                end++;

            var group = clusters[start..(end + 1)];

            // Intensity-weighted average neutral mass
            double totalCorrected = group.Sum(c => c.CorrectedIons);
            double repMass = totalCorrected > 0
                ? group.Sum(c => c.NeutralMass * c.CorrectedIons) / totalCorrected
                : group.Average(c => c.NeutralMass);

            string chargeStates = string.Join("|", group.Select(c => c.Charge).Order().Distinct());

            // Majority-vote hypothesis; tie-break toward "mie" over "ambiguous" over "overlap"
            string hypothesis = group
                .GroupBy(c => c.Hypothesis)
                .OrderByDescending(g => g.Count())
                .ThenBy(g => HypothesisPriority(g.Key))
                .First().Key;

            summaries.Add(new ProteoformSummary(
                NeutralMass:            repMass,
                ChargeStates:           chargeStates,
                ClusterCount:           group.Count,
                TotalObservedIonCount:  group.Sum(c => c.ObservedIons),
                TotalCorrectedIonCount: totalCorrected,
                Hypothesis:             hypothesis,
                MaxConfidence:          group.Max(c => c.Confidence)));

            start = end + 1;
        }

        return summaries;
    }

    private static int HypothesisPriority(string h) => h switch
    {
        "mie"       => 0,
        "ambiguous" => 1,
        "overlap"   => 2,
        _           => 3
    };
}
