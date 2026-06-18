namespace BrainMie.Core;

using System.Text;
using BrainMie.Core.Data;

/// <summary>
/// Produces a self-contained HTML file with one SVG chart per rescued cluster,
/// showing the isotopic envelope before and after MIE peak reconstruction.
/// </summary>
public static class SpectralPlotter
{
    // SVG layout constants
    private const int SvgW    = 400;
    private const int SvgH    = 160;
    private const int Left    = 36;
    private const int Top     = 18;
    private const int Right   = 8;
    private const int Bottom  = 28;
    private const int ChartW  = SvgW - Left - Right;
    private const int ChartH  = SvgH - Top - Bottom;

    /// <summary>
    /// Write an HTML spectral-check report to <paramref name="outputPath"/>.
    /// Only clusters that contain at least one reconstructed peak are included,
    /// sorted by corrected ion count descending and capped at 200 charts.
    /// </summary>
    public static void GenerateHtml(List<ProcessingRow> rows, string outputPath, int nPeaks = 25)
    {
        var clusterGroups = rows
            .GroupBy(r => (r.NeutralMass, r.Charge))
            .Where(g => g.Any(r => r.IsSuppressed))
            .OrderByDescending(g => g.First().CorrectedIonCount)
            .Take(200)
            .ToList();

        var html = new StringBuilder();
        html.AppendLine("<!DOCTYPE html>");
        html.AppendLine("<html><head><meta charset='utf-8'>");
        html.AppendLine("<title>BrainMie — Spectral Check</title>");
        html.AppendLine(Css());
        html.AppendLine("</head><body>");
        html.AppendLine("<h1>BrainMie — Spectral Check</h1>");

        if (clusterGroups.Count == 0)
        {
            html.AppendLine("<p>No MIE-suppressed clusters found in this dataset.</p>");
        }
        else
        {
            html.AppendLine($"<p class='sub'>Showing {clusterGroups.Count} suppression-corrected cluster(s), " +
                            "ranked by corrected ion count. " +
                            "<span class='obs'>&#9646; Observed</span> &nbsp; " +
                            "<span class='est'>&#9646; Suppression correction (stacked)</span> &nbsp; " +
                            "<span class='theo'>--- Theoretical averagine</span></p>");
            html.AppendLine("<div class='grid'>");
            foreach (var g in clusterGroups)
                html.AppendLine(BuildCard(g.Key.NeutralMass, g.Key.Charge, g.ToList(), nPeaks));
            html.AppendLine("</div>");
        }

        html.AppendLine("</body></html>");
        File.WriteAllText(outputPath, html.ToString());
    }

    // -------------------------------------------------------------------------

    private static string BuildCard(
        double neutralMass, int charge, List<ProcessingRow> rows, int nPeaks)
    {
        var observed  = rows.OrderBy(r => r.IsotopicIndex).ToList();
        if (observed.Count == 0) return "";

        var first      = observed.First();
        string hypo    = first.Hypothesis;
        double conf    = first.Confidence;
        double ions    = first.CorrectedIonCount;

        int minIdx = observed.Min(r => r.IsotopicIndex);
        int maxIdx = observed.Max(r => r.IsotopicIndex);
        int nCols  = maxIdx - minIdx + 1;
        if (nCols < 1) return "";

        // Theoretical envelope aligned to observed peaks
        var theoretical = IsotopeDistribution.ComputeEnvelope(neutralMass, charge, nPeaks);
        var obsPeaks    = observed
            .Select(r => new ObservedPeak(r.Mz ?? 0, r.Intensity ?? 0, r.Charge, r.IsotopicIndex))
            .ToList();
        var (offset, scale, _) = Classifier.FitAlignment(obsPeaks, theoretical);

        // Normalisation: suppressed peaks reach CorrectedIntensity, normal peaks reach Intensity.
        double norm = observed.Max(r =>
            r.IsSuppressed ? (r.CorrectedIntensity ?? r.Intensity ?? 0) : (r.Intensity ?? 0));
        if (norm == 0) return "";

        // Helpers — SVG coordinate space
        double ColX(int idx) => Left + (idx - minIdx + 0.175) * ((double)ChartW / nCols);
        double BarW()        => (double)ChartW / nCols * 0.65;
        double ValY(double v) => Top + ChartH * (1.0 - Math.Min(v, 1.0));
        double ValH(double v) => ChartH * Math.Min(Math.Max(v, 0), 1.0);

        var svg = new StringBuilder();
        svg.AppendLine($"<svg viewBox='0 0 {SvgW} {SvgH}' width='{SvgW}' height='{SvgH}'>");

        // Axes
        svg.AppendLine($"<line x1='{Left}' y1='{Top}' x2='{Left}' y2='{Top + ChartH}' stroke='#ccc' stroke-width='1'/>");
        svg.AppendLine($"<line x1='{Left}' y1='{Top + ChartH}' x2='{Left + ChartW}' y2='{Top + ChartH}' stroke='#ccc' stroke-width='1'/>");
        svg.AppendLine($"<text x='{Left - 4}' y='{Top + 4}' text-anchor='end' font-size='9' fill='#999'>1.0</text>");
        svg.AppendLine($"<text x='{Left - 4}' y='{Top + ChartH + 1}' text-anchor='end' font-size='9' fill='#999'>0</text>");

        // Theoretical line (dashed gray)
        var theoPoints = new List<string>();
        for (int idx = minIdx; idx <= maxIdx; idx++)
        {
            int ti = idx + offset;
            if (ti >= 0 && ti < theoretical.Count)
            {
                // Scale theoretical so a perfect-fit cluster would overlay the bars exactly
                double normTheo = scale * theoretical[ti].Intensity / norm;
                double cx = Left + (idx - minIdx + 0.5) * ((double)ChartW / nCols);
                theoPoints.Add($"{cx:F1},{ValY(normTheo):F1}");
            }
        }
        if (theoPoints.Count >= 2)
            svg.AppendLine($"<polyline points='{string.Join(" ", theoPoints)}' " +
                           "fill='none' stroke='#b0b0b0' stroke-width='1.5' stroke-dasharray='4,3'/>");

        // Blue bars: all observed peaks (suppressed or not)
        foreach (var r in observed)
        {
            double v = (r.Intensity ?? 0) / norm;
            double x = ColX(r.IsotopicIndex);
            svg.AppendLine($"<rect x='{x:F1}' y='{ValY(v):F1}' width='{BarW():F1}' height='{ValH(v):F1}' fill='#4472C4' opacity='0.9'/>");
            svg.AppendLine($"<text x='{x + BarW() / 2:F1}' y='{Top + ChartH + 14}' text-anchor='middle' font-size='9' fill='#777'>{r.IsotopicIndex}</text>");
        }

        // Orange stacked segment: suppressed peaks — draws from top of blue bar up to corrected height
        foreach (var r in observed.Where(r => r.IsSuppressed))
        {
            double obsV  = (r.Intensity ?? 0) / norm;
            double corV  = (r.CorrectedIntensity ?? r.Intensity ?? 0) / norm;
            double delta = corV - obsV;
            if (delta < 0.5 / ChartH) continue;  // sub-pixel, skip
            double x = ColX(r.IsotopicIndex);
            // Draw segment from obsV to corV (i.e. y from ValY(corV) down to ValY(obsV))
            svg.AppendLine($"<rect x='{x:F1}' y='{ValY(corV):F1}' width='{BarW():F1}' height='{ValH(delta):F1}' fill='#ED7D31' opacity='0.9'/>");
        }

        svg.AppendLine("</svg>");

        string label = HypoLabel(hypo);
        string title = $"{neutralMass:F2} Da &nbsp;&nbsp; z={charge} &nbsp;&nbsp; " +
                       $"<span class='badge {hypo}'>{label}</span> &nbsp;&nbsp; " +
                       $"conf: {conf:F2} &nbsp;&nbsp; {ions:N0} ions";

        return $"<div class='card'><p class='ctitle'>{title}</p>{svg}</div>";
    }

    private static string HypoLabel(string h) => h switch
    {
        "mie"       => "MIE",
        "ambiguous" => "Ambiguous",
        "overlap"   => "Overlap",
        _           => h
    };

    private static string Css() => @"<style>
body { font-family: Arial, sans-serif; background: #efefef; margin: 24px; color: #333; }
h1 { margin-bottom: 4px; font-size: 20px; }
p.sub { color: #666; margin-top: 0; font-size: 13px; }
span.obs  { color: #4472C4; font-weight: bold; }
span.est  { color: #ED7D31; font-weight: bold; }
span.theo { color: #999; }
.grid { display: flex; flex-wrap: wrap; gap: 14px; margin-top: 16px; }
.card { background: white; border-radius: 8px; padding: 10px 14px 6px;
        box-shadow: 0 1px 4px rgba(0,0,0,0.13); }
p.ctitle { margin: 0 0 6px 0; font-size: 11.5px; color: #444; }
.badge { display: inline-block; padding: 1px 6px; border-radius: 3px;
         font-size: 10px; font-weight: bold; color: white; }
.badge.mie       { background: #4472C4; }
.badge.ambiguous { background: #FFC000; color: #333; }
.badge.overlap   { background: #C00000; }
</style>";
}
