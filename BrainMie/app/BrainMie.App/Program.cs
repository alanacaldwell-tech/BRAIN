using BrainMie.Core;
using BrainMie.Core.Data;

// ---- Argument parsing -------------------------------------------------------
if (args.Length > 0 && args[0] is "-h" or "--help")
{
    PrintHelp();
    return 0;
}

string inputPath;
int    maxGap       = 3;
double massMatchTol = 0.2;
double binWidth     = 0.02;
double massTol      = 1.0;
double minR2        = 0.5;

if (args.Length >= 1)
{
    // Command-line mode: BrainMie.exe input.dmt [options]
    inputPath = args[0];
    for (int i = 1; i < args.Length - 1; i++)
    {
        switch (args[i].ToLowerInvariant())
        {
            case "--max-gap":   maxGap       = int.Parse(args[++i]);    break;
            case "--match-tol": massMatchTol = double.Parse(args[++i]); break;
            case "--bin-width": binWidth     = double.Parse(args[++i]); break;
            case "--mass-tol":  massTol      = double.Parse(args[++i]); break;
            case "--min-r2":    minR2        = double.Parse(args[++i]); break;
        }
    }
}
else
{
    // Interactive mode: run from Visual Studio with no arguments.
    Console.WriteLine("BrainMie — Multi-Ion Event correction");
    Console.WriteLine("======================================");
    Console.WriteLine();

    inputPath = PromptFilePath();

    string? raw;
    raw = Prompt($"Max consecutive missing isotope peaks [{maxGap}]");
    if (!string.IsNullOrWhiteSpace(raw)) maxGap = int.Parse(raw);

    raw = Prompt($"Mass match tolerance in Da [{massMatchTol}]");
    if (!string.IsNullOrWhiteSpace(raw)) massMatchTol = double.Parse(raw);

    raw = Prompt($"Proteoform grouping tolerance in Da [{massTol}]");
    if (!string.IsNullOrWhiteSpace(raw)) massTol = double.Parse(raw);

    Console.WriteLine();
}

// Resolve to an absolute path so the user can see exactly where we're looking.
inputPath = Path.GetFullPath(inputPath.Trim('"').Trim());

if (!File.Exists(inputPath))
{
    Console.Error.WriteLine($"Error: file not found: {inputPath}");
    PauseIfInteractive();
    return 1;
}

// Write all outputs to a BrainMie\ subdirectory so they never collide with
// files from other analysis tools (e.g. ProteoformAnalyzer) in the same folder.
string outputDir  = Path.Combine(Path.GetDirectoryName(inputPath)!, "BrainMie");
Directory.CreateDirectory(outputDir);
string baseName   = Path.GetFileNameWithoutExtension(inputPath);

string outputPath     = Path.Combine(outputDir, baseName + ".corrected.csv");
string proteoformPath = Path.Combine(outputDir, baseName + ".proteoforms.csv");
string spectralPath   = Path.Combine(outputDir, baseName + ".spectral_check.html");

// ---- Run pipeline -----------------------------------------------------------
Console.WriteLine($"Input:      {inputPath}");
Console.WriteLine($"Parameters: max gap: {maxGap}  |  match tol: {massMatchTol} Da  |  mass tol: {massTol} Da  |  min R²: {minR2}");

List<ProcessingRow> rows;
List<(double Mz, int Charge, double? Intensity)> ions;
try
{
    (rows, ions) = Pipeline.ProcessDmt(
        inputPath,
        binWidth:     binWidth,
        maxGap:       maxGap,
        massMatchTol: massMatchTol,
        minFitR2:     minR2);
}
catch (Exception ex)
{
    Console.Error.WriteLine($"Error: {ex.Message}");
    PauseIfInteractive();
    return 1;
}

// ---- Write per-peak CSV output ----------------------------------------------
static string F(double? v) => v.HasValue ? v.Value.ToString("G6") : "";
static string S(bool    b) => b ? "true" : "false";
static string Q(string  s) => $"\"{s.Replace("\"", "\"\"")}\"";

const string CsvHeader =
    "IsotopeMass,Intensity,CentroidNeutralMass,IsotopicIndex,IsEstimated,IsSuppressed," +
    "EstimatedIntensity,CorrectedIntensity,Uncertainty,Confidence,Hypothesis," +
    "FitRSquared,GapAtApex,ClusterIonCount,CorrectedIonCount,Notes";

static void WriteRow(StreamWriter w, ProcessingRow row)
{
    w.WriteLine(
        $"{row.IsotopeMass:G6},{F(row.Intensity)},{row.NeutralMass:G6}," +
        $"{row.IsotopicIndex},{S(row.IsEstimated)},{S(row.IsSuppressed)}," +
        $"{F(row.EstimatedIntensity)},{F(row.CorrectedIntensity)}," +
        $"{F(row.Uncertainty)},{row.Confidence:G4},{row.Hypothesis}," +
        $"{row.FitRSquared:G4},{S(row.GapAtApex)}," +
        $"{row.ClusterIonCount:G6},{row.CorrectedIonCount:G6},{Q(row.Notes)}");
}

using (var w = new StreamWriter(outputPath))
{
    w.WriteLine(CsvHeader);
    foreach (var row in rows)
        WriteRow(w, row);
}

// ---- Write proteoform summary CSV ------------------------------------------
var allProteoforms = ProteoformAggregator.Aggregate(rows, ions, massTol);

// Noise filter: drop proteoforms below 0.5% of the largest, or 100 ions — whichever is higher.
double maxIons      = allProteoforms.Count > 0 ? allProteoforms.Max(pf => pf.TotalCorrectedIonCount) : 0;
double noiseFloor   = Math.Max(maxIons * 0.005, 100.0);
var    proteoforms  = allProteoforms.Where(pf => pf.TotalCorrectedIonCount >= noiseFloor).ToList();
int    noisyRemoved = allProteoforms.Count - proteoforms.Count;

const string PfHeader =
    "CentroidNeutralMass,ChargeStates,ClusterCount,TotalObservedIonCount,TotalCorrectedIonCount," +
    "Hypothesis,MaxConfidence";

using (var w = new StreamWriter(proteoformPath))
{
    w.WriteLine(PfHeader);
    foreach (var pf in proteoforms)
    {
        w.WriteLine(
            $"{pf.NeutralMass:G8},{pf.ChargeStates},{pf.ClusterCount}," +
            $"{pf.TotalObservedIonCount:G6},{pf.TotalCorrectedIonCount:G6}," +
            $"{pf.Hypothesis},{pf.MaxConfidence:G4}");
    }
}

// ---- Write spectral check HTML ----------------------------------------------
SpectralPlotter.GenerateHtml(rows, spectralPath);

// ---- Summary ----------------------------------------------------------------
int totalRows  = rows.Count;
int suppressed = rows.Count(r => r.IsSuppressed);
int mie        = rows.Count(r => r.Hypothesis == "mie");
int overlap    = rows.Count(r => r.Hypothesis == "overlap");
int ambiguous  = rows.Count(r => r.Hypothesis == "ambiguous");

Console.WriteLine($"Output:       {outputPath}  ({totalRows} peaks)");
Console.WriteLine($"Proteoforms:  {proteoformPath}  ({proteoforms.Count} entries, noise floor ≥ {noiseFloor:G3} ions, {noisyRemoved} removed)");
Console.WriteLine($"Spectral:     {spectralPath}");
Console.WriteLine($"              {mie} MIE  |  {overlap} overlap  |  {ambiguous} ambiguous  |  {suppressed} suppression-corrected");

PauseIfInteractive();
return 0;

// ---- Helpers ----------------------------------------------------------------
static string PromptFilePath()
{
    while (true)
    {
        Console.Write("  Path to .dmt file: ");
        string raw = Console.ReadLine()?.Trim().Trim('"') ?? "";
        if (string.IsNullOrWhiteSpace(raw)) continue;

        string full = Path.GetFullPath(raw);
        if (File.Exists(full)) return full;

        Console.WriteLine($"  File not found: {full}");
        Console.WriteLine("  (tip: you can drag-and-drop the file onto this window to paste its path)");
    }
}

static string? Prompt(string label)
{
    Console.Write($"  {label}: ");
    return Console.ReadLine()?.Trim();
}

static void PauseIfInteractive()
{
    // Keep the window open when launched directly from Visual Studio.
    if (!Console.IsInputRedirected)
    {
        Console.WriteLine();
        Console.Write("Press any key to exit...");
        Console.ReadKey(intercept: true);
    }
}

static void PrintHelp() => Console.WriteLine("""
    BrainMie — Multi-Ion Event correction for .dmt mass spectrometry files

    Usage:
      BrainMie <input.dmt> [options]
      BrainMie                         (interactive prompts)

    Options:
      --max-gap N       Max consecutive missing isotope peaks per cluster (default 3)
      --match-tol F     Mass match tolerance in Da for isotope peak assignment (default 0.2)
      --bin-width F     Mass bin width in Da for ion-event aggregation (default 0.02)
      --mass-tol F      Neutral-mass tolerance in Da for proteoform grouping (default 1.0)
      --min-r2 F        Minimum averagine fit R² to keep a cluster (default 0.5)
      -h, --help        Show this help

    Output files written to a BrainMie\ subfolder next to the input:
      <input>.corrected.csv        All observed peaks; CorrectedIntensity filled for suppressed peaks
      <input>.proteoforms.csv      One row per proteoform, ion counts summed across charge states
      <input>.spectral_check.html  SVG charts showing suppression correction per cluster
    """);
