using BrainMie.Core;
using BrainMie.Core.Data;

// ---- Argument parsing -------------------------------------------------------
if (args.Length > 0 && args[0] is "-h" or "--help")
{
    PrintHelp();
    return 0;
}

string inputPath;
int?   chargeMinOverride = null;
int?   chargeMaxOverride = null;
int    maxGap     = 3;
double ppm        = 10.0;
double binWidth   = 0.02;
double massTol    = 1.0;

if (args.Length >= 1)
{
    // Command-line mode: BrainMie.exe input.dmt [options]
    inputPath = args[0];
    for (int i = 1; i < args.Length - 1; i++)
    {
        switch (args[i].ToLowerInvariant())
        {
            case "--charge-min": chargeMinOverride = int.Parse(args[++i]);    break;
            case "--charge-max": chargeMaxOverride = int.Parse(args[++i]);    break;
            case "--max-gap":    maxGap    = int.Parse(args[++i]);    break;
            case "--ppm":        ppm       = double.Parse(args[++i]); break;
            case "--bin-width":  binWidth  = double.Parse(args[++i]); break;
            case "--mass-tol":   massTol   = double.Parse(args[++i]); break;
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

    string raw;
    raw = Prompt($"Max consecutive missing isotope peaks [{maxGap}]");
    if (!string.IsNullOrWhiteSpace(raw)) maxGap = int.Parse(raw);

    raw = Prompt($"m/z tolerance in ppm (Orbitrap: 5–10, Q-TOF: 10–20) [{ppm}]");
    if (!string.IsNullOrWhiteSpace(raw)) ppm = double.Parse(raw);

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

string outputPath      = Path.ChangeExtension(inputPath, ".corrected.csv");
string proteoformPath  = Path.ChangeExtension(inputPath, ".proteoforms.csv");
string spectralPath    = Path.ChangeExtension(inputPath, ".spectral_check.html");

// ---- Auto-detect charge range from file -------------------------------------
var (detectedMin, detectedMax) = DmtReader.ReadChargeRange(inputPath);
int chargeMin = chargeMinOverride ?? detectedMin;
int chargeMax = chargeMaxOverride ?? detectedMax;

// ---- Run pipeline -----------------------------------------------------------
Console.WriteLine($"Input:   {inputPath}");
Console.WriteLine($"Charges: {chargeMin}–{chargeMax}{(chargeMinOverride is null && chargeMaxOverride is null ? " (auto)" : " (override)")}  |  max gap: {maxGap}  |  {ppm} ppm  |  mass tol: {massTol} Da");

List<ProcessingRow> rows;
try
{
    rows = Pipeline.ProcessDmt(
        inputPath,
        chargeRange:  (chargeMin, chargeMax),
        binWidth:     binWidth,
        maxGap:       maxGap,
        ppmTolerance: ppm);
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
    "Mz,Charge,Intensity,NeutralMass,IsotopicIndex,IsEstimated,IsSuppressed," +
    "EstimatedIntensity,CorrectedIntensity,Uncertainty,Confidence,Hypothesis," +
    "FitRSquared,GapAtApex,ClusterIonCount,CorrectedIonCount,Notes";

static void WriteRow(StreamWriter w, ProcessingRow row)
{
    w.WriteLine(
        $"{F(row.Mz)},{row.Charge},{F(row.Intensity)},{row.NeutralMass:G6}," +
        $"{row.IsotopicIndex},{S(row.IsEstimated)},{S(row.IsSuppressed)}," +
        $"{F(row.EstimatedIntensity)},{F(row.CorrectedIntensity)}," +
        $"{F(row.Uncertainty)},{row.Confidence:G4},{row.Hypothesis}," +
        $"{row.FitRSquared:G4},{S(row.GapAtApex)}," +
        $"{row.ClusterIonCount:G6},{row.CorrectedIonCount:G6},{Q(row.Notes)}");
}

// All peaks: CorrectedIntensity populated for suppressed peaks, null otherwise.
using (var w = new StreamWriter(outputPath))
{
    w.WriteLine(CsvHeader);
    foreach (var row in rows)
        WriteRow(w, row);
}

// ---- Write proteoform summary CSV ------------------------------------------
var allProteoforms = ProteoformAggregator.Aggregate(rows, massTol);

// Noise filter: drop proteoforms below 0.5% of the largest, or 100 ions — whichever is higher.
double maxIons       = allProteoforms.Count > 0 ? allProteoforms.Max(pf => pf.TotalCorrectedIonCount) : 0;
double noiseFloor    = Math.Max(maxIons * 0.005, 100.0);
var    proteoforms   = allProteoforms.Where(pf => pf.TotalCorrectedIonCount >= noiseFloor).ToList();
int    noisyRemoved  = allProteoforms.Count - proteoforms.Count;

const string PfHeader =
    "NeutralMass,ChargeStates,ClusterCount,TotalObservedIonCount,TotalCorrectedIonCount," +
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
      --charge-min N    Override minimum charge state (default: auto-detected from file)
      --charge-max N    Override maximum charge state (default: auto-detected from file)
      --max-gap N       Max consecutive missing isotope peaks per cluster (default 3)
      --ppm N           m/z matching tolerance in ppm; Orbitrap: 5–10, Q-TOF: 10–20 (default 10)
      --bin-width F     m/z bin width in Da for ion-event aggregation (default 0.02)
      --mass-tol F      Neutral-mass tolerance in Da for proteoform grouping (default 1.0)
      -h, --help        Show this help

    Output files written next to the input:
      <input>.corrected.csv    All observed peaks; CorrectedIntensity column filled for suppressed peaks
      <input>.proteoforms.csv  One row per proteoform, ion counts summed across charge states
      <input>.spectral_check.html  SVG charts showing suppression correction per cluster
    """);
