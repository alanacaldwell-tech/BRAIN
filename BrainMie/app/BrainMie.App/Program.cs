using BrainMie.Core;
using BrainMie.Core.Data;

// ---- Argument parsing -------------------------------------------------------
if (args.Length > 0 && args[0] is "-h" or "--help")
{
    PrintHelp();
    return 0;
}

string inputPath;
int    chargeMin  = 1;
int    chargeMax  = 10;
int    maxGap     = 3;
double ppm        = 20.0;
double binWidth   = 0.005;
double massTol    = 1.0;

if (args.Length >= 1)
{
    // Command-line mode: BrainMie.exe input.dmt [options]
    inputPath = args[0];
    for (int i = 1; i < args.Length - 1; i++)
    {
        switch (args[i].ToLowerInvariant())
        {
            case "--charge-min": chargeMin = int.Parse(args[++i]);    break;
            case "--charge-max": chargeMax = int.Parse(args[++i]);    break;
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
    raw = Prompt($"Minimum charge state [{chargeMin}]");
    if (!string.IsNullOrWhiteSpace(raw)) chargeMin = int.Parse(raw);

    raw = Prompt($"Maximum charge state [{chargeMax}]");
    if (!string.IsNullOrWhiteSpace(raw)) chargeMax = int.Parse(raw);

    raw = Prompt($"Max consecutive missing isotope peaks [{maxGap}]");
    if (!string.IsNullOrWhiteSpace(raw)) maxGap = int.Parse(raw);

    raw = Prompt($"m/z tolerance in ppm [{ppm}]");
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
string unrescuedPath   = Path.ChangeExtension(inputPath, ".unrescued.csv");
string proteoformPath  = Path.ChangeExtension(inputPath, ".proteoforms.csv");

// ---- Run pipeline -----------------------------------------------------------
Console.WriteLine($"Input:   {inputPath}");
Console.WriteLine($"Charges: {chargeMin}–{chargeMax}  |  max gap: {maxGap}  |  {ppm} ppm  |  mass tol: {massTol} Da");

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
    "Mz,Charge,Intensity,NeutralMass,IsotopicIndex,IsEstimated," +
    "EstimatedIntensity,Uncertainty,Confidence,Hypothesis,FitRSquared,GapAtApex," +
    "ClusterIonCount,CorrectedIonCount,Notes";

static void WriteRow(StreamWriter w, ProcessingRow row)
{
    w.WriteLine(
        $"{F(row.Mz)},{row.Charge},{F(row.Intensity)},{row.NeutralMass:G6}," +
        $"{row.IsotopicIndex},{S(row.IsEstimated)},{F(row.EstimatedIntensity)}," +
        $"{F(row.Uncertainty)},{row.Confidence:G4},{row.Hypothesis}," +
        $"{row.FitRSquared:G4},{S(row.GapAtApex)}," +
        $"{row.ClusterIonCount:G6},{row.CorrectedIonCount:G6},{Q(row.Notes)}");
}

// Unrescued: observed peaks only (no estimated peaks injected)
using (var w = new StreamWriter(unrescuedPath))
{
    w.WriteLine(CsvHeader);
    foreach (var row in rows.Where(r => !r.IsEstimated))
        WriteRow(w, row);
}

// Rescued: observed peaks + reconstructed missing peaks
using (var w = new StreamWriter(outputPath))
{
    w.WriteLine(CsvHeader);
    foreach (var row in rows)
        WriteRow(w, row);
}

// ---- Write proteoform summary CSV ------------------------------------------
var proteoforms = ProteoformAggregator.Aggregate(rows, massTol);

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

// ---- Summary ----------------------------------------------------------------
int totalRows = rows.Count;
int estimated = rows.Count(r => r.IsEstimated);
int mie       = rows.Count(r => r.Hypothesis == "mie"       && !r.IsEstimated);
int overlap   = rows.Count(r => r.Hypothesis == "overlap"   && !r.IsEstimated);
int ambiguous = rows.Count(r => r.Hypothesis == "ambiguous" && !r.IsEstimated);

Console.WriteLine($"Unrescued:    {unrescuedPath}  ({totalRows - estimated} peaks)");
Console.WriteLine($"Rescued:      {outputPath}  ({totalRows} peaks, {estimated} estimated)");
Console.WriteLine($"Proteoforms:  {proteoformPath}  ({proteoforms.Count} entries, {massTol} Da grouping)");
Console.WriteLine($"              {mie} MIE  |  {overlap} overlap  |  {ambiguous} ambiguous");

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
      --charge-min N    Minimum charge state to process (default 1)
      --charge-max N    Maximum charge state to process (default 10)
      --max-gap N       Max consecutive missing isotope peaks per cluster (default 3)
      --ppm N           m/z matching tolerance in ppm (default 20)
      --bin-width F     m/z bin width in Da for ion-event aggregation (default 0.005)
      --mass-tol F      Neutral-mass tolerance in Da for proteoform grouping (default 1.0)
      -h, --help        Show this help

    Output files written next to the input:
      <input>.unrescued.csv    Observed peaks only, before MIE correction
      <input>.corrected.csv    Observed + reconstructed MIE gap peaks
      <input>.proteoforms.csv  One row per proteoform, ion counts summed across charge states
    """);
