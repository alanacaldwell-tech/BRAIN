using BrainMie.Core;
using BrainMie.Core.Data;

// ---- Argument parsing -------------------------------------------------------
if (args.Length < 1 || args[0] is "-h" or "--help")
{
    Console.Error.WriteLine("""
        BrainMie — Multi-Ion Event correction for .dmt mass spectrometry files

        Usage:
          BrainMie <input.dmt> [options]

        Options:
          --charge-min N    Minimum charge state to process (default 1)
          --charge-max N    Maximum charge state to process (default 10)
          --max-gap N       Max consecutive missing isotope peaks per cluster (default 3)
          --ppm N           m/z matching tolerance in ppm (default 20)
          --bin-width F     m/z bin width in Da for ion-event aggregation (default 0.005)
          --output PATH     Output CSV path (default: <input>.corrected.csv)
          -h, --help        Show this help
        """);
    return args.Length < 1 ? 1 : 0;
}

string inputPath  = args[0];
int    chargeMin  = 1;
int    chargeMax  = 10;
int    maxGap     = 3;
double ppm        = 20.0;
double binWidth   = 0.005;
string outputPath = Path.ChangeExtension(inputPath, ".corrected.csv");

for (int i = 1; i < args.Length - 1; i++)
{
    switch (args[i].ToLowerInvariant())
    {
        case "--charge-min": chargeMin  = int.Parse(args[++i]);    break;
        case "--charge-max": chargeMax  = int.Parse(args[++i]);    break;
        case "--max-gap":    maxGap     = int.Parse(args[++i]);    break;
        case "--ppm":        ppm        = double.Parse(args[++i]); break;
        case "--bin-width":  binWidth   = double.Parse(args[++i]); break;
        case "--output":     outputPath = args[++i];               break;
    }
}

// ---- Run pipeline -----------------------------------------------------------
Console.WriteLine($"Input:   {inputPath}");
Console.WriteLine($"Charges: {chargeMin}–{chargeMax}  |  max gap: {maxGap}  |  {ppm} ppm");

List<ProcessingRow> rows;
try
{
    rows = Pipeline.ProcessDmt(
        inputPath,
        chargeRange:      (chargeMin, chargeMax),
        binWidth:         binWidth,
        maxGap:           maxGap,
        ppmTolerance:     ppm);
}
catch (Exception ex)
{
    Console.Error.WriteLine($"Error: {ex.Message}");
    return 1;
}

// ---- Write CSV output -------------------------------------------------------
using var writer = new StreamWriter(outputPath);

writer.WriteLine(
    "Mz,Charge,Intensity,NeutralMass,IsotopicIndex,IsEstimated," +
    "EstimatedIntensity,Uncertainty,Confidence,Hypothesis,FitRSquared,GapAtApex,Notes");

static string F(double? v)  => v.HasValue ? v.Value.ToString("G6") : "";
static string S(bool    b)  => b ? "true" : "false";
static string Q(string  s)  => $"\"{s.Replace("\"", "\"\"")}\"";

foreach (var row in rows)
{
    writer.WriteLine(
        $"{F(row.Mz)},{row.Charge},{F(row.Intensity)},{row.NeutralMass:G6}," +
        $"{row.IsotopicIndex},{S(row.IsEstimated)},{F(row.EstimatedIntensity)}," +
        $"{F(row.Uncertainty)},{row.Confidence:G4},{row.Hypothesis}," +
        $"{row.FitRSquared:G4},{S(row.GapAtApex)},{Q(row.Notes)}");
}

// ---- Summary ----------------------------------------------------------------
int totalRows  = rows.Count;
int estimated  = rows.Count(r => r.IsEstimated);
int mie        = rows.Count(r => r.Hypothesis == "mie" && !r.IsEstimated);
int overlap    = rows.Count(r => r.Hypothesis == "overlap" && !r.IsEstimated);
int ambiguous  = rows.Count(r => r.Hypothesis == "ambiguous" && !r.IsEstimated);

Console.WriteLine($"Output:  {outputPath}");
Console.WriteLine($"Rows:    {totalRows} total  ({estimated} estimated/rescued peaks)");
Console.WriteLine($"         {mie} MIE peaks  |  {overlap} overlap peaks  |  {ambiguous} ambiguous");

return 0;
