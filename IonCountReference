using System.Text.RegularExpressions;
using ProteoformAnalyzer;

Console.WriteLine("╔══════════════════════════════════════════════╗");
Console.WriteLine("║       Proteoform Analyzer  v2.0              ║");
Console.WriteLine("║  UniProt · PRIDE · PTMeXchange/EBI           ║");
Console.WriteLine("╚══════════════════════════════════════════════╝");
Console.WriteLine();

using var http = new HttpClient();
http.DefaultRequestHeaders.Add("User-Agent", "ProteoformAnalyzer/2.0 (research tool)");
http.DefaultRequestHeaders.Add("Accept", "application/json");
http.Timeout = TimeSpan.FromSeconds(30);

// ── Mode selection ────────────────────────────────────────────────────────
Console.WriteLine("Select mode:");
Console.WriteLine("  1  Interactive  — enter one protein sequence or UniProt ID");
Console.WriteLine("  2  Batch CSV    — provide a CSV file containing multiple proteins");
Console.Write("> ");
string? modeInput = Console.ReadLine()?.Trim();
Console.WriteLine();

if (modeInput == "2")
{
    // ── Batch CSV options ─────────────────────────────────────────────────
    string batchCsvPath;
    while (true)
    {
        Console.Write("Path to input CSV file: ");
        batchCsvPath = Console.ReadLine()?.Trim().Trim('"') ?? "";
        if (File.Exists(batchCsvPath)) break;
        Console.WriteLine($"  File not found: {batchCsvPath}");
    }

    Console.Write("Include N- and C-terminal truncations? (y/n, default y): ");
    bool batchTrunc = (Console.ReadLine()?.Trim().ToLower() ?? "y") != "n";

    Console.Write("Estimate false-discovery rate with decoys? (y/n, default y): ");
    bool batchEstimateFdr = (Console.ReadLine()?.Trim().ToLower() ?? "y") != "n";

    int batchMaxOccupancy = 12;
    Console.Write("Max simultaneous modifications per type (default 12): ");
    string? batchMaxOccStr = Console.ReadLine()?.Trim();
    if (!string.IsNullOrEmpty(batchMaxOccStr) &&
        int.TryParse(batchMaxOccStr, out int bmo) && bmo >= 1)
        batchMaxOccupancy = bmo;

    double batchMatchTol = 2.0;
    Console.Write("Match tolerance in Da — max distance from peak centroid to database mass (default 2.0): ");
    string? batchMatchTolStr = Console.ReadLine()?.Trim();
    if (!string.IsNullOrEmpty(batchMatchTolStr) &&
        double.TryParse(batchMatchTolStr, System.Globalization.NumberStyles.Any,
            System.Globalization.CultureInfo.InvariantCulture, out double bmt) && bmt > 0)
        batchMatchTol = bmt;

    double batchIonWindow = 5.0;
    Console.Write("Ion-counting window in Da — signal summed within ±window of database mass (default 5.0): ");
    string? batchIonWinStr = Console.ReadLine()?.Trim();
    if (!string.IsNullOrEmpty(batchIonWinStr) &&
        double.TryParse(batchIonWinStr, System.Globalization.NumberStyles.Any,
            System.Globalization.CultureInfo.InvariantCulture, out double biw) && biw > 0)
        batchIonWindow = biw;

    Console.Write("Path to folder containing .dmt files (press Enter to skip): ");
    string? batchDmtFolder = Console.ReadLine()?.Trim().Trim('"');
    if (!string.IsNullOrEmpty(batchDmtFolder) && !Directory.Exists(batchDmtFolder))
    {
        Console.WriteLine("  Folder not found — skipping ion counting.");
        batchDmtFolder = null;
    }

    Console.Write("Check for contaminating proteins? (y/n, default n): ");
    List<string>? batchContaminants = null;
    if ((Console.ReadLine()?.Trim().ToLower() ?? "n") == "y")
    {
        batchContaminants = new List<string>();
        while (true)
        {
            Console.Write("  Contaminant UniProt accession (press Enter to finish): ");
            string? cId = Console.ReadLine()?.Trim().ToUpper();
            if (string.IsNullOrEmpty(cId)) break;
            batchContaminants.Add(cId);
        }
    }

    Console.WriteLine();
    await CsvBatchMode.RunAsync(batchCsvPath, http, batchTrunc, batchMatchTol, batchIonWindow, batchDmtFolder, batchContaminants, batchMaxOccupancy, batchEstimateFdr);
    return;
}

// ── 1. Sequence or UniProt ID (interactive mode) ──────────────────────────
string sequence = "";
string? uniprotId = null;

Console.WriteLine("Enter a protein sequence (single-letter codes) OR a UniProt accession (e.g. P04637).");
while (true)
{
    Console.Write("> ");
    string? raw = Console.ReadLine()?.Trim();
    if (string.IsNullOrEmpty(raw)) continue;

    if (Regex.IsMatch(raw, @"^[A-Z][0-9][A-Z0-9]{3}[0-9]$", RegexOptions.IgnoreCase) ||
        Regex.IsMatch(raw, @"^[OPQ][0-9][A-Z0-9]{3}[0-9]$", RegexOptions.IgnoreCase))
    {
        // Looks like a UniProt accession
        uniprotId = raw.ToUpper();
        break;
    }

    string candidate = raw.ToUpper().Replace(" ", "").Replace("\t", "");
    var invalid = candidate.Where(c => !AminoAcidData.ResidueFormulas.ContainsKey(c)).Distinct().ToList();
    if (invalid.Count > 0)
    {
        Console.WriteLine($"  Unknown character(s): {string.Join(", ", invalid)}");
        Console.WriteLine("  Enter a valid amino acid sequence or a UniProt accession.");
        continue;
    }

    sequence = candidate;
    break;
}

// ── 2. Fetch from UniProt / databases if accession was given ──────────────
var allPtms = new List<PtmAnnotation>();

if (uniprotId is not null)
{
    Console.WriteLine();
    var uniprotClient = new UniProtClient(http);
    (sequence, var uniprotPtms) = await uniprotClient.FetchAsync(uniprotId);

    if (string.IsNullOrEmpty(sequence))
    {
        Console.WriteLine("Could not retrieve sequence from UniProt. Exiting.");
        return;
    }

    allPtms.AddRange(uniprotPtms);

    var prideClient = new PrideClient(http);
    var pridePtms = await prideClient.FetchAsync(uniprotId, sequence);
    allPtms.AddRange(pridePtms);

    var ptmExClient = new PtmExchangeClient(http);
    var ptmExPtms = await ptmExClient.FetchAsync(uniprotId, sequence);
    allPtms.AddRange(ptmExPtms);
}
else
{
    Console.WriteLine($"Sequence provided directly ({sequence.Length} aa). Skipping database queries.");
}

// ── 3. Summary ────────────────────────────────────────────────────────────
var intactFormula = AminoAcidData.GetFormula(sequence);
double intactMass = AminoAcidData.AverageMass(intactFormula);

Console.WriteLine();
Console.WriteLine($"Sequence  : {(sequence.Length <= 60 ? sequence : sequence[..57] + "...")}");
Console.WriteLine($"Length    : {sequence.Length} aa");
Console.WriteLine($"Formula   : {intactFormula}");
Console.WriteLine($"Avg mass  : {intactMass:F4} Da  (isotope envelope centroid)");
Console.WriteLine($"DB PTMs   : {allPtms.Count} annotations collected");

// ── 4. Custom modifications ───────────────────────────────────────────────
Console.WriteLine();
while (true)
{
    Console.Write("Add a custom modification? (y/n): ");
    if ((Console.ReadLine()?.Trim().ToLower() ?? "n") != "y") break;

    string modName;
    while (true)
    {
        Console.Write("  Modification name: ");
        modName = Console.ReadLine()?.Trim() ?? "";
        if (!string.IsNullOrEmpty(modName)) break;
    }

    double delta = 0;
    while (true)
    {
        Console.Write("  Mass change in Da (e.g. +79.966 or -18.011): ");
        string? d = Console.ReadLine()?.Trim();
        if (double.TryParse(d, System.Globalization.NumberStyles.Any,
                System.Globalization.CultureInfo.InvariantCulture, out delta)) break;
        Console.WriteLine("  Please enter a valid decimal number.");
    }

    int pos = 0;
    Console.Write("  Sequence position (1-based, 0 = whole-protein / unknown): ");
    int.TryParse(Console.ReadLine()?.Trim(), out pos);

    char? residue = pos > 0 && pos <= sequence.Length ? char.ToUpper(sequence[pos - 1]) : null;

    allPtms.Add(new PtmAnnotation
    {
        ModificationName = modName,
        Position = pos,
        Residue = residue,
        MassDelta = delta,
        Source = "Custom"
    });

    Console.WriteLine($"  Added: {modName} @ pos {(pos > 0 ? pos.ToString() : "N/A")} ({(delta >= 0 ? "+" : "")}{delta:F4} Da)");
}

// ── 5. Options ────────────────────────────────────────────────────────────
Console.WriteLine();
Console.Write("Include N- and C-terminal truncations? (y/n, default y): ");
bool includeTrunc = (Console.ReadLine()?.Trim().ToLower() ?? "y") != "n";

int maxOccupancy = 12;
Console.Write("Max simultaneous modifications per type (default 12, e.g. limits hyperphosphorylated proteins): ");
string? maxOccStr = Console.ReadLine()?.Trim();
if (!string.IsNullOrEmpty(maxOccStr) &&
    int.TryParse(maxOccStr, out int mo) && mo >= 1)
    maxOccupancy = mo;

// ── 6. Build proteoform database ──────────────────────────────────────────
Console.WriteLine();
Console.WriteLine($"Building proteoform database (truncations: {includeTrunc}, max occupancy per PTM type: {maxOccupancy})...");
var proteoforms = ProteoformBuilder.Build(sequence, allPtms, includeTrunc, tolerance: 5.0,
                                         proteinLabel: uniprotId ?? "Target",
                                         maxOccupancyPerFamily: maxOccupancy);
Console.WriteLine($"Generated {proteoforms.Count} database entries.");

// ── 6b. Contaminant proteins ──────────────────────────────────────────────
// If the user wants to check for contaminating proteins, fetch their databases
// and merge them into the main list.  Target entries are labelled so results
// from each protein are distinguishable in the output.
Console.WriteLine();
Console.Write("Check for contaminating proteins? (y/n, default n): ");
bool searchContaminants = (Console.ReadLine()?.Trim().ToLower() ?? "n") == "y";
if (searchContaminants)
{
    var contUniProtClient = new UniProtClient(http);
    var contPrideClient   = new PrideClient(http);
    var contPtmExClient   = new PtmExchangeClient(http);

    while (true)
    {
        Console.Write("  Contaminant UniProt accession (press Enter to finish): ");
        string? contId = Console.ReadLine()?.Trim().ToUpper();
        if (string.IsNullOrEmpty(contId)) break;

        Console.Write($"  Fetching {contId} ... ");
        var (contSeq, contUniProtPtms) = await contUniProtClient.FetchAsync(contId);
        if (string.IsNullOrEmpty(contSeq)) { Console.WriteLine("not found — skipping."); continue; }

        var contPridePtms  = await contPrideClient.FetchAsync(contId, contSeq);
        var contPtmExPtms  = await contPtmExClient.FetchAsync(contId, contSeq);
        var contAllPtms    = contUniProtPtms.Concat(contPridePtms).Concat(contPtmExPtms).ToList();

        var contEntries = ProteoformBuilder.Build(contSeq, contAllPtms, includeTrunc, tolerance: 5.0,
                                                   proteinLabel: contId,
                                                   maxOccupancyPerFamily: maxOccupancy);
        proteoforms.AddRange(contEntries);
        Console.WriteLine($"done — {contEntries.Count} entries added ({contAllPtms.Count} PTM annotations).");
    }
}

// ── 7. .dmt spectrum analysis ─────────────────────────────────────────────
Console.WriteLine();
Console.Write("Path to folder containing .dmt files (press Enter to skip): ");
string? dmtFolder = Console.ReadLine()?.Trim().Trim('"');

double matchTol = 2.0;
double ionWindow = 5.0;
List<string> dmtFileNames = new();
List<AnalysisResult> analysisResults = new();

if (!string.IsNullOrEmpty(dmtFolder))
{
    if (!Directory.Exists(dmtFolder))
    {
        Console.WriteLine("  Folder not found — skipping spectrum analysis.");
    }
    else
    {
        Console.Write("Match tolerance in Da — max distance from peak centroid to database mass (default 2.0): ");
        string? matchTolStr = Console.ReadLine()?.Trim();
        if (!string.IsNullOrEmpty(matchTolStr) &&
            double.TryParse(matchTolStr, System.Globalization.NumberStyles.Any,
                System.Globalization.CultureInfo.InvariantCulture, out double mt) && mt > 0)
            matchTol = mt;

        Console.Write("Ion-counting window in Da — signal summed within ±window of database mass (default 5.0): ");
        string? ionWinStr = Console.ReadLine()?.Trim();
        if (!string.IsNullOrEmpty(ionWinStr) &&
            double.TryParse(ionWinStr, System.Globalization.NumberStyles.Any,
                System.Globalization.CultureInfo.InvariantCulture, out double iw) && iw > 0)
            ionWindow = iw;

        Console.Write("Estimate false-discovery rate with decoys? (y/n, default y): ");
        bool estimateFdr = (Console.ReadLine()?.Trim().ToLower() ?? "y") != "n";

        // Add decoy proteoforms so chance matches can be quantified (target–decoy FDR).
        var searchDb = proteoforms;
        if (estimateFdr)
        {
            var decoys = DecoyGenerator.Generate(proteoforms);
            searchDb = proteoforms.Concat(decoys).ToList();
            Console.WriteLine($"  Added {decoys.Count} decoy proteoforms for FDR estimation.");
        }

        RunManifest.RecordParam("matchTolerance", matchTol.ToString("F2"));
        RunManifest.RecordParam("ionWindow", ionWindow.ToString("F2"));
        RunManifest.RecordParam("truncations", includeTrunc.ToString());
        RunManifest.RecordParam("maxOccupancyPerFamily", maxOccupancy.ToString());
        RunManifest.RecordParam("estimateFdr", estimateFdr.ToString());

        // Unmatched peaks are only reported when the user is screening for contaminants.
        var (results, fileNames) = SpectrumBatch.AnalyzeFolder(
            dmtFolder, searchDb, matchTol, ionWindow,
            includeUnmatchedPeaks: searchContaminants);

        if (estimateFdr)
        {
            var fdr = FdrEstimator.Estimate(results);
            Console.WriteLine($"  {fdr.Describe()}");
            RunManifest.Record($"FDR: {fdr.Describe()}");
            results = results.Where(r => !r.DatabaseEntry.IsDecoy).ToList();
        }

        analysisResults = results;
        dmtFileNames    = fileNames;
    }
}

// ── 8. Export CSV ─────────────────────────────────────────────────────────
string defaultCsv = uniprotId is not null ? $"{uniprotId}_proteoforms.csv" : "proteoforms.csv";

string csvPath;
while (true)
{
    Console.Write($"\nOutput CSV path (press Enter for default: {defaultCsv}): ");
    string? raw = Console.ReadLine()?.Trim().Trim('"');
    csvPath = string.IsNullOrEmpty(raw) ? defaultCsv : raw;

    if (!Path.GetFileName(csvPath).Contains('.'))
        csvPath += ".csv";

    string? dir = Path.GetDirectoryName(csvPath);
    if (!string.IsNullOrEmpty(dir) && !Directory.Exists(dir))
    {
        Console.WriteLine($"  Directory not found: {dir}");
        continue;
    }

    try
    {
        if (analysisResults.Count > 0)
            CsvExporter.ExportResults(analysisResults, dmtFileNames, csvPath);
        else
            CsvExporter.ExportDatabase(proteoforms, csvPath);

        Console.WriteLine($"Saved: {Path.GetFullPath(csvPath)}");
        RunManifest.Write(csvPath);
        Console.WriteLine($"Run manifest: {Path.GetFullPath(csvPath)}.manifest.txt");
        break;
    }
    catch (Exception ex)
    {
        Console.WriteLine($"  Could not save: {ex.Message}");
    }
}

// ── 9. Console preview ────────────────────────────────────────────────────
Console.WriteLine();
if (analysisResults.Count > 0)
{
    Console.WriteLine($"{"Modification",-50} {"Pred. Mass",14}  {"Expt. Centroid",15}  {"Total Ions",12}");
    Console.WriteLine(new string('─', 97));
    foreach (var r in analysisResults.Take(15))
    {
        long total = r.IonCountsPerFile.Values.Sum();
        Console.WriteLine($"{r.DatabaseEntry.ModificationName,-50} " +
                          $"{r.DatabaseEntry.CentroidMass,14:F4}  " +
                          $"{r.ExperimentalCentroid,15:F4}  {total,12:N0}");
    }
    if (analysisResults.Count > 15)
        Console.WriteLine($"  ... and {analysisResults.Count - 15} more (see CSV).");
}
else
{
    Console.WriteLine($"{"Modification",-55} {"Pred. Mass (Da)",16}");
    Console.WriteLine(new string('─', 73));
    foreach (var pf in proteoforms.Take(15))
        Console.WriteLine($"{pf.ModificationName,-55} {pf.CentroidMass,16:F4}");
    if (proteoforms.Count > 15)
        Console.WriteLine($"  ... and {proteoforms.Count - 15} more (see CSV).");
}

Console.WriteLine("\nDone.");
