namespace BrainMie.Core;

using Microsoft.Data.Sqlite;

/// <summary>
/// Reads ion detection events from the Ion table of a .dmt SQLite file.
/// </summary>
public static class DmtReader
{
    private static SqliteConnection OpenReadOnly(string filePath)
    {
        // Use the SQLite URI format with immutable=1 so the library never tries to
        // create or lock WAL/journal files alongside the database. This avoids
        // SQLITE_CANTOPEN (error 14) when the file is on a restricted path or is
        // still held open by the acquisition software.
        string uri = "file:" + filePath.Replace('\\', '/').Replace(" ", "%20") + "?mode=ro&immutable=1";
        var conn = new SqliteConnection($"Data Source={uri}");
        conn.Open();
        return conn;
    }

    /// <summary>
    /// Open <paramref name="filePath"/> and return every row from the Ion table.
    /// </summary>
    /// <returns>
    /// A list of (Mz, Charge, Intensity?) tuples. Intensity is <see langword="null"/>
    /// when the Ion table has no Intensity column; in that case the caller should bin
    /// the raw m/z values to derive peak heights by counting.
    /// </returns>
    /// <exception cref="InvalidOperationException">
    /// Thrown when the Ion table is missing the required Mz or Charge column.
    /// </exception>
    public static List<(double Mz, int Charge, double? Intensity)> ReadIons(string filePath)
    {
        var rows = new List<(double Mz, int Charge, double? Intensity)>();

        using var conn = OpenReadOnly(filePath);

        var columnNames = GetColumnNames(conn, "Ion");

        foreach (string required in new[] { "Mz", "Charge" })
        {
            if (!columnNames.Contains(required, StringComparer.OrdinalIgnoreCase))
                throw new InvalidOperationException(
                    $"Ion table is missing required column '{required}'. " +
                    $"Found: {string.Join(", ", columnNames)}");
        }

        bool hasIntensity = columnNames.Contains("Intensity", StringComparer.OrdinalIgnoreCase);
        string query = hasIntensity
            ? "SELECT Mz, Charge, Intensity FROM Ion"
            : "SELECT Mz, Charge FROM Ion";

        using var cmd    = new SqliteCommand(query, conn);
        using var reader = cmd.ExecuteReader();

        while (reader.Read())
        {
            double  mz        = reader.GetDouble(0);
            int     charge    = reader.GetInt32(1);
            double? intensity = hasIntensity ? reader.GetDouble(2) : null;
            rows.Add((mz, charge, intensity));
        }

        return rows;
    }

    private static List<string> GetColumnNames(SqliteConnection conn, string table)
    {
        var names = new List<string>();
        using var cmd    = new SqliteCommand($"PRAGMA table_info({table})", conn);
        using var reader = cmd.ExecuteReader();
        while (reader.Read())
            names.Add(reader.GetString(1));
        return names;
    }
}
