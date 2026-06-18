namespace BrainMie.Core;

using Microsoft.Data.Sqlite;

/// <summary>
/// Reads ion detection events from the Ion table of a .dmt SQLite file.
/// </summary>
public static class DmtReader
{
    private static readonly string[] OptionalColumns = ["Intensity", "RetentionTime", "ScanNumber"];

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

        var builder = new SqliteConnectionStringBuilder
        {
            DataSource = filePath,
            Mode = SqliteOpenMode.ReadOnly,
        };

        using var conn = new SqliteConnection(builder.ToString());
        conn.Open();

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

    /// <summary>
    /// Return the (Min, Max) distinct charge states present in the Ion table.
    /// </summary>
    public static (int Min, int Max) ReadChargeRange(string filePath)
    {
        var builder = new SqliteConnectionStringBuilder
        {
            DataSource = filePath,
            Mode = SqliteOpenMode.ReadOnly,
        };

        using var conn = new SqliteConnection(builder.ToString());
        conn.Open();

        using var cmd    = new SqliteCommand("SELECT MIN(Charge), MAX(Charge) FROM Ion", conn);
        using var reader = cmd.ExecuteReader();
        if (reader.Read() && !reader.IsDBNull(0))
            return (reader.GetInt32(0), reader.GetInt32(1));

        return (1, 10);
    }

    private static IReadOnlyList<string> GetColumnNames(SqliteConnection conn, string table)
    {
        var names = new List<string>();
        using var cmd    = new SqliteCommand($"PRAGMA table_info({table})", conn);
        using var reader = cmd.ExecuteReader();
        while (reader.Read())
            names.Add(reader.GetString(1));
        return names;
    }
}
