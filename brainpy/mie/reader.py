import sqlite3
import pandas as pd


# Columns to use if present beyond the required Mz / Charge pair.
_OPTIONAL_COLUMNS = ("Intensity", "RetentionTime", "ScanNumber")


def read_ions(filepath: str) -> pd.DataFrame:
    """
    Read the Ion table from a .dmt SQLite file.

    Always returns a DataFrame with at least Mz (float) and Charge (int).
    Optional columns (Intensity, RetentionTime, ScanNumber) are included when
    the table contains them.

    Each row represents one ion detection event.  If the file has no Intensity
    column the caller is responsible for deriving peak heights by counting ions
    that fall in the same m/z bin.
    """
    conn = sqlite3.connect(filepath)
    try:
        cursor = conn.cursor()
        cursor.execute("PRAGMA table_info(Ion)")
        available = {row[1] for row in cursor.fetchall()}

        missing = {"Mz", "Charge"} - available
        if missing:
            raise ValueError(
                f"Ion table is missing required columns: {missing}. "
                f"Found: {sorted(available)}"
            )

        select = ["Mz", "Charge"] + [c for c in _OPTIONAL_COLUMNS if c in available]
        df = pd.read_sql_query(f"SELECT {', '.join(select)} FROM Ion", conn)
    finally:
        conn.close()

    df["Mz"] = df["Mz"].astype(float)
    df["Charge"] = df["Charge"].astype(int)
    return df
