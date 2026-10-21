"""Loading and cleaning a USGS earthquake catalog (CSV export)."""
from __future__ import annotations

import pandas as pd

ESSENTIAL = ["time", "latitude", "longitude", "depth", "mag"]
OPTIONAL = ["place", "type", "magType", "status", "id"]


def load_catalog(path, year_min=None, year_max=None):
    """Read a USGS CSV and return (clean_df, cleaning_report).

    Every filtering step is logged so the cleaning is auditable: the report
    lists how many rows survive each step.
    """
    header = pd.read_csv(path, nrows=0).columns
    missing = [c for c in ESSENTIAL if c not in header]
    if missing:
        raise ValueError(f"Catalog is missing required columns: {missing}")

    df = pd.read_csv(path, usecols=[c for c in ESSENTIAL + OPTIONAL if c in header])
    steps = [("Rows in raw file", len(df))]

    df["time"] = pd.to_datetime(df["time"], errors="coerce", utc=True)
    df = df.dropna(subset=ESSENTIAL)
    steps.append(("Drop rows missing time, location, depth or magnitude", len(df)))

    if "type" in df.columns:
        df = df[df["type"].eq("earthquake")]
        steps.append(("Keep type == 'earthquake' (drop blasts, quarries, etc.)", len(df)))

    if "id" in df.columns:
        df = df.drop_duplicates("id")
        steps.append(("Remove duplicate event IDs", len(df)))

    df = df[df["latitude"].between(-90, 90) & df["longitude"].between(-180, 180)]
    steps.append(("Keep valid coordinates", len(df)))

    if year_min is not None:
        df = df[df["time"].dt.year >= year_min]
    if year_max is not None:
        df = df[df["time"].dt.year <= year_max]
    if year_min is not None or year_max is not None:
        steps.append((f"Restrict to years {year_min or '-'}-{year_max or '-'}", len(df)))

    df = df.rename(columns={
        "time": "Time", "latitude": "Latitude", "longitude": "Longitude",
        "depth": "Depth_km", "mag": "Magnitude", "place": "Location",
    })
    # USGS reports a few events above sea level (negative depth); clamp to the surface.
    df["Depth_km"] = df["Depth_km"].clip(lower=0)
    df["Year"] = df["Time"].dt.year
    if "Location" in df.columns:
        df["Region"] = (df["Location"].fillna("Unknown").str.split(",").str[-1]
                        .str.strip().str.replace(r"\s+region$", "", regex=True, case=False))
    else:
        df["Region"] = "Unknown"

    df = df.sort_values("Time").reset_index(drop=True)

    report = pd.DataFrame(steps, columns=["step", "rows_remaining"])
    report["rows_removed"] = (-report["rows_remaining"].diff()).fillna(0).astype(int)
    return df, report
