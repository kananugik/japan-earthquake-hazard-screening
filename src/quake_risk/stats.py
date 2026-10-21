"""Catalog statistics: completeness, Gutenberg-Richter b-value, declustering."""
from __future__ import annotations

import numpy as np
import pandas as pd

EARTH_R_KM = 6371.0088


def haversine_km(lat1, lon1, lat2, lon2):
    """Great-circle distance in km (inputs in degrees, broadcastable)."""
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dphi = p2 - p1
    dlmb = np.radians(lon2) - np.radians(lon1)
    a = np.sin(dphi / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dlmb / 2) ** 2
    return 2 * EARTH_R_KM * np.arcsin(np.sqrt(a))


def magnitude_of_completeness(mags, bin_width=0.1):
    """Mc by the maximum-curvature method (+0.2 correction, Woessner & Wiemer 2005).

    Mc is the magnitude above which the catalog is assumed to record every event.
    """
    binned = np.round(np.round(np.asarray(mags) / bin_width) * bin_width, 2)
    values, counts = np.unique(binned, return_counts=True)
    return float(np.round(values[counts.argmax()] + 0.2, 1))


def gutenberg_richter(mags, mc, dm=0.1):
    """Aki-Utsu maximum-likelihood b-value with binning correction.

    log10 N(>=M) = a - b*M. Returns b, its standard error (b/sqrt(n)), a and n.
    """
    m = np.asarray(mags)
    m = m[m >= mc - 1e-9]
    n = len(m)
    if n < 2:
        raise ValueError("Not enough events above Mc to estimate b")
    b = np.log10(np.e) / (m.mean() - (mc - dm / 2))
    return {"b": float(b), "b_se": float(b / np.sqrt(n)),
            "a": float(np.log10(n) + b * mc), "n": int(n), "mc": float(mc)}


def b_value_by_year(df, mc):
    """b-value per calendar year: a stability check on the fit."""
    rows = []
    for year, g in df[df["Magnitude"] >= mc].groupby("Year"):
        if len(g) >= 50:
            fit = gutenberg_richter(g["Magnitude"], mc)
            rows.append({"Year": year, "b": fit["b"], "b_se": fit["b_se"], "n": fit["n"]})
    return pd.DataFrame(rows)


def decluster_gardner_knopoff(df):
    """Flag mainshocks using Gardner & Knopoff (1974) space-time windows.

    Earthquakes cluster (aftershock sequences), but rate calculations assume
    independent events. Starting from the largest event, any smaller event inside
    its distance/time window is labelled a dependent event. Returns a boolean
    Series (True = mainshock) aligned with df.index.
    """
    d = df.sort_values("Magnitude", ascending=False)
    mag = d["Magnitude"].to_numpy()
    lat = d["Latitude"].to_numpy()
    lon = d["Longitude"].to_numpy()
    t = ((d["Time"] - d["Time"].min()).dt.total_seconds() / 86400.0).to_numpy()

    dist_win = 10 ** (0.1238 * mag + 0.983)                       # km
    time_win = np.where(mag >= 6.5, 10 ** (0.032 * mag + 2.7389),
                        10 ** (0.5409 * mag - 0.547))             # days

    n = len(d)
    dependent = np.zeros(n, dtype=bool)
    idx = np.arange(n)
    for i in range(n):
        if dependent[i]:
            continue
        cand = np.where((np.abs(t - t[i]) <= time_win[i]) & ~dependent & (idx > i))[0]
        if cand.size:
            close = haversine_km(lat[i], lon[i], lat[cand], lon[cand]) <= dist_win[i]
            dependent[cand[close]] = True

    return pd.Series(~dependent, index=d.index).reindex(df.index)


def flag_mainshocks(df, mc):
    """Boolean Series over all of df: True for declustered mainshocks with M >= mc.

    Events below Mc are never mainshocks here (they are excluded from rate
    calculations), so they come back False.
    """
    above = df[df["Magnitude"] >= mc]
    flags = decluster_gardner_knopoff(above)
    return flags.reindex(df.index).fillna(False).astype(bool)


def annual_summary(df, mc):
    """Yearly event counts (raw and mainshocks-only) above Mc, plus mean magnitude."""
    above = df[df["Magnitude"] >= mc]
    out = above.groupby("Year").agg(events=("Magnitude", "size"),
                                    mean_mag=("Magnitude", "mean"),
                                    max_mag=("Magnitude", "max"))
    if "is_mainshock" in above.columns:
        out["mainshocks"] = above[above["is_mainshock"]].groupby("Year").size()
    return out.fillna(0).reset_index()
