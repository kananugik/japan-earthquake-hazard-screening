"""Gridded seismic hazard screening.

Method (a transparent, simplified take on classical probabilistic hazard):
  1. Keep declustered mainshocks above the completeness magnitude Mc.
  2. Count them per grid cell and smooth with a Gaussian kernel (earthquake
     locations are uncertain and sparse data should inform neighbours).
  3. Convert to an annual rate of M >= Mc, then extrapolate to the target
     magnitude with the Gutenberg-Richter law: rate(M>=Mt) = rate(M>=Mc)*10^(-b(Mt-Mc)).
  4. Assuming a Poisson process, P(at least one M>=Mt in T years) = 1 - exp(-rate*T).

This ranks locations by relative earthquake activity. It is NOT a full PSHA: it has
no ground-motion model, site effects or fault-specific recurrence.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import geopandas as gpd
from scipy.ndimage import gaussian_filter
from shapely.geometry import box

TIERS = ["Low", "Moderate", "High", "Very high"]
SHALLOW_KM = 70


def make_edges(df, cell_deg):
    lat0 = np.floor(df["Latitude"].min() / cell_deg) * cell_deg
    lat1 = np.ceil(df["Latitude"].max() / cell_deg) * cell_deg + cell_deg
    lon0 = np.floor(df["Longitude"].min() / cell_deg) * cell_deg
    lon1 = np.ceil(df["Longitude"].max() / cell_deg) * cell_deg + cell_deg
    return (np.round(np.arange(lat0, lat1 + 1e-9, cell_deg), 6),
            np.round(np.arange(lon0, lon1 + 1e-9, cell_deg), 6))


def rate_surface(events, lat_edges, lon_edges, years, sigma_cells=1.0):
    """Smoothed annual rate (events/yr per cell) of the given events."""
    counts, _, _ = np.histogram2d(events["Latitude"], events["Longitude"],
                                  bins=[lat_edges, lon_edges])
    smooth = gaussian_filter(counts, sigma=sigma_cells, mode="constant") if sigma_cells > 0 else counts
    if smooth.sum() > 0:                      # keep total event count unchanged
        smooth = smooth * counts.sum() / smooth.sum()
    return counts, smooth / years


def _years(df):
    return (df["Time"].max() - df["Time"].min()).total_seconds() / (365.25 * 86400)


def hazard_grid(df, mc, b, cell_deg=1.0, sigma_cells=1.0, target_mag=6.0, horizon_years=10):
    """Return (GeoDataFrame of cells, metadata dict). df needs an 'is_mainshock' column."""
    years = _years(df)
    lat_e, lon_e = make_edges(df, cell_deg)
    above = df[df["Magnitude"] >= mc]
    main = above[above["is_mainshock"]]

    counts_main, rate_mc = rate_surface(main, lat_e, lon_e, years, sigma_cells)
    rate_t = rate_mc * 10 ** (-b * (target_mag - mc))
    p = 1 - np.exp(-rate_t * horizon_years)

    i = np.clip(np.digitize(above["Latitude"], lat_e) - 1, 0, len(lat_e) - 2)
    j = np.clip(np.digitize(above["Longitude"], lon_e) - 1, 0, len(lon_e) - 2)
    stats = (above.assign(i=i, j=j).groupby(["i", "j"])
             .agg(n_events=("Magnitude", "size"), max_mag=("Magnitude", "max"),
                  mean_depth_km=("Depth_km", "mean"),
                  shallow_share=("Depth_km", lambda s: (s < SHALLOW_KM).mean())))

    ii, jj = np.indices(p.shape)
    cells = pd.DataFrame({"i": ii.ravel(), "j": jj.ravel(), "p_target": p.ravel(),
                          "rate_mc_per_yr": rate_mc.ravel(),
                          "n_mainshocks": counts_main.ravel().astype(int)})
    cells = cells.merge(stats.reset_index(), on=["i", "j"], how="left")
    cells["n_events"] = cells["n_events"].fillna(0).astype(int)
    cells = cells[(cells["p_target"] >= 1e-4) | (cells["n_events"] > 0)].copy()

    pct = cells["p_target"].rank(pct=True)
    cells["tier"] = pd.cut(pct, [0, .5, .8, .95, 1.0], labels=TIERS, include_lowest=True)
    cells["geometry"] = [box(lon_e[j_], lat_e[i_], lon_e[j_ + 1], lat_e[i_ + 1])
                         for i_, j_ in zip(cells["i"], cells["j"])]
    cells["lat_c"] = [(lat_e[i_] + lat_e[i_ + 1]) / 2 for i_ in cells["i"]]
    cells["lon_c"] = [(lon_e[j_] + lon_e[j_ + 1]) / 2 for j_ in cells["j"]]
    grid = gpd.GeoDataFrame(cells.drop(columns=["i", "j"]), geometry="geometry", crs=4326)

    expected = float(rate_t.sum() * years)
    observed = int((main["Magnitude"] >= target_mag).sum())
    meta = {"years": years, "mc": mc, "b": b, "cell_deg": cell_deg, "sigma_cells": sigma_cells,
            "target_mag": target_mag, "horizon_years": horizon_years,
            "n_mainshocks_above_mc": int(len(main)),
            "expected_events_target_in_catalog": expected,
            "observed_mainshocks_target_in_catalog": observed}
    return grid, meta


def backtest(df, mc, split_year, cell_deg=1.0, sigma_cells=1.0, test_mag=5.5):
    """Out-of-sample check: fit on years <= split_year, test on later events.

    Cells are ranked by predicted rate; the curve shows what share of later
    M >= test_mag mainshocks fall inside the top x% of cells. The diagonal is what
    random guessing would give.
    """
    lat_e, lon_e = make_edges(df, cell_deg)
    train = df[(df["Year"] <= split_year) & (df["Magnitude"] >= mc) & df["is_mainshock"]]
    test = df[(df["Year"] > split_year) & (df["Magnitude"] >= test_mag) & df["is_mainshock"]]
    if len(test) < 10 and test_mag > mc:                       # too few: relax once
        test_mag = max(mc, test_mag - 0.5)
        test = df[(df["Year"] > split_year) & (df["Magnitude"] >= test_mag) & df["is_mainshock"]]
    train_years = _years(df[df["Year"] <= split_year])
    _, rate = rate_surface(train, lat_e, lon_e, train_years, sigma_cells)
    test_counts, _, _ = np.histogram2d(test["Latitude"], test["Longitude"], bins=[lat_e, lon_e])

    order = np.argsort(-rate.ravel(), kind="stable")
    hits = np.cumsum(test_counts.ravel()[order]) / max(test_counts.sum(), 1)
    area = np.arange(1, len(order) + 1) / len(order)
    k = np.searchsorted(area, 0.2)
    return {"area_fraction": area, "hit_fraction": hits, "n_test": int(len(test)),
            "test_mag": float(test_mag), "split_year": int(split_year),
            "hit_at_20pct_area": float(hits[k]),
            "gain_at_20pct_area": float(hits[k] / 0.2)}


def regional_summary(df, mc, top=15):
    """Per-region (country/area label from USGS 'place') activity table."""
    above = df[df["Magnitude"] >= mc]
    out = (above.groupby("Region")
           .agg(events=("Magnitude", "size"), max_mag=("Magnitude", "max"),
                m5_plus=("Magnitude", lambda s: int((s >= 5).sum())),
                mean_depth_km=("Depth_km", "mean"),
                shallow_share=("Depth_km", lambda s: (s < SHALLOW_KM).mean()))
           .sort_values("events", ascending=False).head(top).round(2))
    return out.reset_index()
