"""Unit tests: the numerical pieces are checked against known answers."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from quake_risk import data, hazard, stats  # noqa: E402


def test_haversine_one_degree_latitude():
    assert abs(stats.haversine_km(0, 0, 1, 0) - 111.19) < 0.1


def test_b_value_recovers_known_b():
    rng = np.random.default_rng(1)
    b_true, mc = 1.0, 4.5
    raw = (mc - 0.05) + rng.exponential(1 / (b_true * np.log(10)), 40000)
    mags = np.round(raw, 1)
    fit = stats.gutenberg_richter(mags, mc)
    assert abs(fit["b"] - b_true) < 0.05


def test_completeness_finds_roll_off():
    rng = np.random.default_rng(2)
    mags = np.round(4.0 + rng.exponential(0.43, 30000), 1)
    # simulate under-reporting below 4.5
    keep = (mags >= 4.5) | (rng.random(len(mags)) < 0.3)
    assert 4.4 <= stats.magnitude_of_completeness(mags[keep]) <= 4.9


def test_declustering_flags_aftershocks_only():
    t0 = pd.Timestamp("2020-01-01", tz="UTC")
    df = pd.DataFrame({
        "Time": [t0, t0 + pd.Timedelta(hours=3), t0 + pd.Timedelta(days=1), t0 + pd.Timedelta(days=400)],
        "Latitude": [10.0, 10.05, 10.1, 30.0],
        "Longitude": [100.0, 100.0, 100.1, 120.0],
        "Magnitude": [6.5, 4.6, 4.8, 4.9],
    })
    flags = stats.decluster_gardner_knopoff(df)
    assert flags.tolist() == [True, False, False, True]


def test_load_catalog_logs_every_cleaning_step(tmp_path):
    csv = tmp_path / "q.csv"
    pd.DataFrame({
        "time": ["2020-01-01T00:00:00Z"] * 4 + ["2020-01-02T00:00:00Z"],
        "latitude": [1, 2, 3, 3, 4], "longitude": [1, 2, 3, 3, 4], "depth": [5, 5, -2, -2, 5],
        "mag": [4.5, 4.6, 4.7, 4.7, np.nan], "type": ["earthquake", "quarry blast", "earthquake",
                                                       "earthquake", "earthquake"],
        "id": ["a", "b", "c", "c", "d"], "place": ["x, Japan"] * 5,
    }).to_csv(csv, index=False)
    df, report = data.load_catalog(csv)
    assert len(df) == 2 and df["Depth_km"].min() >= 0 and set(df["Region"]) == {"Japan"}
    assert report["rows_remaining"].iloc[0] == 5 and report["rows_remaining"].iloc[-1] == 2


def test_hazard_probabilities_are_valid_and_back_test_beats_random():
    rng = np.random.default_rng(3)
    n = 4000
    t = pd.Timestamp("2010-01-01", tz="UTC") + pd.to_timedelta(rng.random(n) * 3650, unit="D")
    # two clusters of activity, so a model trained early should beat a random guess later
    centers = np.where(rng.random(n) < 0.5, 0, 1)
    df = pd.DataFrame({"Time": t, "Year": t.year,
                       "Latitude": 35 + centers * 5 + rng.normal(0, 0.5, n),
                       "Longitude": 140 + centers * 5 + rng.normal(0, 0.5, n),
                       "Depth_km": rng.uniform(5, 50, n),
                       "Magnitude": np.round(4.5 - np.log10(1 - rng.random(n)), 1)})
    df["is_mainshock"] = True
    grid, meta = hazard.hazard_grid(df, 4.5, 1.0)
    assert grid["p_target"].between(0, 1).all()
    bt = hazard.backtest(df, 4.5, 2015, test_mag=5.0)
    assert bt["gain_at_20pct_area"] > 2
