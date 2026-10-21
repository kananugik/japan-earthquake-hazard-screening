"""One-call pipeline: catalog CSV in, tables + figures + interactive map out."""
from __future__ import annotations

import json
from pathlib import Path

from . import data, hazard, stats, tectonics, viz


DEFAULT_PLATES = Path(__file__).resolve().parents[2] / "data" / "plate_boundaries.json"


def run(csv_path, out_dir="outputs", plates_path=DEFAULT_PLATES,
        year_min=None, year_max=None, cell_deg=1.0, sigma_cells=1.0,
        target_mag=6.0, horizon_years=10, split_year=None):
    out = Path(out_dir)
    fig_dir = out / "figures"
    out.mkdir(parents=True, exist_ok=True)

    df, report = data.load_catalog(csv_path, year_min, year_max)
    report.to_csv(out / "cleaning_report.csv", index=False)

    mc = stats.magnitude_of_completeness(df["Magnitude"])
    fit = stats.gutenberg_richter(df["Magnitude"], mc)                       # full catalog (for the G-R figure)
    df["is_mainshock"] = stats.flag_mainshocks(df, mc)
    main_mags = df.loc[(df["Magnitude"] >= mc) & df["is_mainshock"], "Magnitude"]
    fit_main = stats.gutenberg_richter(main_mags, mc)                        # declustered: used for hazard rates

    annual = stats.annual_summary(df, mc)
    bty = stats.b_value_by_year(df, mc)
    grid, meta = hazard.hazard_grid(df, mc, fit_main["b"], cell_deg, sigma_cells, target_mag, horizon_years)

    plates = tectonics.load_plates(plates_path)
    above = df[df["Magnitude"] >= mc].copy()
    above["dist_to_boundary_km"] = tectonics.distance_to_boundary_km(above["Latitude"], above["Longitude"], plates)

    split = split_year or int(df["Year"].min() + (df["Year"].max() - df["Year"].min()) * 0.6)
    bt = hazard.backtest(df, mc, split, cell_deg, sigma_cells)
    regions = hazard.regional_summary(df, mc)

    grid.to_file(out / "hazard_grid.geojson", driver="GeoJSON")
    grid.drop(columns="geometry").to_csv(out / "hazard_grid.csv", index=False)
    regions.to_csv(out / "regional_summary.csv", index=False)
    annual.to_csv(out / "annual_summary.csv", index=False)

    viz.save(viz.fig_gutenberg_richter(df, fit), fig_dir / "01_gutenberg_richter.png")
    viz.save(viz.fig_annual(annual), fig_dir / "02_annual_activity.png")
    viz.save(viz.fig_depth(df, mc), fig_dir / "03_depth_magnitude.png")
    viz.save(viz.fig_plate_distance(above["dist_to_boundary_km"].to_numpy(), above["Magnitude"].to_numpy()),
             fig_dir / "04_plate_boundary_distance.png")
    viz.save(viz.fig_hazard_map(grid, meta, plates, df), fig_dir / "05_hazard_map.png")
    viz.save(viz.fig_backtest(bt), fig_dir / "06_backtest.png")
    if len(bty):
        viz.save(viz.fig_b_by_year(bty, fit["b"]), fig_dir / "07_b_value_by_year.png")

    m, map_info = viz.build_map(df, grid, meta, plates_path, mc)
    m.save(out / "earthquake_risk_map.html")

    big = above[above["Magnitude"] >= 5]
    summary = {
        "events_clean": int(len(df)), "period": [str(df["Time"].min().date()), str(df["Time"].max().date())],
        "mc": mc, "b_full_catalog": round(fit["b"], 3), "b_se_full": round(fit["b_se"], 3), "n_above_mc": fit["n"],
        "b_mainshocks": round(fit_main["b"], 3), "b_se_mainshocks": round(fit_main["b_se"], 3),
        "mainshocks_above_mc": int(above["is_mainshock"].sum()),
        "pct_aftershocks_removed": round(100 * (1 - above["is_mainshock"].mean()), 1),
        "share_M5plus_within_100km_of_boundary": round(float((big["dist_to_boundary_km"] <= 100).mean()), 3) if len(big) else None,
        "hazard_meta": {k: (round(v, 3) if isinstance(v, float) else v) for k, v in meta.items()},
        "backtest": {k: (round(v, 3) if isinstance(v, float) else v) for k, v in bt.items()
                     if not hasattr(v, "__len__")},
        "map": map_info,
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    return {"df": df, "report": report, "fit": fit, "fit_main": fit_main, "grid": grid, "meta": meta, "bt": bt,
            "annual": annual, "regions": regions, "summary": summary}


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("csv", nargs="?", default="data/earthquakes.csv")
    ap.add_argument("--out", default="outputs")
    ap.add_argument("--cell-deg", type=float, default=1.0)
    ap.add_argument("--target-mag", type=float, default=6.0)
    ap.add_argument("--horizon", type=int, default=10)
    a = ap.parse_args()
    res = run(a.csv, a.out, cell_deg=a.cell_deg, target_mag=a.target_mag, horizon_years=a.horizon)
    print(json.dumps(res["summary"], indent=2))
