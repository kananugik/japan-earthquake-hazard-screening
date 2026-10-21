"""Static (matplotlib) figures and one lightweight interactive (Folium) map."""
from __future__ import annotations

from pathlib import Path

import folium
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
from folium.plugins import Fullscreen, HeatMap, MiniMap, TimestampedGeoJson

INK, MUTED, ACCENT = "#1f2933", "#6b7785", "#c2410c"
DEPTH_COLORS = {"Shallow (<70 km)": "#fde725", "Intermediate (70-300 km)": "#21918c",
                "Deep (>300 km)": "#440154"}
HAZARD_CMAP = "YlOrRd"


def apply_style():
    mpl.rcParams.update({
        "figure.dpi": 110, "savefig.dpi": 200, "savefig.bbox": "tight",
        "font.family": "DejaVu Sans", "font.size": 10,
        "axes.edgecolor": MUTED, "axes.labelcolor": INK, "axes.titlecolor": INK,
        "axes.titlesize": 12, "axes.titleweight": "bold", "axes.titlelocation": "left",
        "axes.spines.top": False, "axes.spines.right": False,
        "xtick.color": MUTED, "ytick.color": MUTED, "text.color": INK,
        "axes.grid": True, "grid.color": "#e5e7eb", "grid.linewidth": 0.7,
        "axes.axisbelow": True, "legend.frameon": False,
    })


def depth_class(depth):
    return np.where(depth < 70, "Shallow (<70 km)",
                    np.where(depth <= 300, "Intermediate (70-300 km)", "Deep (>300 km)"))


def save(fig, path):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path)
    return path


# ---------------------------------------------------------------- static figures
def fig_gutenberg_richter(df, fit):
    apply_style()
    mags = np.sort(df["Magnitude"].to_numpy())
    grid = np.round(np.arange(mags.min(), mags.max() + 0.05, 0.1), 1)
    cum = np.array([(mags >= m - 1e-9).sum() for m in grid])
    fig, ax = plt.subplots(figsize=(7, 4.6))
    ax.scatter(grid, cum, s=14, color=MUTED, label="Observed, N(M ≥ m)", zorder=3)
    x = np.array([fit["mc"], grid.max()])
    ax.plot(x, 10 ** (fit["a"] - fit["b"] * x), color=ACCENT, lw=2,
            label=f"Fit: b = {fit['b']:.2f} ± {fit['b_se']:.2f}")
    ax.axvline(fit["mc"], color=INK, ls="--", lw=1)
    ax.text(fit["mc"] + 0.04, 0.04, f"Mc = {fit['mc']:.1f}  (catalog complete above this)",
            fontsize=9, va="bottom", transform=ax.get_xaxis_transform())
    ax.set_yscale("log")
    ax.set_xlabel("Magnitude")
    ax.set_ylabel("Cumulative number of events")
    ax.set_title("Frequency\u2013magnitude distribution (Gutenberg-Richter)")
    ax.legend()
    return fig


def fig_annual(annual):
    apply_style()
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.bar(annual["Year"], annual["events"], color="#d1d5db", label="All events")
    if "mainshocks" in annual:
        ax.bar(annual["Year"], annual["mainshocks"], color=ACCENT, label="Mainshocks (declustered)")
    ax.set_xlabel("Year")
    ax.set_ylabel("Events per year (M ≥ Mc)")
    ax.set_title("Events per year: all vs declustered mainshocks")
    ax.set_xticks(annual["Year"])
    ax.tick_params(axis="x", rotation=45)
    ax.legend()
    return fig


def fig_depth(df, mc):
    apply_style()
    d = df[df["Magnitude"] >= mc]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.3), gridspec_kw={"width_ratios": [1, 1.2]})
    ax = axes[0]
    ax.hist(d["Depth_km"], bins=40, color="#9ca3af", edgecolor="white")
    ax.set_yscale("log")
    for v in (70, 300):
        ax.axvline(v, color=INK, ls="--", lw=1)
        ax.text(v + 6, 0.95, f"{v} km", rotation=90, va="top", fontsize=8, color=MUTED,
                transform=ax.get_xaxis_transform())
    ax.set_xlabel("Depth (km)")
    ax.set_ylabel("Events (log scale)")
    ax.set_title("Depth distribution (M \u2265 Mc)")
    ax = axes[1]
    hb = ax.hexbin(d["Depth_km"], d["Magnitude"], gridsize=35, bins="log", cmap="viridis_r",
                   mincnt=1)
    ax.set_xlabel("Depth (km)")
    ax.set_ylabel("Magnitude")
    ax.set_title("Magnitude vs depth")
    ax.grid(False)
    fig.colorbar(hb, ax=ax, label="Events (log)")
    return fig


def fig_b_by_year(bty, overall_b):
    apply_style()
    fig, ax = plt.subplots(figsize=(7, 3.6))
    ax.errorbar(bty["Year"], bty["b"], yerr=bty["b_se"], fmt="o", color=INK, capsize=3)
    ax.axhline(overall_b, color=ACCENT, ls="--", label=f"Overall b = {overall_b:.2f}")
    ax.set_xlabel("Year")
    ax.set_ylabel("b-value")
    ax.set_title("b-value by year (\u00b11 standard error)")
    ax.legend()
    return fig


def fig_plate_distance(dist_km, mags, min_mag=5.0):
    apply_style()
    big = dist_km[mags >= min_mag]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    axes[0].hist(dist_km[dist_km <= 600], bins=40, color="#9ca3af", edgecolor="white")
    axes[0].set_xlabel("Distance to nearest plate boundary (km)")
    axes[0].set_ylabel("Events")
    axes[0].set_title("Distance from events to nearest plate boundary")
    x = np.sort(big)
    y = np.arange(1, len(x) + 1) / len(x)
    axes[1].plot(x, y, color=ACCENT, lw=2)
    within = float((big <= 100).mean()) if len(big) else float("nan")
    axes[1].axvline(100, color=INK, ls="--", lw=1)
    axes[1].set_xlim(0, 600)
    axes[1].set_xlabel("Distance to nearest plate boundary (km)")
    axes[1].set_ylabel(f"Share of M ≥ {min_mag:g} events")
    axes[1].set_title(f"{within:.0%} of M ≥ {min_mag:g} events within 100 km")
    return fig


def fig_backtest(bt):
    apply_style()
    fig, ax = plt.subplots(figsize=(5.6, 5.2))
    ax.plot(bt["area_fraction"], bt["hit_fraction"], color=ACCENT, lw=2.4,
            label="Model trained on early years")
    ax.plot([0, 1], [0, 1], color=MUTED, ls="--", label="Random guess")
    ax.set_xlabel("Share of map area flagged (highest-rate cells first)")
    ax.set_ylabel(f"Share of later M ≥ {bt['test_mag']:g} events captured")
    ax.set_title(f"Back-test: train ≤ {bt['split_year']}, test after (n = {bt['n_test']})")
    ax.text(0.45, 0.12, f"Top 20% of area captures\n{bt['hit_at_20pct_area']:.0%} of events "
            f"({bt['gain_at_20pct_area']:.1f}× random)", fontsize=10)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.02)
    ax.set_aspect("equal")
    ax.legend(loc="lower right")
    return fig


def fig_hazard_map(grid, meta, plates, events, big_mag=6.0):
    apply_style()
    fig, ax = plt.subplots(figsize=(9, 7))
    vmax = float(grid["p_target"].quantile(0.99))
    norm = mcolors.PowerNorm(0.6, vmin=0, vmax=max(vmax, 1e-3))
    shown = grid[grid["p_target"] >= 0.002]
    shown.plot(column="p_target", cmap=HAZARD_CMAP, norm=norm, ax=ax, linewidth=0, alpha=0.95)
    b = grid.total_bounds
    plates.plot(ax=ax, color=INK, linewidth=0.8, alpha=0.7)
    big = events[events["Magnitude"] >= big_mag]
    ax.scatter(big["Longitude"], big["Latitude"], s=(big["Magnitude"] - 4) ** 3 * 4, facecolor="none",
               edgecolor=INK, linewidth=0.9, label=f"Observed M ≥ {big_mag:g}")
    ax.set_xlim(b[0], b[2])
    ax.set_ylim(b[1], b[3])
    ax.set_aspect(1 / np.cos(np.radians((b[1] + b[3]) / 2)))
    ax.grid(False)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_title(f"Probability of an M ≥ {meta['target_mag']:g} earthquake within "
                 f"{meta['horizon_years']} years, per {meta['cell_deg']:g}° cell")
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    cax = make_axes_locatable(ax).append_axes("right", size="3%", pad=0.12)
    sm = mpl.cm.ScalarMappable(norm=norm, cmap=HAZARD_CMAP)
    cb = fig.colorbar(sm, cax=cax)
    cb.set_label("Probability")
    cb.ax.yaxis.set_major_formatter(mpl.ticker.PercentFormatter(1.0, decimals=0))
    ax.legend(loc="lower left")
    return fig


# ------------------------------------------------------------- interactive map
def _cap_by_magnitude(df, start_mag, max_n):
    m = start_mag
    while (df["Magnitude"] >= m).sum() > max_n and m < df["Magnitude"].max():
        m = round(m + 0.1, 1)
    return m


def build_map(df, grid, meta, plates_path, mc, sig_mag=5.5, slider_mag=5.0,
              max_markers=400, max_slider=1500):
    """Interactive map kept deliberately light (a few thousand features, not tens of thousands).

    Layers: hazard grid (choropleth), plate boundaries, mainshock density heatmap,
    significant events coloured by depth, and a time slider for M >= slider_mag.
    """
    lat_c, lon_c = df["Latitude"].median(), df["Longitude"].median()
    m = folium.Map(location=[lat_c, lon_c], zoom_start=5, tiles=None,
                   control_scale=True, prefer_canvas=True)
    folium.TileLayer("https://{s}.basemaps.cartocdn.com/light_all/{z}/{x}/{y}{r}.png", name="Light basemap",
                     attr="&copy; OpenStreetMap contributors &copy; CARTO", subdomains="abcd",
                     max_zoom=19).add_to(m)
    b = grid.total_bounds
    m.fit_bounds([[b[1], b[0]], [b[3], b[2]]])

    # 1) hazard grid
    vmax = float(grid["p_target"].quantile(0.99))
    norm = mcolors.PowerNorm(0.6, vmin=0, vmax=max(vmax, 1e-3))
    cmap = mpl.colormaps[HAZARD_CMAP]
    shown = grid[grid["p_target"] >= 0.002].copy()
    shown["p_pct"] = (shown["p_target"] * 100).round(1)
    shown["mean_depth_km"] = shown["mean_depth_km"].round(0)
    shown["max_mag"] = shown["max_mag"].round(1)
    shown["shallow_pct"] = (shown["shallow_share"] * 100).round(0)
    shown["tier"] = shown["tier"].astype(str)
    cols = ["p_target", "p_pct", "tier", "n_events", "max_mag", "mean_depth_km", "shallow_pct", "geometry"]
    layer = folium.FeatureGroup(name=f"Hazard: P(M≥{meta['target_mag']:g} in "
                                     f"{meta['horizon_years']} yr)", show=True)
    folium.GeoJson(
        shown[cols].to_json(),
        style_function=lambda f: {"fillColor": mcolors.to_hex(cmap(norm(f["properties"]["p_target"]))),
                                  "color": "none", "fillOpacity": 0.72},
        tooltip=folium.GeoJsonTooltip(
            fields=["p_pct", "tier", "n_events", "max_mag", "mean_depth_km", "shallow_pct"],
            aliases=["Probability (%)", "Relative tier", f"Events ≥ {mc:g}", "Largest event (M)",
                     "Mean depth (km)", "Shallow share (%)"], sticky=True),
    ).add_to(layer)
    layer.add_to(m)

    # 2) plate boundaries
    pb = folium.FeatureGroup(name="Plate boundaries", show=True)
    folium.GeoJson(str(plates_path), style_function=lambda f: {"color": "#1f2933", "weight": 1.6,
                                                              "opacity": 0.75}).add_to(pb)
    pb.add_to(m)

    # 3) density heatmap (declustered mainshocks above Mc)
    main = df[(df["Magnitude"] >= mc) & df["is_mainshock"]]
    hm = folium.FeatureGroup(name="Mainshock density heatmap", show=False)
    HeatMap(main[["Latitude", "Longitude"]].values.tolist(), radius=10, blur=14, min_opacity=0.25,
            gradient={0.2: "#fee8c8", 0.5: "#fdbb84", 0.8: "#e34a33", 1.0: "#7f0000"}).add_to(hm)
    hm.add_to(m)

    # 4) significant events, coloured by depth
    sm_mag = _cap_by_magnitude(df, sig_mag, max_markers)
    sig = df[df["Magnitude"] >= sm_mag]
    ev = folium.FeatureGroup(name=f"Events M≥{sm_mag:g} (colour = depth)", show=True)
    for r in sig.itertuples():
        color = DEPTH_COLORS[str(depth_class(np.array([r.Depth_km]))[0])]
        folium.CircleMarker(
            [r.Latitude, r.Longitude], radius=2 + (r.Magnitude - 4) * 2.6, color="#111827", weight=0.8,
            fill=True, fill_color=color, fill_opacity=0.9,
            popup=folium.Popup(f"<b>M {r.Magnitude:.1f}</b><br>Depth: {r.Depth_km:.0f} km<br>"
                               f"{r.Time:%Y-%m-%d %H:%M} UTC<br>{getattr(r, 'Location', '')}", max_width=260),
        ).add_to(ev)
    ev.add_to(m)

    # 5) time slider (capped so the page stays light)
    sl_mag = _cap_by_magnitude(df, slider_mag, max_slider)
    sl = df[df["Magnitude"] >= sl_mag]
    feats = []
    for r in sl.itertuples():
        feats.append({
            "type": "Feature",
            "geometry": {"type": "Point", "coordinates": [r.Longitude, r.Latitude]},
            "properties": {
                "time": f"{r.Time:%Y-%m-%dT%H:%M:%S}",
                "icon": "circle",
                "iconstyle": {"fillColor": DEPTH_COLORS[str(depth_class(np.array([r.Depth_km]))[0])],
                              "fillOpacity": 0.85, "stroke": "true", "color": "#111827", "weight": 1,
                              "radius": float(3 + (r.Magnitude - 4) * 2.4)},
                "popup": f"M {r.Magnitude:.1f} | {r.Depth_km:.0f} km | {r.Time:%Y-%m-%d}",
            }})
    TimestampedGeoJson({"type": "FeatureCollection", "features": feats}, period="P1M", duration="P1M",
                       transition_time=150, auto_play=False, loop=False, add_last_point=False,
                       date_options="YYYY-MM", time_slider_drag_update=True).add_to(m)

    # legend
    grad = ", ".join(mcolors.to_hex(cmap(norm(v))) for v in np.linspace(0, vmax, 6))
    depth_rows = "".join(f'<div><span style="display:inline-block;width:11px;height:11px;border-radius:50%;'
                         f'background:{c};border:1px solid #111827;margin-right:6px"></span>{k}</div>'
                         for k, c in DEPTH_COLORS.items())
    legend = f"""
    <div style="position:fixed;top:95px;left:12px;z-index:9999;background:#fff;padding:10px 12px;
                border:1px solid #cbd5e1;border-radius:8px;font:12px/1.5 Arial,sans-serif;width:215px;
                box-shadow:0 2px 8px rgba(0,0,0,.15)">
      <b>P(M≥{meta['target_mag']:g} within {meta['horizon_years']} yr)</b>
      <div style="height:10px;margin:5px 0 2px;background:linear-gradient(to right,{grad});border-radius:2px"></div>
      <div style="display:flex;justify-content:space-between"><span>0%</span><span>{vmax:.0%}+</span></div>
      <hr style="margin:8px 0"><b>Event depth</b>{depth_rows}
      <div style="margin-top:4px;color:#6b7785">Circle size = magnitude</div>
    </div>"""
    m.get_root().html.add_child(folium.Element(legend))
    MiniMap(toggle_display=True, position="bottomright").add_to(m)
    Fullscreen(position="topright").add_to(m)
    folium.LayerControl(collapsed=False).add_to(m)
    return m, {"marker_min_mag": sm_mag, "n_markers": int(len(sig)),
               "slider_min_mag": sl_mag, "n_slider": int(len(sl))}
