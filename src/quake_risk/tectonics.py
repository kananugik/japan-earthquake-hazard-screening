"""Distance from events to tectonic plate boundaries (PB2002, Bird 2003)."""
from __future__ import annotations

import geopandas as gpd
import numpy as np
import shapely
from scipy.spatial import cKDTree

from .stats import EARTH_R_KM


def load_plates(path):
    return gpd.read_file(path).to_crs(4326)


def _to_xyz(lat, lon):
    la, lo = np.radians(lat), np.radians(lon)
    return np.column_stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)])


def distance_to_boundary_km(lats, lons, plates, step_deg=0.05):
    """Great-circle distance (km) from each point to the nearest plate boundary.

    Boundaries are densified to ~5 km vertex spacing and searched with a KD-tree
    on the unit sphere, so the result is geodesic (not distorted by map projection).
    """
    dense = shapely.segmentize(plates.geometry.values, step_deg)
    coords = shapely.get_coordinates(dense)                      # (lon, lat)
    tree = cKDTree(_to_xyz(coords[:, 1], coords[:, 0]))
    chord, _ = tree.query(_to_xyz(np.asarray(lats), np.asarray(lons)))
    return 2 * EARTH_R_KM * np.arcsin(np.clip(chord / 2, 0, 1))
