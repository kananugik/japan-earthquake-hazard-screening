"""Synthetic USGS-style catalog used ONLY to smoke-test the pipeline code.
Not real data; never use its outputs for analysis."""
import sys
import numpy as np
import pandas as pd

rng = np.random.default_rng(7)
N = 9000
b, mmin = 1.0, 4.3
mag = np.round(mmin - np.log10(1 - rng.random(N)) / b, 1)
mag = np.clip(mag, mmin, 8.3)
# subduction-like arc: trench line, depth increases westward
t = rng.random(N)
lat = 31 + 14 * t + rng.normal(0, 0.6, N)
trench_lon = 141 + 3 * t
dip = np.abs(rng.normal(0, 3.0, N))
lon = trench_lon - dip * 0.8 + rng.normal(0, 0.3, N)
depth = np.clip(10 + dip * 45 + rng.normal(0, 8, N), 1, 600)
times = pd.to_datetime("2015-01-01", utc=True) + pd.to_timedelta(rng.random(N) * 3652, unit="D")
df = pd.DataFrame({"time": times, "latitude": lat, "longitude": lon, "depth": depth, "mag": mag})
# aftershock swarm after one big event
main = df.sort_values("mag").iloc[-1]
k = 600
sw = pd.DataFrame({"time": main.time + pd.to_timedelta(rng.exponential(5, k), unit="D"),
                   "latitude": main.latitude + rng.normal(0, 0.4, k),
                   "longitude": main.longitude + rng.normal(0, 0.4, k),
                   "depth": np.clip(main.depth + rng.normal(0, 5, k), 1, 600),
                   "mag": np.clip(np.round(mmin - np.log10(1 - rng.random(k)) / 1.1, 1), mmin, 6.8)})
df = pd.concat([df, sw], ignore_index=True)
df["magType"] = "mb"
df["place"] = rng.choice(["Izu Islands, Japan region", "86 km WNW of Amagi, Japan",
                          "26 km E of Oma, Japan", "off the east coast of Honshu, Japan"], len(df))
df["type"] = "earthquake"
df["status"] = "reviewed"
df["id"] = [f"us{i:07d}" for i in range(len(df))]
df = df.sort_values("time", ascending=False).reset_index(drop=True)
# messy rows to exercise cleaning
df.loc[5:12, "type"] = "quarry blast"
df.loc[20:24, "mag"] = np.nan
df = pd.concat([df, df.iloc[30:40]], ignore_index=True)       # duplicate ids
df.loc[50, "depth"] = -1.2
df["time"] = df["time"].dt.strftime("%Y-%m-%dT%H:%M:%S.%f").str[:-3] + "Z"
df.to_csv(sys.argv[1], index=False)
print(len(df), "rows ->", sys.argv[1])
