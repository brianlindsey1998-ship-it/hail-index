# hail-index

Radar hail index for the Forensic Toolbox **Radar Hail Estimate** (Weather tab → Estimated).

- **Source:** NOAA Multi-Radar/Multi-Sensor System (MRMS), MESH_Max_1440min (daily maximum estimated hail size), NOAA Open Data on AWS — https://registry.opendata.aws/noaa-mrms-pds. Archive begins 2020-10-14.
- **Storm day:** 12:00 UTC to 12:00 UTC (6:00 a.m. to 6:00 a.m. CST).
- **Nightly job** (`.github/workflows/index.yml`, 14:23 UTC): rebuilds the current year and publishes it on the `data` branch (a single commit replaced each night). Manual run with `years=all` rebuilds everything.
- **Format:** `data/y{YEAR}/{floor(lat*2)}_{floor(lon*2)}.bin.gz` — 0.5° tiles, every 0.01° cell with MESH ≥ 0.50 in. See `scripts/build_index.py` for the byte layout. `data/manifest.json` gives the last storm day processed.
- **Validation:** `validation/` — point scans for the Cotality (Denver NC) and GAF/PSAI (Land O' Lakes FL) comparison addresses, and a check against measured hail in NWS Local Storm Reports.
