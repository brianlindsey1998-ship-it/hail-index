#!/usr/bin/env python3
"""Write data/manifest.json from each year's meta.json."""
import glob, json, os, sys, datetime as dt
root = sys.argv[1]
years = {}
for m in sorted(glob.glob(os.path.join(root, "y*/meta.json"))):
    j = json.load(open(m))
    years[str(j["year"])] = {k: j[k] for k in ("first", "last_processed", "days", "missing", "tiles", "bytes")}
last = max((v["last_processed"] for v in years.values() if v["last_processed"]), default=None)
json.dump({
    "format": "HIX1", "tile_deg": 0.5, "cell_deg": 0.01, "thresh_in": 0.5, "value_unit": "0.5 mm",
    "storm_day": "12:00 UTC to 12:00 UTC (06:00 CST to 06:00 CST)",
    "source": "NOAA MRMS MESH_Max_1440min_00.50, NOAA Open Data on AWS (s3://noaa-mrms-pds)",
    "archive_start": "2020-10-14", "last_storm_day": last,
    "built_utc": dt.datetime.utcnow().strftime("%Y-%m-%dT%H:%MZ"), "years": years,
}, open(os.path.join(root, "manifest.json"), "w"), indent=1)
print("manifest: last storm day", last)
