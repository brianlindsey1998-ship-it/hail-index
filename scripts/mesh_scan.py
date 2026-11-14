#!/usr/bin/env python3
"""Point scan of NOAA MRMS MESH_Max_1440min (daily max estimated hail size).

For each storm day D (06:00 CST D -> 06:00 CST D+1 = 12:00 UTC D -> 12:00 UTC D+1)
the 1440-minute rolling-max file valid at 12:00 UTC on D+1 is read, and for each
site the value at the site (pixel containing it) and the max within 1, 3 and 5
miles are recorded.  Sizes are reported in inches (MESH is stored in mm).

Source: NOAA MRMS via NOAA Open Data on AWS, s3://noaa-mrms-pds (anonymous).
"""
import argparse, datetime as dt, gzip, json, math, os, sys, tempfile, time
import xml.etree.ElementTree as ET
from concurrent.futures import ProcessPoolExecutor, as_completed
import urllib.request

import numpy as np

BUCKET = "https://noaa-mrms-pds.s3.amazonaws.com"
PRODUCT = "MESH_Max_1440min_00.50"
MI_KM = 1.609344
R_KM = 6371.0088


def http(url, tries=4):
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return r.read()
        except Exception as e:
            if i == tries - 1:
                raise
            time.sleep(2 * (i + 1))


def list_keys(prefix):
    keys, token = [], None
    while True:
        url = f"{BUCKET}/?list-type=2&prefix={prefix}"
        if token:
            url += "&continuation-token=" + urllib.request.quote(token)
        root = ET.fromstring(http(url))
        ns = {"s": root.tag.split("}")[0].strip("{")}
        keys += [c.find("s:Key", ns).text for c in root.findall("s:Contents", ns)]
        if root.findtext("s:IsTruncated", namespaces=ns) == "true":
            token = root.findtext("s:NextContinuationToken", namespaces=ns)
        else:
            return keys


def pick_key(storm_day):
    """File valid at (or first after) 12:00 UTC the next calendar day."""
    nd = storm_day + dt.timedelta(days=1)
    keys = list_keys(f"CONUS/{PRODUCT}/{nd:%Y%m%d}/")
    best = None
    for k in keys:
        t = k.rsplit("_", 1)[-1].split(".")[0]  # YYYYMMDD-HHMMSS
        try:
            ts = dt.datetime.strptime(t, "%Y%m%d-%H%M%S")
        except ValueError:
            continue
        target = dt.datetime(nd.year, nd.month, nd.day, 12, 0, 0)
        d = (ts - target).total_seconds()
        if -120 <= d <= 3600 and (best is None or abs(d) < abs(best[0])):
            best = (d, k)
    return best[1] if best else None


def hav(lat1, lon1, lat2, lon2):
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dphi = p2 - p1
    dl = np.radians(lon2 - lon1)
    a = np.sin(dphi / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2
    return 2 * R_KM * np.arcsin(np.sqrt(a)) / MI_KM


def scan_day(day_iso, sites):
    import pygrib
    day = dt.date.fromisoformat(day_iso)
    key = pick_key(day)
    if not key:
        return {"date": day_iso, "missing": True}
    raw = gzip.decompress(http(f"{BUCKET}/{key}"))
    with tempfile.NamedTemporaryFile(suffix=".grib2", delete=False) as f:
        f.write(raw)
        path = f.name
    try:
        grbs = pygrib.open(path)
        g = grbs.message(1)
        vals = g.values
        if np.ma.isMaskedArray(vals):
            vals = vals.filled(0)
        lat0 = g["latitudeOfFirstGridPointInDegrees"]
        lon0 = g["longitudeOfFirstGridPointInDegrees"]
        dlat = g["jDirectionIncrementInDegrees"]
        dlon = g["iDirectionIncrementInDegrees"]
        jpos = g["jScansPositively"]
        grbs.close()
    finally:
        os.unlink(path)
    if lon0 > 180:
        lon0 -= 360
    ny, nx = vals.shape
    out = {"date": day_iso, "key": key, "sites": {}}
    for s in sites:
        lat, lon = s["lat"], s["lon"]
        fj = (lat - lat0) / dlat if jpos else (lat0 - lat) / dlat
        fi = (lon - lon0) / dlon
        j, i = int(round(fj)), int(round(fi))
        rj = int(math.ceil(5.5 * MI_KM / 111.0 / dlat)) + 1
        ri = int(math.ceil(5.5 * MI_KM / (111.0 * math.cos(math.radians(lat))) / dlon)) + 1
        j0, j1 = max(0, j - rj), min(ny, j + rj + 1)
        i0, i1 = max(0, i - ri), min(nx, i + ri + 1)
        win = np.asarray(vals[j0:j1, i0:i1], dtype=float)
        win[win < 0] = 0  # -1 / -3 = no coverage / missing
        jj, ii = np.mgrid[j0:j1, i0:i1]
        plat = (lat0 + jj * dlat) if jpos else (lat0 - jj * dlat)
        plon = lon0 + ii * dlon
        dist = hav(lat, lon, plat, plon)
        at = float(np.asarray(vals[j, i], dtype=float)) if 0 <= j < ny and 0 <= i < nx else 0.0
        at = max(at, 0.0)
        r = {"at": round(at / 25.4, 3)}
        for m in (1, 3, 5):
            sel = win[dist <= m]
            r[f"w{m}"] = round(float(sel.max()) / 25.4, 3) if sel.size else 0.0
        out["sites"][s["id"]] = r
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--start", required=True)
    ap.add_argument("--end", required=True)
    ap.add_argument("--sites", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()
    sites = json.load(open(a.sites))
    d0, d1 = dt.date.fromisoformat(a.start), dt.date.fromisoformat(a.end)
    days = [(d0 + dt.timedelta(n)).isoformat() for n in range((d1 - d0).days + 1)]
    results, errors = [], []
    with ProcessPoolExecutor(a.workers) as ex:
        futs = {ex.submit(scan_day, d, sites): d for d in days}
        for n, f in enumerate(as_completed(futs), 1):
            d = futs[f]
            try:
                results.append(f.result())
            except Exception as e:
                errors.append({"date": d, "error": repr(e)})
            if n % 25 == 0:
                print(f"{n}/{len(days)}", flush=True)
    results.sort(key=lambda r: r["date"])
    json.dump({"product": PRODUCT, "results": results, "errors": errors}, open(a.out, "w"))
    print(f"done: {len(results)} days, {len(errors)} errors, "
          f"{sum(1 for r in results if r.get('missing'))} missing")


if __name__ == "__main__":
    main()
