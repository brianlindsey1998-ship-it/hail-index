#!/usr/bin/env python3
"""Build one year of the radar hail index from NOAA MRMS MESH_Max_1440min.

Storm day D = 12:00 UTC D -> 12:00 UTC D+1 (06:00 CST convention); its value is
the 1440-minute max file valid at 12:00 UTC on D+1.  Every 0.01-degree cell at
or above THRESH_MM is kept, binned into 0.5-degree tiles, one gzip file per tile
per year:

  header : b"HIX1", uint16 year, uint16 nblocks          (little-endian)
  block  : uint16 day-of-year, uint16 n, n x uint16 cell, n x uint8 value
  cell   : row*50 + col inside the tile; cell centre = tile_south + (row+0.5)*0.01,
           tile_west + (col+0.5)*0.01
  value  : MESH in 0.5 mm steps (value/2 = mm), 255 = 127.5 mm or more
  tile   : file name "{floor(lat*2)}_{floor(lon*2)}.bin.gz"

Source: NOAA MRMS via NOAA Open Data on AWS (s3://noaa-mrms-pds), anonymous.
"""
import argparse, datetime as dt, gzip, json, os, struct, sys, tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from mesh_scan import BUCKET, http, pick_key  # noqa: E402

THRESH_MM = 12.7          # 0.50 in
ARCHIVE_START = dt.date(2020, 10, 14)


def day_cells(day_iso):
    """Return (tile_lat_idx, tile_lon_idx, cell, val) arrays for one storm day."""
    import pygrib
    day = dt.date.fromisoformat(day_iso)
    key = pick_key(day)
    if not key:
        return day_iso, None, None
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
    jj, ii = np.nonzero(vals >= THRESH_MM)
    v = vals[jj, ii]
    lat = (lat0 + jj * dlat) if jpos else (lat0 - jj * dlat)
    lon = lon0 + ii * dlon
    # integer hundredths of a degree of the cell's south-west corner
    la_c = np.floor(lat * 100 + 1e-6).astype(np.int64)
    lo_c = np.floor(lon * 100 + 1e-6).astype(np.int64)
    tla, tlo = np.floor_divide(la_c, 50), np.floor_divide(lo_c, 50)
    cell = ((la_c - tla * 50) * 50 + (lo_c - tlo * 50)).astype(np.uint16)
    q = np.minimum(255, np.round(v * 2)).astype(np.uint8)
    return day_iso, key, (tla.astype(np.int16), tlo.astype(np.int16), cell, q)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--year", type=int, required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()
    y = a.year
    first = max(dt.date(y, 1, 1), ARCHIVE_START)
    last = min(dt.date(y, 12, 31), dt.datetime.utcnow().date() - dt.timedelta(days=1))
    days = [(first + dt.timedelta(n)).isoformat() for n in range((last - first).days + 1)]
    tiles = {}            # (tla, tlo) -> list of (doy, cells, vals)
    done, missing, errors = [], [], []
    ncell = 0
    with ProcessPoolExecutor(a.workers) as ex:
        futs = {ex.submit(day_cells, d): d for d in days}
        for n, f in enumerate(as_completed(futs), 1):
            d = futs[f]
            try:
                _, key, arrs = f.result()
            except Exception as e:  # retried by the next nightly run
                errors.append({"date": d, "error": repr(e)})
                continue
            if key is None:
                missing.append(d)
                continue
            done.append(d)
            tla, tlo, cell, q = arrs
            ncell += len(cell)
            doy = dt.date.fromisoformat(d).timetuple().tm_yday
            if len(cell):
                order = np.lexsort((cell, tlo, tla))
                tla, tlo, cell, q = tla[order], tlo[order], cell[order], q[order]
                kk = tla.astype(np.int64) * 100000 + tlo
                cuts = np.flatnonzero(np.diff(kk)) + 1
                for s, e in zip(np.r_[0, cuts], np.r_[cuts, len(kk)]):
                    tiles.setdefault((int(tla[s]), int(tlo[s])), []).append((doy, cell[s:e], q[s:e]))
            if n % 50 == 0:
                print(f"{y}: {n}/{len(days)} days, {ncell} cells, {len(tiles)} tiles", flush=True)
    os.makedirs(a.out, exist_ok=True)
    total = 0
    for (tla, tlo), blocks in tiles.items():
        blocks.sort(key=lambda b: b[0])
        buf = bytearray(b"HIX1" + struct.pack("<HH", y, len(blocks)))
        for doy, c, q in blocks:
            buf += struct.pack("<HH", doy, len(c)) + c.astype("<u2").tobytes() + q.tobytes()
        data = gzip.compress(bytes(buf), 9, mtime=0)
        total += len(data)
        with open(os.path.join(a.out, f"{tla}_{tlo}.bin.gz"), "wb") as fh:
            fh.write(data)
    meta = {"year": y, "first": first.isoformat(), "last_processed": max(done) if done else None,
            "days": len(done), "missing": sorted(missing), "errors": errors,
            "tiles": len(tiles), "cells": ncell, "bytes": total, "thresh_mm": THRESH_MM}
    json.dump(meta, open(os.path.join(a.out, "meta.json"), "w"), indent=1)
    print(json.dumps({k: v for k, v in meta.items() if k not in ("missing", "errors")}),
          "missing", len(missing), "errors", len(errors))


if __name__ == "__main__":
    main()
