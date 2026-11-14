#!/usr/bin/env python3
"""Ground-truth check: MESH vs measured hail in NWS Local Storm Reports.

Pulls hail LSRs (IEM archive) for a year, samples reports with a stated size,
and records the MESH storm-day value at the report point and the max within
1 mile.  Storm day = 12:00 UTC -> 12:00 UTC (06:00 CST convention).
"""
import argparse, csv, datetime as dt, io, json, random, sys, os
from concurrent.futures import ProcessPoolExecutor, as_completed
import urllib.request

sys.path.insert(0, os.path.dirname(__file__))
from mesh_scan import scan_day  # noqa: E402

IEM = "https://mesonet.agron.iastate.edu/cgi-bin/request/gis/lsr.py"


def fetch(year, start=None):
    sts = start or f"{year}-01-01T00:00Z"
    url = f"{IEM}?sts={sts}&ets={year + 1}-01-01T00:00Z&fmt=csv"
    with urllib.request.urlopen(url, timeout=300) as r:
        txt = r.read().decode("utf-8", "replace")
    return list(csv.DictReader(io.StringIO(txt)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--year", type=int, required=True)
    ap.add_argument("--n", type=int, default=150)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    start = "2020-10-14T12:00Z" if a.year == 2020 else None
    rows = fetch(a.year, start)
    hail = []
    for r in rows:
        if (r.get("TYPECODE") or "").strip() != "H":
            continue
        try:
            mag = float(r["MAG"]); lat = float(r["LAT"]); lon = float(r["LON"])
        except (ValueError, KeyError):
            continue
        if mag <= 0:
            continue
        v = dt.datetime.strptime(r["VALID"][:12], "%Y%m%d%H%M")
        day = (v - dt.timedelta(hours=12)).date()
        hail.append({"day": day.isoformat(), "valid": r["VALID"], "size": mag, "lat": lat,
                     "lon": lon, "state": r.get("STATE"), "source": r.get("SOURCE"),
                     "city": r.get("CITY")})
    print(len(rows), "LSR rows,", len(hail), "hail with size", flush=True)
    random.seed(a.year)
    sample = random.sample(hail, min(a.n, len(hail)))
    by_day = {}
    for k, h in enumerate(sample):
        h["id"] = f"r{k}"
        by_day.setdefault(h["day"], []).append(h)
    out = []
    with ProcessPoolExecutor(4) as ex:
        futs = {ex.submit(scan_day, d, [{"id": h["id"], "lat": h["lat"], "lon": h["lon"]} for h in hs]): d
                for d, hs in by_day.items()}
        for f in as_completed(futs):
            d = futs[f]
            try:
                res = f.result()
            except Exception as e:
                print("err", d, e); continue
            for h in by_day[d]:
                m = (res.get("sites") or {}).get(h["id"])
                if m:
                    out.append({**h, **{"mesh_" + k: v for k, v in m.items()}})
    json.dump(out, open(a.out, "w"))
    print(len(out), "reports checked")


if __name__ == "__main__":
    main()
