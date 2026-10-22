#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""地震活动性复核：b 值、a 值与目标震级概率。

用法：
    python3 src/seismic.py --catalog data/catalog.json --out out/report.json
"""
import argparse
import datetime as dt
import json
import math
import os


def load(path):
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def analyse(data):
    cat = data["catalog"]
    events = cat["events"]
    mc = float(cat["mc"])
    target = float(data["forecast"]["target_magnitude"])
    horizon = float(data["forecast"]["horizon_years"])

    mags = [float(e["mag"]) for e in events]
    times = [dt.datetime.fromisoformat(e["time_utc"]) for e in events]

    span_days = (max(times) - dt.datetime(2011, 1, 1)).total_seconds() / 86400.0
    years = span_days / 365.25

    n = len(mags)
    mean_m = sum(mags) / n
    b = math.log10(math.e) / (mean_m - mc)
    sigma = b / math.sqrt(n)
    a = math.log10(n) + b * mc
    rate = 10.0 ** (a - b * target)
    prob = 1.0 - math.exp(-rate * horizon)
    return {"n": n, "years": years, "b": b, "sigma": sigma, "a": a, "rate": rate, "prob": prob}


def main():
    ap = argparse.ArgumentParser(description="地震活动性 b 值复核")
    ap.add_argument("--catalog", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    data = load(args.catalog)
    res = analyse(data)
    cat = data["catalog"]

    report = {
        "n_total": len(cat["events"]),
        "n_selected": res["n"],
        "mc": float(cat["mc"]),
        "magnitude_bin": float(cat["magnitude_bin"]),
        "observation_years": round(res["years"], 2),
        "b_value": round(res["b"], 3),
        "b_sigma": round(res["sigma"], 3),
        "a_value": round(res["a"], 3),
        "annual_rate_target": round(res["rate"], 4),
        "probability_horizon": round(res["prob"], 4),
        "meta": {"method": "aki-utsu",
                 "target_magnitude": float(data["forecast"]["target_magnitude"]),
                 "horizon_years": float(data["forecast"]["horizon_years"])},
    }

    out_dir = os.path.dirname(os.path.abspath(args.out))
    if out_dir and not os.path.isdir(out_dir):
        os.makedirs(out_dir, exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(report, f, ensure_ascii=False, indent=2)
        f.write("\n")


if __name__ == "__main__":
    main()
