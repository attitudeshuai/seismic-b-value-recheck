#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""地震活动性复核：Gutenberg-Richter 的 b 值、a 值与目标震级概率。

命令行：
    python3 src/seismic.py --catalog data/catalog.json --out out/report.json

口径：
  · 只统计 M ≥ mc 的事件（目录里 mc 以下的小震要剔掉，mc 由目录元数据给出）；
  · b 值用 Aki-Utsu 极大似然：b = log10(e) / (M̄ - M0)，其中 M0 = mc - ΔM/2；
    震级是分档记录的（ΔM 是档宽），M0 取档中心而不是下边界，这一步漏了会系统性偏大；
  · b 值标准差按 b/√N；
  · a 值按年化：log10 λ(M) = a - b·M，用 N/T（T 是观测年数）与 mc 定出
    a = log10(N/T) + b·mc；
  · 目标震级年发生率 λ = 10^(a - b·M_target)，T 年内至少一次的概率按泊松过程 1 - exp(-λT)；
  · 分段：在 mc 以上的事件序列里找分段时刻，使两段各自的 Aki-Utsu 对数似然之和最大
    （单段对数似然 lnL = N·(ln β̂ - 1)，β̂ = 1/(M̄ - M0)），两段各至少 60 条；
  · 输出先写临时文件再原子替换；非法输入以退出码 2 结束且不落文件。
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import math
import os
import sys
import tempfile

EXIT_INVALID = 2


def _fail(msg):
    sys.stderr.write("非法输入：%s\n" % msg)
    raise SystemExit(EXIT_INVALID)


def _parse_time(text, where):
    if not isinstance(text, str) or not text.strip():
        _fail("%s 时间缺失" % where)
    try:
        return dt.datetime.fromisoformat(text.replace("Z", "+00:00"))
    except ValueError:
        _fail("%s 时间格式无法解析：%r" % (where, text))


def load_catalog(path):
    try:
        with open(path, encoding="utf-8") as f:
            raw = json.load(f)
    except FileNotFoundError:
        _fail("读不到目录 %s" % path)
    except ValueError as e:
        _fail("目录不是合法 JSON：%s" % e)
    if not isinstance(raw, dict):
        _fail("目录顶层要是对象")
    cat = raw.get("catalog")
    if not isinstance(cat, dict):
        _fail("缺 catalog")
    mc = cat.get("mc")
    width = cat.get("magnitude_bin")
    if not isinstance(mc, (int, float)) or isinstance(mc, bool) or float(mc) <= 0:
        _fail("catalog.mc 要是正数")
    if not isinstance(width, (int, float)) or isinstance(width, bool) or not (0 < float(width) <= 1.0):
        _fail("catalog.magnitude_bin 要在 0 与 1 之间")
    events = cat.get("events")
    if not isinstance(events, list) or len(events) < 100:
        _fail("events 至少要 100 条")

    forecast = raw.get("forecast") or {}
    target = forecast.get("target_magnitude")
    horizon = forecast.get("horizon_years")
    if not isinstance(target, (int, float)) or isinstance(target, bool):
        _fail("forecast.target_magnitude 要是数值")
    if float(target) < float(mc):
        _fail("forecast.target_magnitude 不能小于 mc")
    if not isinstance(horizon, (int, float)) or isinstance(horizon, bool) or float(horizon) <= 0:
        _fail("forecast.horizon_years 要是正数")

    events_sorted, mags, ids = [], [], set()
    for i, ev in enumerate(events):
        if not isinstance(ev, dict):
            _fail("events[%d] 不是对象" % i)
        rid = ev.get("id")
        if not isinstance(rid, str) or not rid.strip():
            _fail("events[%d].id 不合法" % i)
        if rid in ids:
            _fail("events[%d].id 重复：%s" % (i, rid))
        ids.add(rid)
        mag = ev.get("mag")
        if not isinstance(mag, (int, float)) or isinstance(mag, bool) or not math.isfinite(float(mag)):
            _fail("events[%d].mag 要是有限数值" % i)
        t = _parse_time(ev.get("time_utc"), "events[%d]" % i)
        events_sorted.append(t)
        mags.append(float(mag))
    return {"times": events_sorted, "mags": mags, "mc": float(mc), "width": float(width),
            "target": float(target), "horizon": float(horizon)}


def analyse(scan):
    mc, width = scan["mc"], scan["width"]
    sel_m = [m for m in scan["mags"] if m >= mc]
    if len(sel_m) < 50:
        _fail("mc 以上的事件只有 %d 条，样本不够估 b 值" % len(sel_m))

    sel_t = [t for t, m in zip(scan["times"], scan["mags"]) if m >= mc]
    span_days = (max(sel_t) - min(sel_t)).total_seconds() / 86400.0
    years = span_days / 365.25
    if years <= 0:
        _fail("观测时长算出来不是正数")

    n = len(sel_m)
    mean_m = sum(sel_m) / n
    m0 = mc - width / 2.0
    if mean_m - m0 <= 0:
        _fail("平均震级不大于 mc 的档中心，估不出 b 值")
    b = math.log10(math.e) / (mean_m - m0)
    sigma = b / math.sqrt(n)
    a = math.log10(n / years) + b * mc
    rate = 10.0 ** (a - b * scan["target"])
    prob = 1.0 - math.exp(-rate * scan["horizon"])
    return {"n_selected": n, "years": years, "b": b, "sigma": sigma, "a": a,
            "rate": rate, "prob": prob, "mean_m": mean_m}


MIN_SEGMENT_EVENTS = 60


def _segment_lnl(mags, m0):
    n = len(mags)
    mean_m = sum(mags) / n
    beta = 1.0 / (mean_m - m0)
    return n * (math.log(beta) - 1.0)


def best_split(times, mags, m0):
    """返回 (最优分段下标, 两段对数似然之和, 单段对数似然)；样本不够时返回 None。"""
    n = len(mags)
    if n < 2 * MIN_SEGMENT_EVENTS:
        return None
    base = _segment_lnl(mags, m0)
    best = None
    for i in range(MIN_SEGMENT_EVENTS, n - MIN_SEGMENT_EVENTS + 1):
        val = _segment_lnl(mags[:i], m0) + _segment_lnl(mags[i:], m0)
        if best is None or val > best[1] + 1e-12:
            best = (i, val, base)
    return best


def build_report(scan, res):
    return {
        "n_total": len(scan["mags"]),
        "n_selected": res["n_selected"],
        "mc": scan["mc"],
        "magnitude_bin": scan["width"],
        "observation_years": round(res["years"], 2),
        "b_value": round(res["b"], 3),
        "b_sigma": round(res["sigma"], 3),
        "a_value": round(res["a"], 3),
        "annual_rate_target": round(res["rate"], 4),
        "probability_horizon": round(res["prob"], 4),
        "meta": {"method": "aki-utsu", "target_magnitude": scan["target"],
                 "horizon_years": scan["horizon"]},
    }


def atomic_write(path, obj):
    out_dir = os.path.dirname(os.path.abspath(path))
    if out_dir and not os.path.isdir(out_dir):
        os.makedirs(out_dir, exist_ok=True)
    body = json.dumps(obj, ensure_ascii=False, indent=2) + "\n"
    fd, tmp = tempfile.mkstemp(prefix=".seis-", suffix=".tmp", dir=out_dir or ".")
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="\n") as f:
            f.write(body)
        os.replace(tmp, path)
    except BaseException:
        try:
            os.unlink(tmp)
        except OSError:
            pass
        raise


def main(argv=None):
    ap = argparse.ArgumentParser(description="地震活动性 b 值复核")
    ap.add_argument("--catalog", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)

    scan = load_catalog(args.catalog)
    res = analyse(scan)
    report = build_report(scan, res)

    idx = [i for i, m in enumerate(scan["mags"]) if m >= scan["mc"]]
    sel_times = [scan["times"][i] for i in idx]
    sel_mags = [scan["mags"][i] for i in idx]
    order = sorted(range(len(sel_times)), key=lambda k: sel_times[k])
    sel_times = [sel_times[k] for k in order]
    sel_mags = [sel_mags[k] for k in order]
    m0 = scan["mc"] - scan["width"] / 2.0
    split = best_split(sel_times, sel_mags, m0)
    if split is None:
        report["segmentation"] = None
    else:
        i, total, base = split
        left, right = sel_mags[:i], sel_mags[i:]
        span = (sel_times[i] - sel_times[0]).total_seconds() / 86400.0 / 365.25
        report["segmentation"] = {
            "split_years": round(span, 2),
            "b_before": round(math.log10(math.e) / (sum(left) / len(left) - m0), 3),
            "b_after": round(math.log10(math.e) / (sum(right) / len(right) - m0), 3),
            "delta_lnl": round(total - base, 3),
        }
    atomic_write(args.out, report)
    return 0


if __name__ == "__main__":
    sys.exit(main())
