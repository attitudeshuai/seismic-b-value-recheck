#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""判分输出：10 条行为测试，按重要程度加权，权重和 = 1。

每条测试自己造一份地震目录，并按 Gutenberg-Richter 与 Aki-Utsu 的定义式独立复算
b 值、a 值、年发生率与泊松概率，再驱动被测程序的命令行入口比对报告。
期望值不读被测代码，也不读参考解的输出。
"""
import datetime as dt
import json
import math
import os
import subprocess
import sys
import tempfile

APP = os.environ.get("CT_APP_DIR", "/app")
ENTRY = os.path.join(APP, "src", "seismic.py")

WEIGHTS = {
    "test_b_value_matches_aki_utsu_on_base_catalog": 0.12,
    "test_b_value_second_catalog_with_other_parameters": 0.10,
    "test_magnitude_bin_correction_matters": 0.09,
    "test_a_value_and_annual_rate_match": 0.07,
    "test_target_probability_matches_poisson": 0.06,
    "test_catalog_truncated_at_completeness_magnitude": 0.07,
    "test_segment_split_matches_mle_search": 0.09,
    "test_segmentation_null_when_sample_small": 0.08,
    "test_zero_time_span_rejected": 0.07,
    "test_invalid_catalog_rejected_without_output": 0.09,
    "test_report_schema_and_rounding": 0.06,
    "test_rerun_produces_identical_bytes": 0.05,
    "test_failed_run_keeps_previous_report": 0.05,
}


# ---------------------------------------------------------------- 造目录（独立实现）
class _Lcg(object):
    def __init__(self, seed):
        self.state = seed & 0xFFFFFFFF

    def next_u(self):
        self.state = (1103515245 * self.state + 12345) & 0x7FFFFFFF
        return (self.state % 1000003) / 1000003.0


def make_catalog(b_value=1.05, mc=2.5, width=0.1, n=800, n_below=0, years=10.0,
                 seed=101, target=6.5, horizon=30.0, start_year=2012,
                 b_after=None, split_frac=0.55):
    rng = _Lcg(seed)
    raw = []
    for _ in range(n):
        raw.append((rng.next_u() * years, rng.next_u(), rng.next_u(), rng.next_u(), rng.next_u()))
    raw.sort(key=lambda r: r[0])
    events = []
    for pos, (off, u_mag, u_lon, u_lat, u_depth) in enumerate(raw):
        frac = pos / float(n)
        b_here = b_value if (b_after is None or frac < split_frac) else b_after
        mag = mc + (-math.log(1.0 - u_mag * 0.999999)) / (b_here * math.log(10.0))
        mag = math.floor(mag / width) * width
        events.append({
            "id": "a%05d" % (pos + 1),
            "time_utc": "%04d-%02d-%02dT%02d:00:00Z" % (start_year + int(off), 1 + int(u_lon * 12) % 12,
                                                       1 + int(u_lat * 28) % 28, int(u_depth * 24) % 24),
            "lon": 100.0, "lat": 25.0, "depth_km": 10.0,
            "mag": round(mag, 1),
        })
    for i in range(n_below):
        events.append({
            "id": "b%05d" % (i + 1),
            "time_utc": "%04d-%02d-%02dT%02d:00:00Z" % (start_year + int(rng.next_u() * years),
                                                       1 + int(rng.next_u() * 11),
                                                       1 + int(rng.next_u() * 27), int(rng.next_u() * 23)),
            "lon": 100.0, "lat": 25.0, "depth_km": 10.0,
            "mag": round(mc - width * (1 + int(rng.next_u() * 5)), 1),
        })
    events.sort(key=lambda e: e["time_utc"])
    return {"catalog": {"network": "t", "mc": mc, "magnitude_bin": width, "events": events},
            "forecast": {"target_magnitude": target, "horizon_years": horizon}}


def expect(data):
    """按定义式独立复算：只留 M>=mc，b 用档中心修正，a 年化，概率按泊松。"""
    cat = data["catalog"]
    mc, width = float(cat["mc"]), float(cat["magnitude_bin"])
    ev = [e for e in cat["events"] if float(e["mag"]) >= mc]
    mags = [float(e["mag"]) for e in ev]
    times = [dt.datetime.fromisoformat(e["time_utc"].replace("Z", "+00:00")) for e in ev]
    years = (max(times) - min(times)).total_seconds() / 86400.0 / 365.25
    n = len(mags)
    mean_m = sum(mags) / n
    b = math.log10(math.e) / (mean_m - (mc - width / 2.0))
    a = math.log10(n / years) + b * mc
    target = float(data["forecast"]["target_magnitude"])
    horizon = float(data["forecast"]["horizon_years"])
    rate = 10.0 ** (a - b * target)
    prob = 1.0 - math.exp(-rate * horizon)
    return {"n_selected": n, "years": years, "b": b, "sigma": b / math.sqrt(n),
            "a": a, "rate": rate, "prob": prob}


# ---------------------------------------------------------------- 驱动被测程序
def _run(cat_path, out_path, timeout=240):
    cmd = [sys.executable, ENTRY, "--catalog", cat_path, "--out", out_path]
    return subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8",
                          errors="replace", timeout=timeout)


def run_case(data):
    tmp = tempfile.mkdtemp(prefix="eq-case-")
    cpath = os.path.join(tmp, "catalog.json")
    opath = os.path.join(tmp, "report.json")
    with open(cpath, "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False)
    proc = _run(cpath, opath)
    report = None
    if os.path.exists(opath):
        try:
            with open(opath, encoding="utf-8") as f:
                report = json.load(f)
        except ValueError:
            report = None
    return report, {"cat": cpath, "out": opath, "proc": proc, "dir": tmp}


def _close(got, want, rel):
    try:
        v = float(got)
    except (TypeError, ValueError):
        return False
    if not math.isfinite(v):
        return False
    return abs(v - want) <= rel * abs(want)


# ---------------------------------------------------------------- 测试
def test_b_value_matches_aki_utsu_on_base_catalog():
    """题面原文：「b 值用 Aki-Utsu 极大似然，M0 取 mc 减半个档宽」。

    800 条目录（mc=2.5、档宽 0.1）：b 值与独立复算相差不超过 0.2%。
    """
    data = make_catalog(1.05, 2.5, 0.1, n=800, seed=101)
    report, _ = run_case(data)
    if not isinstance(report, dict):
        return False
    exp = expect(data)
    return _close(report.get("b_value"), exp["b"], 0.002) and report.get("n_selected") == exp["n_selected"]


def test_b_value_second_catalog_with_other_parameters():
    """题面原文：「复检时会喂别的目录，样本数、档宽与目标震级都不同」。

    换 b≈0.85、mc=3.0、档宽 0.05、目标 6.0 的目录后仍要吻合。
    """
    data = make_catalog(0.85, 3.0, 0.05, n=900, seed=211, target=6.0, horizon=20.0)
    report, _ = run_case(data)
    if not isinstance(report, dict):
        return False
    exp = expect(data)
    return _close(report.get("b_value"), exp["b"], 0.002) and _close(report.get("b_sigma"), exp["sigma"], 0.005)


def test_magnitude_bin_correction_matters():
    """题面原文：「M0 = mc - ΔM/2，用档中心而不是下边界」。

    档宽 0.2 的目录里，漏掉半个档宽的修正会让 b 偏低约 8%，本用例卡 0.3%。
    """
    data = make_catalog(1.1, 2.4, 0.2, n=900, seed=307)
    report, _ = run_case(data)
    if not isinstance(report, dict):
        return False
    exp = expect(data)
    return _close(report.get("b_value"), exp["b"], 0.003)


def test_a_value_and_annual_rate_match():
    """题面原文：「a 值按年化给出，log10 λ(M) = a - b·M」。

    a 值相差不超过 0.2%，目标震级年发生率相差不超过 1%。
    """
    data = make_catalog(1.0, 2.5, 0.1, n=1000, seed=409)
    report, _ = run_case(data)
    if not isinstance(report, dict):
        return False
    exp = expect(data)
    return (_close(report.get("a_value"), exp["a"], 0.002)
            and _close(report.get("annual_rate_target"), exp["rate"], 0.01))


def test_target_probability_matches_poisson():
    """题面原文：「T 年内至少一次的概率按 1 - exp(-λT)」。

    概率与独立复算相差不超过 1%。
    """
    data = make_catalog(1.15, 2.5, 0.1, n=900, seed=503, target=6.2, horizon=25.0)
    report, _ = run_case(data)
    if not isinstance(report, dict):
        return False
    exp = expect(data)
    return _close(report.get("probability_horizon"), exp["prob"], 0.01)


def test_catalog_truncated_at_completeness_magnitude():
    """题面原文：「只统计 M 大于等于 mc 的事件，mc 以下的小震要剔掉」。

    目录里混进 300 条 mc 以下的小震：n_selected 只数 mc 以上的，b 值也不能被带偏。
    """
    data = make_catalog(1.05, 2.5, 0.1, n=700, n_below=300, seed=601)
    report, _ = run_case(data)
    if not isinstance(report, dict):
        return False
    exp = expect(data)
    return (report.get("n_selected") == exp["n_selected"]
            and report.get("n_total") == len(data["catalog"]["events"])
            and _close(report.get("b_value"), exp["b"], 0.002))


def test_segmentation_null_when_sample_small():
    """题面原文：「样本不足 120 条时 segmentation 给 null」。

    110 条同一量级的目录：报告仍要出全其余字段，分段一栏必须是 null（不是空对象、不是省略）。
    """
    data = make_catalog(1.0, 2.5, 0.1, n=110, seed=1311)
    report, _ = run_case(data)
    if not isinstance(report, dict):
        return False
    if report.get("n_selected") != 110:
        return False
    if "segmentation" not in report:
        return False
    return report["segmentation"] is None and math.isfinite(float(report["b_value"]))


def test_zero_time_span_rejected():
    """题面原文：「时间窗不是正数、或观测时长算不出正数时以退出码 2 结束」。

    所有事件同一时刻：观测时长为零，必须拒收、不落文件、stderr 非空。
    """
    data = make_catalog(1.0, 2.5, 0.1, n=150, seed=1409)
    for e in data["catalog"]["events"]:
        e["time_utc"] = "2015-06-01T00:00:00Z"
    tmp = tempfile.mkdtemp(prefix="eq-span-")
    cpath = os.path.join(tmp, "catalog.json")
    opath = os.path.join(tmp, "report.json")
    with open(cpath, "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False)
    proc = _run(cpath, opath)
    return (proc.returncode == 2 and not os.path.exists(opath)
            and bool((proc.stderr or "").strip()))


def test_invalid_catalog_rejected_without_output():
    """题面原文：「目录不合法时以退出码 2 结束，不产出报告」。

    四种坏目录：事件不足 100 条、id 重复、震级不是数值、目标震级低于 mc。
    """
    good = make_catalog(1.0, 2.5, 0.1, n=200, seed=701)
    few = json.loads(json.dumps(good))
    few["catalog"]["events"] = few["catalog"]["events"][:60]
    dup = json.loads(json.dumps(good))
    dup["catalog"]["events"][1]["id"] = dup["catalog"]["events"][0]["id"]
    badmag = json.loads(json.dumps(good))
    badmag["catalog"]["events"][3]["mag"] = "M4.5"
    lowtarget = json.loads(json.dumps(good))
    lowtarget["forecast"]["target_magnitude"] = 1.0
    for data in (few, dup, badmag, lowtarget):
        tmp = tempfile.mkdtemp(prefix="eq-bad-")
        cpath = os.path.join(tmp, "catalog.json")
        opath = os.path.join(tmp, "report.json")
        with open(cpath, "w", encoding="utf-8") as f:
            json.dump(data, f, ensure_ascii=False)
        proc = _run(cpath, opath)
        if proc.returncode != 2 or os.path.exists(opath):
            return False
        if not (proc.stderr or "").strip():
            return False
        if [n for n in os.listdir(tmp) if n.endswith(".tmp")]:
            return False
    return True


def test_report_schema_and_rounding():
    """题面原文：「报告结构固定，b/a/sigma 保留三位小数，年发生率与概率保留四位」。

    校验顶层字段、位数与有限性，文件里不出现 NaN 或 Infinity。
    """
    data = make_catalog(1.0, 2.5, 0.1, n=300, seed=809)
    report, ctx = run_case(data)
    if not isinstance(report, dict):
        return False
    need = {"n_total", "n_selected", "mc", "magnitude_bin", "observation_years", "b_value",
            "b_sigma", "a_value", "annual_rate_target", "probability_horizon", "meta",
            "segmentation"}
    if not need.issubset(set(report)):
        return False
    for key in ("b_value", "b_sigma", "a_value"):
        v = float(report[key])
        if not math.isfinite(v) or abs(v - round(v, 3)) > 1e-12:
            return False
    for key in ("annual_rate_target", "probability_horizon"):
        v = float(report[key])
        if not math.isfinite(v) or abs(v - round(v, 4)) > 1e-12:
            return False
    seg = report.get("segmentation")
    if not isinstance(seg, dict):
        return False
    if abs(float(seg["split_years"]) - round(float(seg["split_years"]), 2)) > 1e-12:
        return False
    for key in ("b_before", "b_after", "delta_lnl"):
        v = float(seg[key])
        if not math.isfinite(v) or abs(v - round(v, 3)) > 1e-12:
            return False
    raw = open(ctx["out"], encoding="utf-8").read()
    return "NaN" not in raw and "Infinity" not in raw


def test_rerun_produces_identical_bytes():
    """题面原文：「同一份目录重复运行，报告逐字节一致」。

    连跑两次，两个报告文件字节完全相同。
    """
    data = make_catalog(1.05, 2.5, 0.1, n=500, seed=907)
    first, ctx1 = run_case(data)
    if first is None:
        return False
    body1 = open(ctx1["out"], "rb").read()
    second, ctx2 = run_case(data)
    if second is None:
        return False
    return body1 == open(ctx2["out"], "rb").read()


def test_failed_run_keeps_previous_report():
    """题面原文：「运行失败时不能覆盖上一次已经写好的报告」。

    先跑一次合格目录留下报告，再用 id 重复的坏目录跑：退出码 2，旧报告字节不变。
    """
    good = make_catalog(1.0, 2.5, 0.1, n=300, seed=1009)
    tmp = tempfile.mkdtemp(prefix="eq-keep-")
    cpath = os.path.join(tmp, "catalog.json")
    opath = os.path.join(tmp, "report.json")
    with open(cpath, "w", encoding="utf-8") as f:
        json.dump(good, f, ensure_ascii=False)
    first = _run(cpath, opath)
    if first.returncode != 0 or not os.path.exists(opath):
        return False
    before = open(opath, "rb").read()
    broken = json.loads(json.dumps(good))
    broken["catalog"]["events"][1]["id"] = broken["catalog"]["events"][0]["id"]
    with open(cpath, "w", encoding="utf-8") as f:
        json.dump(broken, f, ensure_ascii=False)
    second = _run(cpath, opath)
    if second.returncode != 2:
        return False
    return open(opath, "rb").read() == before


MIN_SEGMENT_EVENTS = 60


def _seg_lnl(mags, m0):
    n = len(mags)
    beta = 1.0 / (sum(mags) / n - m0)
    return n * (math.log(beta) - 1.0)


def segment_expected(data):
    """独立复算最优分段：两段对数似然之和最大，两段各至少 60 条。"""
    cat = data["catalog"]
    mc, width = float(cat["mc"]), float(cat["magnitude_bin"])
    ev = [e for e in cat["events"] if float(e["mag"]) >= mc]
    ev.sort(key=lambda e: e["time_utc"])
    mags = [float(e["mag"]) for e in ev]
    times = [dt.datetime.fromisoformat(e["time_utc"].replace("Z", "+00:00")) for e in ev]
    m0 = mc - width / 2.0
    n = len(mags)
    if n < 2 * MIN_SEGMENT_EVENTS:
        return None
    best = None
    for i in range(MIN_SEGMENT_EVENTS, n - MIN_SEGMENT_EVENTS + 1):
        val = _seg_lnl(mags[:i], m0) + _seg_lnl(mags[i:], m0)
        if best is None or val > best[1] + 1e-12:
            best = (i, val)
    i, val = best
    return {"split_years": (times[i] - times[0]).total_seconds() / 86400.0 / 365.25,
            "b_before": math.log10(math.e) / (sum(mags[:i]) / i - m0),
            "b_after": math.log10(math.e) / (sum(mags[i:]) / (n - i) - m0),
            "delta_lnl": val - _seg_lnl(mags, m0)}


def test_segment_split_matches_mle_search():
    """题面原文：「找分段时刻使两段对数似然之和最大，两段各至少 60 条」。

    目录在观测期中段有一次 b 值抬升：最优分段时刻要在独立复算的值附近 0.15 年以内，
    两段 b 值相差不超过 0.5%，对数似然增量相差不超过 0.05。
    """
    data = make_catalog(0.9, 2.5, 0.1, n=900, seed=1201, b_after=1.4, split_frac=0.5)
    report, _ = run_case(data)
    if not isinstance(report, dict):
        return False
    exp = segment_expected(data)
    seg = report.get("segmentation")
    if not isinstance(seg, dict) or exp is None:
        return False
    try:
        span = float(seg["split_years"])
        before = float(seg["b_before"])
        after = float(seg["b_after"])
        delta = float(seg["delta_lnl"])
    except (KeyError, TypeError, ValueError):
        return False
    return (abs(span - exp["split_years"]) <= 0.15
            and abs(before - exp["b_before"]) <= 0.005 * exp["b_before"]
            and abs(after - exp["b_after"]) <= 0.005 * exp["b_after"]
            and abs(delta - exp["delta_lnl"]) <= 0.05)


def main():
    total = 0.0
    for name, weight in WEIGHTS.items():
        fn = globals().get(name)
        try:
            passed = bool(fn and fn())
        except Exception as e:  # noqa: BLE001 —— 任何异常都按不通过计
            print("%s: 异常 %r" % (name, e), file=sys.stderr)
            passed = False
        total += weight * passed
        print("%s: %s (权重 %.2f)" % (name, "PASS" if passed else "FAIL", weight), file=sys.stderr)
    print("%.6f" % total)


if __name__ == "__main__":
    main()
