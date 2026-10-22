#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""生成地震目录（Gutenberg-Richter 分布）+ 复核输入。

数据口径：业务模型 + 确定性生成规则 —— 震级服从 G-R（b 已知），时间在观测窗口内均匀分布，
震级按 0.1 分档四舍五入（真实目录就是这样报的），并混入少量 mc 以下的小震当干扰。
固定伪随机序列，不用真随机，整份目录可逐字节复现。
"""
import json
import math
import sys


class Lcg(object):
    def __init__(self, seed):
        self.state = seed & 0xFFFFFFFF

    def next_u(self):
        self.state = (1103515245 * self.state + 12345) & 0x7FFFFFFF
        return (self.state % 1000003) / 1000003.0


def make_catalog(b_value=1.05, mc=2.5, bin_width=0.1, rate_above_mc=220.0, years=12.0,
                 n_below=40, seed=20261023, target_magnitude=6.5, horizon_years=30.0,
                 start_year=2011):
    rng = Lcg(seed)
    n = int(round(rate_above_mc * years))
    beta = b_value * math.log(10.0)
    events = []
    t0 = start_year
    for i in range(n):
        u = rng.next_u()
        # 逆变换：M = mc + (-ln(1-u))/beta
        mag = mc + (-math.log(1.0 - u * 0.999999)) / beta
        mag = math.floor(mag / bin_width) * bin_width   # 震级按档向下取整入库
        years_off = rng.next_u() * years
        months = int(years_off * 12) % 12
        day = int(years_off * 365) % 28 + 1
        hour = int(years_off * 24) % 24
        events.append({
            "id": "e%05d" % (i + 1),
            "time_utc": "%04d-%02d-%02dT%02d:00:00Z" % (t0 + int(years_off), months + 1, day, hour),
            "lon": round(100.0 + rng.next_u() * 4.0, 2),
            "lat": round(24.0 + rng.next_u() * 3.0, 2),
            "depth_km": round(5.0 + rng.next_u() * 25.0, 1),
            "mag": round(mag, 1),
        })
    for i in range(n_below):
        mag = mc - bin_width * (1 + int(rng.next_u() * 6))
        events.append({
            "id": "s%05d" % (i + 1),
            "time_utc": "%04d-%02d-%02dT%02d:00:00Z" % (t0 + int(rng.next_u() * years), 1 + int(rng.next_u() * 11),
                                                       1 + int(rng.next_u() * 27), int(rng.next_u() * 23)),
            "lon": round(100.0 + rng.next_u() * 4.0, 2),
            "lat": round(24.0 + rng.next_u() * 3.0, 2),
            "depth_km": round(5.0 + rng.next_u() * 25.0, 1),
            "mag": round(mag, 1),
        })
    events.sort(key=lambda e: e["time_utc"])
    return {
        "catalog": {
            "network": "区域台网",
            "mc": mc,
            "magnitude_bin": bin_width,
            "events": events,
        },
        "forecast": {"target_magnitude": target_magnitude, "horizon_years": horizon_years},
    }


if __name__ == "__main__":
    out = sys.argv[1] if len(sys.argv) > 1 else "catalog.json"
    data = make_catalog()
    with open(out, "w", encoding="utf-8", newline="\n") as f:
        json.dump(data, f, ensure_ascii=False, separators=(",", ":"))
        f.write("\n")
    ev = data["catalog"]["events"]
    print("wrote %s  events=%d  M>=mc: %d" % (out, len(ev), sum(1 for e in ev if e["mag"] >= data["catalog"]["mc"])))
