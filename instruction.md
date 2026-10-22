# 区域台网地震活动性复核：b 值、a 值与目标震级概率

区域台网的年度会商要给出活动性参数：这套目录的 b 值（震级-频度关系的斜率）、年化 a 值，以及未来若干年目标震级的发震概率。目录由台网按统一规则汇总，震级按 0.1 档向下取整入库，另有一列完备震级 mc 写在目录元数据里。

上一班把复核程序 `src/seismic.py` 接上今年的目录时，它跑出一行报错就退了；拿去年一份已经人工算过的目录回测，程序给的 b 值比会商稿低了一截，目标震级的年发生率也偏高一个量级。

## 现状

- 程序骨架已经在了，命令行入口固定：

  ```
  python3 src/seismic.py --catalog 目录JSON路径 --out 报告JSON路径
  ```

- 目录结构：`catalog.mc` 是完备震级，`catalog.magnitude_bin` 是震级档宽，`catalog.events` 每条有 `id` / `time_utc` / `lon` / `lat` / `depth_km` / `mag`。目录里会混有 mc 以下的小震，评估只统计 mc 以上的事件。
- `forecast.target_magnitude` 与 `forecast.horizon_years` 是这次会商要报的目标震级与时间窗。

## 任务

修好 `src/seismic.py`：

1. 只统计 `mag >= mc` 的事件；观测时长按这些事件的最早与最晚时刻计算。
2. b 值用 Aki-Utsu 极大似然：`b = log10(e) / (M̄ - M0)`，`M0 = mc - ΔM/2`（震级是分档记录的，M0 取档中心）。
3. b 值标准差按 `b / √N`。
4. a 值按年化给出：`a = log10(N/T) + b·mc`，其中 `T` 是观测年数。
5. 目标震级年发生率 `λ = 10^(a - b·M_target)`；时间窗内至少一次的概率按泊松过程 `1 - exp(-λT_h)`。

## 输出契约

报告是 UTF-8 编码的 JSON 对象：

```json
{
  "n_total": 2680,
  "n_selected": 2640,
  "mc": 2.5,
  "magnitude_bin": 0.1,
  "observation_years": 11.99,
  "b_value": 1.026,
  "b_sigma": 0.02,
  "a_value": 4.908,
  "annual_rate_target": 0.0173,
  "probability_horizon": 0.4046,
  "meta": {"method": "aki-utsu", "target_magnitude": 6.5, "horizon_years": 30.0}
}
```

- `n_total` 是目录里的记录条数（含 mc 以下），`n_selected` 是参与统计的条数。
- `observation_years` 保留两位小数；`b_value` / `b_sigma` / `a_value` 三位；`annual_rate_target` / `probability_horizon` 四位。
- `meta` 里照抄输入的目标震级与时间窗；数值必须是有限值，不许出现 `NaN` 或 `Infinity`。

## 运行与容错

- 正常跑完退出码 0；目录不合法时以**退出码 2** 结束、往 stderr 写原因、**不产出报告**。不合法包括：事件少于 100 条、id 重复或缺失、震级不是数值、时间无法解析、目标震级低于 mc、时间窗不是正数。
- 同一份目录重复运行，报告逐字节一致。
- 运行失败时不能覆盖上一次已经写好的报告。

## 约束

- 容器里只有 Python 标准库，不联网装包；只改 `src/seismic.py`，命令行入口与 JSON 字段名不要动。
- 报告只能由当前目录算出来，不得写死样例数值；复检会喂别的台网目录，样本量、档宽、mc 与目标震级都不同。
