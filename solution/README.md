# 难点与测试对应

## 难点 1：震级分档的半个档宽

推理轴：领域公式

落在推理路径上的哪里：Aki-Utsu 的 `b = log10(e)/(M̄ - M0)` 里，`M0` 取的是档中心 `mc - ΔM/2`，不是 `mc`。震级是分档记录的，数据里的 `mc` 落在档的下边界，直接用 `mc` 会把分母放大半个档宽，b 值系统性偏低（档宽 0.2 时约低 8%）。原实现就是这么写的，所以回测时比会商稿低一截。

对应测试：test_magnitude_bin_correction_matters、test_b_value_matches_aki_utsu_on_base_catalog

## 难点 2：目录截断与观测时长

推理轴：数据口径

落在推理路径上的哪里：评估只统计 `mag >= mc` 的事件，mc 以下的小震既不能进 b 值的样本，也不能参与观测时长（最短时刻被小震拉早，`T` 会变大，年化 a 值与年发生率跟着偏）。原实现用了全目录。

对应测试：test_catalog_truncated_at_completeness_magnitude、test_a_value_and_annual_rate_match

## 难点 3：a 值年化与泊松概率

推理轴：业务规则

落在推理路径上的哪里：`a = log10(N/T) + b·mc` 是年化的，`T` 用观测年数（不是天数、不是年数加一）；目标震级年发生率 `10^(a - b·M_target)`，时间窗概率 `1 - exp(-λT)`，两个指数都要用对。原实现既没年化，又把目标震级的公式写成了 `10^(a - b·mc)`。

对应测试：test_target_probability_matches_poisson、test_a_value_and_annual_rate_match

## 跨运行状态（重复运行与不污染）

同一份目录重复运行，报告逐字节一致；运行失败时不能覆盖上一次的报告：先写临时文件再原子替换，非法目录以退出码 2 收场且不落文件。对应测试：test_rerun_produces_identical_bytes、test_failed_run_keeps_previous_report。

第二套参数（b≈0.85、mc=3.0、档宽 0.05）另有 test_b_value_second_catalog_with_other_parameters 单独覆盖，防止把某一份目录的口径写死。

## 权重口径

10 条测试：基准 b 值 0.16、第二套参数 0.13、档宽修正 0.12、a 值与年发生率 0.10、泊松概率 0.09、目录截断 0.10、非法目录 0.12、输出契约 0.08、重复运行 0.05、失败不覆盖 0.05，和为 1。

## 定标依据（自检用）

判分侧独立按定义式复算：截断后取平均震级、b 值用档中心、a 值年化、概率按泊松，与被测程序走的是不同的代码路径。本地实测参考解 1.000000、缺陷版 0.000000；缺陷版在解析带时区的时间戳后与朴素时间相减，直接 `TypeError` 退出，所以空跑拿 0。
