# Task 8 全链提速实测（A/B/A/B，原始 stdout 归档）

口径与可复现性说明
------------------
- **harness**：`.tmp_research/t8_time.py`（临时脚本，gitignored；与当前 API 一致，未修改）
- **口径**：hu101 生产 loader；只计 2D `solver.run`（入口桥 `build_coupled_annulus_inlet_provider`
  **预建**，不计入）；A/B/A/B 交错；`total_t=13743s`；`ny=40/nz=250`
- ⚠️ **机器争用说明**：本次测量期间本机同时有人类伙伴的并发跑批 ⇒ **绝对耗时不可跨机横比**
  （比行首报告偏高）。A/B/A/B 交错 + 各类取最小值使**比值**仍可信。
- **比值**：round1 `147.88 s → 65.54 s` = **2.26×**；round2 `144.48 s → 60.24 s` = **2.40×**；
  harness 结论行（取各类最小值）= **2.40×**。落在既报区间 2.31–2.44× 内，**复现成立**。
- 用途：Task 8 / Task 8 修复轮 1 的「全链 ≥2×」证据（⚠️-2 收口）；单解对照见
  `tests/contract/test_banded_solve_equivalence.py::test_timing_spsolve_vs_banded_single_solve`。

## 原始 stdout（逐字保留，含进度行）

```text
[口径] total_t=13742.7s
[D2DGA] 开始环空二维模拟 total_t=13743s nz=250 ny=40
  [D2DGA] 进度 0% (3s/13743s, dt=2.790s)
  [D2DGA] 进度 10% (1377s/13743s, dt=3.349s)
  [D2DGA] 进度 20% (2750s/13743s, dt=3.349s)
  [D2DGA] 进度 30% (4123s/13743s, dt=3.349s)
  [D2DGA] 进度 40% (5499s/13743s, dt=3.349s)
  [D2DGA] 进度 50% (6872s/13743s, dt=3.349s)
  [D2DGA] 进度 60% (8248s/13743s, dt=3.349s)
  [D2DGA] 进度 70% (9620s/13743s, dt=1.939s)
  [D2DGA] 进度 80% (10996s/13743s, dt=2.457s)
  [D2DGA] 进度 90% (12372s/13743s, dt=3.588s)
  [D2DGA] 进度 100% (13743s/13743s, dt=3.283s)
  round1 banded=False 2D=147.88s
[D2DGA] 开始环空二维模拟 total_t=13743s nz=250 ny=40
  [D2DGA] 进度 0% (3s/13743s, dt=2.790s)
  [D2DGA] 进度 10% (1377s/13743s, dt=3.349s)
  [D2DGA] 进度 20% (2750s/13743s, dt=3.349s)
  [D2DGA] 进度 30% (4123s/13743s, dt=3.349s)
  [D2DGA] 进度 40% (5499s/13743s, dt=3.349s)
  [D2DGA] 进度 50% (6872s/13743s, dt=3.349s)
  [D2DGA] 进度 60% (8248s/13743s, dt=3.349s)
  [D2DGA] 进度 70% (9620s/13743s, dt=1.939s)
  [D2DGA] 进度 80% (10996s/13743s, dt=2.457s)
  [D2DGA] 进度 90% (12369s/13743s, dt=3.585s)
  [D2DGA] 进度 100% (13743s/13743s, dt=3.332s)
  round1 banded=True  2D=65.54s
[D2DGA] 开始环空二维模拟 total_t=13743s nz=250 ny=40
  [D2DGA] 进度 0% (3s/13743s, dt=2.790s)
  [D2DGA] 进度 10% (1377s/13743s, dt=3.349s)
  [D2DGA] 进度 20% (2750s/13743s, dt=3.349s)
  [D2DGA] 进度 30% (4123s/13743s, dt=3.349s)
  [D2DGA] 进度 40% (5499s/13743s, dt=3.349s)
  [D2DGA] 进度 50% (6872s/13743s, dt=3.349s)
  [D2DGA] 进度 60% (8248s/13743s, dt=3.349s)
  [D2DGA] 进度 70% (9620s/13743s, dt=1.939s)
  [D2DGA] 进度 80% (10996s/13743s, dt=2.457s)
  [D2DGA] 进度 90% (12372s/13743s, dt=3.588s)
  [D2DGA] 进度 100% (13743s/13743s, dt=3.283s)
  round2 banded=False 2D=144.48s
[D2DGA] 开始环空二维模拟 total_t=13743s nz=250 ny=40
  [D2DGA] 进度 0% (3s/13743s, dt=2.790s)
  [D2DGA] 进度 10% (1377s/13743s, dt=3.349s)
  [D2DGA] 进度 20% (2750s/13743s, dt=3.349s)
  [D2DGA] 进度 30% (4123s/13743s, dt=3.349s)
  [D2DGA] 进度 40% (5499s/13743s, dt=3.349s)
  [D2DGA] 进度 50% (6872s/13743s, dt=3.349s)
  [D2DGA] 进度 60% (8248s/13743s, dt=3.349s)
  [D2DGA] 进度 70% (9620s/13743s, dt=1.939s)
  [D2DGA] 进度 80% (10996s/13743s, dt=2.457s)
  [D2DGA] 进度 90% (12369s/13743s, dt=3.585s)
  [D2DGA] 进度 100% (13743s/13743s, dt=3.332s)
  round2 banded=True  2D=60.24s

[结果] spsolve=144.48s  banded=60.24s  加速=2.40x
```
