# runner 链 T-on 试点（D8-A，2026-10-06）

口径：cemdisp/runners/<井>.py 构造的忠实镜像（1D T1 生产开关、F2 停算、2D 纯默认 nz=250）；Ton_static = enable_temperature_rheology=True + temperature_mode=static（无瞬态表井取静温剖面）。
⚠️ 本试点按 Phase 0.0 密度「就近取」口径（Q16）执行。

| 井 | zero η_E | zero η_N | Ton_static η_E | Ton_static η_N | ΔTon−zero η_N/pp | 权威 η_N | zero−权威 η_N/pp | 判定 |
|---|---|---|---|---|---|---|---|---|
| 呼101 | 0.9735 | 0.8306 | 0.9621 | 0.7803 | -5.030 | 0.8306 | +0.000 | PASS(|Δ|≤0.05pp) |
| 呼1-004 | 0.9997 | 0.9983 | 0.9998 | 0.9989 | +0.061 | 0.9983 | +0.000 | PASS(|Δ|≤0.05pp) |
| 呼1-003 | 0.9986 | 0.9932 | 0.9989 | 0.9945 | +0.134 | 0.9932 | +0.000 | PASS(|Δ|≤0.05pp) |