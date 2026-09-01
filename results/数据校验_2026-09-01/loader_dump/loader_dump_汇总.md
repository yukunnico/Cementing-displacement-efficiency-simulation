# loader_dump 校验汇总（2026-09-01，脚本自动生成）

数据源：cemdisp loader 默认入口 × 现场提取层CSV（2026-08-29校准轮）+ 01_总表对照表。

## 检查条目统计

| 井 | ✓ | ✗ | ⚠ | - |
|---|---|---|---|---|
| hu101 | 61 | 17 | 23 | 8 |
| hu102 | 35 | 3 | 18 | 8 |
| hu103 | 75 | 10 | 26 | 9 |
| hu103[variant] | 27 | 0 | 6 | 0 |
| hu1 | 28 | 0 | 21 | 0 |
| hu2 | 78 | 3 | 30 | 2 |
| ht1_001 | 72 | 12 | 43 | 13 |
| ht1_003 | 87 | 4 | 32 | 2 |
| ht1_003[variant] | 25 | 0 | 4 | 0 |
| ht1_004 | 57 | 22 | 38 | 17 |
| ht1_004[variant] | 28 | 0 | 3 | 0 |

## hu101（呼101）

- 模型域 5400.0–7868.0 m；鞋 7868.0 m；悬挂器 5407.46 m；liner OD/ID 139.7/91.7328 mm
- cbl_pass_rate = 0.6277
- 日程总注入 221.4 m³，9 步
- 环空体积重算 65.5 m³；浆体 70.0 m³；库存比 1.07
- 管内容积重算 16.3 m³；顶替流体合计 101.4 m³

### hu101 需关注条目（✗/⚠，共 40 条）

| 类别 | 条目 | loader | CSV | 状态 | 备注 |
|---|---|---|---|---|---|
| fluid | displacement/中置液[第2对2/4] 密度kg/m3 | 1960 | 1850 | ✗ | CSV g/cm3×1000；容差10kg/m3 |
| fluid | displacement/井浆[第3对3/4] 密度kg/m3 | 1960 | 2000 | ✗ | CSV g/cm3×1000；容差10kg/m3 |
| fluid | displacement/井浆[第3对3/4] PV Pa·s | 0.058 | 0.03 | ✗ | CSV mPa·s/1000；容差2% |
| fluid | role=lead[第1对1/2] n | 0.844 | 0.719 | ✗ | 容差2% |
| fluid | role=lead[第1对1/2] K Pa·s^n | 0.381 | 0.815 | ✗ | 容差3% |
| fluid | role=tail[第1对1/2] n | 0.83 | 0.722 | ✗ | 容差2% |
| fluid | role=tail[第1对1/2] K Pa·s^n | 0.352 | 0.684 | ✗ | 容差3% |
| schedule | role=displacement[1] 体积m3(注后置液(管内)) | 2 | 26 | ✗ | 容差0.2m³/1% |
| schedule | role=displacement[1] 排量m3/min | 0.6 | 1.5 | ✗ | 容差2% |
| schedule | role=displacement[1] 开始min(推算vs现场) | 115.8 | 119 | ⚠ | loader时间为体积/排量推算；容差3min |
| schedule | role=displacement[1] 结束min(推算vs现场) | 119.2 | 137 | ⚠ | 容差3min |
| schedule | role=displacement[2] 体积m3(注轻泥浆) | 26 | 10 | ✗ | 容差0.2m³/1% |
| schedule | role=displacement[2] 排量m3/min | 1.5 | 1 | ✗ | 容差2% |
| schedule | role=displacement[2] 开始min(推算vs现场) | 119.2 | 137 | ⚠ | loader时间为体积/排量推算；容差3min |
| schedule | role=displacement[2] 结束min(推算vs现场) | 136.5 | 147 | ⚠ | 容差3min |
| schedule | role=displacement[3] 体积m3(注中置液) | 10 | 63.4 | ✗ | 容差0.2m³/1% |
| schedule | role=displacement[3] 开始min(推算vs现场) | 136.5 | 147 | ⚠ | loader时间为体积/排量推算；容差3min |
| schedule | role=displacement[3] 结束min(推算vs现场) | 146.5 | 217 | ⚠ | 容差3min |
| schedule | role=displacement[4] 体积m3(井浆快替) | 40 | 85 | ✗ | 容差0.2m³/1% |
| schedule | role=displacement[4] 排量m3/min | 1 | 2.4 | ✗ | 容差2% |
| schedule | role=displacement[4] 开始min(推算vs现场) | 146.5 | 162 | ⚠ | loader时间为体积/排量推算；容差3min |
| schedule | role=displacement[4] 结束min(推算vs现场) | 186.5 | 200 | ⚠ | 容差3min |
| schedule | role=displacement[5] 体积m3(井浆慢替) | 23.4 | 85 | ✗ | 容差0.2m³/1% |
| schedule | role=displacement[5] 排量m3/min | 0.55 | 2.4 | ✗ | 容差2% |
| schedule | role=displacement[5] 开始min(推算vs现场) | 186.5 | 162 | ⚠ | loader时间为体积/排量推算；容差3min |
| schedule | role=displacement[5] 结束min(推算vs现场) | 229 | 200 | ⚠ | 容差3min |
| schedule | role=spacer[1] 排量m3/min | 1 | 1.2 | ✗ | 容差2% |
| schedule | role=spacer[1] 开始min(推算vs现场) | 20.8 | 0 | ⚠ | loader时间为体积/排量推算；容差3min |
| schedule | role=spacer[1] 结束min(推算vs现场) | 45.8 | 21 | ⚠ | 容差3min |
| schedule | role=wash | 1步 | 0步 | ⚠ | 单侧无该角色步骤（口径/序列差异） |
| schedule | 日程总时长min(推算vs现场末步) | 229 | 217 | ⚠ | 现场末步含试压/停泵等非注入步，差异仅提示 |
| geometry | well_geometry.csv 同item多值行 | — | float_collar_md,liner_bottom_md,liner_inner_diameter,liner_outer_diameter,liner_top_md | ⚠ | CSV长表混入回接段同名item，比对取首现（尾管段）值 |
| geometry | 井径剖面点数(CSV实测 vs loader模型域) | 110 | 109 | ⚠ | loader仅保留模型段内点并可能做等效换算 |
| geometry | 井径max mm(CSV vs loader) | 258.2 | 265.2 | ⚠ | 上段等效换算会改变数值 |
| cbl | CBL测量段 top/bottom | 窗['5700-7810'] | 5390–7810 | ⚠ | loader评价窗常为可评价子段（扣除双层套管/悬空段），见notes |
| cbl | 地层目标[CBL评价重点-裸眼段] 5700-7868 | 未覆盖 | target_intervals.csv | ⚠ |  |
| 内部 | 环空体积重算m3(hole_profile×OD积分) | 65.53 | — | ⚠ | 重算口径：hole_diameter_profile×(双径)liner OD；技套重叠段按外推井径计 |
| 内部 | 水泥浆总体积/环空体积(η_E库存比上限) | 1.068 | — | ⚠ | >1 表示水泥浆量超过环空容量（返高以上/漏失/计量口径） |
| 内部 | 管内容积重算m3 vs 替浆总体积m3 | 16.3 vs 101.4 | — | ⚠ | 替浆总量应≈碰压时管内容积；差异大提示 liner_id/pipe_id 口径为等效或代理 |
| sanity | 密度范围 领浆 | 2100 | 水泥浆1.85–1.95g/cc | ⚠ | 超出常规水泥浆范围：加重/低密配方需与配方记录核对 |

## hu102（呼102）

- 模型域 6823.1–7735.0 m；鞋 7735.0 m；悬挂器 6823.1 m；liner OD/ID 139.7/108.1 mm
- cbl_pass_rate = 0.6665
- 日程总注入 180.0 m³，9 步
- 环空体积重算 12.4 m³；浆体 17.0 m³；库存比 1.37
- 管内容积重算 8.4 m³；顶替流体合计 74.0 m³

### hu102 需关注条目（✗/⚠，共 21 条）

| 类别 | 条目 | loader | CSV | 状态 | 备注 |
|---|---|---|---|---|---|
| fluid | role=? | 0条 | 70条 | ⚠ | loader或CSV侧该角色无对应条目（角色映射或口径差异） |
| fluid | role=displacement | 1条 | 0条 | ⚠ | loader或CSV侧该角色无对应条目（角色映射或口径差异） |
| fluid | role=flusher | 1条 | 0条 | ⚠ | loader或CSV侧该角色无对应条目（角色映射或口径差异） |
| fluid | role=lead | 1条 | 0条 | ⚠ | loader或CSV侧该角色无对应条目（角色映射或口径差异） |
| fluid | role=mud | 1条 | 0条 | ⚠ | loader或CSV侧该角色无对应条目（角色映射或口径差异） |
| fluid | role=other | 2条 | 0条 | ⚠ | loader或CSV侧该角色无对应条目（角色映射或口径差异） |
| fluid | role=spacer | 1条 | 0条 | ⚠ | loader或CSV侧该角色无对应条目（角色映射或口径差异） |
| fluid | role=tail | 1条 | 0条 | ⚠ | loader或CSV侧该角色无对应条目（角色映射或口径差异） |
| fluid | role=wash | 2条 | 0条 | ⚠ | loader或CSV侧该角色无对应条目（角色映射或口径差异） |
| schedule | role=displacement[1] 体积m3(替浆（钻井液）) | 74 | 7 | ✗ | 容差0.2m³/1% |
| schedule | role=mud[1] 体积m3(循环排混浆) | 41 | 72 | ✗ | 容差0.2m³/1% |
| schedule | role=other | 2步 | 0步 | ⚠ | 单侧无该角色步骤（口径/序列差异） |
| schedule | role=保护液 | 0步 | 1步 | ⚠ | 单侧无该角色步骤（口径/序列差异） |
| geometry | 悬挂器深度hanger_md_m | 6823.1 | 7665.52 | ✗ | well_geometry.csv；容差0.5m |
| geometry | 井径剖面点数(CSV实测 vs loader模型域) | 63 | 62 | ⚠ | loader仅保留模型段内点并可能做等效换算 |
| cbl | 地层目标[裸眼段] 6720-7735 | 未覆盖 | target_intervals.csv | ⚠ |  |
| 内部 | 环空体积重算m3(hole_profile×OD积分) | 12.38 | — | ⚠ | 重算口径：hole_diameter_profile×(双径)liner OD；技套重叠段按外推井径计 |
| 内部 | 水泥浆总体积/环空体积(η_E库存比上限) | 1.374 | — | ⚠ | >1 表示水泥浆量超过环空容量（返高以上/漏失/计量口径） |
| 内部 | 管内容积重算m3 vs 替浆总体积m3 | 8.4 vs 74.0 | — | ⚠ | 替浆总量应≈碰压时管内容积；差异大提示 liner_id/pipe_id 口径为等效或代理 |
| sanity | 密度范围 领浆 | 2100 | 水泥浆1.85–1.95g/cc | ⚠ | 超出常规水泥浆范围：加重/低密配方需与配方记录核对 |
| sanity | 密度范围 尾管水泥浆 | 2100 | 水泥浆1.85–1.95g/cc | ⚠ | 超出常规水泥浆范围：加重/低密配方需与配方记录核对 |

## hu103（呼103）

- 模型域 5536.662–7770.0 m；鞋 7770.0 m；悬挂器 5545.97 m；liner OD/ID 139.7/107.94 mm；双径上段 168.3/138.9 mm 底 7330.694 m
- cbl_pass_rate = 0.1206
- 日程总注入 236.2 m³，10 步
- 环空体积重算 75.0 m³；浆体 83.0 m³；库存比 1.11
- 管内容积重算 20.4 m³；顶替流体合计 83.2 m³

### hu103 需关注条目（✗/⚠，共 36 条）

| 类别 | 条目 | loader | CSV | 状态 | 备注 |
|---|---|---|---|---|---|
| fluid | displacement/中置液 密度kg/m3 | 1950 | 2100 | ✗ | CSV g/cm3×1000；容差10kg/m3 |
| fluid | displacement/替浆液 密度kg/m3 | 1980 | 1950 | ✗ | CSV g/cm3×1000；容差10kg/m3 |
| fluid | role=intermediate | 1条 | 0条 | ⚠ | loader或CSV侧该角色无对应条目（角色映射或口径差异） |
| fluid | role=intermediate_cement | 0条 | 1条 | ⚠ | loader或CSV侧该角色无对应条目（角色映射或口径差异） |
| fluid | mud/钻井液[第1对1/6] 流变模型 | bingham | power_law | ⚠ |  |
| fluid | mud/平衡液[第2对2/6] 密度kg/m3 | 1880 | 1980 | ✗ | CSV g/cm3×1000；容差10kg/m3 |
| fluid | mud/平衡液[第2对2/6] PV Pa·s | 0.025 | 0.088 | ✗ | CSV mPa·s/1000；容差2% |
| fluid | mud/平衡液[第2对2/6] YP Pa | 1.5 | 15.5 | ✗ | 容差2% |
| fluid | role=other | 1条 | 0条 | ⚠ | loader或CSV侧该角色无对应条目（角色映射或口径差异） |
| fluid | spacer/隔离液(实际) 密度kg/m3 | 1920 | 2000 | ✗ | CSV g/cm3×1000；容差10kg/m3 |
| fluid | 隔离液1 n(六速反算参考) | 0.54 | 0.791 | ⚠ | 140℃ 六速θ600/θ300=128.0/74.0；反算仅供温度口径参考 |
| fluid | 隔离液1 K(六速反算参考) | 2.12 | 0.273 | ⚠ | 反算 K=0.511·θ300/(1.7023·300)^n |
| fluid | 隔离液2 n(六速反算参考) | 0.54 | 0.791 | ⚠ | 140℃ 六速θ600/θ300=128.0/74.0；反算仅供温度口径参考 |
| fluid | 隔离液2 K(六速反算参考) | 2.12 | 0.273 | ⚠ | 反算 K=0.511·θ300/(1.7023·300)^n |
| fluid | 隔离液(实际) n(六速反算参考) | 0.54 | 0.791 | ⚠ | 140℃ 六速θ600/θ300=128.0/74.0；反算仅供温度口径参考 |
| fluid | 隔离液(实际) K(六速反算参考) | 2.12 | 0.273 | ⚠ | 反算 K=0.511·θ300/(1.7023·300)^n |
| schedule | role=displacement[1] 体积m3(替钻井液(一段)) | 23 | 7 | ✗ | 容差0.2m³/1% |
| schedule | role=displacement[2] 体积m3(替中置液) | 15 | 23 | ✗ | 容差0.2m³/1% |
| schedule | role=displacement[3] 体积m3(替钻井液(二段)) | 45.2 | 15 | ✗ | 容差0.2m³/1% |
| schedule | role=intermediate | 1步 | 0步 | ⚠ | 单侧无该角色步骤（口径/序列差异） |
| schedule | role=intermediate_cement | 0步 | 1步 | ⚠ | 单侧无该角色步骤（口径/序列差异） |
| schedule | role=other | 1步 | 0步 | ⚠ | 单侧无该角色步骤（口径/序列差异） |
| geometry | well_geometry.csv 同item多值行 | — | borehole_nominal_diameter,drill_bit_size,float_collar_md,liner_inner_diameter,liner_outer_diameter | ⚠ | CSV长表混入回接段同名item，比对取首现（尾管段）值 |
| geometry | 井径剖面点数(CSV实测 vs loader模型域) | 70 | 69 | ⚠ | loader仅保留模型段内点并可能做等效换算 |
| cbl | CBL测量段 top/bottom | 窗['7338-7712', '5540-7331'] | 7338–7712 | ⚠ | loader评价窗常为可评价子段（扣除双层套管/悬空段），见notes |
| cbl | 地层目标[渗漏点] 7652.6-7652.6 | 未覆盖 | target_intervals.csv | ⚠ |  |
| cbl | 地层目标[领浆封固段] 5350-6000 | 未覆盖 | target_intervals.csv | ⚠ |  |
| cbl | 地层目标[中间浆封固段] 6000-7000 | 未覆盖 | target_intervals.csv | ⚠ |  |
| cbl | 地层目标[尾浆封固段] 7000-7770 | 未覆盖 | target_intervals.csv | ⚠ |  |
| 总表 | 总表[program/tail_volume_m3] | 26 | 19 | ✗ | 含下塞1.18m3; intermediate slurry remains proxy |
| 内部 | 环空体积重算m3(hole_profile×OD积分) | 75.01 | — | ⚠ | 重算口径：hole_diameter_profile×(双径)liner OD；技套重叠段按外推井径计 |
| 内部 | 水泥浆总体积/环空体积(η_E库存比上限) | 1.107 | — | ⚠ | >1 表示水泥浆量超过环空容量（返高以上/漏失/计量口径） |
| 内部 | 管内容积重算m3 vs 替浆总体积m3 | 20.4 vs 83.2 | — | ⚠ | 替浆总量应≈碰压时管内容积；差异大提示 liner_id/pipe_id 口径为等效或代理 |
| sanity | 密度范围 领浆 | 2050 | 水泥浆1.85–1.95g/cc | ⚠ | 超出常规水泥浆范围：加重/低密配方需与配方记录核对 |
| sanity | 密度范围 中间浆 | 2050 | 水泥浆1.85–1.95g/cc | ⚠ | 超出常规水泥浆范围：加重/低密配方需与配方记录核对 |
| sanity | 密度范围 尾浆 | 2050 | 水泥浆1.85–1.95g/cc | ⚠ | 超出常规水泥浆范围：加重/低密配方需与配方记录核对 |

## hu103[variant]（呼103）

- 模型域 5536.662–7770.0 m；鞋 7770.0 m；悬挂器 5545.97 m；liner OD/ID 139.7/107.94 mm；双径上段 168.3/138.9 mm 底 7330.694 m
- cbl_pass_rate = 0.1206
- 日程总注入 259.9 m³，7 步
- 环空体积重算 75.0 m³；浆体 83.0 m³；库存比 1.11
- 管内容积重算 20.4 m³；顶替流体合计 90.9 m³

### hu103[variant] 需关注条目（✗/⚠，共 6 条）

| 类别 | 条目 | loader | CSV | 状态 | 备注 |
|---|---|---|---|---|---|
| 内部 | 环空体积重算m3(hole_profile×OD积分) | 75.01 | — | ⚠ | 重算口径：hole_diameter_profile×(双径)liner OD；技套重叠段按外推井径计 |
| 内部 | 水泥浆总体积/环空体积(η_E库存比上限) | 1.107 | — | ⚠ | >1 表示水泥浆量超过环空容量（返高以上/漏失/计量口径） |
| 内部 | 管内容积重算m3 vs 替浆总体积m3 | 20.4 vs 90.9 | — | ⚠ | 替浆总量应≈碰压时管内容积；差异大提示 liner_id/pipe_id 口径为等效或代理 |
| sanity | 密度范围 领浆 | 2050 | 水泥浆1.85–1.95g/cc | ⚠ | 超出常规水泥浆范围：加重/低密配方需与配方记录核对 |
| sanity | 密度范围 中间浆 | 2050 | 水泥浆1.85–1.95g/cc | ⚠ | 超出常规水泥浆范围：加重/低密配方需与配方记录核对 |
| sanity | 密度范围 尾浆 | 2050 | 水泥浆1.85–1.95g/cc | ⚠ | 超出常规水泥浆范围：加重/低密配方需与配方记录核对 |

## hu1（呼探1井）

- 模型域 3523.27–7601.0 m；鞋 7601.0 m；悬挂器 3523.27 m；liner OD/ID 139.7/111.16 mm
- 日程总注入 177.0 m³，7 步
- 环空体积重算 59.9 m³；浆体 60.0 m³；库存比 1.00
- 管内容积重算 41.3 m³；顶替流体合计 85.0 m³

### hu1 需关注条目（✗/⚠，共 21 条）

| 类别 | 条目 | loader | CSV | 状态 | 备注 |
|---|---|---|---|---|---|
| fluid | role=? | 0条 | 28条 | ⚠ | loader或CSV侧该角色无对应条目（角色映射或口径差异） |
| fluid | role=displacement | 2条 | 0条 | ⚠ | loader或CSV侧该角色无对应条目（角色映射或口径差异） |
| fluid | role=lead | 1条 | 0条 | ⚠ | loader或CSV侧该角色无对应条目（角色映射或口径差异） |
| fluid | role=mud | 1条 | 0条 | ⚠ | loader或CSV侧该角色无对应条目（角色映射或口径差异） |
| fluid | role=other | 1条 | 0条 | ⚠ | loader或CSV侧该角色无对应条目（角色映射或口径差异） |
| fluid | role=spacer | 1条 | 0条 | ⚠ | loader或CSV侧该角色无对应条目（角色映射或口径差异） |
| fluid | role=tail | 1条 | 0条 | ⚠ | loader或CSV侧该角色无对应条目（角色映射或口径差异） |
| fluid | role=wash | 1条 | 0条 | ⚠ | loader或CSV侧该角色无对应条目（角色映射或口径差异） |
| schedule | role=displacement | 2步 | 0步 | ⚠ | 单侧无该角色步骤（口径/序列差异） |
| schedule | role=flusher | 0步 | 1步 | ⚠ | 单侧无该角色步骤（口径/序列差异） |
| schedule | role=mud | 0步 | 1步 | ⚠ | 单侧无该角色步骤（口径/序列差异） |
| schedule | role=other | 1步 | 0步 | ⚠ | 单侧无该角色步骤（口径/序列差异） |
| schedule | role=wash | 1步 | 0步 | ⚠ | 单侧无该角色步骤（口径/序列差异） |
| schedule | 日程总时长min(推算vs现场末步) | 349.1 | 240 | ⚠ | 现场末步含试压/停泵等非注入步，差异仅提示 |
| geometry | well_geometry.csv 同item多值行 | — |  | ⚠ | CSV长表混入回接段同名item，比对取首现（尾管段）值 |
| geometry | 井径剖面点数(CSV实测 vs loader模型域) | 101 | 96 | ⚠ | loader仅保留模型段内点并可能做等效换算 |
| cbl | 地层目标[365.1mm技套固井段] 1000-3036 | 未覆盖 | target_intervals.csv | ⚠ |  |
| 内部 | 环空体积重算m3(hole_profile×OD积分) | 59.9 | — | ⚠ | 重算口径：hole_diameter_profile×(双径)liner OD；技套重叠段按外推井径计 |
| 内部 | 水泥浆总体积/环空体积(η_E库存比上限) | 1.002 | — | ⚠ | >1 表示水泥浆量超过环空容量（返高以上/漏失/计量口径） |
| 内部 | 管内容积重算m3 vs 替浆总体积m3 | 41.3 vs 85.0 | — | ⚠ | 替浆总量应≈碰压时管内容积；差异大提示 liner_id/pipe_id 口径为等效或代理 |
| sanity | 密度范围 领浆 | 2100 | 水泥浆1.85–1.95g/cc | ⚠ | 超出常规水泥浆范围：加重/低密配方需与配方记录核对 |

## hu2（呼探1-002井（HT1-002））

- 模型域 5292.5–7554.0 m；鞋 7554.0 m；悬挂器 5292.5 m；liner OD/ID 139.7/107.94 mm
- 日程总注入 152.0 m³，10 步
- 环空体积重算 33.5 m³；浆体 38.0 m³；库存比 1.13
- 管内容积重算 20.7 m³；顶替流体合计 74.0 m³

### hu2 需关注条目（✗/⚠，共 33 条）

| 类别 | 条目 | loader | CSV | 状态 | 备注 |
|---|---|---|---|---|---|
| fluid | role=displacement | 3条 | 0条 | ⚠ | loader或CSV侧该角色无对应条目（角色映射或口径差异） |
| fluid | role=intermediate | 1条 | 0条 | ⚠ | loader或CSV侧该角色无对应条目（角色映射或口径差异） |
| fluid | role=intermediate_cement | 0条 | 2条 | ⚠ | loader或CSV侧该角色无对应条目（角色映射或口径差异） |
| fluid | role=lead[第1对1/2] 流变模型 | power_law | herschel_bulkley | ⚠ |  |
| fluid | role=mud 流变模型 | bingham | herschel_bulkley | ⚠ |  |
| fluid | role=other | 1条 | 0条 | ⚠ | loader或CSV侧该角色无对应条目（角色映射或口径差异） |
| fluid | role=plug | 0条 | 1条 | ⚠ | loader或CSV侧该角色无对应条目（角色映射或口径差异） |
| fluid | role=spacer[第1对1/4] 流变模型 | power_law | herschel_bulkley | ⚠ |  |
| fluid | role=tail[第1对1/2] 流变模型 | power_law | herschel_bulkley | ⚠ |  |
| fluid | 隔离液 n(六速反算参考) | 0.545 | 0.842 | ⚠ | 94℃ 六速θ600/θ300=95.0/53.0；反算仅供温度口径参考 |
| fluid | 隔离液 K(六速反算参考) | 1.338 | 0.142 | ⚠ | 反算 K=0.511·θ300/(1.7023·300)^n |
| schedule | role=displacement[1] 开始min(推算vs现场) | 97.5 | 101 | ⚠ | loader时间为体积/排量推算；容差3min |
| schedule | role=displacement[1] 结束min(推算vs现场) | 112.5 | 117 | ⚠ | 容差3min |
| schedule | role=displacement[2] 开始min(推算vs现场) | 112.5 | 137 | ⚠ | loader时间为体积/排量推算；容差3min |
| schedule | role=displacement[2] 结束min(推算vs现场) | 131.2 | 157 | ⚠ | 容差3min |
| schedule | role=displacement[3] 体积m3(井浆快替) | 30 | 5 | ✗ | 容差0.2m³/1% |
| schedule | role=displacement[3] 排量m3/min | 0.8 | 0.5 | ✗ | 容差2% |
| schedule | role=displacement[3] 开始min(推算vs现场) | 131.2 | 157 | ⚠ | loader时间为体积/排量推算；容差3min |
| schedule | role=displacement[4] 体积m3(井浆慢替) | 17 | 18 | ✗ | 容差0.2m³/1% |
| schedule | role=intermediate | 1步 | 0步 | ⚠ | 单侧无该角色步骤（口径/序列差异） |
| schedule | role=intermediate_cement | 0步 | 1步 | ⚠ | 单侧无该角色步骤（口径/序列差异） |
| schedule | role=other | 1步 | 0步 | ⚠ | 单侧无该角色步骤（口径/序列差异） |
| schedule | role=plug | 0步 | 1步 | ⚠ | 单侧无该角色步骤（口径/序列差异） |
| schedule | 日程总时长min(推算vs现场末步) | 225.4 | 252 | ⚠ | 现场末步含试压/停泵等非注入步，差异仅提示 |
| geometry | 井径剖面点数(CSV实测 vs loader模型域) | 65 | 64 | ⚠ | loader仅保留模型段内点并可能做等效换算 |
| cbl | 地层目标[侏罗系喀拉扎组(J3k2)] 6826.9-7559 | 未覆盖 | target_intervals.csv | ⚠ |  |
| cbl | 地层目标[连木沁组(K1l)] 5292.5-5960 | 未覆盖 | target_intervals.csv | ⚠ |  |
| cbl | 地层目标[领浆封固段] 5068-5796 | 未覆盖 | target_intervals.csv | ⚠ |  |
| cbl | 地层目标[尾浆封固段] 6739-7554 | 未覆盖 | target_intervals.csv | ⚠ |  |
| 内部 | 环空体积重算m3(hole_profile×OD积分) | 33.52 | — | ⚠ | 重算口径：hole_diameter_profile×(双径)liner OD；技套重叠段按外推井径计 |
| 内部 | 水泥浆总体积/环空体积(η_E库存比上限) | 1.134 | — | ⚠ | >1 表示水泥浆量超过环空容量（返高以上/漏失/计量口径） |
| 内部 | 管内容积重算m3 vs 替浆总体积m3 | 20.7 vs 74.0 | — | ⚠ | 替浆总量应≈碰压时管内容积；差异大提示 liner_id/pipe_id 口径为等效或代理 |
| sanity | 密度范围 领浆 | 2100 | 水泥浆1.85–1.95g/cc | ⚠ | 超出常规水泥浆范围：加重/低密配方需与配方记录核对 |

## ht1_001（呼探1-001井（HT1-001））

- 模型域 5469.711–7746.0 m；鞋 7746.0 m；悬挂器 5469.71 m；liner OD/ID 139.7/107.94 mm
- 日程总注入 225.1 m³，11 步
- 环空体积重算 72.8 m³；浆体 71.4 m³；库存比 0.98
- 管内容积重算 20.8 m³；顶替流体合计 91.7 m³

### ht1_001 需关注条目（✗/⚠，共 55 条）

| 类别 | 条目 | loader | CSV | 状态 | 备注 |
|---|---|---|---|---|---|
| fluid | displacement/替浆液[第1对1/2] 密度kg/m3 | 1920 | 1980 | ✗ | CSV g/cm3×1000；容差10kg/m3 |
| fluid | displacement/替钻井液[第2对2/2] 密度kg/m3 | 1920 | 1020 | ✗ | CSV g/cm3×1000；容差10kg/m3 |
| fluid | displacement/中置液[第3对2/2] 密度kg/m3 | 1980 | 1020 | ✗ | CSV g/cm3×1000；容差10kg/m3 |
| fluid | displacement/井浆[第5对2/2] 密度kg/m3 | 1920 | 1020 | ✗ | CSV g/cm3×1000；容差10kg/m3 |
| fluid | role=intermediate | 1条 | 0条 | ⚠ | loader或CSV侧该角色无对应条目（角色映射或口径差异） |
| fluid | role=intermediate_cement | 0条 | 1条 | ⚠ | loader或CSV侧该角色无对应条目（角色映射或口径差异） |
| fluid | role=mud[第1对1/2] 流变模型 | bingham | power_law | ⚠ |  |
| fluid | role=other | 1条 | 0条 | ⚠ | loader或CSV侧该角色无对应条目（角色映射或口径差异） |
| fluid | 隔离液 n(六速反算参考) | 0.545 | 0.511 | ⚠ | 128->93℃ 六速θ600/θ300=114.0/80.0；反算仅供温度口径参考 |
| fluid | 隔离液 K(六速反算参考) | 1.338 | 1.689 | ⚠ | 反算 K=0.511·θ300/(1.7023·300)^n |
| fluid | 领浆 n(六速反算参考) | 0.811 | 0.163 | ⚠ | 128->93℃ 六速θ600/θ300=300.0/268.0；反算仅供温度口径参考 |
| fluid | 领浆 K(六速反算参考) | 0.876 | 49.643 | ⚠ | 反算 K=0.511·θ300/(1.7023·300)^n |
| fluid | 中间浆 n(六速反算参考) | 0.871 | 0.441 | ⚠ | 128->93℃ 六速θ600/θ300=300.0/221.0；反算仅供温度口径参考 |
| fluid | 中间浆 K(六速反算参考) | 0.504 | 7.223 | ⚠ | 反算 K=0.511·θ300/(1.7023·300)^n |
| fluid | 尾浆 n(六速反算参考) | 0.886 | 0.461 | ⚠ | 128->93℃ 六速θ600/θ300=300.0/218.0；反算仅供温度口径参考 |
| fluid | 尾浆 K(六速反算参考) | 0.453 | 6.301 | ⚠ | 反算 K=0.511·θ300/(1.7023·300)^n |
| schedule | role=displacement[1] 体积m3(替钻井液(快)) | 25 | 2 | ✗ | 容差0.2m³/1% |
| schedule | role=displacement[1] 排量m3/min | 1.4 | 0.6 | ✗ | 容差2% |
| schedule | role=displacement[1] 开始min(推算vs现场) | 123.3 | 165 | ⚠ | loader时间为体积/排量推算；容差3min |
| schedule | role=displacement[1] 结束min(推算vs现场) | 141.2 | 168 | ⚠ | 容差3min |
| schedule | role=displacement[2] 体积m3(替保护液/中置液) | 10 | 25 | ✗ | 容差0.2m³/1% |
| schedule | role=displacement[2] 排量m3/min | 1 | 1.4 | ✗ | 容差2% |
| schedule | role=displacement[2] 开始min(推算vs现场) | 141.2 | 168 | ⚠ | loader时间为体积/排量推算；容差3min |
| schedule | role=displacement[2] 结束min(推算vs现场) | 151.2 | 187 | ⚠ | 容差3min |
| schedule | role=displacement[3] 体积m3(替基液) | 3 | 10 | ✗ | 容差0.2m³/1% |
| schedule | role=displacement[3] 开始min(推算vs现场) | 151.2 | 187 | ⚠ | loader时间为体积/排量推算；容差3min |
| schedule | role=displacement[3] 结束min(推算vs现场) | 154.2 | 198 | ⚠ | 容差3min |
| schedule | role=displacement[4] 体积m3(井浆快替) | 35 | 3 | ✗ | 容差0.2m³/1% |
| schedule | role=displacement[4] 开始min(推算vs现场) | 154.2 | 198 | ⚠ | loader时间为体积/排量推算；容差3min |
| schedule | role=displacement[4] 结束min(推算vs现场) | 189.2 | 200 | ⚠ | 容差3min |
| schedule | role=displacement[5] 体积m3(井浆慢替) | 18.7 | 53.7 | ✗ | 容差0.2m³/1% |
| schedule | role=displacement[5] 排量m3/min | 0.6 | 0.9 | ✗ | 容差2% |
| schedule | role=displacement[5] 开始min(推算vs现场) | 189.2 | 200 | ⚠ | loader时间为体积/排量推算；容差3min |
| schedule | role=displacement[5] 结束min(推算vs现场) | 220.3 | 259 | ⚠ | 容差3min |
| schedule | role=intermediate | 1步 | 0步 | ⚠ | 单侧无该角色步骤（口径/序列差异） |
| schedule | role=intermediate_cement | 0步 | 1步 | ⚠ | 单侧无该角色步骤（口径/序列差异） |
| schedule | role=lead[1] 开始min(推算vs现场) | 48.6 | 70 | ⚠ | loader时间为体积/排量推算；容差3min |
| schedule | role=lead[1] 结束min(推算vs现场) | 69.2 | 93 | ⚠ | 容差3min |
| schedule | role=other | 1步 | 0步 | ⚠ | 单侧无该角色步骤（口径/序列差异） |
| schedule | role=spacer[1] 开始min(推算vs现场) | 28.6 | 50 | ⚠ | loader时间为体积/排量推算；容差3min |
| schedule | role=spacer[1] 结束min(推算vs现场) | 48.6 | 70 | ⚠ | 容差3min |
| schedule | role=tail[1] 开始min(推算vs现场) | 97.9 | 125 | ⚠ | loader时间为体积/排量推算；容差3min |
| schedule | role=tail[1] 结束min(推算vs现场) | 120 | 156 | ⚠ | 容差3min |
| schedule | role=wash[1] 开始min(推算vs现场) | 0 | 10 | ⚠ | loader时间为体积/排量推算；容差3min |
| schedule | role=wash[1] 结束min(推算vs现场) | 28.6 | 45 | ⚠ | 容差3min |
| schedule | 日程总时长min(推算vs现场末步) | 220.3 | 259 | ⚠ | 现场末步含试压/停泵等非注入步，差异仅提示 |
| geometry | 井径剖面点数(CSV实测 vs loader模型域) | 70 | 69 | ⚠ | loader仅保留模型段内点并可能做等效换算 |
| cbl | 地层目标[裸眼模拟段] 5670-7746 | 未覆盖 | target_intervals.csv | ⚠ |  |
| cbl | 地层目标[重叠段] 5460.16-5670 | 未覆盖 | target_intervals.csv | ⚠ |  |
| cbl | 地层目标[领浆封固段] 5185.7-5900 | 未覆盖 | target_intervals.csv | ⚠ |  |
| cbl | 地层目标[中间浆封固段] 5900-7000 | 未覆盖 | target_intervals.csv | ⚠ |  |
| 内部 | 环空体积重算m3(hole_profile×OD积分) | 72.76 | — | ⚠ | 重算口径：hole_diameter_profile×(双径)liner OD；技套重叠段按外推井径计 |
| 内部 | 水泥浆总体积/环空体积(η_E库存比上限) | 0.981 | — | ⚠ | >1 表示水泥浆量超过环空容量（返高以上/漏失/计量口径） |
| 内部 | 管内容积重算m3 vs 替浆总体积m3 | 20.8 vs 91.7 | — | ⚠ | 替浆总量应≈碰压时管内容积；差异大提示 liner_id/pipe_id 口径为等效或代理 |
| sanity | 密度范围 领浆 | 2050 | 水泥浆1.85–1.95g/cc | ⚠ | 超出常规水泥浆范围：加重/低密配方需与配方记录核对 |

## ht1_003（呼1-003井（HT1-003））

- 模型域 5307.539–7618.0 m；鞋 7618.0 m；悬挂器 5307.54 m；liner OD/ID 139.7/107.94 mm
- cbl_pass_rate = 0.787
- 日程总注入 214.9 m³，7 步
- 环空体积重算 70.0 m³；浆体 67.0 m³；库存比 0.96
- 管内容积重算 31.8 m³；顶替流体合计 91.9 m³

### ht1_003 需关注条目（✗/⚠，共 36 条）

| 类别 | 条目 | loader | CSV | 状态 | 备注 |
|---|---|---|---|---|---|
| fluid | displacement/替钻井液[第3对3/3] YP Pa | 10 | 62 | ✗ | 容差2% |
| fluid | displacement/井浆[第4对3/3] YP Pa | 9.3 | 62 | ✗ | 容差2% |
| fluid | role=mud[第1对1/3] n | 0.631 | 0.82 | ✗ | 容差2% |
| fluid | role=mud[第1对1/3] K Pa·s^n | 0.751 | 0.21 | ✗ | 容差3% |
| fluid | role=other | 1条 | 0条 | ⚠ | loader或CSV侧该角色无对应条目（角色映射或口径差异） |
| fluid | role=plug | 0条 | 1条 | ⚠ | loader或CSV侧该角色无对应条目（角色映射或口径差异） |
| fluid | 隔离液1 n(六速反算参考) | 0.668 | 0.433 | ⚠ | 93℃ 六速θ600/θ300=239.0/177.0；反算仅供温度口径参考 |
| fluid | 隔离液1 K(六速反算参考) | 1.245 | 6.068 | ⚠ | 反算 K=0.511·θ300/(1.7023·300)^n |
| fluid | 隔离液2 n(六速反算参考) | 0.668 | 0.433 | ⚠ | 93℃ 六速θ600/θ300=239.0/177.0；反算仅供温度口径参考 |
| fluid | 隔离液2 K(六速反算参考) | 1.245 | 6.068 | ⚠ | 反算 K=0.511·θ300/(1.7023·300)^n |
| fluid | 领浆 n(六速反算参考) | 0.597 | 0.775 | ⚠ | 93℃ 六速θ600/θ300=207.0/121.0；反算仅供温度口径参考 |
| fluid | 领浆 K(六速反算参考) | 1.622 | 0.494 | ⚠ | 反算 K=0.511·θ300/(1.7023·300)^n |
| fluid | 尾浆 n(六速反算参考) | 0.585 | 0.732 | ⚠ | 93℃ 六速θ600/θ300=196.0/118.0；反算仅供温度口径参考 |
| fluid | 尾浆 K(六速反算参考) | 1.673 | 0.628 | ⚠ | 反算 K=0.511·θ300/(1.7023·300)^n |
| schedule | role=displacement[1] 开始min(推算vs现场) | 106.7 | 22.83 | ⚠ | loader时间为体积/排量推算；容差3min |
| schedule | role=displacement[1] 结束min(推算vs现场) | 172.3 | 24.17 | ⚠ | 容差3min |
| schedule | role=lead[1] 开始min(推算vs现场) | 49.3 | 21.58 | ⚠ | loader时间为体积/排量推算；容差3min |
| schedule | role=lead[1] 结束min(推算vs现场) | 81.8 | 22.27 | ⚠ | 容差3min |
| schedule | role=other | 1步 | 0步 | ⚠ | 单侧无该角色步骤（口径/序列差异） |
| schedule | role=spacer[1] 结束min(推算vs现场) | 39.3 | 21.25 | ⚠ | 容差3min |
| schedule | role=spacer[2] 开始min(推算vs现场) | 39.3 | 21.25 | ⚠ | loader时间为体积/排量推算；容差3min |
| schedule | role=spacer[2] 结束min(推算vs现场) | 49.3 | 21.58 | ⚠ | 容差3min |
| schedule | role=tail[1] 开始min(推算vs现场) | 81.8 | 22.27 | ⚠ | loader时间为体积/排量推算；容差3min |
| schedule | role=tail[1] 结束min(推算vs现场) | 105.2 | 22.83 | ⚠ | 容差3min |
| schedule | role=wash[1] 开始min(推算vs现场) | 0 | 18.6 | ⚠ | loader时间为体积/排量推算；容差3min |
| schedule | 日程总时长min(推算vs现场末步) | 172.3 | 27.17 | ⚠ | 现场末步含试压/停泵等非注入步，差异仅提示 |
| geometry | well_geometry.csv 同item多值行 | — | borehole_nominal_diameter,drill_bit_size,float_collar_md,liner_inner_diameter,liner_outer_diameter | ⚠ | CSV长表混入回接段同名item，比对取首现（尾管段）值 |
| geometry | 井径剖面点数(CSV实测 vs loader模型域) | 69 | 68 | ⚠ | loader仅保留模型段内点并可能做等效换算 |
| cbl | CBL测量段 top/bottom | 窗['5568-7618'] | 5307.54–7514.21 | ⚠ | loader评价窗常为可评价子段（扣除双层套管/悬空段），见notes |
| cbl | 地层目标[水层-K1q] 6860-6880 | 未覆盖 | target_intervals.csv | ⚠ |  |
| cbl | 地层目标[干层-K1q] 7110-7114 | 未覆盖 | target_intervals.csv | ⚠ |  |
| cbl | 地层目标[主力油层段] 6860-7522 | 未覆盖 | target_intervals.csv | ⚠ |  |
| 内部 | 环空体积重算m3(hole_profile×OD积分) | 69.99 | — | ⚠ | 重算口径：hole_diameter_profile×(双径)liner OD；技套重叠段按外推井径计 |
| 内部 | 水泥浆总体积/环空体积(η_E库存比上限) | 0.957 | — | ⚠ | >1 表示水泥浆量超过环空容量（返高以上/漏失/计量口径） |
| 内部 | 管内容积重算m3 vs 替浆总体积m3 | 31.8 vs 91.9 | — | ⚠ | 替浆总量应≈碰压时管内容积；差异大提示 liner_id/pipe_id 口径为等效或代理 |
| sanity | 密度范围 领浆 | 2050 | 水泥浆1.85–1.95g/cc | ⚠ | 超出常规水泥浆范围：加重/低密配方需与配方记录核对 |

## ht1_003[variant]（呼1-003井（HT1-003））

- 模型域 5307.539–7618.0 m；鞋 7618.0 m；悬挂器 5307.54 m；liner OD/ID 139.7/107.94 mm
- cbl_pass_rate = 0.787
- 日程总注入 212.9 m³，12 步
- 环空体积重算 70.0 m³；浆体 67.0 m³；库存比 0.96
- 管内容积重算 31.8 m³；顶替流体合计 89.9 m³

### ht1_003[variant] 需关注条目（✗/⚠，共 4 条）

| 类别 | 条目 | loader | CSV | 状态 | 备注 |
|---|---|---|---|---|---|
| 内部 | 环空体积重算m3(hole_profile×OD积分) | 69.99 | — | ⚠ | 重算口径：hole_diameter_profile×(双径)liner OD；技套重叠段按外推井径计 |
| 内部 | 水泥浆总体积/环空体积(η_E库存比上限) | 0.957 | — | ⚠ | >1 表示水泥浆量超过环空容量（返高以上/漏失/计量口径） |
| 内部 | 管内容积重算m3 vs 替浆总体积m3 | 31.8 vs 89.9 | — | ⚠ | 替浆总量应≈碰压时管内容积；差异大提示 liner_id/pipe_id 口径为等效或代理 |
| sanity | 密度范围 领浆 | 2050 | 水泥浆1.85–1.95g/cc | ⚠ | 超出常规水泥浆范围：加重/低密配方需与配方记录核对 |

## ht1_004（呼1-004井（HT1-004））

- 模型域 5243.207–7660.0 m；鞋 7660.0 m；悬挂器 5243.21 m；liner OD/ID 139.7/107.94 mm；双径上段 168.3/138.9 mm 底 7376.656 m
- cbl_pass_rate = 0.003
- 日程总注入 226.4 m³，14 步
- 环空体积重算 59.8 m³；浆体 76.0 m³；库存比 1.27
- 管内容积重算 34.9 m³；顶替流体合计 97.4 m³

### ht1_004 需关注条目（✗/⚠，共 60 条）

| 类别 | 条目 | loader | CSV | 状态 | 备注 |
|---|---|---|---|---|---|
| fluid | displacement/替钻井液[第2对2/3] 密度kg/m3 | 1900 | 1020 | ✗ | CSV g/cm3×1000；容差10kg/m3 |
| fluid | displacement/基液[第3对3/3] 密度kg/m3 | 1020 | 1910 | ✗ | CSV g/cm3×1000；容差10kg/m3 |
| fluid | displacement/基液[第3对3/3] YP Pa | 9 | 50 | ✗ | 容差2% |
| fluid | displacement/井浆[第4对3/3] YP Pa | 9.5 | 50 | ✗ | 容差2% |
| fluid | role=flusher | 1条 | 0条 | ⚠ | loader或CSV侧该角色无对应条目（角色映射或口径差异） |
| fluid | role=lead[第1对1/2] 流变模型 | bingham | power_law | ⚠ |  |
| fluid | role=other | 1条 | 0条 | ⚠ | loader或CSV侧该角色无对应条目（角色映射或口径差异） |
| fluid | role=plug | 0条 | 1条 | ⚠ | loader或CSV侧该角色无对应条目（角色映射或口径差异） |
| fluid | spacer/隔离液1 流变模型 | bingham | power_law | ⚠ |  |
| fluid | spacer/隔离液2 流变模型 | bingham | power_law | ⚠ |  |
| fluid | role=tail[第1对1/2] 流变模型 | bingham | power_law | ⚠ |  |
| fluid | role=wash PV Pa·s | 0.058 | 0.03 | ✗ | CSV mPa·s/1000；容差2% |
| fluid | role=wash YP Pa | 9.8 | 7 | ✗ | 容差2% |
| schedule | role=displacement[1] 体积m3(替钻井液(快)) | 29 | 97.1 | ✗ | 容差0.2m³/1% |
| schedule | role=displacement[1] 开始min(推算vs现场) | 103.6 | 170 | ⚠ | loader时间为体积/排量推算；容差3min |
| schedule | role=displacement[1] 结束min(推算vs现场) | 122.9 | 260 | ⚠ | 容差3min |
| schedule | role=displacement[2] 体积m3(替保护液) | 14 | 2 | ✗ | 容差0.2m³/1% |
| schedule | role=displacement[2] 排量m3/min | 1.4 | 1.2 | ✗ | 容差2% |
| schedule | role=displacement[3] 体积m3(替基液) | 1 | 29 | ✗ | 容差0.2m³/1% |
| schedule | role=displacement[3] 排量m3/min | 1.4 | 1.5 | ✗ | 容差2% |
| schedule | role=displacement[4] 排量m3/min | 1.15 | 1.5 | ✗ | 容差2% |
| schedule | role=displacement[5] 体积m3(井浆替入2) | 10 | 30 | ✗ | 容差0.2m³/1% |
| schedule | role=displacement[6] 体积m3(井浆替入3) | 10 | 22 | ✗ | 容差0.2m³/1% |
| schedule | role=displacement[7] 体积m3(井浆替入4) | 10 | 200 | ✗ | 容差0.2m³/1% |
| schedule | role=displacement[7] 排量m3/min | 0.85 | 2.4 | ✗ | 容差2% |
| schedule | role=displacement[7] 开始min(推算vs现场) | 165.9 | 390 | ⚠ | loader时间为体积/排量推算；容差3min |
| schedule | role=displacement[7] 结束min(推算vs现场) | 177.6 | 470 | ⚠ | 容差3min |
| schedule | role=displacement[8] 体积m3(井浆替入5) | 9.4 | 2 | ✗ | 容差0.2m³/1% |
| schedule | role=displacement[8] 排量m3/min | 0.75 | 1.2 | ✗ | 容差2% |
| schedule | role=lead[1] 开始min(推算vs现场) | 39.5 | 73 | ⚠ | loader时间为体积/排量推算；容差3min |
| schedule | role=lead[1] 结束min(推算vs现场) | 79.5 | 140 | ⚠ | 容差3min |
| schedule | role=other | 1步 | 0步 | ⚠ | 单侧无该角色步骤（口径/序列差异） |
| schedule | role=spacer[1] 排量m3/min | 1.2 | 1 | ✗ | 容差2% |
| schedule | role=spacer[1] 开始min(推算vs现场) | 17.9 | 25 | ⚠ | loader时间为体积/排量推算；容差3min |
| schedule | role=spacer[1] 结束min(推算vs现场) | 31.2 | 60 | ⚠ | 容差3min |
| schedule | role=spacer[2] 排量m3/min | 1.2 | 1 | ✗ | 容差2% |
| schedule | role=spacer[2] 开始min(推算vs现场) | 31.2 | 60 | ⚠ | loader时间为体积/排量推算；容差3min |
| schedule | role=spacer[2] 结束min(推算vs现场) | 39.5 | 73 | ⚠ | 容差3min |
| schedule | role=tail[1] 排量m3/min | 1.25 | 1 | ✗ | 容差2% |
| schedule | role=tail[1] 开始min(推算vs现场) | 79.5 | 140 | ⚠ | loader时间为体积/排量推算；容差3min |
| schedule | role=tail[1] 结束min(推算vs现场) | 101.9 | 170 | ⚠ | 容差3min |
| schedule | role=wash[1] 排量m3/min | 1.4 | 1.2 | ✗ | 容差2% |
| schedule | role=wash[1] 开始min(推算vs现场) | 0 | 5 | ⚠ | loader时间为体积/排量推算；容差3min |
| schedule | role=wash[1] 结束min(推算vs现场) | 17.9 | 25 | ⚠ | 容差3min |
| schedule | 日程总时长min(推算vs现场末步) | 190.2 | 470 | ⚠ | 现场末步含试压/停泵等非注入步，差异仅提示 |
| geometry | well_geometry.csv 同item多值行 | — | borehole_nominal_diameter,drill_bit_size,float_collar_md,liner_inner_diameter,liner_outer_diameter | ⚠ | CSV长表混入回接段同名item，比对取首现（尾管段）值 |
| cbl | CBL测量段 top/bottom | 窗['5245-7581'] | 5245–7581 | ⚠ | loader评价窗常为可评价子段（扣除双层套管/悬空段），见notes |
| cbl | 地层目标[干层-K1q] 6908-6912 | 未覆盖 | target_intervals.csv | ⚠ |  |
| cbl | 地层目标[干层-K1q] 6938-6943 | 未覆盖 | target_intervals.csv | ⚠ |  |
| cbl | 地层目标[干层-K1q] 6971-6975 | 未覆盖 | target_intervals.csv | ⚠ |  |
| cbl | 地层目标[干层-K1q] 7049-7055 | 未覆盖 | target_intervals.csv | ⚠ |  |
| cbl | 地层目标[差气层-J3k2] 7601-7603 | 未覆盖 | target_intervals.csv | ⚠ |  |
| cbl | 地层目标[差气层-J3k2] 7610-7613 | 未覆盖 | target_intervals.csv | ⚠ |  |
| cbl | 地层目标[差气层-J3k2] 7617-7621 | 未覆盖 | target_intervals.csv | ⚠ |  |
| cbl | 地层目标[主力气层段] 7495-7550 | 未覆盖 | target_intervals.csv | ⚠ |  |
| cbl | 地层目标[目的井段] 7495-7550 | 未覆盖 | target_intervals.csv | ⚠ |  |
| cbl | 地层目标[油气水层段] 7482-7560 | 未覆盖 | target_intervals.csv | ⚠ |  |
| 内部 | 环空体积重算m3(hole_profile×OD积分) | 59.8 | — | ⚠ | 重算口径：hole_diameter_profile×(双径)liner OD；技套重叠段按外推井径计 |
| 内部 | 水泥浆总体积/环空体积(η_E库存比上限) | 1.271 | — | ⚠ | >1 表示水泥浆量超过环空容量（返高以上/漏失/计量口径） |
| 内部 | 管内容积重算m3 vs 替浆总体积m3 | 34.9 vs 97.4 | — | ⚠ | 替浆总量应≈碰压时管内容积；差异大提示 liner_id/pipe_id 口径为等效或代理 |

## ht1_004[variant]（呼1-004井（HT1-004））

- 模型域 5243.207–7660.0 m；鞋 7660.0 m；悬挂器 5243.21 m；liner OD/ID 139.7/107.94 mm；双径上段 168.3/138.9 mm 底 7376.656 m
- cbl_pass_rate = 0.003
- 日程总注入 224.0 m³，10 步
- 环空体积重算 59.8 m³；浆体 76.0 m³；库存比 1.27
- 管内容积重算 34.9 m³；顶替流体合计 95.0 m³

### ht1_004[variant] 需关注条目（✗/⚠，共 3 条）

| 类别 | 条目 | loader | CSV | 状态 | 备注 |
|---|---|---|---|---|---|
| 内部 | 环空体积重算m3(hole_profile×OD积分) | 59.8 | — | ⚠ | 重算口径：hole_diameter_profile×(双径)liner OD；技套重叠段按外推井径计 |
| 内部 | 水泥浆总体积/环空体积(η_E库存比上限) | 1.271 | — | ⚠ | >1 表示水泥浆量超过环空容量（返高以上/漏失/计量口径） |
| 内部 | 管内容积重算m3 vs 替浆总体积m3 | 34.9 vs 95.0 | — | ⚠ | 替浆总量应≈碰压时管内容积；差异大提示 liner_id/pipe_id 口径为等效或代理 |
