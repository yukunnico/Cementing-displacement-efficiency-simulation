# -*- coding: utf-8 -*-
"""分析 20312.xlsm 固井传感器时间序列：识别各浆体阶段的累计量/密度/排量锚点。"""
import pandas as pd

f = r"D:\users\desktop\research\控压固井项目\0708\2\203\2031\20312.xlsm"
df = pd.read_excel(f, sheet_name="数据表", header=0)
df.columns = [str(c).strip() for c in df.columns]
print("列名:", list(df.columns)[:12])
tcol = df.columns[0]
# 数值列
num = df.select_dtypes("number")
print("时间范围:", df[tcol].min(), "->", df[tcol].max(), " 行数:", len(df))
tot_col = [c for c in num.columns if "TotalAll" in c][0]
rate_col = [c for c in num.columns if "RateAll" in c][0]
den_col = [c for c in num.columns if "MassFlowDensity" in c][0]

df = df.sort_values(tcol).reset_index(drop=True)
df["_t"] = pd.to_datetime(df[tcol])
tot = pd.to_numeric(df[tot_col], errors="coerce")
rate = pd.to_numeric(df[rate_col], errors="coerce")
den = pd.to_numeric(df[den_col], errors="coerce")

print("TotalAll 量纲探查: max=%.3f  末值=%.3f" % (tot.max(), tot.iloc[-1]))
print("RateAll: max=%.4f" % rate.max())
print("Density: min=%.4f max=%.4f" % (den.min(), den.max()))

# 找泵启停事件：rate 从 <0.05 到 >0.05 的跳变
act = rate > 0.05
starts = df["_t"][(~act.shift(fill_value=False)) & act]
stops = df["_t"][act & (~act.shift(-1, fill_value=False))]
print("\n泵启时刻:", [str(s) for s in starts.head(30).tolist()])
print("泵停时刻:", [str(s) for s in stops.head(30).tolist()])

# 每 15 分钟快照
snap = df.set_index("_t")[[tot_col, rate_col, den_col]].resample("10min").last().dropna(how="all")
print("\n每10分钟快照(TotalAll/Rate/Density):")
print(snap.to_string(float_format=lambda x: "%.3f" % x))
