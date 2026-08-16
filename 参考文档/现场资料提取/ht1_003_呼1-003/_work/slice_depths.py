# -*- coding: utf-8 -*-
"""计算每个切片的深度覆盖范围"""
PX2PT = 0.6
SLICE_H_PT = 3600.0
OVERLAP_PT = 120.0
STEP_PT = SLICE_H_PT - OVERLAP_PT
PH = 59069.0
k = 0.07057; b = 3436.84
def y2d(y): return k*y + b
y0 = 0.0; i = 0
while y0 < PH - 1:
    y1 = min(y0 + SLICE_H_PT, PH)
    print(f"slice_{i:02d}: y_pt[{y0:.0f},{y1:.0f}] -> depth {y2d(y0):.0f}-{y2d(y1):.0f}m")
    if y1 >= PH - 1: break
    y0 += STEP_PT; i += 1
