"""对比新旧套管段模型的顶替效率差异。

本脚本比较两种配置的模拟结果：
- 旧模型：无重力修正、无轴向弥散（enable_gravity=False, enable_axial_dispersion=False）
- 新模型：启用重力修正和轴向弥散（enable_gravity=True, enable_axial_dispersion=True）
"""

from cemdisp.data.loaders import load_hu101_tailpipe
from cemdisp.transport1d import CasingFlowSolver
from cemdisp.models2d import AnnulusD2DGASolver
from cemdisp.models2d.boundary_bridge import build_coupled_annulus_inlet_provider


def run_model(enable_gravity, enable_dispersion, label):
    """运行模型并返回结果。"""
    well_spec, fluids, schedule, _ = load_hu101_tailpipe()

    # 套管段模拟
    casing_solver = CasingFlowSolver(
        enable_gravity=enable_gravity,
        enable_axial_dispersion=enable_dispersion,
    )
    casing_result = casing_solver.run(well_spec, fluids, schedule)

    # 获取水泥到达时间
    annulus_stop_time = casing_result.cement_end_time_s

    # 构建耦合边界
    coupled_provider = build_coupled_annulus_inlet_provider(
        casing_result, casing_solver, fluids
    )

    # 环空2D模拟
    solver = AnnulusD2DGASolver(total_t=annulus_stop_time)
    result = solver.run(well_spec, fluids, coupled_provider)

    # 获取最终效率
    final_efficiency = result.metrics.iloc[-1]['effective_efficiency']

    print(f"\n{label}:")
    print(f"  水泥到达鞋口时间: {annulus_stop_time:.1f} s ({annulus_stop_time/60:.1f} min)")
    print(f"  最终有效顶替效率: {final_efficiency:.4f}")

    return casing_result, result


# 运行旧模型（无重力、无弥散）
old_casing, old_result = run_model(
    enable_gravity=False,
    enable_dispersion=False,
    label="旧模型（无重力、无弥散）"
)

# 运行新模型（启用重力和弥散）
new_casing, new_result = run_model(
    enable_gravity=True,
    enable_dispersion=True,
    label="新模型（启用重力和弥散）"
)

# 计算差异
old_eff = old_result.metrics.iloc[-1]['effective_efficiency']
new_eff = new_result.metrics.iloc[-1]['effective_efficiency']
diff = new_eff - old_eff
relative_diff = diff / old_eff * 100 if old_eff > 0 else float('inf')

print("\n" + "="*50)
print("对比结果:")
print(f"  旧模型效率: {old_eff:.4f}")
print(f"  新模型效率: {new_eff:.4f}")
print(f"  绝对差异: {diff:+.4f}")
print(f"  相对差异: {relative_diff:+.2f}%")

# 对比水泥到达时间
old_arrival = old_casing.cement_end_time_s
new_arrival = new_casing.cement_end_time_s
time_diff = new_arrival - old_arrival

print(f"\n水泥到达时间:")
print(f"  旧模型: {old_arrival:.1f} s")
print(f"  新模型: {new_arrival:.1f} s")
print(f"  差异: {time_diff:+.1f} s ({time_diff/60:+.1f} min)")

print("\n物理解释:")
if time_diff > 0:
    print("  - 水泥到达时间延迟，可能由重力沉降导致密度大的水泥浆加速前缘推进后又被弥散效应平滑")
elif time_diff < 0:
    print("  - 水泥到达时间提前，重力沉降效应使高密度水泥浆在管内加速下沉")
else:
    print("  - 水泥到达时间相同，重力和弥散效应相互抵消")

if diff > 0:
    print("  - 新模型效率更高，弥散效应减少了尖锐前缘的窜槽风险")
elif diff < 0:
    print("  - 新模型效率降低，重力沉降导致水泥浆分布不均")
else:
    print("  - 新旧模型效率相同")
