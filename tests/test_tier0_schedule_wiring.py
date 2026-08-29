"""D1 接线测试：AnnulusD2DGASolver.run() 透传 schedule → T0-6 停泵衰减诊断。

背景：run() 末尾内嵌 compute_all_tier0_diagnostics 原先不传 schedule，
T0-6 shutdown_decay 恒跳过（notes 记 "未提供 schedule"）。
接线后：
- 不传 schedule：保持既有行为（shutdown_decay=None + notes 优雅降级）——向后兼容；
- 传入含 SHUTDOWN 事件的 toy schedule：shutdown_decay 有值（dict，非 None）。
"""

from __future__ import annotations

from cemdisp.data.fluid_spec import FluidSpec, FluidRole, RheologyModel
from cemdisp.data.pumping_schedule import (
    PumpingSchedule,
    PumpingScheduleStep,
    PumpingStageEvent,
)
from cemdisp.data.well_spec import DepthValuePoint, EvaluationWindow, WellSpec
from cemdisp.models2d.annulus_d2dga import AnnulusD2DGASolver


def _toy_well() -> WellSpec:
    """最小井身结构（垂直井 β=0 → S_buoy=0，判据自动满足；e≈0.17）。"""
    pts = lambda d, v: DepthValuePoint(depth_md_m=d, value=v)  # noqa: E731
    return WellSpec(
        well_name="toy_schedule_wiring",
        top_md_m=1000.0,
        bottom_md_m=1100.0,
        shoe_md_m=1100.0,
        hanger_md_m=1000.0,
        casing_id_mm=200.0,
        liner_od_mm=139.7,
        liner_id_mm=108.0,
        hole_diameter_profile=[pts(1000.0, 215.9), pts(1100.0, 215.9)],
        inclination_profile=[pts(1000.0, 0.0), pts(1100.0, 0.0)],
        standoff_profile=[pts(1000.0, 0.83), pts(1100.0, 0.83)],
        evaluation_windows=[EvaluationWindow(name="w", top_md_m=1000.0,
                                             bottom_md_m=1100.0, window_type="full")],
    )


def _toy_fluids() -> tuple:
    """Bingham 屈服性流体组合（mud tau_y=8.5 Pa，保证 T0-6 判据可满足）。"""
    mud = FluidSpec(name="mud", role=FluidRole.MUD, density_kg_m3=1900.0,
                    rheology_model=RheologyModel.BINGHAM,
                    plastic_viscosity_pa_s=0.053, yield_stress_pa=8.5)
    cement = FluidSpec(name="cement", role=FluidRole.LEAD, density_kg_m3=2160.0,
                       rheology_model=RheologyModel.BINGHAM,
                       plastic_viscosity_pa_s=0.180, yield_stress_pa=14.0)
    return (mud, cement)


def _toy_schedule_with_shutdown() -> PumpingSchedule:
    """含 SHUTDOWN 事件的两段式 toy 泵注程序：注入 → 停泵。"""
    return PumpingSchedule(
        steps=(
            PumpingScheduleStep(
                step_name="注入水泥浆",
                fluid_name="cement",
                volume_m3=0.8,
                rate_m3_min=0.7,
            ),
            PumpingScheduleStep(
                step_name="停泵观察",
                fluid_name="cement",
                volume_m3=0.0,
                rate_m3_min=0.0,
                event_tag=PumpingStageEvent.SHUTDOWN,
            ),
        ),
        notes=("toy schedule：注入后停泵，供 T0-6 诊断触发。",),
    )


class TestTier0ScheduleWiring:
    """run() 的 schedule 参数接线回归测试。"""

    def _run(self, schedule):
        """共用 toy 求解流程（nz=20 快速档），返回 summary。"""
        from cemdisp.models2d.boundary_bridge import AnnulusInletState

        def _inlet(t: float) -> AnnulusInletState:
            return AnnulusInletState(
                time_s=t, flow_rate_m3_s=0.02, stage_name="pump",
                phase_fractions=(("cement", 1.0), ("lead", 1.0)),
            )

        solver = AnnulusD2DGASolver(dt=4.0, nz=20, ny=8, total_t=40.0)
        return solver.run(
            _toy_well(), _toy_fluids(), _inlet, schedule=schedule
        ).summary

    def test_no_schedule_keeps_backward_compat(self):
        """不传 schedule：shutdown_decay=None 且 notes 记跳过（向后兼容）。"""
        summary = self._run(schedule=None)
        tier0 = summary["tier0_diagnostics"]
        assert tier0["shutdown_decay"] is None
        assert any("未提供 schedule" in n for n in tier0["notes"])
        assert "tier0_diagnostics_error" not in summary

    def test_with_shutdown_schedule_shutdown_decay_computed(self):
        """传含 SHUTDOWN 的 schedule：shutdown_decay 有值（dict，非 None）。"""
        summary = self._run(schedule=_toy_schedule_with_shutdown())
        tier0 = summary["tier0_diagnostics"]
        sd = tier0["shutdown_decay"]
        assert isinstance(sd, dict), f"shutdown_decay 应有值，实际：{sd}"
        assert set(sd) == {"condition_satisfied", "freeze_time_s",
                           "tau_y_min", "physical_interpretation"}
        # toy 井垂直（β=0，S_buoy=0）+ Bingham 屈服流体 → 式 3.35 判据必满足
        assert sd["condition_satisfied"] is True
        assert sd["freeze_time_s"] is not None and sd["freeze_time_s"] > 0
        # 诊断层 notes 不再出现 "未提供 schedule" 跳过说明
        assert not any("未提供 schedule" in n for n in tier0["notes"])
        assert "tier0_diagnostics_error" not in summary
