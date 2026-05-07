
from __future__ import annotations
from cemdisp.data.fluid_spec import FluidRole, FluidSpec, RheologyModel
from cemdisp.data.pumping_schedule import PumpingSchedule, PumpingScheduleStep
from cemdisp.data.well_spec import DepthValuePoint, EvaluationWindow, WellSpec


def load_hu101_tailpipe():
    """呼101尾管段输入。单井硬编码数据集中在加载器，保持项目结构一致。"""
    well = WellSpec(
        well_name="呼101",
        top_md_m=5402.85,
        bottom_md_m=7868.00,
        shoe_md_m=7868.00,
        hanger_md_m=5407.46,
        liner_od_mm=139.70,
        liner_id_mm=108.10,
        casing_lag_volume_m3=52.0,
        upper_liner_od_mm=168.30,
        lower_liner_od_mm=139.70,
        upper_lower_transition_md_m=6796.0,
        hole_diameter_profile=(
            DepthValuePoint(5402.85, 273.10),
            DepthValuePoint(6100.00, 273.10),
            DepthValuePoint(6796.00, 273.10),
            DepthValuePoint(7200.00, 241.30),
            DepthValuePoint(7600.00, 241.30),
            DepthValuePoint(7868.00, 241.30),
        ),
        inclination_profile=(
            DepthValuePoint(5402.85, 0.6), DepthValuePoint(6100.0, 1.1),
            DepthValuePoint(6796.0, 1.4), DepthValuePoint(7400.0, 1.7),
            DepthValuePoint(7868.0, 1.9),
        ),
        standoff_profile=(
            DepthValuePoint(5402.85, 0.45), DepthValuePoint(6100.0, 0.38),
            DepthValuePoint(6796.0, 0.44), DepthValuePoint(7200.0, 0.48),
            DepthValuePoint(7600.0, 0.42), DepthValuePoint(7868.0, 0.46),
        ),
        evaluation_windows=(
            EvaluationWindow("全井段", 5402.85, 7868.00, "full"),
            EvaluationWindow("CBL评价段", 5700.00, 7810.00, "cbl"),
            EvaluationWindow("目的层段", 7492.00, 7735.00, "target"),
        ),
        notes="呼101现场尾管固井验证数据；环空入口由套管段1D前沿追踪给出。",
    )
    fluids = (
        FluidSpec("钻井液", FluidRole.MUD, 1960.0, RheologyModel.BINGHAM, plastic_viscosity_pa_s=0.058, yield_stress_pa=5.0),
        FluidSpec("平衡液", FluidRole.WASH, 1850.0, RheologyModel.BINGHAM, plastic_viscosity_pa_s=0.030, yield_stress_pa=3.0),
        FluidSpec("驱油隔离液", FluidRole.SPACER, 2000.0, RheologyModel.BINGHAM, plastic_viscosity_pa_s=0.030, yield_stress_pa=5.0),
        FluidSpec("领浆", FluidRole.LEAD, 2100.0, RheologyModel.POWER_LAW, power_law_n=0.719, consistency_k=0.815),
        FluidSpec("尾浆", FluidRole.TAIL, 1900.0, RheologyModel.POWER_LAW, power_law_n=0.722, consistency_k=0.684),
        FluidSpec("后置液", FluidRole.DISPLACEMENT, 1850.0, RheologyModel.BINGHAM, plastic_viscosity_pa_s=0.025, yield_stress_pa=3.0),
        FluidSpec("轻泥浆", FluidRole.DISPLACEMENT, 1850.0, RheologyModel.BINGHAM, plastic_viscosity_pa_s=0.035, yield_stress_pa=4.0),
        FluidSpec("中置液", FluidRole.DISPLACEMENT, 1900.0, RheologyModel.BINGHAM, plastic_viscosity_pa_s=0.040, yield_stress_pa=4.5),
        FluidSpec("井浆", FluidRole.DISPLACEMENT, 1960.0, RheologyModel.BINGHAM, plastic_viscosity_pa_s=0.058, yield_stress_pa=5.0),
    )
    schedule = PumpingSchedule((
        PumpingScheduleStep("注平衡液", "平衡液", 25.0, 1.20),
        PumpingScheduleStep("注驱油隔离液", "驱油隔离液", 25.0, 1.20),
        PumpingScheduleStep("注领浆", "领浆", 47.0, 1.20),
        PumpingScheduleStep("注尾浆", "尾浆", 23.0, 1.20),
        PumpingScheduleStep("注后置液", "后置液", 2.0, 0.60),
        PumpingScheduleStep("注轻泥浆", "轻泥浆", 26.0, 1.50),
        PumpingScheduleStep("注中置液", "中置液", 10.0, 1.20),
        PumpingScheduleStep("井浆快替", "井浆", 40.0, 1.00),
        PumpingScheduleStep("井浆慢替", "井浆", 23.4, 0.55),
    ))
    validation = {"cbl_quality_pass_rate": 0.6277, "cbl_interval_m": [5700.0, 7810.0]}
    return well, fluids, schedule, validation
