"""
WellSpec 扩展功能测试

测试新增功能：
- 双径向井（dual-diameter）标识字段：upper_section_bottom_md_m、upper_liner_od_mm、upper_liner_id_mm
- 衬管壁厚字段：liner_wall_thickness_mm
- 鞋口延迟体积（shoe-lag）字段：shoe_lag_volume_m3
- WellSpec.is_dual_diameter 属性
- 全有或全无（all-present / all-absent）校验逻辑

验证点：
1. 全部上段字段缺失时 is_dual_diameter == False
2. 全部上段字段完整时 is_dual_diameter == True
3. 上段字段部分提供时触发 ValueError
4. shoe_lag_volume_m3 为正有限值
5. liner_wall_thickness_mm 为正有限值（当提供时）
6. 向后兼容：不含新字段的现有实例仍可正常构造
"""

import unittest
from pathlib import Path
from cemdisp.data.well_spec import WellSpec, DepthValuePoint, EvaluationWindow
from cemdisp.data.fluid_spec import FluidRole, FluidSpec


class TestWellSpecDualDiameterFields(unittest.TestCase):
    """测试双径向井可选字段的添加与 is_dual_diameter 逻辑。"""

    def test_dual_diameter_false_when_no_upper_fields(self):
        """
        上段字段全部缺失时 is_dual_diameter 为 False。
        """
        spec = WellSpec(
            well_name="呼101",
            top_md_m=0.1,
            bottom_md_m=3000.0,
            shoe_md_m=3000.0,
        )
        self.assertFalse(spec.is_dual_diameter)

    def test_dual_diameter_true_when_all_upper_fields_present(self):
        """
        上段字段全部提供时 is_dual_diameter 为 True。
        """
        spec = WellSpec(
            well_name="呼101",
            top_md_m=0.1,
            bottom_md_m=3000.0,
            shoe_md_m=3000.0,
            upper_section_bottom_md_m=1500.0,
            upper_liner_od_mm=244.5,
            upper_liner_id_mm=224.5,
        )
        self.assertTrue(spec.is_dual_diameter)

    def test_dual_diameter_raises_when_upper_fields_partial(self):
        """
        上段字段只提供部分时触发 ValueError（全有或全无约束）。
        """
        cases = [
            # 只提供 upper_section_bottom_md_m
            {
                "upper_section_bottom_md_m": 1500.0,
            },
# 只提供 upper_liner_od_mm
            {
                "upper_liner_od_mm": 244.5,
            },
            # 只提供 upper_liner_id_mm
            {
                "upper_liner_id_mm": 224.5,
            },
            # 只提供 upper_section_bottom_md_m 和 upper_liner_od_mm（缺少 upper_liner_id_mm）
            {
                "upper_section_bottom_md_m": 1500.0,
                "upper_liner_od_mm": 244.5,
            },
            # upper_section_bottom_md_m + upper_liner_id_mm（缺少 upper_liner_od_mm）
            {
                "upper_section_bottom_md_m": 1500.0,
                "upper_liner_id_mm": 224.5,
            },
            # upper_liner_od_mm + upper_liner_id_mm（缺少 upper_section_bottom_md_m）
            {
                "upper_liner_od_mm": 244.5,
                "upper_liner_id_mm": 224.5,
            },
        ]
        for kwargs in cases:
            with self.subTest(kwargs=kwargs):
                with self.assertRaises(ValueError) as ctx:
                    WellSpec(
                        well_name="呼101",
                        top_md_m=0.1,
                        bottom_md_m=3000.0,
                        shoe_md_m=3000.0,
                        **kwargs,
                    )
                self.assertIn("全有或全无", str(ctx.exception))

    def test_upper_section_bottom_md_m_must_be_positive(self):
        """
        upper_section_bottom_md_m 必须为正数。
        """
        with self.assertRaises(ValueError):
            WellSpec(
                well_name="呼101",
                top_md_m=0.1,
                bottom_md_m=3000.0,
                shoe_md_m=3000.0,
                upper_section_bottom_md_m=-100.0,
                upper_liner_od_mm=244.5,
                upper_liner_id_mm=224.5,
            )

    def test_upper_section_bottom_md_m_must_be_within_section(self):
        """
        upper_section_bottom_md_m 必须在井段范围内。
        """
        with self.assertRaises(ValueError):
            WellSpec(
                well_name="呼101",
                top_md_m=0.1,
                bottom_md_m=3000.0,
                shoe_md_m=3000.0,
                upper_section_bottom_md_m=4000.0,  # > bottom_md_m
                upper_liner_od_mm=244.5,
                upper_liner_id_mm=224.5,
            )

    def test_upper_liner_od_mm_must_be_positive(self):
        """
        upper_liner_od_mm 必须为正数。
        """
        with self.assertRaises(ValueError):
            WellSpec(
                well_name="呼101",
                top_md_m=0.1,
                bottom_md_m=3000.0,
                shoe_md_m=3000.0,
                upper_section_bottom_md_m=1500.0,
                upper_liner_od_mm=-1.0,
                upper_liner_id_mm=224.5,
            )

    def test_upper_liner_id_mm_must_be_positive(self):
        """
        upper_liner_id_mm 必须为正数。
        """
        with self.assertRaises(ValueError):
            WellSpec(
                well_name="呼101",
                top_md_m=0.1,
                bottom_md_m=3000.0,
                shoe_md_m=3000.0,
                upper_section_bottom_md_m=1500.0,
                upper_liner_od_mm=244.5,
                upper_liner_id_mm=-1.0,
            )


class TestWellSpecShoeLagVolume(unittest.TestCase):
    """测试 shoe_lag_volume_m3 字段的校验。"""

    def test_shoe_lag_volume_optional(self):
        """
        shoe_lag_volume_m3 未提供时为 None（向后兼容）。
        """
        spec = WellSpec(
            well_name="呼101",
            top_md_m=0.1,
            bottom_md_m=3000.0,
            shoe_md_m=3000.0,
        )
        self.assertIsNone(spec.shoe_lag_volume_m3)

    def test_shoe_lag_volume_must_be_positive(self):
        """
        shoe_lag_volume_m3 提供时必须大于零。
        """
        with self.assertRaises(ValueError):
            WellSpec(
                well_name="呼101",
                top_md_m=0.1,
                bottom_md_m=3000.0,
                shoe_md_m=3000.0,
                shoe_lag_volume_m3=0.0,
            )
        with self.assertRaises(ValueError):
            WellSpec(
                well_name="呼101",
                top_md_m=0.1,
                bottom_md_m=3000.0,
                shoe_md_m=3000.0,
                shoe_lag_volume_m3=-5.0,
            )

    def test_shoe_lag_volume_must_be_finite(self):
        """
        shoe_lag_volume_m3 提供时必须为有限值。
        """
        with self.assertRaises(ValueError):
            WellSpec(
                well_name="呼101",
                top_md_m=0.1,
                bottom_md_m=3000.0,
                shoe_md_m=3000.0,
                shoe_lag_volume_m3=float("inf"),
            )


class TestWellSpecLinerWallThickness(unittest.TestCase):
    """测试 liner_wall_thickness_mm 字段的校验。"""

    def test_liner_wall_thickness_optional(self):
        """
        liner_wall_thickness_mm 未提供时为 None。
        """
        spec = WellSpec(
            well_name="呼101",
            top_md_m=0.1,
            bottom_md_m=3000.0,
            shoe_md_m=3000.0,
        )
        self.assertIsNone(spec.liner_wall_thickness_mm)

    def test_liner_wall_thickness_must_be_positive(self):
        """
        liner_wall_thickness_mm 提供时必须大于零。
        """
        with self.assertRaises(ValueError):
            WellSpec(
                well_name="呼101",
                top_md_m=0.1,
                bottom_md_m=3000.0,
                shoe_md_m=3000.0,
                liner_wall_thickness_mm=0.0,
            )
        with self.assertRaises(ValueError):
            WellSpec(
                well_name="呼101",
                top_md_m=0.1,
                bottom_md_m=3000.0,
                shoe_md_m=3000.0,
                liner_wall_thickness_mm=-2.0,
            )


class TestWellSpecPositionalBackwardCompatibility(unittest.TestCase):
    """验证旧代码使用位置参数构造 WellSpec 时，剖面/窗口/reference_root/notes 仍正确落地。"""

    def test_positional_profile_and_windows_land_correctly(self):
        """
        通过位置参数传入 hole_diameter_profile、evaluation_windows，
        验证它们仍落入正确的现有字段，而非被新字段吸收。
        """
        pt = DepthValuePoint(depth_md_m=500.0, value=215.0)
        win = EvaluationWindow(name="目的层", top_md_m=800.0, bottom_md_m=1200.0)
        spec = WellSpec(
            "呼102",          # pos 0 well_name
            0.1,              # pos 1 top_md_m
            2500.0,           # pos 2 bottom_md_m
            2500.0,           # pos 3 shoe_md_m
            None,             # pos 4 hanger_md_m
            None,             # pos 5 casing_id_mm
            None,             # pos 6 liner_od_mm
            None,             # pos 7 liner_id_mm
            (pt,),            # pos 8 hole_diameter_profile
            (),               # pos 9 liner_od_profile
            (),               # pos 10 pipe_id_profile
            (),               # pos 11 inclination_profile
            (),               # pos 12 standoff_profile
            (win,),           # pos 13 evaluation_windows
        )
        self.assertEqual(len(spec.hole_diameter_profile), 1)
        self.assertEqual(spec.hole_diameter_profile[0].depth_md_m, 500.0)
        self.assertEqual(len(spec.evaluation_windows), 1)
        self.assertEqual(spec.evaluation_windows[0].name, "目的层")

    def test_positional_notes_land_correctly(self):
        """
        通过位置参数传入 notes，验证仍落入正确的现有字段。
        """
        spec = WellSpec(
            "呼102",          # pos 0 well_name
            0.1,              # pos 1 top_md_m
            2500.0,           # pos 2 bottom_md_m
            2500.0,           # pos 3 shoe_md_m
            None, None, None, None,  # pos 4-7 hanger/casing/liner_od/liner_id
            (),               # pos 8 hole_diameter_profile
            (),               # pos 9 liner_od_profile
            (),               # pos 10 pipe_id_profile
            (),               # pos 11 inclination_profile
            (),               # pos 12 standoff_profile
            (),               # pos 13 evaluation_windows
            None,             # pos 14 reference_root
            ("备注1", "备注2"),  # pos 15 notes
        )
        self.assertEqual(spec.notes, ("备注1", "备注2"))

    def test_positional_reference_root_land_correctly(self):
        """
        通过位置参数传入 reference_root（Path），验证仍落入正确的现有字段。
        """
        ref = Path("C:/Users/30525/Documents")
        spec = WellSpec(
            "呼102",          # pos 0 well_name
            0.1,              # pos 1 top_md_m
            2500.0,           # pos 2 bottom_md_m
            2500.0,           # pos 3 shoe_md_m
            None, None, None, None,  # pos 4-7 hanger/casing/liner_od/liner_id
            (),               # pos 8 hole_diameter_profile
            (),               # pos 9 liner_od_profile
            (),               # pos 10 pipe_id_profile
            (),               # pos 11 inclination_profile
            (),               # pos 12 standoff_profile
            (),               # pos 13 evaluation_windows
            ref,              # pos 14 reference_root
            ("note",),        # pos 15 notes
        )
        self.assertEqual(spec.reference_root, ref)

    def test_minimal_spec_still_works(self):
        """
        只含必选字段的实例仍可正常构造。
        """
        spec = WellSpec(
            well_name="呼102",
            top_md_m=0.1,
            bottom_md_m=2500.0,
            shoe_md_m=2500.0,
        )
        self.assertEqual(spec.well_name, "呼102")
        self.assertFalse(spec.is_dual_diameter)
        self.assertIsNone(spec.shoe_lag_volume_m3)
        self.assertIsNone(spec.liner_wall_thickness_mm)

    def test_existing_optional_fields_still_work(self):
        """
        原有可选字段（hanger_md_m、casing_id_mm、liner_od_mm、liner_id_mm）仍可用。
        """
        spec = WellSpec(
            well_name="呼102",
            top_md_m=0.1,
            bottom_md_m=2500.0,
            shoe_md_m=2500.0,
            hanger_md_m=100.0,
            casing_id_mm=250.0,
            liner_od_mm=244.5,
            liner_id_mm=224.5,
        )
        self.assertEqual(spec.hanger_md_m, 100.0)
        self.assertEqual(spec.casing_id_mm, 250.0)
        self.assertEqual(spec.liner_od_mm, 244.5)
        self.assertEqual(spec.liner_id_mm, 224.5)


class TestFluidRoleFlusher(unittest.TestCase):
    """测试 FluidRole.FLUSHER 枚举值与 FluidSpec 构造。"""

    def test_flusher_enum_value(self):
        """FluidRole.FLUSHER.value 应为 'flusher'。"""
        self.assertEqual(FluidRole.FLUSHER.value, "flusher")

    def test_fluidspec_with_flusher_role_constructs(self):
        """FluidSpec(role=FluidRole.FLUSHER) 构造不应报错。"""
        spec = FluidSpec(
            role=FluidRole.FLUSHER,
            name="f",
            density_kg_m3=1000,
            plastic_viscosity_pa_s=0.1,
        )
        self.assertEqual(spec.role, FluidRole.FLUSHER)
        self.assertEqual(spec.name, "f")
        self.assertEqual(spec.density_kg_m3, 1000)

    def test_flusher_role_is_str_enum(self):
        """FluidRole.FLUSHER 是字符串枚举。"""
        self.assertEqual(FluidRole.FLUSHER, "flusher")


class TestFlusherLoaderMapping(unittest.TestCase):
    """验证 hu102 / ht1_004 加载器按 mud-spacer-flusher-cement 序列映射 FLUSHER。"""

    def _step_fluid_roles(self, schedule, fluids):
        """根据 schedule step 的 fluid_name 返回对应 FluidRole 列表。"""
        role_by_name = {fluid.name: fluid.role for fluid in fluids}
        return [role_by_name.get(step.fluid_name) for step in schedule.steps]

    def test_hu102_strict_mode_no_flusher_step(self):
        """呼102严格现场模式（默认）schedule 中不含 FLUSHER step，但流体定义保留。"""
        from cemdisp.data.loaders.hu102_loader import load_hu102_tailpipe
        _, fluids, schedule, _ = load_hu102_tailpipe()
        roles = [fluid.role for fluid in fluids]
        self.assertIn(FluidRole.FLUSHER, roles)
        step_roles = self._step_fluid_roles(schedule, fluids)
        self.assertNotIn(FluidRole.FLUSHER, step_roles)
        # 现有"冲洗液"为套管清洗液，保留 WASH
        wash = next((f for f in fluids if f.name == "冲洗液"), None)
        self.assertIsNotNone(wash)
        self.assertEqual(wash.role, FluidRole.WASH)

    def test_hu102_include_wash_spacer_has_flusher_step(self):
        """呼102 include_wash_spacer=True 时，冲洗液仍为 WASH，且 schedule 含 FLUSHER。"""
        from cemdisp.data.loaders.hu102_loader import load_hu102_tailpipe
        _, fluids, schedule, _ = load_hu102_tailpipe(include_wash_spacer=True)
        roles = [fluid.role for fluid in fluids]
        self.assertIn(FluidRole.FLUSHER, roles)
        self.assertIn(FluidRole.WASH, roles)
        step_roles = self._step_fluid_roles(schedule, fluids)
        self.assertIn(FluidRole.FLUSHER, step_roles)
        wash = next((f for f in fluids if f.name == "冲洗液"), None)
        self.assertIsNotNone(wash)
        self.assertEqual(wash.role, FluidRole.WASH)

    def test_ht1_004_strict_mode_no_flusher_step(self):
        """呼1-004 严格现场模式（默认）schedule 中不含 FLUSHER step，先导浆保留 WASH。"""
        from cemdisp.data.loaders.ht1_004_loader import load_ht1_004_tailpipe
        _, fluids, schedule, _ = load_ht1_004_tailpipe()
        roles = [fluid.role for fluid in fluids]
        self.assertIn(FluidRole.FLUSHER, roles)
        step_roles = self._step_fluid_roles(schedule, fluids)
        self.assertNotIn(FluidRole.FLUSHER, step_roles)
        lead_mud = next((f for f in fluids if f.name == "先导浆"), None)
        self.assertIsNotNone(lead_mud)
        self.assertEqual(lead_mud.role, FluidRole.WASH)

    def test_ht1_004_include_wash_spacer_has_flusher_step(self):
        """呼1-004 include_wash_spacer=True 时，先导浆保留 WASH，且 schedule 含 FLUSHER。"""
        from cemdisp.data.loaders.ht1_004_loader import load_ht1_004_tailpipe
        _, fluids, schedule, _ = load_ht1_004_tailpipe(include_wash_spacer=True)
        roles = [fluid.role for fluid in fluids]
        self.assertIn(FluidRole.FLUSHER, roles)
        step_roles = self._step_fluid_roles(schedule, fluids)
        self.assertIn(FluidRole.FLUSHER, step_roles)
        lead_mud = next((f for f in fluids if f.name == "先导浆"), None)
        self.assertIsNotNone(lead_mud)
        self.assertEqual(lead_mud.role, FluidRole.WASH)


if __name__ == "__main__":
    unittest.main()