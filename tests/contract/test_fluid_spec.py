"""FluidSpec / FluidRole 基础测试

覆盖 fluid_spec 模块的核心契约：FluidRole 枚举、FluidSpec 构造、
FLUSHER 角色标识等。
"""

import unittest
from cemdisp.data.fluid_spec import FluidRole, FluidSpec


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


if __name__ == "__main__":
    unittest.main()
