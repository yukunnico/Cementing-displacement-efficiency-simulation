"""多井来源口径聚合模型测试。

这些测试先固定 Task 4 的对外行为：在不丢失既有流体来源登记的前提下，
新增每口井的几何、施工程序和同步边界来源口径摘要。
"""

import unittest


class TestMultiwellProvenanceRegistry(unittest.TestCase):
    """测试六井聚合来源口径登记。"""

    def test_registry_contains_all_six_wells_with_four_sections(self):
        """聚合登记覆盖六口井，且每口井都有 fluid/geometry/program/sync 四段。"""
        from cemdisp.data.provenance import WELL_PROVENANCE

        expected_wells = {"呼101", "呼102", "呼103", "呼探1", "呼探1-002", "呼探1-001"}

        self.assertEqual(set(WELL_PROVENANCE), expected_wells)
        for well_name, provenance in WELL_PROVENANCE.items():
            with self.subTest(well_name=well_name):
                self.assertEqual(provenance.well_name, well_name)
                self.assertTrue(provenance.fluid)
                self.assertIn("geometry", provenance.sections)
                self.assertIn("program", provenance.sections)
                self.assertIn("sync", provenance.sections)

    def test_fluid_registry_is_preserved_in_aggregate_model(self):
        """聚合模型中的流体登记不重命名、不裁剪既有覆盖范围。"""
        from cemdisp.data.fluid_provenance import WELL_FLUID_PROVENANCE as legacy_registry
        from cemdisp.data.provenance import WELL_FLUID_PROVENANCE, WELL_PROVENANCE

        self.assertIs(WELL_FLUID_PROVENANCE, legacy_registry)
        self.assertEqual(set(WELL_PROVENANCE), set(legacy_registry))
        for well_name, fluid_registry in legacy_registry.items():
            with self.subTest(well_name=well_name):
                self.assertEqual(WELL_PROVENANCE[well_name].fluid, fluid_registry)


class TestFluidProvenanceCompatibility(unittest.TestCase):
    """测试旧 fluid_provenance API 仍可使用。"""

    def test_legacy_module_reexports_new_provenance_objects(self):
        """兼容层导出的类、登记表和摘要函数与新模块保持同一对象。"""
        from cemdisp.data import provenance
        from cemdisp.data import fluid_provenance

        self.assertIs(fluid_provenance.FluidProvenance, provenance.FluidProvenance)
        self.assertIs(fluid_provenance.WELL_FLUID_PROVENANCE, provenance.WELL_FLUID_PROVENANCE)
        self.assertIs(
            fluid_provenance.build_injected_fluid_provenance_summary,
            provenance.build_injected_fluid_provenance_summary,
        )
        self.assertIs(
            fluid_provenance.format_injected_fluid_provenance_markdown,
            provenance.format_injected_fluid_provenance_markdown,
        )


if __name__ == "__main__":
    _ = unittest.main()
