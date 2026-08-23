"""
窄边屈服死区正则化模型测试。

验证 Papanastasiou 正则化在环空二维求解器中的行为：
- 高剪切区：正则化不改变正常流动
- 低剪切区：黏度显著增大，形成死区
- 零屈服应力流体：行为退化到原模型
"""

import unittest
import numpy as np

from cemdisp.data.fluid_spec import FluidSpec, FluidRole, RheologyModel
from cemdisp.models2d.annulus_d2dga import AnnulusD2DGASolver


class TestYieldRegularization(unittest.TestCase):
    """测试 Papanastasiou 正则化黏度行为。"""

    def _make_solver(self, M: float = 100.0) -> AnnulusD2DGASolver:
        """创建带指定正则化参数的求解器。"""
        return AnnulusD2DGASolver(
            dt=4.0,
            nz=10,
            ny=5,
            total_t=100.0,
            yield_regularization_M=M,
        )

    def _make_bingham_fluid(self, name: str, tau_y: float) -> FluidSpec:
        """创建 Bingham 流体。"""
        return FluidSpec(
            name=name,
            role=FluidRole.MUD,
            density_kg_m3=1200.0,
            rheology_model=RheologyModel.BINGHAM,
            plastic_viscosity_pa_s=0.05,
            yield_stress_pa=tau_y,
        )

    def _make_newtonian_fluid(self, name: str) -> FluidSpec:
        """创建牛顿流体。"""
        return FluidSpec(
            name=name,
            role=FluidRole.MUD,
            density_kg_m3=1000.0,
            rheology_model=RheologyModel.NEWTONIAN,
            plastic_viscosity_pa_s=0.01,
        )

    def test_high_shear_regularity_negligible(self) -> None:
        """高剪切 (γ̇ > 100) 时正则化项替代 tau_y/gamma，μ_reg ≈ μ_app。"""
        solver = self._make_solver(M=100.0)
        fluid = self._make_bingham_fluid("mud", tau_y=10.0)
        gamma = np.array([200.0, 500.0, 1000.0])
        mu_app = solver._apparent_viscosity(fluid, gamma)

        # 模拟 _compute_velocity 中的正则化计算（修正版）
        gamma_safe = np.maximum(gamma, 1.0e-8)
        regularization_factor = (1.0 - np.exp(-solver.yield_regularization_M * gamma)) / gamma_safe
        # 减去纯屈服应力贡献，再加正则化项
        mu_shear = np.maximum(mu_app - fluid.yield_stress_pa / gamma_safe, 1.0e-6)
        mu_reg = mu_shear + fluid.yield_stress_pa * regularization_factor

        # 高剪切下正则化应近似等于原表观黏度（差异 < 5%）
        np.testing.assert_allclose(mu_reg, mu_app, rtol=0.05)

    def test_low_shear_regularity_significant(self) -> None:
        """低剪切 (γ̇ < 0.01) 时正则化项显著增大，μ_reg >> μ_app。"""
        solver = self._make_solver(M=100.0)
        fluid = self._make_bingham_fluid("mud", tau_y=10.0)
        gamma = np.array([0.001, 0.005, 0.01])
        mu_app = solver._apparent_viscosity(fluid, gamma)

        # 模拟修正后的正则化计算
        gamma_safe = np.maximum(gamma, 1.0e-8)
        regularization_factor = (1.0 - np.exp(-solver.yield_regularization_M * gamma)) / gamma_safe
        mu_shear = np.maximum(mu_app - fluid.yield_stress_pa / gamma_safe, 1.0e-6)
        mu_reg = mu_shear + fluid.yield_stress_pa * regularization_factor

        # 低剪切下正则化黏度应显著大于表观黏度（> 5 倍）
        self.assertTrue(np.all(mu_reg > 5.0 * mu_app))

    def test_zero_yield_stress_no_effect(self) -> None:
        """τ_y = 0 的流体（牛顿/幂律）不受正则化影响，μ_reg ≈ μ_app。"""
        solver = self._make_solver(M=100.0)
        fluid = self._make_newtonian_fluid("mud")
        gamma = np.array([0.001, 1.0, 100.0, 1000.0])
        mu_app = solver._apparent_viscosity(fluid, gamma)

        regularization_factor = (1.0 - np.exp(-solver.yield_regularization_M * gamma)) / np.maximum(gamma, 1.0e-8)
        tau_y = solver._fluid_yield_stress(fluid)
        mu_reg = mu_app + tau_y * regularization_factor

        # 牛顿流体 tau_y = 0，正则化不应改变黏度
        np.testing.assert_allclose(mu_reg, mu_app, rtol=1e-12)

    def test_fluid_yield_stress_extraction(self) -> None:
        """_fluid_yield_stress 正确提取各流变模型的屈服应力。"""
        bingham = self._make_bingham_fluid("bingham", tau_y=5.0)
        newtonian = self._make_newtonian_fluid("newtonian")

        self.assertEqual(AnnulusD2DGASolver._fluid_yield_stress(bingham), 5.0)
        self.assertEqual(AnnulusD2DGASolver._fluid_yield_stress(newtonian), 0.0)

    def test_mixed_yield_stress_computation(self) -> None:
        """混合屈服应力按体积分数加权正确计算。"""
        solver = self._make_solver()
        mud = self._make_bingham_fluid("mud", tau_y=8.0)
        spacer = FluidSpec(
            name="spacer",
            role=FluidRole.SPACER,
            density_kg_m3=1100.0,
            rheology_model=RheologyModel.NEWTONIAN,
            plastic_viscosity_pa_s=0.02,
        )

        # 创建简单的二维场
        ny, nz = 5, 10
        lead = np.zeros((ny, nz))
        tail = np.zeros((ny, nz))
        spacer_field = np.full((ny, nz), 0.3)
        w_prev = np.ones((ny, nz)) * 0.5

        # 构造最小几何参数
        geom = {
            "b": np.ones((ny, nz)) * 0.02,
            "effective_b": np.ones((ny, nz)) * 0.02,
        }

        mu, rho, mud_frac, tau_y, m_field, eta1, eta2, _n_mix, _kappa_mix = solver._compute_props(
            lead, tail, spacer_field, np.zeros_like(lead), w_prev, geom,
            mud, None, None, spacer,
        )

        # 混合屈服应力 = 0.7 * 8.0 + 0.3 * 0.0 = 5.6
        expected_tau_y = 0.7 * 8.0 + 0.3 * 0.0
        np.testing.assert_allclose(
            np.mean(tau_y), expected_tau_y, rtol=1e-6,
        )

    def test_yield_regularization_M_parameter(self) -> None:
        """__init__ 正确保存 yield_regularization_M 参数。"""
        solver1 = AnnulusD2DGASolver(yield_regularization_M=50.0)
        solver2 = AnnulusD2DGASolver(yield_regularization_M=200.0)
        self.assertEqual(solver1.yield_regularization_M, 50.0)
        self.assertEqual(solver2.yield_regularization_M, 200.0)

    def test_default_M_is_100(self) -> None:
        """默认正则化参数 M = 100.0。"""
        solver = AnnulusD2DGASolver()
        self.assertEqual(solver.yield_regularization_M, 100.0)


class TestNarrowSideDeadzoneEffect(unittest.TestCase):
    """测试窄边死区对速度场的影响。"""

    def test_narrow_side_velocity_reduced(self) -> None:
        """有屈服应力时，窄边速度应显著低于宽边速度。"""
        solver = AnnulusD2DGASolver(
            dt=4.0,
            nz=20,
            ny=10,
            total_t=100.0,
            yield_regularization_M=100.0,
        )

        # Bingham 泥浆有屈服应力
        mud = FluidSpec(
            name="mud",
            role=FluidRole.MUD,
            density_kg_m3=1200.0,
            rheology_model=RheologyModel.BINGHAM,
            plastic_viscosity_pa_s=0.05,
            yield_stress_pa=15.0,
        )
        lead = FluidSpec(
            name="lead",
            role=FluidRole.LEAD,
            density_kg_m3=1900.0,
            rheology_model=RheologyModel.NEWTONIAN,
            plastic_viscosity_pa_s=0.3,
        )

        ny, nz = 10, 20
        lead_field = np.zeros((ny, nz))
        tail_field = np.zeros((ny, nz))
        spacer_field = np.zeros((ny, nz))
        w_prev = np.full((ny, nz), 0.3)

        # 构造偏心环空几何：窄边间隙小，宽边间隙大
        y = np.linspace(0.0, 1.0, ny)
        phi = y / y[-1]
        e = 0.4
        half_gap_mean = 0.015
        b = np.zeros((ny, nz))
        for j in range(nz):
            b[:, j] = 2.0 * half_gap_mean * (1.0 + e * np.cos(np.pi * phi))

        geom = {
            "y": y,
            "phi": phi,
            "b": b,
            "effective_b": b.copy(),
            "e": np.full(nz, e),
            "s": np.linspace(0.0, 100.0, nz),
        }

        w, v, mu_reg, rho, mud_frac, Re, mu_turbulent, m_field, _tau_y, _eta2, _n_mix, _kappa_mix = solver._compute_velocity(
            lead_field, tail_field, spacer_field, np.zeros_like(lead_field), geom,
            q_m3s=0.01, w_prev=w_prev,
            mud_fluid=mud, lead_fluid=lead, tail_fluid=None, spacer_fluid=None,
        )

        # 宽边 (phi=0) 速度应大于窄边 (phi=-1) 速度
        wide_velocity = np.mean(w[0, :])
        narrow_velocity = np.mean(w[-1, :])
        self.assertGreater(
            wide_velocity, narrow_velocity * 1.5,
            f"宽边速度 ({wide_velocity:.4f}) 应显著大于窄边速度 ({narrow_velocity:.4f})",
        )

    def test_compute_props_returns_nine_values(self) -> None:
        """_compute_props 返回九元组 (mu, rho, mud, tau_y, m_field, eta1, eta2, n_mix, kappa_mix)。"""
        solver = AnnulusD2DGASolver()
        mud = FluidSpec(
            name="mud",
            role=FluidRole.MUD,
            density_kg_m3=1200.0,
            rheology_model=RheologyModel.NEWTONIAN,
            plastic_viscosity_pa_s=0.01,
        )

        ny, nz = 5, 10
        lead = np.zeros((ny, nz))
        tail = np.zeros((ny, nz))
        spacer = np.zeros((ny, nz))
        w_prev = np.ones((ny, nz)) * 0.5
        geom = {"b": np.ones((ny, nz)) * 0.02, "effective_b": np.ones((ny, nz)) * 0.02}

        result = solver._compute_props(
            lead, tail, spacer, np.zeros_like(lead), w_prev, geom,
            mud, None, None, None,
        )
        self.assertEqual(len(result), 9)
        mu, rho, mud_frac, tau_y, m_field, eta1, eta2, n_mix, kappa_mix = result
        self.assertEqual(mu.shape, (ny, nz))
        self.assertEqual(tau_y.shape, (ny, nz))
        self.assertEqual(eta1.shape, (ny, nz))
        self.assertEqual(eta2.shape, (ny, nz))
        self.assertEqual(n_mix.shape, (ny, nz))
        self.assertEqual(kappa_mix.shape, (ny, nz))


if __name__ == "__main__":
    unittest.main()
