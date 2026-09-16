"""水泥屈服应力（τ_y）场外常数汇总 —— Phase A Task 0（两路线并行取证版）。

本模块是**唯一供 Task 6（`annulus_d2dga` 接线）消费的进口面**：把 8 口井
各水泥相（领浆 / 中间浆 / 尾浆）的 τ_y 与它的"读数 + 出处 + 温度"冻结记录
放在一起，并给出可复算的拟合函数，使每一条常数都能被契约测试从读数重算。

两路线并行取证（用户裁定：**两条路线都出证据、不做取舍**）
----------------------------------------------------------
- **路线 A `report_given`**：化验报告**原文直接给出**的屈服应力。
  取值面 = CSV `source_file` + `source_location` 指明的报告表。
- **路线 B `fitted_new`**：从六速读数**新推导**（HB 三参数取 τ_y；θ₆₀ 可用时
  另算 Bingham `YP = 2θ₃₀ − θ₆₀₀`）。附带现场工具箱口径的旁证
  （Bingham 线性最小二乘截距）。

两条路线**并列记录、分列标注来源性质**：不合并、不取平均、不互相填补。
某井某相某路线取不到时，登记为该路线的**显式缺失**（`MISSING` 哨兵），
**不是静默 `None`**，也**不许跨路线借值**。

路线 A 结论（本模块取证的事实，不含"哪条更权威"的判断）
------------------------------------------------------
8 口井的化验报告"流变性能"表（HT1-001 表7 / HT1-002 表8 / HT1-003 表7 /
HT1-004 表7 / hu102 20234.doc / hu103 203111.docx 表7 / hu1 204131.doc /
hu101 2011121.doc·2011122.pdf）**只给 `n` 与 `K`（"流变模式：幂律模式"），
没有任何一条给出水泥浆屈服应力**。设计书、施工记录表、作业史里"屈服值"
只出现在**钻井液 / 前置液**行（如 HT1-004 施工记录表"结束前钻井液性能
屈服值 10Pa"、hu103 20214.doc"平衡液 屈服值 1.5Pa"），**水泥浆行只有密度与
抗压强度**。⇒ 路线 A 在 8 口井 × 21 条"井×相"记录上**覆盖率为 0**。

另有一条**不可归属**的旁证：`0708/2/206/2061/20611/206114/2061144.xls`
（固井工程计算工具箱）的 `流变参数计算`/`流变曲线1` 两表**确实**算过水泥浆的
"屈服值，Pa"与"赫巴模型 Txs/YP/n/k"，但该工作簿**不含任何井号/井段标识**，
其读数与 8 口井任一浆体都不对应，故**不构成路线 A 的 `report_given` 值**；
它只被用来锚定路线 B 的**计算口径**（见下）。

路线 B 口径与换算常数
---------------------
Fann 35 换算（全篇统一；另一常见口径 0.4788 不得混用）：
- 剪切速率  γ̇ [s⁻¹] = 1.703 × RPM
- 剪切应力  τ [Pa]  = 0.511 × θ（θ = 刻度盘读数，单位"格"）

- `fann_herschel_bulkley`：τ = τ_y + K·γ̇ⁿ 的三参数最小二乘（τ 残差），
  约束 τ_y ∈ [0, min τ_i]、K > 0、n ∈ [0.05, 1.6]。**本数据上不可辨识**
  （见 §可辨识性），故其值是"新推导值"、带明确方法学代价。
- `fann_bingham_two_point`：API 口径 `YP = 0.511·(2θ₃₀₀ − θ₆₀)`，
  **仅 θ₆₀ 是真实读数时可用**；本语料只有 ht1_003 两相满足。
- `fann_bingham_yield_stress`：Bingham 两参数**线性最小二乘截距**——
  现场固井工程计算表 `2061144.xls` 的"屈服值，Pa"就是这个量，本函数对其
  `流变曲线1` 三组读数**逐位复现**（见契约测试）。⚠️ 它**不是**路线 A 的值
  （工作簿不可归属），只作路线 B 的口径旁证。

路线 B 与报告给定 n/K 的口径自洽性
----------------------------------
**不自洽，已知且不得掩盖**：报告给定 `n/K` 是幂律两参数拟合的产物，
屈服应力需要第三个参数；本模块的 τ_y 来自另一个模型（HB 或 Bingham），
**与报告给定 `n/K` 不是同一次拟合的结果**。R8 的口径裁定因此是：
**保留 loader 既有幂律 `n/K` 不动、只新增 τ_y**——不改任何 `FluidSpec`、
不重拟 `n/K`，2026-08-29 校准基线逐位不变。

可辨识性（路线 B 的代价，供裁定参考）
------------------------------------
5 点 Fann 读数跨 100:1 剪切速率，HB 三参数**不可辨识**：
- 绝对残差（τ）与相对残差（ln τ）两种度量给出的 τ_y 相差 0–3 Pa；
- 严格三点恰定解的 τ_y 在 −70…+50 Pa 之间跳变；
- 21 条记录里有 **12 条**落在 τ_y = 0 边界（数据支持纯幂律）：
  hu102 lead/tail、hu103 lead、hu1 lead、hu2 lead/intermediate/tail、
  ht1_001 lead/intermediate/tail、ht1_004 lead/tail；
- 读数 ±0.5 格量化扰动下，HB 的 τ_y 标准差均值 0.075 Pa、最大 0.157 Pa。

未经外部验证
------------
全部 τ_y **无外部（文献/第三方）验证**，也未与任何独立测量（静胶凝强度、
vane 法屈服应力等）对照。Phase A 的 A/B 归因必须标注其来源性质。
**不得引 Zhang & Frigaard 2022/2023 作为 Herschel-Bulkley 或斜井依据**
（该两篇自述仅适用竖直井 + 牛顿流体）。

R7 说明（本模块的边界）
-----------------------
本模块**只新增常数**。它不改任何水泥相 `FluidSpec` 的 `rheology_model`，
也不给它填 `yield_stress_pa`：`annulus_d2dga.py` 的既有屈服门消费
`FluidSpec.yield_stress_pa`（`_fluid_yield_stress` → 相体积加权 τ_y 场），
一旦改动就会在默认参数下改变求解结果，破坏 Phase A"默认关逐位 = HEAD"的硬约束。
本模块的值只供 Task 6 的 `hb_fix_cement_tau_y`（默认 False）消费，
且**哪一套被生产消费由 controller 裁定**——本模块两套并列暴露、不预设。
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Final, Mapping, Union

import numpy as np
from scipy.optimize import minimize_scalar

from cemdisp.data.loaders.ht1_001_loader import (
    HT1_001_INTERMEDIATE_YIELD_STRESS_FITTED_PA,
    HT1_001_LEAD_YIELD_STRESS_FITTED_PA,
    HT1_001_TAIL_YIELD_STRESS_FITTED_PA,
)
from cemdisp.data.loaders.ht1_003_loader import (
    HT1_003_LEAD_YIELD_STRESS_FITTED_PA,
    HT1_003_TAIL_YIELD_STRESS_FITTED_PA,
)
from cemdisp.data.loaders.ht1_004_loader import (
    HT1_004_LEAD_YIELD_STRESS_FITTED_PA,
    HT1_004_TAIL_YIELD_STRESS_FITTED_PA,
)
from cemdisp.data.loaders.hu1_loader import (
    HU1_LEAD_YIELD_STRESS_FITTED_PA,
    HU1_TAIL_YIELD_STRESS_FITTED_PA,
)
from cemdisp.data.loaders.hu101_loader import (
    HU101_LEAD_RECHECK_YIELD_STRESS_FITTED_PA,
    HU101_LEAD_YIELD_STRESS_FITTED_PA,
    HU101_TAIL_RECHECK_YIELD_STRESS_FITTED_PA,
    HU101_TAIL_YIELD_STRESS_FITTED_PA,
)
from cemdisp.data.loaders.hu102_loader import (
    HU102_LEAD_YIELD_STRESS_FITTED_PA,
    HU102_TAIL_YIELD_STRESS_FITTED_PA,
)
from cemdisp.data.loaders.hu103_loader import (
    HU103_INTERMEDIATE_YIELD_STRESS_FITTED_PA,
    HU103_LEAD_YIELD_STRESS_FITTED_PA,
    HU103_TAIL_YIELD_STRESS_FITTED_PA,
)
from cemdisp.data.loaders.hu2_loader import (
    HU2_INTERMEDIATE_YIELD_STRESS_FITTED_PA,
    HU2_LEAD_YIELD_STRESS_FITTED_PA,
    HU2_TAIL_YIELD_STRESS_FITTED_PA,
)

# --------------------------------------------------------------------------
# 缺项哨兵：取不到的值一律显式登记为 MISSING，**不许静默 None**
# --------------------------------------------------------------------------
MISSING: Final[str] = "MISSING"

# 两路线标识（供 Task 6/7 显式选择，语料里不预设默认路线）
ROUTE_REPORT_GIVEN: Final[str] = "report_given"
ROUTE_FITTED: Final[str] = "fitted_new"
ROUTES: Final[tuple[str, str]] = (ROUTE_REPORT_GIVEN, ROUTE_FITTED)

# --------------------------------------------------------------------------
# Fann 35 换算（全篇统一口径；另一常见口径 0.4788 不得混用）
# --------------------------------------------------------------------------
FANN_STRESS_PA_PER_DEGREE: Final[float] = 0.511
FANN_SHEAR_RATE_PER_RPM: Final[float] = 1.703

_HB_N_BOUNDS: Final[tuple[float, float]] = (0.05, 1.60)
_OPT_XATOL: Final[float] = 1e-12
# loader 常数按 4 位小数取整；容差取取整半量，既有物理意义又足以锁住读数改动。
_FITTED_TOLERANCE_PA: Final[float] = 5e-5

Readings = tuple[tuple[int, float], ...]
MaybeFloat = Union[float, str]


# --------------------------------------------------------------------------
# 拟合函数（确定性；契约测试用同一函数从冻结读数复算）
# --------------------------------------------------------------------------
def _to_si(readings: Readings) -> tuple[np.ndarray, np.ndarray]:
    """把 (RPM, 刻度读数) 转成 (γ̇ [s⁻¹], τ [Pa])。"""
    gamma = np.array([FANN_SHEAR_RATE_PER_RPM * rpm for rpm, _ in readings], dtype=float)
    tau = np.array([FANN_STRESS_PA_PER_DEGREE * theta for _, theta in readings], dtype=float)
    return gamma, tau


def fann_herschel_bulkley(readings: Readings) -> tuple[float, float, float]:
    """HB 三参数最小二乘（τ 残差），约束 τ_y∈[0, min τ]、K>0、n∈[0.05, 1.6]。

    ⚠️ 本数据上**不可辨识**（见模块 docstring §可辨识性），故这是"新推导值"。
    外层对 n 有界极小化、内层对 τ_y 有界极小化、K 由线性最小二乘解析给出，
    故结果对初值不敏感、完全确定。

    Returns:
        (tau_y_pa, consistency_k_pa_s_n, power_law_n)
    """
    gamma, tau = _to_si(readings)
    tau_min = float(tau.min())

    def residual(tau_y: float, n: float) -> tuple[float, float]:
        powered = gamma**n
        k = float((tau - tau_y) @ powered / (powered @ powered))
        if k <= 0.0:
            return 1e18, k
        return float(np.sum((tau - tau_y - k * powered) ** 2)), k

    def best_for_n(n: float) -> tuple[float, float]:
        result = minimize_scalar(
            lambda tau_y: residual(tau_y, n)[0],
            bounds=(0.0, tau_min),
            method="bounded",
            options={"xatol": _OPT_XATOL},
        )
        return float(result.fun), float(result.x)

    outer = minimize_scalar(
        lambda n: best_for_n(n)[0],
        bounds=_HB_N_BOUNDS,
        method="bounded",
        options={"xatol": _OPT_XATOL},
    )
    n = float(outer.x)
    _, tau_y = best_for_n(n)
    _, k = residual(tau_y, n)
    return tau_y, k, n


def fann_bingham_yield_stress(readings: Readings) -> tuple[float, float, float]:
    """Bingham 两参数线性最小二乘：τ = τ_y + μ_p·γ̇。

    与现场固井工程计算表 `2061144.xls / 流变曲线1` 的"屈服值，Pa"口径一致
    （该表已逐位复现，见契约测试）。**注意**：该工作簿不可归属到 8 口井，
    故这是**口径旁证**，不是路线 A 的 `report_given` 值。

    Returns:
        (tau_y_pa, plastic_viscosity_pa_s, r_squared)
    """
    gamma, tau = _to_si(readings)
    design = np.vstack([np.ones_like(gamma), gamma]).T
    solution, *_ = np.linalg.lstsq(design, tau, rcond=None)
    residual = tau - design @ solution
    total = float(np.sum((tau - tau.mean()) ** 2))
    r_squared = 1.0 if total == 0.0 else 1.0 - float(np.sum(residual**2)) / total
    return float(solution[0]), float(solution[1]), r_squared


def fann_bingham_two_point(readings: Readings) -> MaybeFloat:
    """API 口径 Bingham 动切力 `YP = 0.511·(2θ₃₀₀ − θ₆₀)` [Pa]。

    **仅 θ₆₀₀ 是真实读数时可用**（超量程 `>300` 与未记录 `/` 都不可用——
    把 300 当读数会直接把 YP 算小一半以上）。本语料只有 ht1_003 两相满足；
    其余返回 `MISSING`，**不许用 2θ₂₀₀−θ₃₀₀ 之类凑一个数顶上**。
    """
    by_rpm = dict(readings)
    if 600 not in by_rpm or 300 not in by_rpm:
        return MISSING
    return FANN_STRESS_PA_PER_DEGREE * (2.0 * by_rpm[300] - by_rpm[600])


def fann_power_law(readings: Readings) -> tuple[float, float, float]:
    """幂律两参数对数-对数最小二乘：τ = K·γ̇ⁿ（**与化验报告给定 n/K 同口径**）。

    回归在 (ln γ̇, ln θ) 空间做，截距 b 满足 τ = 0.511·e^b·γⁿ，
    故 **K = 0.511·exp(b)** [Pa·sⁿ]；返回 (n, K, r_squared)。

    实测锚（2026-09-16）：本口径逐位复现化验报告给定 K（hu101 复检
    0.8147→0.815 等），也**逐位复现**现场固井工程计算表 `2061144.xls / 流变曲线1`
    的全部三列——领浆 0.9744824296、尾浆 0.4169271535、钻井液 0.3975288753
    （契约测试按 rel=1e-9 断言，见
    `test_toolkit_yield_stress_convention_is_bitwise_reproduced`）。

    订正留痕：本篇曾一度称"该表领浆列漏乘 0.511、属工作簿内部不一致"——
    经复核为**手工换算口误**，三列实为全部一致；**不要**据此改动本函数口径
    （改成不带 0.511 会立刻复现不出任何一条化验给定 K，并把测试搞红）。
    """
    gamma, _ = _to_si(readings)
    theta = np.array([value for _, value in readings], dtype=float)
    log_gamma, log_theta = np.log(gamma), np.log(theta)
    slope, intercept = np.polyfit(log_gamma, log_theta, 1)
    n = float(slope)
    consistency_k = FANN_STRESS_PA_PER_DEGREE * float(np.exp(intercept))
    predicted = n * log_gamma + float(intercept)
    total = float(np.sum((log_theta - log_theta.mean()) ** 2))
    r_squared = 1.0 if total == 0.0 else 1.0 - float(np.sum((log_theta - predicted) ** 2)) / total
    return n, consistency_k, r_squared


# --------------------------------------------------------------------------
# 逐井逐浆体的冻结读数 + 出处（读数改动即测试变红）
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class CementYieldStress:
    """单个"井 × 水泥相"的两路线 τ_y 记录（读数 + 出处 + 两套值）。"""

    well_key: str
    """井键：hu101 / hu102 / hu103 / hu1 / hu2 / ht1_001 / ht1_003 / ht1_004。"""

    phase: str
    """水泥相：lead / intermediate / tail（hu101 另有 lead_recheck / tail_recheck）。"""

    fluid_name: str
    """loader `FluidSpec` 里的流体名（与该井 loader 实际构造一致）。"""

    temperature_c: str
    """源报告记录的试验温度（保持原文写法，可能是"128→93"这类区间）。"""

    readings: Readings
    """冻结的 Fann 六速读数 ((RPM, 刻度读数), ...)；超量程/未记录的速度点已剔除。"""

    source_file: str
    """出处文件；CSV 有该行时取 `rheometer_readings.csv` 的 `source_file` 原字符串。"""

    source_location: str
    """出处小节；CSV 有该行时取 `rheometer_readings.csv` 的 `source_location` 原字符串。"""

    in_rheometer_csv: bool
    """读数是否能在 `参考文档/现场资料提取/<井>/rheometer_readings.csv` 中直接查到。"""

    note: str
    """补充说明（读取口径、双报告并存、超量程记录方式等）。"""

    report_given_pa: MaybeFloat
    """**路线 A**：报告原文直接给出的屈服应力 [Pa]；未给出时为 `MISSING`。"""

    report_given_note: str
    """路线 A 的取证说明：查了哪张表、为什么没有。"""

    fitted_hb_tau_y_pa: float
    """**路线 B**：HB 三参数拟合 τ_y [Pa]（loader 常数 `..._FITTED_PA` 的同一口径）。"""

    fitted_hb_consistency_k_pa_s_n: float
    """路线 B：HB 拟合的 K [Pa·sⁿ]。"""

    fitted_hb_power_law_n: float
    """路线 B：HB 拟合的 n。"""

    fitted_bingham_two_point_pa: MaybeFloat
    """路线 B 旁枝：API 口径 `0.511·(2θ₀₀ − θ₆₀)`；θ₆₀₀ 不可用时 `MISSING`。"""

    fitted_bingham_ls_intercept_pa: float
    """路线 B 旁证：Bingham 线性最小二乘截距（现场工程计算表口径）。"""

    fitted_bingham_ls_r_squared: float
    """该截距拟合的决定系数 R²。"""


# 路线 A 的逐井取证说明（**同一结论、逐井留痕**：报告只给 n/K）
_A_NO_CEMENT_YS = (
    "该井化验报告『流变性能』表（{table}）对水泥浆只给 `n` 与 `K`"
    "（『流变模式：幂律模式』），无 YP / 屈服值 / PV-YP 对；"
    "设计书与施工记录表里的『屈服值』全部落在钻井液/前置液行，水泥浆行只有密度与抗压强度。"
)

_SOURCE_HU101_PRIMARY = "2/201/2011/20111/201112/2011122.pdf"
_SOURCE_HU101_RECHECK = "2/201/2011/20111/201112/2011121.doc"
_SOURCE_HU102 = "2021/20234.doc"
_SOURCE_HU103 = "2\\203\\2031\\20311\\203111.docx"
_SOURCE_HU1 = "2/204/2041/20413/204131.doc"
_SOURCE_HT1_002 = "化验报告 Table7"
_SOURCE_HT1_001 = "化验报告/Table 8"
_SOURCE_HT1_003 = "化验报告/HT1-003 油层尾管 化验报告.docx 表7"
_SOURCE_HT1_004 = "化验报告/HT1-004 油层尾管 化验报告.docx 表7"


def _record(
    *,
    well_key: str,
    phase: str,
    fluid_name: str,
    temperature_c: str,
    readings: Readings,
    source_file: str,
    source_location: str,
    in_rheometer_csv: bool,
    report_given: MaybeFloat,
    report_given_note: str,
    fitted: float,
    note: str = "",
) -> CementYieldStress:
    """按冻结读数复算路线 B 各量，并与 loader 常数做一致性自检。"""
    hb_tau_y, hb_k, hb_n = fann_herschel_bulkley(readings)
    ls_tau_y, _, ls_r_squared = fann_bingham_yield_stress(readings)
    if abs(hb_tau_y - fitted) > _FITTED_TOLERANCE_PA:
        raise AssertionError(
            f"{well_key}/{phase} 路线 B 采用值与读数复算不一致："
            f"loader={fitted!r} 复算={hb_tau_y!r}"
        )
    return CementYieldStress(
        well_key=well_key,
        phase=phase,
        fluid_name=fluid_name,
        temperature_c=temperature_c,
        readings=readings,
        source_file=source_file,
        source_location=source_location,
        in_rheometer_csv=in_rheometer_csv,
        note=note,
        report_given_pa=report_given,
        report_given_note=report_given_note,
        fitted_hb_tau_y_pa=hb_tau_y,
        fitted_hb_consistency_k_pa_s_n=hb_k,
        fitted_hb_power_law_n=hb_n,
        fitted_bingham_two_point_pa=fann_bingham_two_point(readings),
        fitted_bingham_ls_intercept_pa=ls_tau_y,
        fitted_bingham_ls_r_squared=ls_r_squared,
    )


CEMENT_YIELD_STRESS: Final[Mapping[str, Mapping[str, CementYieldStress]]] = {
    "hu101": {
        "lead": _record(
            well_key="hu101", phase="lead", fluid_name="领浆", temperature_c="89",
            readings=((300, 144), (200, 104), (100, 57), (6, 8), (3, 4)),
            source_file=_SOURCE_HU101_PRIMARY,
            source_location="检测分析结果·领浆·流变性能（委托书 W301-22094）",
            in_rheometer_csv=False,
            report_given=MISSING,
            report_given_note=(
                _A_NO_CEMENT_YS.format(table="2011122.pdf 第2页『领浆·流变性能』")
                + " 取证受限说明：2011122.pdf 为**扫描件、无文字层**"
                  "（`_extracted/2011122_text.txt` 仅 34 字节），其流变读数由页面图像读取；"
                  "表内『流变模式 幂律模式 / 流变参数 n=0.844 K=0.381』即全部流变结论。"
            ),
            fitted=HU101_LEAD_YIELD_STRESS_FITTED_PA,
            note="主检报告（loader 生产口径 n/K=0.844/0.381 的同一份报告，89℃）；θ600 未测。",
        ),
        "tail": _record(
            well_key="hu101", phase="tail", fluid_name="尾浆", temperature_c="89",
            readings=((300, 122), (200, 82), (100, 49), (6, 7), (3, 4)),
            source_file=_SOURCE_HU101_PRIMARY,
            source_location="检测分析结果·尾浆·流变性能（委托书 W301-22094）",
            in_rheometer_csv=False,
            report_given=MISSING,
            report_given_note=(
                _A_NO_CEMENT_YS.format(table="2011122.pdf 第3页『尾浆·流变性能』")
                + " 取证受限说明同领浆：扫描件无文字层。"
            ),
            fitted=HU101_TAIL_YIELD_STRESS_FITTED_PA,
            note="主检报告（loader 生产口径 n/K=0.830/0.352 的同一份报告，89℃）；θ600 未测。",
        ),
        "lead_recheck": _record(
            well_key="hu101", phase="lead_recheck", fluid_name="领浆（复检口径）",
            temperature_c="93",
            readings=((300, 149), (200, 106), (100, 60), (6, 9), (3, 5)),
            source_file=_SOURCE_HU101_RECHECK, source_location="领浆流变",
            in_rheometer_csv=True,
            report_given=MISSING,
            report_given_note=_A_NO_CEMENT_YS.format(table="2011121.doc『领浆·流变性能』"),
            fitted=HU101_LEAD_RECHECK_YIELD_STRESS_FITTED_PA,
            note="复检报告口径（93℃），对应 loader 的 `HU101_LEAD_RECHECK_POWER_LAW_N/K`。",
        ),
        "tail_recheck": _record(
            well_key="hu101", phase="tail_recheck", fluid_name="尾浆（复检口径）",
            temperature_c="93",
            readings=((300, 126), (200, 87), (100, 53), (6, 8), (3, 4)),
            source_file=_SOURCE_HU101_RECHECK, source_location="尾浆流变",
            in_rheometer_csv=True,
            report_given=MISSING,
            report_given_note=_A_NO_CEMENT_YS.format(table="2011121.doc『尾浆·流变性能』"),
            fitted=HU101_TAIL_RECHECK_YIELD_STRESS_FITTED_PA,
            note="复检报告口径（93℃），对应 loader 的 `HU101_TAIL_RECHECK_*`。",
        ),
    },
    "hu102": {
        "lead": _record(
            well_key="hu102", phase="lead", fluid_name="领浆", temperature_c="93",
            readings=((300, 163), (200, 128), (100, 106), (6, 10), (3, 6)),
            source_file=_SOURCE_HU102, source_location="slurry rheology",
            in_rheometer_csv=True,
            report_given=MISSING,
            report_given_note=_A_NO_CEMENT_YS.format(table="20234.doc『slurry rheology』"),
            fitted=HU102_LEAD_YIELD_STRESS_FITTED_PA,
            note="θ600 缺失（CSV 原注）；领/尾浆同体系，读数与 n/K 相同。",
        ),
        "tail": _record(
            well_key="hu102", phase="tail", fluid_name="尾管水泥浆", temperature_c="93",
            readings=((300, 163), (200, 128), (100, 106), (6, 10), (3, 6)),
            source_file=_SOURCE_HU102, source_location="slurry rheology",
            in_rheometer_csv=True,
            report_given=MISSING,
            report_given_note=_A_NO_CEMENT_YS.format(table="20234.doc『slurry rheology』"),
            fitted=HU102_TAIL_YIELD_STRESS_FITTED_PA,
            note="与领浆同读数（20234.doc 记同体系）；θ600 缺失。",
        ),
    },
    "hu103": {
        "lead": _record(
            well_key="hu103", phase="lead", fluid_name="领浆", temperature_c="140→93",
            readings=((300, 224), (200, 158), (100, 93), (6, 9), (3, 5)),
            source_file=_SOURCE_HU103, source_location="Table7",
            in_rheometer_csv=True,
            report_given=MISSING,
            report_given_note=_A_NO_CEMENT_YS.format(table="203111.docx 表7"),
            fitted=HU103_LEAD_YIELD_STRESS_FITTED_PA,
            note="θ600 原文记 `/`（未记录），按缺测剔除，不得当 0 或 300 用。",
        ),
        "intermediate": _record(
            well_key="hu103", phase="intermediate", fluid_name="中间浆",
            temperature_c="140→93",
            readings=((300, 254), (200, 179), (100, 112), (6, 14), (3, 7)),
            source_file=_SOURCE_HU103, source_location="Table7",
            in_rheometer_csv=True,
            report_given=MISSING,
            report_given_note=_A_NO_CEMENT_YS.format(table="203111.docx 表7"),
            fitted=HU103_INTERMEDIATE_YIELD_STRESS_FITTED_PA,
            note="θ600 原文记 `/`（未记录）。",
        ),
        "tail": _record(
            well_key="hu103", phase="tail", fluid_name="尾浆", temperature_c="133→93",
            readings=((300, 265), (200, 186), (100, 115), (6, 15), (3, 7)),
            source_file=_SOURCE_HU103, source_location="Table7",
            in_rheometer_csv=True,
            report_given=MISSING,
            report_given_note=_A_NO_CEMENT_YS.format(table="203111.docx 表7"),
            fitted=HU103_TAIL_YIELD_STRESS_FITTED_PA,
            note="θ600 原文记 `/`；温度口径 133℃ 与领/中间浆 140℃ 不同。",
        ),
    },
    "hu1": {
        "lead": _record(
            well_key="hu1", phase="lead", fluid_name="领浆", temperature_c="93",
            readings=((300, 171), (200, 132), (100, 81), (6, 10), (3, 6)),
            source_file=_SOURCE_HU1,
            source_location="油井水泥浆物理性能试验结果·第1个浆体块·流变性能",
            in_rheometer_csv=True,
            report_given=MISSING,
            report_given_note=(
                _A_NO_CEMENT_YS.format(table="204131.doc 第1个浆体块『流变性能』")
                + " 补录留痕（2026-09-16 用户裁定）：该读数原**不在** hu1 的"
                  " rheometer_readings.csv 内（该 CSV 的 hu1 水泥行只有污染实验/密度高点/"
                  "温度高点三类），已按报告原文回填进该 CSV（`notes` 标『报告原文补录"
                  "（2026-09-16）』），故 `in_rheometer_csv` 已由 False 翻为 True。"
                  "出处为 loader 既有 n/K（0.732/0.933）所引用的同一份化验报告"
                  "204131.doc（签发 2020-10-24）原文浆体块，属同一出处链、非新数据源。"
            ),
            fitted=HU1_LEAD_YIELD_STRESS_FITTED_PA,
            note="读数取自 204131.doc 原文浆体块；θ600 未记录。",
        ),
        "tail": _record(
            well_key="hu1", phase="tail", fluid_name="尾浆", temperature_c="93",
            readings=((300, 117), (200, 86), (100, 52), (6, 9), (3, 5)),
            source_file=_SOURCE_HU1,
            source_location="油井水泥浆物理性能试验结果·第2个浆体块·流变性能",
            in_rheometer_csv=True,
            report_given=MISSING,
            report_given_note=(
                _A_NO_CEMENT_YS.format(table="204131.doc 第2个浆体块『流变性能』")
                + " 补录留痕同领浆（2026-09-16 用户裁定）：已按报告原文回填进"
                  " rheometer_readings.csv，`in_rheometer_csv` 由 False 翻为 True。"
            ),
            fitted=HU1_TAIL_YIELD_STRESS_FITTED_PA,
            note="取自 204131.doc（loader n/K=0.666/0.906 的同一出处）；θ600 未记录。",
        ),
    },
    "hu2": {
        "lead": _record(
            well_key="hu2", phase="lead", fluid_name="领浆", temperature_c="133→93",
            readings=((300, 268), (200, 194), (100, 113), (6, 10), (3, 7)),
            source_file=_SOURCE_HT1_002, source_location="化验报告 Table7",
            in_rheometer_csv=True,
            report_given=MISSING,
            report_given_note=_A_NO_CEMENT_YS.format(table="HT1-002 尾管化验报告 表8"),
            fitted=HU2_LEAD_YIELD_STRESS_FITTED_PA,
            note="θ600 记 `>300`（超量程），按缺测剔除，不得当 300 用。",
        ),
        "intermediate": _record(
            well_key="hu2", phase="intermediate", fluid_name="中间浆", temperature_c="133→93",
            readings=((300, 221), (200, 165), (100, 92), (6, 7), (3, 4)),
            source_file=_SOURCE_HT1_002, source_location="化验报告 Table7",
            in_rheometer_csv=True,
            report_given=MISSING,
            report_given_note=_A_NO_CEMENT_YS.format(table="HT1-002 尾管化验报告 表8"),
            fitted=HU2_INTERMEDIATE_YIELD_STRESS_FITTED_PA,
            note="θ600 记 `>300`（超量程），剔除。",
        ),
        "tail": _record(
            well_key="hu2", phase="tail", fluid_name="尾浆", temperature_c="133→93",
            readings=((300, 218), (200, 163), (100, 88), (6, 6), (3, 4)),
            source_file=_SOURCE_HT1_002, source_location="化验报告 Table7",
            in_rheometer_csv=True,
            report_given=MISSING,
            report_given_note=_A_NO_CEMENT_YS.format(table="HT1-002 尾管化验报告 表8"),
            fitted=HU2_TAIL_YIELD_STRESS_FITTED_PA,
            note="θ600 记 `>300`（超量程），剔除。",
        ),
    },
    "ht1_001": {
        "lead": _record(
            well_key="ht1_001", phase="lead", fluid_name="领浆", temperature_c="128→93",
            readings=((300, 268), (200, 194), (100, 113), (6, 10), (3, 7)),
            source_file=_SOURCE_HT1_001, source_location="128->93C降温测试",
            in_rheometer_csv=True,
            report_given=MISSING,
            report_given_note=_A_NO_CEMENT_YS.format(table="HT1-001 尾管化验报告 表7"),
            fitted=HT1_001_LEAD_YIELD_STRESS_FITTED_PA,
            note=(
                "θ600 化验报告记 `>300`（超量程），CSV 却写作 300——**数据陷阱之一**，"
                "按超量程剔除，不得当 300 用。"
            ),
        ),
        "intermediate": _record(
            well_key="ht1_001", phase="intermediate", fluid_name="中间浆",
            temperature_c="128→93",
            readings=((300, 221), (200, 165), (100, 92), (6, 7), (3, 4)),
            source_file=_SOURCE_HT1_001, source_location="128->93C降温测试",
            in_rheometer_csv=True,
            report_given=MISSING,
            report_given_note=_A_NO_CEMENT_YS.format(table="HT1-001 尾管化验报告 表7"),
            fitted=HT1_001_INTERMEDIATE_YIELD_STRESS_FITTED_PA,
            note="θ600 化验报告记 `>300`，CSV 写作 300；按超量程剔除。",
        ),
        "tail": _record(
            well_key="ht1_001", phase="tail", fluid_name="尾浆", temperature_c="128→93",
            readings=((300, 218), (200, 163), (100, 88), (6, 6), (3, 4)),
            source_file=_SOURCE_HT1_001, source_location="128->93C降温测试",
            in_rheometer_csv=True,
            report_given=MISSING,
            report_given_note=_A_NO_CEMENT_YS.format(table="HT1-001 尾管化验报告 表7"),
            fitted=HT1_001_TAIL_YIELD_STRESS_FITTED_PA,
            note="θ600 化验报告记 `>300`，CSV 写作 300；按超量程剔除。",
        ),
    },
    "ht1_003": {
        "lead": _record(
            well_key="ht1_003", phase="lead", fluid_name="领浆", temperature_c="129→93",
            readings=((600, 207), (300, 121), (200, 93), (100, 60), (6, 15), (3, 9)),
            source_file=_SOURCE_HT1_003, source_location="表7 流变性能 · 领浆",
            in_rheometer_csv=True,
            report_given=MISSING,
            report_given_note=(
                _A_NO_CEMENT_YS.format(table="HT1-003 油层尾管化验报告 表7")
                + " 本表列头为『温度|浆体名称|Ф600|Ф300|Ф200|Ф100|Ф6|Ф3|n|K Pa.sn』——"
                  "**列头里没有 YP/屈服值列**（该井 θ600 是真实读数 207，仍只给 n/K）。"
            ),
            fitted=HT1_003_LEAD_YIELD_STRESS_FITTED_PA,
            note=(
                "本井 θ600 是**真实读数**（207），未超量程，故参与拟合（唯一含 600 的井，"
                "因而也是唯一能算 API 口径 `2θ₃₀₀−θ₆₀` 的井）。"
                "⚠️ 报告同表给定的 n=0.597/K=1.622 与本表读数**不自洽**（Step 1 不可复现样本）。"
            ),
        ),
        "tail": _record(
            well_key="ht1_003", phase="tail", fluid_name="尾浆", temperature_c="129→93",
            readings=((600, 196), (300, 118), (200, 89), (100, 56), (6, 14), (3, 10)),
            source_file=_SOURCE_HT1_003, source_location="表7 流变性能 · 尾浆",
            in_rheometer_csv=True,
            report_given=MISSING,
            report_given_note=_A_NO_CEMENT_YS.format(table="HT1-003 油层尾管化验报告 表7"),
            fitted=HT1_003_TAIL_YIELD_STRESS_FITTED_PA,
            note="θ600=196 为真实读数；同表给定 n=0.585/K=1.673 与读数不自洽（同上）。",
        ),
    },
    "ht1_004": {
        "lead": _record(
            well_key="ht1_004", phase="lead", fluid_name="领浆", temperature_c="132→93",
            readings=((300, 296), (200, 213), (100, 120), (6, 10), (3, 6)),
            source_file=_SOURCE_HT1_004, source_location="表7 流变性能 · 领浆",
            in_rheometer_csv=True,
            report_given=MISSING,
            report_given_note=(
                _A_NO_CEMENT_YS.format(table="HT1-004 油层尾管化验报告 表7")
                + " ⚠️ 本井**生产口径水泥相是 Bingham**（优化参数化 PV0.17/YP13.0 与 "
                  "PV0.18/YP14.0，来源『优化参数.docx 2026-06-11』，非化验报告）——"
                  "那两个 YP **不是**化验报告的 `report_given` 值，而是 MATLAB 反演优化值，"
                  "故同样登记为 `MISSING`。"
            ),
            fitted=HT1_004_LEAD_YIELD_STRESS_FITTED_PA,
            note="θ600 记 `>300`（CSV 原注 Ф600 超量程），剔除。",
        ),
        "tail": _record(
            well_key="ht1_004", phase="tail", fluid_name="尾浆", temperature_c="132→93",
            readings=((300, 293), (200, 211), (100, 120), (6, 10), (3, 5)),
            source_file=_SOURCE_HT1_004, source_location="表7 流变性能 · 尾浆",
            in_rheometer_csv=True,
            report_given=MISSING,
            report_given_note=_A_NO_CEMENT_YS.format(table="HT1-004 油层尾管化验报告 表7"),
            fitted=HT1_004_TAIL_YIELD_STRESS_FITTED_PA,
            note="θ600 记 `>300`，剔除；生产口径同领浆为 Bingham（见上）。",
        ),
    },
}


# --------------------------------------------------------------------------
# 供 Task 6 消费的取用接口（**两路线并列暴露，本模块不预设生产路线**）
# --------------------------------------------------------------------------
def get_cement_yield_stress(well_key: str, phase: str) -> CementYieldStress:
    """按井键 + 相取两路线 τ_y 记录；取不到时抛 KeyError（不返回 None）。

    Args:
        well_key: `hu101` / `hu102` / `hu103` / `hu1` / `hu2` / `ht1_001` /
            `ht1_003` / `ht1_004`。
        phase: `lead` / `intermediate` / `tail`；hu101 另有 `lead_recheck` /
            `tail_recheck`（复检报告口径）。

    Raises:
        KeyError: 井键或相不存在。
    """
    try:
        return CEMENT_YIELD_STRESS[well_key][phase]
    except KeyError as error:
        raise KeyError(
            f"未登记的水泥 τ_y 记录：well={well_key!r} phase={phase!r}；"
            f"已登记井={sorted(CEMENT_YIELD_STRESS)}"
        ) from error


def yield_stress_by_route(well_key: str, phase: str, route: str) -> MaybeFloat:
    """显式按路线取值：`report_given` 或 `fitted_new`。

    取不到该路线的值时返回 `MISSING` 哨兵（**不返回 None、不跨路线借值**）。
    生产消费哪条路线由 controller 裁定，调用方必须显式传 `route`。

    Raises:
        KeyError: 井键或相不存在。
        ValueError: route 不是 `ROUTES` 之一。
    """
    if route not in ROUTES:
        raise ValueError(f"route 必须为 {ROUTES} 之一，得到 {route!r}")
    record = get_cement_yield_stress(well_key, phase)
    if route == ROUTE_REPORT_GIVEN:
        return record.report_given_pa
    return record.fitted_hb_tau_y_pa


def yield_stress_by_role(well_key: str, route: str) -> Mapping[str, MaybeFloat]:
    """返回该井 `FluidRole` 名 → τ_y 的映射（供 Task 6 按 role 取用）。

    只覆盖 loader 实际消费的水泥相；hu101 的复检口径不在其中（仅供敏感性）。
    """
    roles = {"lead": "LEAD", "intermediate": "INTERMEDIATE", "tail": "TAIL"}
    return {
        roles[phase]: yield_stress_by_route(well_key, phase, route)
        for phase in CEMENT_YIELD_STRESS[well_key]
        if phase in roles
    }