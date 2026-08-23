"""M2 局部流态修正的纯函数闭包（Maleki & Frigaard 2017 式58-66）。

所有函数数组友好：标量或 np.ndarray 逐元运算，返回 ndarray。
几何约定：b=effective_b（全间隙，与求解器 γ̇=6|w|/b 自洽）。
"""
import numpy as np


def metzner_reed_re(w, rho, n, kappa, b):
    """式58：Re_p = 6ρw²/(κ γ̇_N^n)，γ̇_N=6|w|/b。"""
    w = np.asarray(w, dtype=float)
    b = np.maximum(np.asarray(b, dtype=float), 1e-9)
    gamma = np.maximum(6.0 * np.abs(w) / b, 1e-6)
    tau_app = np.maximum(np.asarray(kappa, dtype=float) * gamma ** np.asarray(n, dtype=float), 1e-12)
    return 6.0 * np.asarray(rho, dtype=float) * w * w / tau_app


def hedstrom_number(tau_y, rho, n, kappa, b):
    """式59：He = τy [ρ^n b^(2n)/κ²]^(1/(2-n))。"""
    tau_y = np.asarray(tau_y, dtype=float)
    rho = np.asarray(rho, dtype=float)
    n = np.asarray(n, dtype=float)
    kappa = np.maximum(np.asarray(kappa, dtype=float), 1e-12)
    b = np.maximum(np.asarray(b, dtype=float), 1e-9)
    inner = rho ** n * b ** (2.0 * n) / kappa ** 2
    return tau_y * np.power(np.maximum(inner, 1e-300), 1.0 / (2.0 - n))


def friction_laminar(re_p, he, n):
    """式60-61：层流 f=24/Re_p，含塞流核修正 yY=He/Hw。

    无屈服(he=0)时退化为 24/Re_p。塞流核幂次以 Maleki 式60-61 为准，
    实现后用文献算例复核；此处给标准槽流形式 (1-yY)^2。
    """
    re_p = np.maximum(np.asarray(re_p, dtype=float), 1e-9)
    he = np.asarray(he, dtype=float)
    hw = 24.0 / re_p
    yY = np.clip(he / np.maximum(hw, 1e-9), 0.0, 0.95)
    return 24.0 * (1.0 - yY) ** 2 / re_p


def friction_dodge_metzner(re_p, n, f0=0.01, tol=1e-4, it_max=15):
    """式65：1/√f = (4/n^0.75) log10(Re f^(1-n/2)) - 0.4/n^1.2。固定点迭代。"""
    re_p = np.asarray(re_p, dtype=float)
    n = np.asarray(n, dtype=float)
    f = np.full_like(re_p, f0, dtype=float)
    for _ in range(it_max):
        rhs = ((4.0 / n ** 0.75) * np.log10(np.maximum(re_p * f ** (1.0 - n / 2.0), 1e-30))
               - 0.4 / n ** 1.2)
        f_new = 1.0 / np.maximum(rhs, 1e-9) ** 2
        if np.all(np.abs(f_new - f) / np.maximum(f, 1e-12) < tol):
            f = f_new
            break
        f = f_new
    return f


def friction_transition(re_p, re_crit, re_turb, f_lam_cr, f_turb):
    """过渡区 log-Re 空间线性插值（禁硬切换）。"""
    log_re = np.log10(np.maximum(np.asarray(re_p, dtype=float), 1e-9))
    lo = np.log10(np.maximum(np.asarray(re_crit, dtype=float), 1e-9))
    hi = np.log10(np.maximum(np.asarray(re_turb, dtype=float), 1e-9))
    t = np.clip((log_re - lo) / np.maximum(hi - lo, 1e-9), 0.0, 1.0)
    log_f = (np.log10(np.maximum(f_lam_cr, 1e-12))
             + t * (np.log10(np.maximum(f_turb, 1e-12)) - np.log10(np.maximum(f_lam_cr, 1e-12))))
    return 10.0 ** log_f


def drag_weight(re_p, he, n, re_crit, re_turb_ratio=1.8):
    """返回 (R, regime_mask)。R 乘到 pref 上；mask: 0 层流/1 过渡/2 湍流。

    层流 R=1（保持现有 b² 幂律，不降指数）；过渡/湍流 R=clip(f_lam_cr/f_eff,0.3,1.5)。
    总面积归一保证总流量守恒。Re_crit 由调用方按 He 给出（屈服推迟转捩）。
    """
    re_p = np.asarray(re_p, dtype=float)
    he = np.asarray(he, dtype=float)
    n = np.asarray(n, dtype=float)
    re_crit = np.asarray(re_crit, dtype=float)
    re_turb = re_crit * re_turb_ratio
    f_lam_cr = friction_laminar(re_crit, he, n)
    f_turb = friction_dodge_metzner(np.maximum(re_turb, 1e-9), n)
    f_eff = friction_transition(re_p, re_crit, re_turb, f_lam_cr, f_turb)
    R = np.clip(f_lam_cr / np.maximum(f_eff, 1e-9), 0.3, 1.5)
    R = np.where(re_p <= re_crit, 1.0, R)
    mask = np.where(re_p <= re_crit, 0, np.where(re_p >= re_turb, 2, 1)).astype(int)
    return R, mask
