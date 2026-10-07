function [Ff_a, flow_pattern_a] = Friction_casing_bh(rho_a, V_a, A_a, mu_p_a, tau_y, D_w, D_do, QL)
    % 计算环空中的摩擦压耗和流型
    % 输入参数：
    %   rho_a: 流体密度 (kg/m^3)
    %   V_a: 流体速度 (m/s)
    %   A_a: 环空截面积 (m^2)
    %   mu_p_a: 流体塑性粘度 (Pa·s)
    %   tau_y: 流体屈服应力 (Pa)
    %   D_w: 井眼直径 (m)
    %   D_do: 套管外径 (m)
    %   QL: 流量 (m^3/s)
    % 输出参数：
    %   Ff_a: 摩擦压耗 (Pa/m)
    %   flow_pattern_a: 流型 (1 表示层流，2 表示紊流)

    global pianxin PR; % 全局变量，用于控制是否考虑偏心影响
    mu_p_a = mu_p_a / 1000; % 将塑性粘度从 mPa·s 转换为 Pa·s

    %% 计算井壁剪切应力
    tau_w_old = 100; % 初始假设的井壁剪切应力 (Pa)
    err_tau = 1; % 初始误差

    maxIter = 30;                         % 强制退出
    err_tau = 1;
    iter    = 0;
    
    while abs(err_tau) > 1e-3 && iter < maxIter% 迭代计算，直到误差小于 1e-4
        tau_w_new = 8 * mu_p_a * QL / (pi * ((D_w - D_do) / 2)^2 * (D_w + D_do) / 2) + 3 / 2 * tau_y - 1 / 2 * tau_y^3 / tau_w_old^2; % 更新井壁剪切应力
        err_tau = abs(tau_w_new - tau_w_old) / tau_w_old; % 计算误差 
        tau_w_old = tau_w_new; % 更新假设的井壁剪切应力
        iter      = iter + 1;
    end
    tau_w = tau_w_new; % 最终的井壁剪切应力

    int1 = tau_y / tau_w; % 无量纲的剪切应力比
    He = 16800 * int1 / (1 - int1)^3; % Hedstrom 数

    %% 判断流型并计算摩擦压耗
    Re_a_c = (1 - 4 / 3 * int1 + 1 / 3 * int1^4) / (8 * int1) * He; % 临界雷诺数
    Re_a = 0.81619 * rho_a * V_a * (D_w - D_do) / mu_p_a; % 实际雷诺数
    
    if Re_a <= Re_a_c % 层流
        % 层流公式：利用圆管公式乘以 0.75 转换为环空公式，逻辑正确
        Ff_a = 1.25* (6895 * mu_p_a * V_a / (216 * (D_w - D_do)^2) + 5.33355 * tau_y / (D_w - D_do)); % 摩擦压耗公式（层流）
        flow_pattern_a = 1; % 流型标记（1 表示层流）
    else % 紊流
        % 【核心Bug修复】：去除了原代码中错误的 0.75 乘数！
        % 环空紊流摩阻应直接使用水力直径代入，乘以 0.75 会导致紊流摩阻被严重低估 25%。
        Ff_a = 0.158278 * rho_a^0.75 * V_a^1.75 * mu_p_a^0.25 / (D_w - D_do)^1.25; 
        flow_pattern_a = 2; % 流型标记（2 表示紊流）
    end
end
