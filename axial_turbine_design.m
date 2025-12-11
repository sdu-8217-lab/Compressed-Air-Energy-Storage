function turbine_result = axial_turbine_design(T_in, P_in, m, c_in, expansion_ratio, eta_ts, R, k, Cp)
% 单级轴流式透平机热力计算
% 输入参数：
%   T_in: 入口温度 (K)
%   P_in: 入口压力 (Pa)
%   m: 燃气质量流量 (kg/s)
%   c_in: 入口速度 (m/s)
%   expansion_ratio: 膨胀比 (入口压力/出口压力)
%   eta_ts: 等熵效率
%   R: 气体常数 (J/kg·K)
%   k: 绝热指数
%   Cp: 定压比热容 (J/kg·K)
%
% 输出结构体 turbine_result 包含：
%   T_out: 出口温度 (K)
%   P_out: 出口压力 (Pa)
%   deltaT: 温降 (K)
%   power: 输出功率 (W)
%   T_in_stag: 入口滞止温度 (K)
%   P_in_stag: 入口滞止压力 (Pa)
%   H_is: 等熵焓降 (J/kg)
%   work_output: 实际输出功 (J/kg)
%   efficiency: 实际效率

% 1. 计算出口压力
P_out = P_in / expansion_ratio;

% 2. 计算等熵膨胀出口温度 (理想过程)
T_out_s = T_in / (expansion_ratio^((k-1)/k));

% 3. 计算实际出口温度
T_out = T_in - eta_ts * (T_in - T_out_s);

% 4. 计算温降
deltaT = T_in - T_out;

% 5. 计算等熵焓降
H_is = Cp * (T_in - T_out_s);

% 6. 计算实际输出功
work_output = eta_ts * H_is;

% 7. 计算输出功率
power = m * work_output;

% 8. 滞止参数计算
T_in_stag = T_in + c_in^2/(2*Cp);
P_in_stag = P_in * (T_in_stag/T_in)^(k/(k-1));

% 9. 计算实际效率
efficiency = work_output / (Cp * T_in_stag);

% 将结果保存到结构体
turbine_result = struct(...
    'T_in', T_in, ...
    'P_in', P_in, ...
    'm', m, ...
    'c_in', c_in, ...
    'expansion_ratio', expansion_ratio, ...
    'eta_ts', eta_ts, ...
    'R', R, ...
    'k', k, ...
    'Cp', Cp, ...
    'T_out', T_out, ...
    'P_out', P_out, ...
    'deltaT', deltaT, ...
    'power', power, ...
    'T_in_stag', T_in_stag, ...
    'P_in_stag', P_in_stag, ...
    'H_is', H_is, ...
    'work_output', work_output, ...
    'efficiency', efficiency);

% 显示主要计算结果
fprintf('单级轴流式透平机设计结果:\n');
fprintf('出口温度: %.2f K\n', T_out);
fprintf('出口压力: %.2f Pa\n', P_out);
fprintf('温降: %.2f K\n', deltaT);
fprintf('输出功率: %.2f W\n', power);
fprintf('等熵焓降: %.2f J/kg\n', H_is);
fprintf('实际效率: %.2f%%\n', efficiency*100);

end