function system_result = multi_stage_turbine_system(...
    N_stage, ...                     % 级数
    T_in1, P_in1, m, c_in1, ...      % 第一级透平机入口参数
    expansion_ratio_list, eta_ts_list, ...   % 各级透平机膨胀比和效率
    T_hot_in_list, m_dot_hot_list, cp_hot_list, P_hot_list, ... % 各级加热器热流体参数
    U_list, A_list, HX_type_list, ... % 各级加热器参数
    R, k, Cp, ...                    % 气体常数
    pressure_drop_ratio)             % 级间压力损失比例

%% 参数验证
if length(expansion_ratio_list) ~= N_stage || length(eta_ts_list) ~= N_stage || ...
   length(T_hot_in_list) ~= N_stage || length(m_dot_hot_list) ~= N_stage || ...
   length(cp_hot_list) ~= N_stage || length(P_hot_list) ~= N_stage || ...
   length(U_list) ~= N_stage || length(A_list) ~= N_stage || ...
   length(HX_type_list) ~= N_stage
    error('所有级数相关参数列表长度必须等于N_stage');
end

%% 初始化变量
turbine_results = cell(1, N_stage);    % 存储各级透平机结果
hx_T_cold_out = zeros(1, N_stage);     % 各级加热器冷流体（工质）出口温度
hx_T_hot_out = zeros(1, N_stage);      % 各级加热器热流体出口温度
hx_Q = zeros(1, N_stage);              % 各级加热量
hx_effectiveness = zeros(1, N_stage);  % 各级加热器效能
total_power = 0;                       % 总输出功率
total_heat_input = 0;                  % 总加热量
turbine_efficiency = zeros(1, N_stage); % 透平机效率
interstage_pressure_drop = zeros(1, N_stage-1); % 级间压力损失

% 第一级透平机入口参数
T_in = T_in1;
P_in = P_in1;
c_in = c_in1;

%% 多级系统计算
for i = 1:N_stage
    fprintf('\n--- 第 %d 级计算开始 ---\n', i);
    
    % 1. 级间加热器计算（除第一级外）
    if i > 1
        fprintf('加热器%d计算: 工质入口=%.2fK, 热源入口=%.2fK\n', i, T_in, T_hot_in_list(i));
        
        % 调用换热器函数，工质作为冷流体被加热
        [T_hot_out, T_cold_out, Q, effectiveness] = GasLiquidHXDesign(...
            m_dot_hot_list(i), cp_hot_list(i), T_hot_in_list(i), P_hot_list(i), ... % 热流体参数
            m, Cp, T_in, P_in, ...          % 冷流体（工质）参数
            U_list(i), A_list(i), HX_type_list{i});
        
        % 保存加热器结果
        hx_T_hot_out(i) = T_hot_out;
        hx_T_cold_out(i) = T_cold_out;
        hx_Q(i) = Q;
        hx_effectiveness(i) = effectiveness;
        
        total_heat_input = total_heat_input + Q;
        T_in = T_cold_out;  % 加热后工质温度作为透平机入口温度
        
        fprintf('加热器%d输出: 工质出口=%.2fK, 热源出口=%.2fK, 加热量=%.2fkW\n', ...
            i, T_cold_out, T_hot_out, Q/1000);
    else
        % 第一级不使用加热器，直接使用入口参数
        hx_T_cold_out(1) = T_in1;
        hx_T_hot_out(1) = T_hot_in_list(1);
        hx_Q(1) = 0;
        hx_effectiveness(1) = 0;
    end
    
    % 2. 透平机计算
    fprintf('透平机%d计算: 入口温度=%.2fK, 入口压力=%.2fPa\n', i, T_in, P_in);
    turbine_results{i} = axial_turbine_design(...
        T_in, P_in, m, c_in, ...
        expansion_ratio_list(i), eta_ts_list(i), R, k, Cp);
    
    total_power = total_power + turbine_results{i}.power;
    turbine_efficiency(i) = eta_ts_list(i);
    
    fprintf('透平机%d输出: 出口温度=%.2fK, 出口压力=%.2fPa, 输出功率=%.2fkW\n', ...
        i, turbine_results{i}.T_out, turbine_results{i}.P_out, turbine_results{i}.power/1000);
    
    % 3. 考虑级间压力损失并为下一级准备参数
    if i < N_stage
        % 考虑加热器中的压力损失
        P_in = turbine_results{i}.P_out * (1 - pressure_drop_ratio);
        interstage_pressure_drop(i) = turbine_results{i}.P_out * pressure_drop_ratio;
        T_in = turbine_results{i}.T_out;  % 下一级加热器入口温度
        
        fprintf('级间压力损失: %.2fPa\n', interstage_pressure_drop(i));
        fprintf('下一级加热器入口: 温度=%.2fK, 压力=%.2fPa\n', T_in, P_in);
    end
end

%% 计算系统性能指标
% 总膨胀比
total_expansion_ratio = P_in1 / turbine_results{end}.P_out;

% 总温降
total_temperature_drop = T_in1 - turbine_results{end}.T_out;

% 透平机平均效率
avg_turbine_efficiency = mean(turbine_efficiency);

% 加热器平均效能
avg_hx_effectiveness = mean(hx_effectiveness(hx_effectiveness > 0));

% 系统热效率
system_thermal_efficiency = total_power / total_heat_input;

%% 整理输出结果
system_result = struct(...
    'total_power', total_power, ...
    'total_heat_input', total_heat_input, ...
    'final_temperature', turbine_results{end}.T_out, ...
    'final_pressure', turbine_results{end}.P_out, ...
    'total_expansion_ratio', total_expansion_ratio, ...
    'total_temperature_drop', total_temperature_drop, ...
    'avg_turbine_efficiency', avg_turbine_efficiency, ...
    'avg_hx_effectiveness', avg_hx_effectiveness, ...
    'system_thermal_efficiency', system_thermal_efficiency, ...
    'number_of_stages', N_stage, ...
    'interstage_pressure_drop', interstage_pressure_drop, ...
    'turbine_results', {turbine_results}, ...
    'heater_results', struct(...
        'T_cold_out', hx_T_cold_out, ...  % 工质出口温度
        'T_hot_out', hx_T_hot_out, ...    % 热源出口温度
        'Q', hx_Q, ...
        'effectiveness', hx_effectiveness), ...
    'pressure_drop_ratio', pressure_drop_ratio, ...
    'turbine_efficiency', turbine_efficiency, ...
    'hx_effectiveness', hx_effectiveness, ...
    'initial_conditions', struct('T_in1', T_in1, 'P_in1', P_in1, 'm', m, 'c_in1', c_in1));

%% 输出系统性能摘要
fprintf('\n========== 多级轴流透平机系统性能摘要 ==========\n');
fprintf('系统级数: %d\n', N_stage);
fprintf('总输出功率: %.2f kW\n', total_power/1000);
fprintf('总加热量: %.2f kW\n', total_heat_input/1000);
fprintf('总膨胀比: %.4f\n', total_expansion_ratio);
fprintf('总温降: %.2f K\n', total_temperature_drop);
fprintf('最终出口温度: %.2f K\n', system_result.final_temperature);
fprintf('最终出口压力: %.2f MPa\n', system_result.final_pressure/1e6);
fprintf('透平机平均效率: %.4f\n', avg_turbine_efficiency);
fprintf('加热器平均效能: %.4f\n', avg_hx_effectiveness);
fprintf('系统热效率: %.4f\n', system_thermal_efficiency);
fprintf('级间压力损失比例: %.2f%%\n', pressure_drop_ratio*100);
fprintf('==================================================\n');

%% 输出各级加热器温度结果
fprintf('\n===== 各级加热器出口温度 =====\n');
for i = 1:N_stage
    if i == 1
        fprintf('第 1 级: 无加热器（直接入口）\n');
    else
        fprintf('第 %d 级加热器: 工质出口温度 = %.2f K, 热源出口温度 = %.2f K\n', ...
            i, hx_T_cold_out(i), hx_T_hot_out(i));
    end
end

%% 生成详细性能表格
fprintf('\n----------- 各级详细性能 -----------\n');
fprintf('级数 | 透平机功率(kW) | 加热量(kW) | 透平机效率 | 加热器效能 | 工质出口温度(K)\n');
fprintf('-----|----------------|------------|------------|------------|-----------------\n');
for i = 1:N_stage
    if i == 1
        fprintf('%3d  | %13.2f  | %10.2f |   %.4f   |     -      |      %.2f\n', ...
            i, turbine_results{i}.power/1000, hx_Q(i)/1000, ...
            turbine_efficiency(i), hx_T_cold_out(i));
    else
        fprintf('%3d  | %13.2f  | %10.2f |   %.4f   |   %.4f   |      %.2f\n', ...
            i, turbine_results{i}.power/1000, hx_Q(i)/1000, ...
            turbine_efficiency(i), hx_effectiveness(i), hx_T_cold_out(i));
    end
end
fprintf('---------------------------------------------------------------------\n');
end