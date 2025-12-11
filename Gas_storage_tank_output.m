function dX = Gas_storage_tank_output(t, X, param)
% --------------------------------------------------------
% X(1) = T  —— 罐内温度 (K)
% X(2) = m  —— 罐内气体质量 (kg)
%
% param 结构体：
%   param.R       —— 气体常数 (kJ/kg·K 或一致单位)
%   param.V       —— 储气罐体积
%   param.h       —— 换热系数
%   param.A       —— 换热面积
%   param.Tamb    —— 环境温度 (K)
%   param.mdot    —— 出口质量流量 (kg/s)
%
% 输出 dX = [dT/dt; dm/dt]
% --------------------------------------------------------

T = X(1);
m = X(2);

% --------- cp(T) 多项式 ----------------
cp = 0.9705 + 0.00006791*T + 1.658e-7*T^2 - 6.788e-11*T^3;

% --------- cv = cp - R -----------------
cv = cp - param.R;

% --------- 能量方程 --------------------
Q_out = param.mdot * cp * T;                  % 出流带走的能量
Q_loss = param.h * param.A * (T - param.Tamb);% 与环境换热
dTdt = -( Q_out + Q_loss ) / ( m * cv );

% --------- 质量方程 ---------------------
dmdt = -param.mdot;

dX = [dTdt; dmdt];
end
