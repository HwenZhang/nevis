%==========================================================================
% 完整 MATLAB 代码：离散润滑方程（q = -h^3 dh/dx）数值模拟
% 并在每个时间步可视化 h(x,t) 分布
%==========================================================================
clear; clc;

%% 参数设定
L  = 1.0;        % 空间域长度
Nx = 100;        % 网格点数
dx = L / (Nx-1); % 空间步长 （注意：Nx-1 个区间）
x  = linspace(0, L, Nx)';  % 空间坐标列向量

dt   = 1e-6;     % 时间步长
Tend = 5e-2;     % 末时刻
Nt   = round(Tend / dt);

%% 初始条件：中间一点有非零厚度，其余为0
h = zeros(Nx, 1);
% h(round(Nx/2)) = 0.1;
h = (abs(x-0.5*L)<0.25*L).*(0.01.*(1-(x-L/2).^2/(L/4)^2).^2); % 

%% 预分配通量数组
q = zeros(Nx+1, 1);  % q 在 i+1/2 处有 Nx+1 个值（从 0.5,1.5,...,Nx+0.5）

%% 主循环：显式 Euler 时间推进
% D5 = get_D5_matrix(Nx, dx);  % 获取五阶导数矩阵
radius = zeros(Nt, 1); % 存储每个时间步的半径
time = zeros(Nt, 1);   % 存储每个时间步的时间
figure(1);
for n = 1:Nt
    h_old = h;
    
    % 计算每个 i+1/2 处的通量 q_{i+1/2}
    % 索引 i 从 1 到 Nx-1，对应 q(2) 到 q(Nx)
    % 边界处 q(1) 和 q(Nx+1) 保持为 0（固壁或零通量）
    for i = 1:Nx-1
        hL = h_old(i);
        hR = h_old(i+1);
        h3avg = 0.5 * (max(hL,0)^3 + max(hR,0)^3);

        % d5h_i   = D5(i+1,:) * h_old;
        % d5h_ip1 = D5(i+1,:) * h_old;
        % d5h_avg = 0.5 * (d5h_i + d5h_ip1);

        dhdx = (hR - hL) / dx;
        if i>= 3 && i < Nx-2
            d5hdx = (h_old(i-2) - 5*h_old(i-1) + 10*h_old(i) - 10*h_old(i+1) + 5*h_old(i+2) - h_old(i+3)) / dx^5; % 五阶导数近似
            d5hdx = (h(i-2) - 5*h(i-1) + 10*h(i) - 10*h(i+1) + 5*h(i+2) - h(i+3)) / dx^5; % 五阶导数近似
        else
            d5hdx = 0; % 边界处五阶导数为0
        end
        % q(i+1) = - 1000*h3avg * dhdx;
        q(i+1) = - 1*h3avg * d5hdx;
    end
    
    % 计算通量散度：dq/dx 在格点 i 处
    % dq/dx |_i ≈ (q_{i+1/2} - q_{i-1/2}) / dx
    dhdt = zeros(Nx,1);
    for i = 1:Nx
        dhdt(i) = - (q(i+1) - q(i)) / dx;
    end
    
    % 更新 h
    h = h_old + dt * dhdt;
    h(round(Nx/2)) = h(round(Nx/2)) + 0e4*dt; % 在中间点添加一个小的正项，模拟源项

    % 保证厚度非负
    h(h < 0) = 0;

    % 计算半径
    % [~, maxIdx] = max(h);
    % nonzeroIdx = find(h > 1e-4*h(maxIdx));
    % [~, localIdx] = min(h(nonzeroIdx));
    % minidx = nonzeroIdx(localIdx);
    % radius(n) = abs(x(maxIdx) - x(minidx));
    % time(n)= n * dt;

    % 可视化：绘制 h(x) 曲线
    plot(x, h, 'b-', 'LineWidth', 1.5);
    xlim([0.2*L, 0.8*L]);
    ylim([0, max(0.1, 1.1*max(h_old))]);
    ylim([0, 1.1]);
    title(sprintf('时间 t = %.4f', n*dt), 'FontSize', 14);
    xlabel('x', 'FontSize', 12);
    ylabel('h(x,t)', 'FontSize', 12);
    grid on;
    drawnow;
end

%% 结束
disp('模拟完成。');
figure(2);
plot(time, radius, 'r-', 'LineWidth', 1.5);
xlabel('时间 t', 'FontSize', 12);
ylabel('半径 r', 'FontSize', 12);
grid on;

function D5 = get_D5_matrix(N, dx)
    D5 = zeros(N, N);
    coeffs = [1, -5, 10, -10, 5, -1] / dx^5;

    for i = 4:N-3
        D5(i, i-3:i+2) = coeffs;
    end

    % Zero rows at edges to avoid instability
    D5(1:3, :)     = 0;
    D5(N-2:N, :) = 0;
end

figure(3);
index = radius > 0; % 只保留半径大于0的点
radius = radius(index); % 去掉第一个半径点
time = time(index); % 去掉第一个时间点

plot(time, radius, 'r-', 'LineWidth', 1.5);
xlabel('时间 t', 'FontSize', 12);
ylabel('半径 r', 'FontSize', 12);
grid on;
% Fit a power law to the data
log10t = log10(time);
log10R_b = log10(radius);
p = polyfit(log10t, log10R_b, 1);
% Calculate the fitted values
fitted_Rb = 10.^(polyval(p, log10(time)));
% Plot the fitted line
hold on;
plot(time, fitted_Rb, 'k--', 'LineWidth', 2);
% Display the power law exponent
exponent = p(1);
disp(['Fitted power law exponent: ', num2str(exponent)]);
% 添加图例
legend('Blister radius', 'Fitted power law', 'Location', 'NorthWest');
% 设置图形属性
set(gca, 'FontSize', 14);
% 设置标题
title('气泡半径随时间变化', 'FontSize', 16);