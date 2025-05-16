% 参数
L = 5;
nx = 201;
x = linspace(-L, L, nx);
dx = x(2) - x(1);
dt = 1e-4;
Tmax = 0.2;
nt = round(Tmax/dt);
m = 3;

% 初始条件（近似 Dirac）：高斯尖峰
u = exp(-x.^2/0.01);
u = u / trapz(x, u);  % 归一化质量为1

% 时间推进
for n = 1:nt
    u_old = u;

    % 计算 flux = d/dx(u^m)
    u_pad = [u(1), u, u(end)];
    um = u_pad.^m;
    flux = (um(3:end) - um(1:end-2)) / (2*dx);  % 中心差分

    % 计算 d^2/dx^2(u^m)
    dfluxdx = (flux(2:end) - flux(1:end-1)) / dx;

    % 显式更新
    u = u_old + dt * dfluxdx;

    % 可视化
    if mod(n, 100) == 0
        plot(x, u, 'b-', 'LineWidth', 1.5)
        title(['t = ', num2str(n*dt)])
        xlabel('x'); ylabel('u(x,t)');
        drawnow
    end
end
