% profile off           % 关闭并清除之前的 profiling 记录
% profile clear         % 清空所有旧的 profile 数据
% profile on -timer real % 开启 profiler，使用真实时间（wall‑clock）

% n2d_100m3s_kappa1e_10_mu1e1_hbreg5e_3_spinup

% profile off           % 运行结束后关闭 profiler
% profile viewer        % 打开交互式 Profiler 界面
% profsave(profile('info'), 'profiler_results');

disp('Hello from VSCode');
a = 1;
b = 2;
c = a + b;

pause(2);   % 延迟两秒便于观察