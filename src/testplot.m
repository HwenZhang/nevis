% 创建主 Figure
f = figure('Position',[100,100,1200,600]);

% 顶层布局：1 行 4 列，等宽
mainTL = tiledlayout(f, 1, 4, ...
    'TileSpacing','compact', 'Padding','compact');

%% 左侧子布局：6 行 1 列，放在第 1 列
leftTL = tiledlayout(mainTL, 6, 1, ...
    'TileSpacing','compact', 'Padding','compact');
leftTL.Layout.Tile = 1;   % 指定放在主布局的第 1 格

for i = 1:6
    ax = nexttile(leftTL);
    title(ax, sprintf('Left %d', i));
    plot(ax, rand(1,10));   % 示例内容
end

%% 右侧子布局：3 行 2 列，跨占主布局的第 2～4 列
rightTL = tiledlayout(mainTL, 3, 2, ...
    'TileSpacing','compact', 'Padding','compact');
rightTL.Layout.Tile     = 2;      % 从主布局的第 2 格开始
rightTL.Layout.TileSpan = [1, 3]; % 向下跨 1 行、向右跨 3 列

for i = 1:6
    ax = nexttile(rightTL);
    title(ax, sprintf('Right %d', i));
    imagesc(ax, rand(5)); 
    colorbar(ax);             % 示例内容
end
