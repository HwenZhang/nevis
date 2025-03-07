% 创建一个 6 行 2 列的布局
t = tiledlayout(6,2,'TileSpacing','compact','Padding','compact');

% 左侧一整列 (合并6行1列为一个大图)
ax_left = nexttile([6 1]);  % 左侧图合并6个panel（6行1列）
plot(ax_left, rand(10,1), 'LineWidth',2);
title(ax_left, '左侧大图');

% 右侧6个小图
for i = 1:6
    ax_right = nexttile(); % 默认占据1行1列
    plot(ax_right, rand(10,1));
    title(ax_right, ['右侧子图 ',num2str(i)]);
end

% 可选：添加总标题
title(t, '示例：左侧合并为大图，右侧6个子图');
