# GPS Station 位置数据与时间序列提取指南

## 1. GPS Station 位置数据地点

### 位置数据文件
- **`station_names.mat`** - Station名称数据
  - 包含所有GPS station的标记名称（如SQ11, SQ12, MHIH, MLOW等）
  - 在脚本中加载后转换为大写：`station_names = upper(station_names);`

- **`polarstereo_stations_2022_short.mat`** - Station坐标数据（**关键文件**）
  - Polar Stereographic投影系统中的站点位置
  - 脚本加载：`load polarstereo_stations_2022_short.mat`
  - 这是最主要的位置数据来源

### 坐标系统信息
脚本中的坐标参考系统：
```matlab
% XY相对于M1 moulin (Steven et al. 2015)
origin = [68.72, -49.53]; % M1 moulin (lat, lon)
radius = 6378137.0;        % 地球半径 (m)
eccen = 0.08181919;        % 离心率
lat_true = 70;             % 真实纬度
lon_posy = -45;            % 参考经度

% 使用polarstereo_fwd进行坐标转换
[moulin_x, moulin_y] = polarstereo_fwd(origin(1), origin(2), radius, eccen, lat_true, lon_posy);
moulin_x_km = moulin_x./1e3;   % 转换为km
moulin_y_km = moulin_y./1e3;   % 转换为km
```

### Station分组分类
根据脚本中的数据结构，21个station分为4个地理区域：

| 地区代码 | Station索引 | Station名称 | 用途 |
|---------|----------|-----------|------|
| 100s | 4-9 | SQ11-16 | Within-basin (950s区域) |
| 200s | 10-16 | SQ21-27 | Within-basin (1150s区域) |
| 300s | 17-21 | SQ31-37 | Within-basin (1350s区域) |
| Tiepoints | 1-3 | MHIH, MLOW, QIET | Out-of-basin (1050s-1100s区域) |

## 2. 时间序列数据提取方法

### 2.1 主要数据文件

#### paperfig5脚本使用的文件（DOY 150-254）：
```matlab
% 垂直速度与床开放 (c_dot) 时间序列
load('daily_epsilon_zz_2022R_BF2_30min_BF2_UP4_sZERO_clean_w36_t12_260119.mat')  % loose约束, 36小时窗口
load('daily_epsilon_zz_2022R_BF2_30min_BF2_UP4_sZERO_clean_w18_t6_260119.mat')   % tight约束, 18小时窗口

% 应变率 (strain rates) 时间序列
load('daily_strain_rates_2022R_30min_BF2_UP4_sZERO_clean_w36_t12_260119.mat')    % loose约束
load('daily_strain_rates_2022R_30min_BF2_UP4_sZERO_clean_w18_t6_260119.mat')     % tight约束
```

#### suppfig01a/b/c脚本使用的文件（DOY 150-340，整个季节）：
```matlab
% 垂直速度、水平速度、床开放时间序列 (30分钟分辨率)
load('daily_epsilon_zz_2022R_BF2_30min_150_340_w36_t18_260318.mat')  % loose约束，36小时窗口
load('daily_epsilon_zz_2022R_BF2_30min_150_340_w18_t6_260318.mat')   % tight约束，18小时窗口
```

**脚本对应关系：**
- `suppfig01a_horizontal_vels_runoff_2022_260318.m` - 水平速度 $u_s$ vs 流出量
- `suppfig01b_vertical_vels_runoff_2022_260318.m` - 垂直速度 $w_s$ vs 流出量
- `suppfig01c_separation_runoff_2022_260318.m` - 床开放量 $h$ vs 流出量

### 2.2 时间序列数据结构与提取

#### daily_epsilon_zz结构（速度与变形数据）
```matlab
daily_epsilon_zz(i).t22                         % 时间向量 (DOY, day of year)
daily_epsilon_zz(i).c_dot_delta_t               % 日均垂直速度变化
daily_epsilon_zz(i).c_dot_delta_t_cumulative    % 累积床开放量 (m)
daily_epsilon_zz(i).u_s                         % 水平速度 (m/yr)
daily_epsilon_zz(i).w_s                         % 垂直速度 (m/yr) [正值=上升]
daily_epsilon_zz(i).u_s_combo                   % 融合loose/tight的水平速度
daily_epsilon_zz(i).w_s_combo                   % 融合loose/tight的垂直速度
daily_epsilon_zz(i).c_dot_delta_t_combo         % 融合loose/tight的垂直速度变化
daily_epsilon_zz(i).c_dot_delta_t_cumulative_combo  % 融合loose/tight的累积床开放
daily_epsilon_zz(i).epsilon_dot_lon             % 纵向应变率 (1/yr)
daily_epsilon_zz(i).epsilon_dot_trans           % 横向应变率 (1/yr)
daily_epsilon_zz(i).epsilon_dot_zz              % 垂直应变率 (1/yr)
```

**提取方法（单个station）：**
```matlab
% 对于station i (1-21)
time = daily_epsilon_zz(i).t22;

% 水平速度相关
horiz_vel = daily_epsilon_zz(i).u_s;           % 原始36小时窗口
horiz_vel_tight = daily_epsilon_zz(i).u_s_combo;  % 融合loose/tight的版本

% 垂直速度相关
vert_vel = daily_epsilon_zz(i).w_s;            % 原始36小时窗口
vert_vel_tight = daily_epsilon_zz(i).w_s_combo;   % 融合loose/tight的版本

% 床开放相关
bed_opening = daily_epsilon_zz(i).c_dot_delta_t_cumulative;  % 累积值

% 应变率相关
lon_strain = daily_epsilon_zz(i).epsilon_dot_lon;   % 纵向应变率
trans_strain = daily_epsilon_zz(i).epsilon_dot_trans; % 横向应变率

% 示例：SQ13 (索引6) 的完整轨迹
time_SQ13 = daily_epsilon_zz(6).t22;
u_s_SQ13 = daily_epsilon_zz(6).u_s_combo;      % 水平速度
w_s_SQ13 = daily_epsilon_zz(6).w_s_combo;      % 垂直速度
h_SQ13 = daily_epsilon_zz(6).c_dot_delta_t_cumulative_combo;  % 床开放
```

#### daily_strain_rates结构（应变率时间序列）
```matlab
daily_strain_rates.time                    % 时间向量
daily_strain_rates.station_names           % Station名称
daily_strain_rates.lon_yr_100s(i,j,:)     % 100s区域station i-j的纵向应变率
daily_strain_rates.lon_yr_200s(i,j,:)     % 200s区域station i-j的纵向应变率
daily_strain_rates.lon_yr_300s(i,j,:)     % 300s区域station i-j的纵向应变率
daily_strain_rates.delta_lon_yr_100s(i,j,:)  % 对应的误差估计
daily_strain_rates.delta_lon_yr_200s(i,j,:)  % 对应的误差估计
daily_strain_rates.delta_lon_yr_300s(i,j,:)  % 对应的误差估计
```

**提取方法（station对的应变率）：**
```matlab
% 获取时间向量
strain_time = daily_strain_rates_combo.time;

% 提取特定station对的纵向应变率及误差
% 示例：100s区域的SQ11-12对 (索引4,5)
strain_SQ11_12 = squeeze(daily_strain_rates_combo.lon_yr_100s(4,5,:));
error_SQ11_12 = squeeze(daily_strain_rates_combo.delta_lon_yr_100s(4,5,:));

% 基准化处理（从特定DOY开始计算偏差）
ID_t0 = find(strain_time >= 165.01, 1, 'first');  % 找参考时间索引
strain_normalized = strain_SQ11_12 - strain_SQ11_12(ID_t0);

% 用于绘图的误差包络
tt = [strain_time; flipud(strain_time)];
ee = [strain_SQ11_12 + 3.*error_SQ11_12;
      flipud(strain_SQ11_12 - 3.*error_SQ11_12)];
```

### 2.3 数据融合策略（Loose vs Tight约束）

两个脚本都使用了"松散"和"紧密"约束的组合策略：

```matlab
% Loose约束：36小时滑动窗口，需要18小时数据
daily_epsilon_zz_loose = ...  % 更平滑，长期趋势

% Tight约束：18小时滑动窗口，需要6小时数据
daily_epsilon_zz_tight = ...  % 更灵敏，捕捉短期变化

% 在冰湖排放期间使用tight约束（时间更高分辨率）
% 在其他时期使用loose约束（更稳定）
tight_100s = 187:196;           % DOY 187-196使用tight
tight_200s = 206:217;           % DOY 206-217使用tight
tight_300s = 206:212;           % DOY 206-212使用tight
tight_MLOW = 207:217;           % MLOW station的tight期间
tight_MHIH = 207:217;           % MHIH station的tight期间
```

### 2.4 关键变量时间范围
- **paperfig5脚本**: DOY 150-254 (5月30日 - 9月10日，仅夏季融雪期)
- **suppfig01a脚本** (水平速度): DOY 150-340 (5月30日 - 12月6日，整个季节)
- **suppfig01b脚本** (垂直速度): DOY 150-340 (5月30日 - 12月6日，整个季节)
- **suppfig01c脚本** (床开放): DOY 150-340 (5月30日 - 12月6日，整个季节)
- **RACMO时间向量**: 1.5:1:334.5 (全年DOY)

**数据质量注意：**
- SQ13 (索引6): 数据到DOY 197截止（站点翻倒）
- SQ15 (索引8): 数据到DOY 210截止（站点翻倒）
- SQ14 (索引7): 在DOY 210之后使用SQ15进行横向应变率计算

## 3. 实际数据提取示例

### 示例1：提取单个station的完整速度时间序列（suppfig01a/b风格）
```matlab
% 加载数据
load station_names.mat
load('daily_epsilon_zz_2022R_BF2_30min_150_340_w36_t18_260318.mat')
load('../RACMO/racmo_station_2022_index.mat')

% 获取SQ13 (索引6) 的所有可用参数
station_index = 6;
time = daily_epsilon_zz(station_index).t22;
horiz_vel = daily_epsilon_zz(station_index).u_s_combo;    % 水平速度
vert_vel = daily_epsilon_zz(station_index).w_s_combo;     % 垂直速度
bed_sep = daily_epsilon_zz(station_index).c_dot_delta_t_cumulative_combo;  % 床开放

% 从RACMO获取该位置的流出量
runoff_SQ13 = runoff_2022_nevis(:, ID(6,1)) ./ 10;  % 转换为cm

% 绘制水平/垂直速度与流出量的关系
subplot(2,1,1)
yyaxis left
plot(time, horiz_vel, 'LineWidth', 1.3);
ylabel('Horizontal Velocity $u_s$ [m/yr]');
yyaxis right
bar(racmo_time, runoff_SQ13, 'FaceAlpha', 0.3);
ylabel('Runoff [cm w.e.]');

subplot(2,1,2)
yyaxis left
plot(time, vert_vel, 'LineWidth', 1.3);
ylabel('Vertical Velocity $w_s$ [m/yr]');
yyaxis right
bar(racmo_time, runoff_SQ13, 'FaceAlpha', 0.3);
ylabel('Runoff [cm w.e.]');
```

### 示例2：批量提取所有100s区域stations的水平速度（suppfig01a风格）
```matlab
load('daily_epsilon_zz_2022R_BF2_30min_150_340_w36_t18_260318.mat')
load('../RACMO/racmo_station_2022_index.mat')
load station_names.mat

% 100s区域: 索引4-9 (SQ11-16)
figure; hold on;
for i=4:9
    plot(daily_epsilon_zz(i).t22, daily_epsilon_zz(i).u_s_combo, 'LineWidth', 1.3);
end
legend('SQ11','SQ12','SQ13','SQ14','SQ15','SQ16');
ylabel('Horizontal Velocity $u_s$ [m/yr]');
xlabel('Day of Year, 2022');
grid on;

% 添加RACMO平均流出量
yyaxis right
bar(racmo_time, nanmean(runoff_2022_nevis(:,ID(4:9,1)),2)./10, 'FaceAlpha', 0.25);
ylabel('Runoff at 950s [cm w.e.]');
```

### 示例3：提取单个station的床开放演化（suppfig01c风格）
```matlab
load('daily_epsilon_zz_2022R_BF2_30min_150_340_w36_t18_260318.mat')
load('../RACMO/racmo_station_2022_index.mat')

% 比较所有tiepoint stations (MHIH, MLOW, QIET)的床开放
figure; hold on;
plot(daily_epsilon_zz(1).t22, daily_epsilon_zz(1).c_dot_delta_t_cumulative_combo, 'LineWidth', 1.3);  % MHIH
plot(daily_epsilon_zz(2).t22, daily_epsilon_zz(2).c_dot_delta_t_cumulative_combo, 'LineWidth', 1.3);  % MLOW
plot(daily_epsilon_zz(3).t22, daily_epsilon_zz(3).c_dot_delta_t_cumulative_combo, 'LineWidth', 1.3);  % QIET

legend('MHIH', 'MLOW', 'QIET');
ylabel('Bed Separation, $h$ [m]');
xlabel('Day of Year, 2022');
grid on;

% 添加RACMO平均流出量
yyaxis right
bar(racmo_time, nanmean(runoff_2022_nevis(:,ID(1:3,1)),2)./10, 'FaceAlpha', 0.25);
ylabel('Runoff at 1100s [cm w.e.]');
```

### 示例4：提取应变率以及相应的误差包络（paperfig5风格）
```matlab
load('daily_strain_rates_2022R_30min_BF2_UP4_sZERO_clean_w36_t12_260119.mat')

% 提取100s区域的多个station对
strain_time = daily_strain_rates_2023.time;
ID_t0 = find(strain_time >= 165.01, 1, 'first');  % 找参考时间索引

% SQ11-12对 (索引4,5)
strain_SQ11_12 = squeeze(daily_strain_rates_2023.lon_yr_100s(4,5,:));
error_SQ11_12 = squeeze(daily_strain_rates_2023.delta_lon_yr_100s(4,5,:));
strain_normalized = strain_SQ11_12 - strain_SQ11_12(ID_t0);

% 绘制含误差包络的应变率
figure; hold on;
tt = [strain_time; flipud(strain_time)];
ee_upper = strain_normalized + 3.*error_SQ11_12;
ee_lower = strain_normalized - 3.*error_SQ11_12;
ee = [ee_upper; flipud(ee_lower)];
patch(tt, ee, [0.7 0.7 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
plot(strain_time, strain_normalized, 'k-', 'LineWidth', 1.3);
grid on;
ylabel('Longitudinal Strain Rate (normalized)');
xlabel('Day of Year');
```

## 4. 数据目录结构

所有这些`.mat`文件位置：
```
/Users/hwenzhang/matlabprojects/nevis/data/GNSS_2022/
├── station_names.mat
├── polarstereo_stations_2022_short.mat
├── daily_epsilon_zz_2022R_BF2_30min_*.mat (多个版本)
│   ├── *_BF2_UP4_sZERO_clean_w36_t12_260119.mat  (paperfig5, DOY150-254)
│   ├── *_BF2_UP4_sZERO_clean_w18_t6_260119.mat   (paperfig5, DOY150-254)
│   ├── *_150_340_w36_t18_260318.mat              (suppfig01a/b/c, DOY150-340)
│   └── *_150_340_w18_t6_260318.mat               (suppfig01a/b/c, DOY150-340)
├── daily_strain_rates_2022R_30min_*.mat (多个版本)
│   ├── *_BF2_UP4_sZERO_clean_w36_t12_260119.mat  (paperfig5, DOY150-254)
│   └── *_BF2_UP4_sZERO_clean_w18_t6_260119.mat   (paperfig5, DOY150-254)
├── paperfig5_2022_distillations_260116.m
├── suppfig01a_horizontal_vels_runoff_2022_260318.m  ★ 水平速度
├── suppfig01b_vertical_vels_runoff_2022_260318.m    ★ 垂直速度
├── suppfig01c_separation_runoff_2022_260318.m       ★ 床开放
└── ... (其他支持文件)
```

## 5. 关键说明

- **Station索引 (i)**: 1-21，对应不同的GPS接收器
- **Station对索引 (i,j)**: 用于计算相邻station间的应变率
- **DOY (Day of Year)**: 2022年的儒略日(1-365)
- **时间分辨率**: 30分钟（0.5小时）
- **坐标系**: Polar Stereographic (极地立体投影)
- **参考原点**: M1 moulin at [68.72°N, 49.53°W]

## 6. 脚本功能汇总

| 脚本名称 | 数据类型 | 时间范围 | 主要变量 | 关键特点 |
|---------|---------|--------|---------|--------|
| **paperfig5_2022_distillations_260116.m** | 应变率 | DOY150-254 | $\dot{\epsilon}_{lon}$ | 纵向应变率 + 误差包络 |
| **suppfig01a_horizontal_vels_runoff_2022_260318.m** | 水平速度 | DOY150-340 | $u_s$ | 水平速度与流出量关系 |
| **suppfig01b_vertical_vels_runoff_2022_260318.m** | 垂直速度 | DOY150-340 | $w_s$ | 垂直速度与流出量关系 |
| **suppfig01c_separation_runoff_2022_260318.m** | 床开放 | DOY150-340 | $h$ | 累积床开放与流出量关系 |

## 7. Loose vs Tight约束的应用

脚本中的约束切换逻辑：

```matlab
% Loose约束：36小时滑动窗口（平滑，稳定）
w36_loop = 36;      % 小时
t18_data = 18;      % 需要的最少数据小时数

% Tight约束：18小时滑动窗口（灵敏，捕捉快速变化）
w18_loop = 18;      % 小时
t6_data = 6;        % 需要的最少数据小时数

% 应用规则：
% - 平时使用loose (平滑趋势)
% - 冰湖排放期间 (DOY 187-217) 使用tight (高分辨率)
% - 这样结合了两者的优点
```

**各区域的tight约束期间：**
- **100s** (4-9): DOY 187-196
- **200s** (10-16): DOY 206-217
- **300s** (17-21): DOY 206-212
- **MLOW** (2): DOY 207-217
- **MHIH** (1): DOY 207-217
- **QIET** (3): DOY 207-217

## 8. 参考坐标如何转化为 nevis 坐标（基于 MEaSUREs 脚本）

下面是 `measures_for_nevis_BMv5.m` 里用到的完整坐标转换链路。

### 8.1 参考系设置

脚本使用 Greenland 常见的 Polar Stereographic 参数：

```matlab
radius = 6378137.0;   % a, WGS84 半径 (m)
eccen = 0.08181919;   % e, WGS84 离心率
lat_true = 70;        % 标准纬线 (deg)
lon_posy = -45;       % +Y 对应经线 (deg)
```

用 `polarstereo_fwd` 将经纬度转到投影平面坐标 $(x, y)$（单位 m）。

### 8.2 以 M1 moulin 作为 nevis 原点

核心思想是：

1. 先把 M1 的经纬度投影为 $(x_{M1}, y_{M1})$。
2. 任意点投影得到 $(x, y)$ 后，做平移：

$$
x_{nevis} = x - x_{M1}, \qquad y_{nevis} = y - y_{M1}
$$

在脚本中对应：

```matlab
[moulin_x,moulin_y] = polarstereo_fwd(moulin_lat,moulin_lon,radius,eccen,lat_true,lon_posy);
X_rel = X - moulin_x;
Y_rel = Y - moulin_y;
```

所以 nevis 坐标系本质上是“同一投影参数下的局部平移坐标系”，其原点定义为 M1。

### 8.3 Lake 与 BedMachine 数据如何进入 nevis 坐标

- `lake.mat` 的经纬度先通过 `polarstereo_fwd` 投影，再减去 M1 的投影坐标。
- `BM5_lake6.mat` 的 `X, Y` 已是投影平面坐标，直接减去 `moulin_x, moulin_y` 得到 `X_rel, Y_rel`。

这一步之后，所有地形量（`Bed`, `Surface`, `Bederr`）都和 `X_rel, Y_rel` 对齐，可直接作为 nevis 几何输入。

### 8.4 MEaSUREs 速度如何映射到 nevis 网格

在脚本里，`mosaicV1.xx_moulin`、`mosaicV1.yy_moulin` 已经是以 moulin 为原点的相对坐标（与 nevis 坐标同一物理定义）。

接着通过插值映射到 nevis 模型网格节点：

```matlab
F  = scatteredInterpolant(X_rel_mosaic, Y_rel_mosaic, double(mosaic_speed));
Fx = scatteredInterpolant(X_rel_mosaic, Y_rel_mosaic, double(mosaic_vel_x));
Fy = scatteredInterpolant(X_rel_mosaic, Y_rel_mosaic, double(mosaic_vel_y));

mosaic_nevis_noSK   = F(gg.nx.*ps.x, gg.ny.*ps.x);
mosaic_nevis_noSK_x = Fx(gg.nx.*ps.x, gg.ny.*ps.x);
mosaic_nevis_noSK_y = Fy(gg.nx.*ps.x, gg.ny.*ps.x);
```

其中 `gg.nx, gg.ny` 是无量纲网格坐标，乘 `ps.x`（默认 10000）还原为米后参与插值。

### 8.5 nevis 坐标反算回经纬度

若你有某点的 nevis 坐标 $(x_{nevis}, y_{nevis})$，先加回 M1 平移，再用 `polarstereo_inv`：

$$
x = x_{nevis} + x_{M1}, \qquad y = y_{nevis} + y_{M1}
$$

```matlab
x_ps = x_nevis + moulin_x;
y_ps = y_nevis + moulin_y;
[lat, lon] = polarstereo_inv(x_ps, y_ps, radius, eccen, lat_true, lon_posy);
```

### 8.6 使用时的两个注意点

1. **单位一致性**：脚本里大部分是米（m），绘图时常除以 `1e3` 转 km；字段名有时写 `X_km` 但变量实际可能仍是 m，使用前建议核对数值量级。
2. **参数必须一致**：`a, e, lat_true, lon_posy` 在正向和反向转换中必须完全一致，否则会出现系统性位置偏差。
