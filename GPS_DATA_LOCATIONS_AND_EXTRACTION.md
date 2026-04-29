# GPS Station Location Data and Time Series Extraction Guide

## 1. GPS Station Location Data

### Location Data Files

- **`station_names.mat`** - Station name data
  - Contains the marker names for all GPS stations, such as SQ11, SQ12, MHIH, and MLOW.
  - After loading, the script converts names to uppercase: `station_names = upper(station_names);`

- **`polarstereo_stations_2022_short.mat`** - Station coordinate data (**key file**)
  - Contains station locations in the Polar Stereographic projection.
  - Loaded in scripts with: `load polarstereo_stations_2022_short.mat`
  - This is the main source of station location data.

### Coordinate System Information

The coordinate reference system used in the scripts is:

```matlab
% XY relative to the M1 moulin (Stevens et al. 2015)
origin = [68.72, -49.53]; % M1 moulin (lat, lon)
radius = 6378137.0;        % Earth radius (m)
eccen = 0.08181919;        % Eccentricity
lat_true = 70;             % Latitude of true scale
lon_posy = -45;            % Reference longitude

% Convert coordinates with polarstereo_fwd
[moulin_x, moulin_y] = polarstereo_fwd(origin(1), origin(2), radius, eccen, lat_true, lon_posy);
moulin_x_km = moulin_x./1e3;   % Convert to km
moulin_y_km = moulin_y./1e3;   % Convert to km
```

### Station Groups

Based on the script data structure, the 21 stations are divided into four geographic groups:

| Region Code | Station Indices | Station Names | Purpose |
| --- | --- | --- | --- |
| 100s | 4-9 | SQ11-16 | Within-basin, 950s region |
| 200s | 10-16 | SQ21-27 | Within-basin, 1150s region |
| 300s | 17-21 | SQ31-37 | Within-basin, 1350s region |
| Tiepoints | 1-3 | MHIH, MLOW, QIET | Out-of-basin, 1050s-1100s region |

## 2. Time Series Data Extraction

### 2.1 Main Data Files

#### Files used by the `paperfig5` script, DOY 150-254:

```matlab
% Vertical velocity and bed opening (c_dot) time series
load('daily_epsilon_zz_2022R_BF2_30min_BF2_UP4_sZERO_clean_w36_t12_260119.mat')  % loose constraint, 36-hour window
load('daily_epsilon_zz_2022R_BF2_30min_BF2_UP4_sZERO_clean_w18_t6_260119.mat')   % tight constraint, 18-hour window

% Strain-rate time series
load('daily_strain_rates_2022R_30min_BF2_UP4_sZERO_clean_w36_t12_260119.mat')    % loose constraint
load('daily_strain_rates_2022R_30min_BF2_UP4_sZERO_clean_w18_t6_260119.mat')     % tight constraint
```

#### Files used by the `suppfig01a/b/c` scripts, DOY 150-340, full season:

```matlab
% Vertical velocity, horizontal velocity, and bed opening time series at 30-minute resolution
load('daily_epsilon_zz_2022R_BF2_30min_150_340_w36_t18_260318.mat')  % loose constraint, 36-hour window
load('daily_epsilon_zz_2022R_BF2_30min_150_340_w18_t6_260318.mat')   % tight constraint, 18-hour window
```

**Script mapping:**

- `suppfig01a_horizontal_vels_runoff_2022_260318.m` - Horizontal velocity $u_s$ vs runoff
- `suppfig01b_vertical_vels_runoff_2022_260318.m` - Vertical velocity $w_s$ vs runoff
- `suppfig01c_separation_runoff_2022_260318.m` - Bed opening $h$ vs runoff

### 2.2 Time Series Data Structures and Extraction

#### `daily_epsilon_zz` structure, velocity and deformation data

```matlab
daily_epsilon_zz(i).t22                         % Time vector (DOY, day of year)
daily_epsilon_zz(i).c_dot_delta_t               % Daily mean vertical velocity change
daily_epsilon_zz(i).c_dot_delta_t_cumulative    % Cumulative bed opening (m)
daily_epsilon_zz(i).u_s                         % Horizontal velocity (m/yr)
daily_epsilon_zz(i).w_s                         % Vertical velocity (m/yr), positive = uplift
daily_epsilon_zz(i).u_s_combo                   % Horizontal velocity merged from loose/tight constraints
daily_epsilon_zz(i).w_s_combo                   % Vertical velocity merged from loose/tight constraints
daily_epsilon_zz(i).c_dot_delta_t_combo         % Vertical velocity change merged from loose/tight constraints
daily_epsilon_zz(i).c_dot_delta_t_cumulative_combo  % Cumulative bed opening merged from loose/tight constraints
daily_epsilon_zz(i).epsilon_dot_lon             % Longitudinal strain rate (1/yr)
daily_epsilon_zz(i).epsilon_dot_trans           % Transverse strain rate (1/yr)
daily_epsilon_zz(i).epsilon_dot_zz              % Vertical strain rate (1/yr)
```

**Extraction method for one station:**

```matlab
% For station i (1-21)
time = daily_epsilon_zz(i).t22;

% Horizontal velocity
horiz_vel = daily_epsilon_zz(i).u_s;              % Raw 36-hour window
horiz_vel_tight = daily_epsilon_zz(i).u_s_combo;  % Merged loose/tight version

% Vertical velocity
vert_vel = daily_epsilon_zz(i).w_s;               % Raw 36-hour window
vert_vel_tight = daily_epsilon_zz(i).w_s_combo;   % Merged loose/tight version

% Bed opening
bed_opening = daily_epsilon_zz(i).c_dot_delta_t_cumulative;  % Cumulative value

% Strain rates
lon_strain = daily_epsilon_zz(i).epsilon_dot_lon;       % Longitudinal strain rate
trans_strain = daily_epsilon_zz(i).epsilon_dot_trans;   % Transverse strain rate

% Example: full trajectory for SQ13 (index 6)
time_SQ13 = daily_epsilon_zz(6).t22;
u_s_SQ13 = daily_epsilon_zz(6).u_s_combo;      % Horizontal velocity
w_s_SQ13 = daily_epsilon_zz(6).w_s_combo;      % Vertical velocity
h_SQ13 = daily_epsilon_zz(6).c_dot_delta_t_cumulative_combo;  % Bed opening
```

#### `daily_strain_rates` structure, strain-rate time series

```matlab
daily_strain_rates.time                    % Time vector
daily_strain_rates.station_names           % Station names
daily_strain_rates.lon_yr_100s(i,j,:)      % Longitudinal strain rate for station pair i-j in the 100s region
daily_strain_rates.lon_yr_200s(i,j,:)      % Longitudinal strain rate for station pair i-j in the 200s region
daily_strain_rates.lon_yr_300s(i,j,:)      % Longitudinal strain rate for station pair i-j in the 300s region
daily_strain_rates.delta_lon_yr_100s(i,j,:)  % Corresponding uncertainty estimate
daily_strain_rates.delta_lon_yr_200s(i,j,:)  % Corresponding uncertainty estimate
daily_strain_rates.delta_lon_yr_300s(i,j,:)  % Corresponding uncertainty estimate
```

**Extraction method for strain rate between a station pair:**

```matlab
% Get the time vector
strain_time = daily_strain_rates_combo.time;

% Extract longitudinal strain rate and uncertainty for a specific station pair
% Example: SQ11-12 pair in the 100s region (indices 4,5)
strain_SQ11_12 = squeeze(daily_strain_rates_combo.lon_yr_100s(4,5,:));
error_SQ11_12 = squeeze(daily_strain_rates_combo.delta_lon_yr_100s(4,5,:));

% Normalize relative to a reference DOY
ID_t0 = find(strain_time >= 165.01, 1, 'first');  % Reference-time index
strain_normalized = strain_SQ11_12 - strain_SQ11_12(ID_t0);

% Uncertainty envelope for plotting
tt = [strain_time; flipud(strain_time)];
ee = [strain_SQ11_12 + 3.*error_SQ11_12;
      flipud(strain_SQ11_12 - 3.*error_SQ11_12)];
```

### 2.3 Loose vs Tight Constraint Merging Strategy

Both scripts use a merged strategy combining loose and tight constraints:

```matlab
% Loose constraint: 36-hour sliding window, requiring 18 hours of data
daily_epsilon_zz_loose = ...  % Smoother, better for long-term trends

% Tight constraint: 18-hour sliding window, requiring 6 hours of data
daily_epsilon_zz_tight = ...  % More responsive, better for short-term changes

% Use the tight constraint during lake drainage events for higher time resolution.
% Use the loose constraint outside those intervals for stability.
tight_100s = 187:196;           % Use tight constraint for DOY 187-196
tight_200s = 206:217;           % Use tight constraint for DOY 206-217
tight_300s = 206:212;           % Use tight constraint for DOY 206-212
tight_MLOW = 207:217;           % Tight interval for MLOW
tight_MHIH = 207:217;           % Tight interval for MHIH
```

### 2.4 Key Variable Time Ranges

- **`paperfig5` script**: DOY 150-254, May 30 to September 10, summer melt season only
- **`suppfig01a` script**, horizontal velocity: DOY 150-340, May 30 to December 6, full season
- **`suppfig01b` script**, vertical velocity: DOY 150-340, May 30 to December 6, full season
- **`suppfig01c` script**, bed opening: DOY 150-340, May 30 to December 6, full season
- **RACMO time vector**: 1.5:1:334.5, full-year DOY

**Data quality notes:**

- SQ13, index 6: data end at DOY 197 because the station tipped over.
- SQ15, index 8: data end at DOY 210 because the station tipped over.
- SQ14, index 7: after DOY 210, SQ15 is used for the transverse strain-rate calculation.

## 3. Practical Data Extraction Examples

### Example 1: Extract the full velocity time series for a single station, `suppfig01a/b` style

```matlab
% Load data
load station_names.mat
load('daily_epsilon_zz_2022R_BF2_30min_150_340_w36_t18_260318.mat')
load('../RACMO/racmo_station_2022_index.mat')

% Get all available parameters for SQ13 (index 6)
station_index = 6;
time = daily_epsilon_zz(station_index).t22;
horiz_vel = daily_epsilon_zz(station_index).u_s_combo;    % Horizontal velocity
vert_vel = daily_epsilon_zz(station_index).w_s_combo;     % Vertical velocity
bed_sep = daily_epsilon_zz(station_index).c_dot_delta_t_cumulative_combo;  % Bed opening

% Get runoff at this location from RACMO
runoff_SQ13 = runoff_2022_nevis(:, ID(6,1)) ./ 10;  % Convert to cm

% Plot horizontal/vertical velocity against runoff
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

### Example 2: Extract horizontal velocities for all 100s-region stations, `suppfig01a` style

```matlab
load('daily_epsilon_zz_2022R_BF2_30min_150_340_w36_t18_260318.mat')
load('../RACMO/racmo_station_2022_index.mat')
load station_names.mat

% 100s region: indices 4-9 (SQ11-16)
figure; hold on;
for i=4:9
    plot(daily_epsilon_zz(i).t22, daily_epsilon_zz(i).u_s_combo, 'LineWidth', 1.3);
end
legend('SQ11','SQ12','SQ13','SQ14','SQ15','SQ16');
ylabel('Horizontal Velocity $u_s$ [m/yr]');
xlabel('Day of Year, 2022');
grid on;

% Add mean RACMO runoff
yyaxis right
bar(racmo_time, nanmean(runoff_2022_nevis(:,ID(4:9,1)),2)./10, 'FaceAlpha', 0.25);
ylabel('Runoff at 950s [cm w.e.]');
```

### Example 3: Extract bed-opening evolution for individual stations, `suppfig01c` style

```matlab
load('daily_epsilon_zz_2022R_BF2_30min_150_340_w36_t18_260318.mat')
load('../RACMO/racmo_station_2022_index.mat')

% Compare bed opening for all tiepoint stations (MHIH, MLOW, QIET)
figure; hold on;
plot(daily_epsilon_zz(1).t22, daily_epsilon_zz(1).c_dot_delta_t_cumulative_combo, 'LineWidth', 1.3);  % MHIH
plot(daily_epsilon_zz(2).t22, daily_epsilon_zz(2).c_dot_delta_t_cumulative_combo, 'LineWidth', 1.3);  % MLOW
plot(daily_epsilon_zz(3).t22, daily_epsilon_zz(3).c_dot_delta_t_cumulative_combo, 'LineWidth', 1.3);  % QIET

legend('MHIH', 'MLOW', 'QIET');
ylabel('Bed Separation, $h$ [m]');
xlabel('Day of Year, 2022');
grid on;

% Add mean RACMO runoff
yyaxis right
bar(racmo_time, nanmean(runoff_2022_nevis(:,ID(1:3,1)),2)./10, 'FaceAlpha', 0.25);
ylabel('Runoff at 1100s [cm w.e.]');
```

### Example 4: Extract strain rate and the corresponding uncertainty envelope, `paperfig5` style

```matlab
load('daily_strain_rates_2022R_30min_BF2_UP4_sZERO_clean_w36_t12_260119.mat')

% Extract several station pairs in the 100s region
strain_time = daily_strain_rates_2023.time;
ID_t0 = find(strain_time >= 165.01, 1, 'first');  % Reference-time index

% SQ11-12 pair (indices 4,5)
strain_SQ11_12 = squeeze(daily_strain_rates_2023.lon_yr_100s(4,5,:));
error_SQ11_12 = squeeze(daily_strain_rates_2023.delta_lon_yr_100s(4,5,:));
strain_normalized = strain_SQ11_12 - strain_SQ11_12(ID_t0);

% Plot strain rate with an uncertainty envelope
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

## 4. Data Directory Structure

All of these `.mat` files are located under:

```text
/Users/hwenzhang/matlabprojects/nevis/data/GNSS_2022/
|-- station_names.mat
|-- polarstereo_stations_2022_short.mat
|-- daily_epsilon_zz_2022R_BF2_30min_*.mat (multiple versions)
|   |-- *_BF2_UP4_sZERO_clean_w36_t12_260119.mat  (paperfig5, DOY150-254)
|   |-- *_BF2_UP4_sZERO_clean_w18_t6_260119.mat   (paperfig5, DOY150-254)
|   |-- *_150_340_w36_t18_260318.mat              (suppfig01a/b/c, DOY150-340)
|   `-- *_150_340_w18_t6_260318.mat               (suppfig01a/b/c, DOY150-340)
|-- daily_strain_rates_2022R_30min_*.mat (multiple versions)
|   |-- *_BF2_UP4_sZERO_clean_w36_t12_260119.mat  (paperfig5, DOY150-254)
|   `-- *_BF2_UP4_sZERO_clean_w18_t6_260119.mat   (paperfig5, DOY150-254)
|-- paperfig5_2022_distillations_260116.m
|-- suppfig01a_horizontal_vels_runoff_2022_260318.m  * horizontal velocity
|-- suppfig01b_vertical_vels_runoff_2022_260318.m    * vertical velocity
|-- suppfig01c_separation_runoff_2022_260318.m       * bed opening
`-- ... (other support files)
```

## 5. Key Notes

- **Station index (`i`)**: 1-21, corresponding to individual GPS receivers.
- **Station-pair index (`i,j`)**: Used to calculate strain rates between neighboring stations.
- **DOY (Day of Year)**: Julian day of 2022, from 1 to 365.
- **Time resolution**: 30 minutes, or 0.5 hours.
- **Coordinate system**: Polar Stereographic.
- **Reference origin**: M1 moulin at [68.72 deg N, 49.53 deg W].

## 6. Script Function Summary

| Script Name | Data Type | Time Range | Main Variable | Key Feature |
| --- | --- | --- | --- | --- |
| **paperfig5_2022_distillations_260116.m** | Strain rate | DOY150-254 | $\dot{\epsilon}_{lon}$ | Longitudinal strain rate plus uncertainty envelope |
| **suppfig01a_horizontal_vels_runoff_2022_260318.m** | Horizontal velocity | DOY150-340 | $u_s$ | Horizontal velocity versus runoff |
| **suppfig01b_vertical_vels_runoff_2022_260318.m** | Vertical velocity | DOY150-340 | $w_s$ | Vertical velocity versus runoff |
| **suppfig01c_separation_runoff_2022_260318.m** | Bed opening | DOY150-340 | $h$ | Cumulative bed opening versus runoff |

## 7. Loose vs Tight Constraint Usage

The scripts use the following constraint-switching logic:

```matlab
% Loose constraint: 36-hour sliding window (smooth, stable)
w36_loop = 36;      % hours
t18_data = 18;      % minimum required data hours

% Tight constraint: 18-hour sliding window (responsive, captures rapid changes)
w18_loop = 18;      % hours
t6_data = 6;        % minimum required data hours

% Usage rule:
% - Use loose constraints during ordinary periods for smoothed trends.
% - Use tight constraints during lake drainage events (DOY 187-217) for higher resolution.
% - This combines the advantages of both approaches.
```

**Tight-constraint intervals by region:**

- **100s** (4-9): DOY 187-196
- **200s** (10-16): DOY 206-217
- **300s** (17-21): DOY 206-212
- **MLOW** (2): DOY 207-217
- **MHIH** (1): DOY 207-217
- **QIET** (3): DOY 207-217

## 8. Converting Reference Coordinates to nevis Coordinates, Based on the MEaSUREs Script

This section summarizes the full coordinate-conversion chain used in `measures_for_nevis_BMv5.m`.

### 8.1 Reference System Setup

The script uses common Greenland Polar Stereographic parameters:

```matlab
radius = 6378137.0;   % a, WGS84 radius (m)
eccen = 0.08181919;   % e, WGS84 eccentricity
lat_true = 70;        % Standard parallel (deg)
lon_posy = -45;       % Meridian corresponding to +Y (deg)
```

`polarstereo_fwd` converts latitude/longitude to projected-plane coordinates $(x, y)$ in meters.

### 8.2 Using the M1 Moulin as the nevis Origin

The key idea is:

1. Project the M1 latitude/longitude to $(x_{M1}, y_{M1})$.
2. For any projected point $(x, y)$, translate by subtracting the M1 coordinates:

$$
x_{nevis} = x - x_{M1}, \qquad y_{nevis} = y - y_{M1}
$$

In the script:

```matlab
[moulin_x,moulin_y] = polarstereo_fwd(moulin_lat,moulin_lon,radius,eccen,lat_true,lon_posy);
X_rel = X - moulin_x;
Y_rel = Y - moulin_y;
```

Thus, the nevis coordinate system is a locally translated coordinate system using the same projection parameters, with the origin defined at M1.

### 8.3 How Lake and BedMachine Data Enter nevis Coordinates

- Latitude/longitude from `lake.mat` is first projected with `polarstereo_fwd`, then shifted by subtracting the M1 projected coordinates.
- `X, Y` in `BM5_lake6.mat` are already projected-plane coordinates, so `moulin_x, moulin_y` are subtracted directly to obtain `X_rel, Y_rel`.

After this step, all topographic fields (`Bed`, `Surface`, `Bederr`) are aligned with `X_rel, Y_rel` and can be used directly as nevis geometry input.

### 8.4 Mapping MEaSUREs Velocity to the nevis Grid

In the script, `mosaicV1.xx_moulin` and `mosaicV1.yy_moulin` are already relative coordinates with the moulin as the origin, matching the physical definition of nevis coordinates.

The data are then interpolated onto the nevis model grid nodes:

```matlab
F  = scatteredInterpolant(X_rel_mosaic, Y_rel_mosaic, double(mosaic_speed));
Fx = scatteredInterpolant(X_rel_mosaic, Y_rel_mosaic, double(mosaic_vel_x));
Fy = scatteredInterpolant(X_rel_mosaic, Y_rel_mosaic, double(mosaic_vel_y));

mosaic_nevis_noSK   = F(gg.nx.*ps.x, gg.ny.*ps.x);
mosaic_nevis_noSK_x = Fx(gg.nx.*ps.x, gg.ny.*ps.x);
mosaic_nevis_noSK_y = Fy(gg.nx.*ps.x, gg.ny.*ps.x);
```

Here, `gg.nx` and `gg.ny` are nondimensional grid coordinates. Multiplying by `ps.x`, default value 10000, converts them back to meters before interpolation.

### 8.5 Converting nevis Coordinates Back to Latitude/Longitude

If a point has nevis coordinates $(x_{nevis}, y_{nevis})$, first add back the M1 translation and then call `polarstereo_inv`:

$$
x = x_{nevis} + x_{M1}, \qquad y = y_{nevis} + y_{M1}
$$

```matlab
x_ps = x_nevis + moulin_x;
y_ps = y_nevis + moulin_y;
[lat, lon] = polarstereo_inv(x_ps, y_ps, radius, eccen, lat_true, lon_posy);
```

### 8.6 Two Usage Notes

1. **Unit consistency**: Most values in the script are in meters (m). Plotting often divides by `1e3` to convert to km. Some field names include `X_km`, even though the stored values may still be in meters, so check value magnitudes before use.
2. **Projection parameters must match**: `a`, `e`, `lat_true`, and `lon_posy` must be identical in forward and inverse transformations. Otherwise, systematic position offsets will occur.
