from pathlib import Path
import warnings

import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio
from scipy.integrate import cumulative_trapezoid


class NevisGpsClusterComparisonPlotter:
    def __init__(
        self,
        casename,
        repo_root='.',
        clusters_to_plot=None,
        tmin_yr=0.42,
        tmax_yr=0.82,
        right_column_mode='speed',
    ):
        self.casename = casename
        self.repo_root = Path(repo_root).resolve()
        self.results_path = self.repo_root / 'results' / casename / f'{casename}.mat'
        self.station_path = self.repo_root / 'data' / 'station_timeseries_2022.mat'
        self.clusters_to_plot = list(clusters_to_plot) if clusters_to_plot is not None else [1, 2, 3, 4, 5]
        self.tmin_yr = float(tmin_yr)
        self.tmax_yr = float(tmax_yr)
        self.right_column_mode = str(right_column_mode).lower()
        if self.right_column_mode not in {'speed', 'displacement'}:
            raise ValueError("right_column_mode must be 'speed' or 'displacement'.")

        self.cluster_colors = [
            (0.85, 0.33, 0.10),
            (0.13, 0.55, 0.13),
            (0.12, 0.47, 0.71),
            (0.58, 0.00, 0.83),
            (0.72, 0.53, 0.04),
        ]
        self.cluster_labels = ['S1', 'S2', 'S3', 'S4', 'S5']
        self._loaded = False

    def load_data(self):
        raw = sio.loadmat(self.results_path, squeeze_me=True)
        st = sio.loadmat(self.station_path, squeeze_me=True)

        self.ps = raw['ps']
        self.pd_r = raw['pd']
        self.pp = raw['pp']
        self.oo = raw['oo']
        self.gg = raw['gg']
        self.aa = raw['aa']
        self.tt = raw['tt']
        self.stations = st['station_data']

        self.ps_t = float(self.ps['t'])
        self.pd_td = float(self.pd_r['td'])
        self.pd_ty = float(self.pd_r['ty'])
        self.ps_h = float(self.ps['h'])
        self.ps_hb = float(self.ps['hb'])
        self.ps_m = float(self.ps['m'])
        self.ps_Q = float(self.ps['Q'])
        self.ps_Q0 = float(self.ps['Q0'])
        self.ps_ub = float(self.ps['u_b'])
        self.ps_x = float(self.ps['x'])

        self.tmin = self.tmin_yr * 365.0
        self.tmax = self.tmax_yr * 365.0

        self.t_nondim = np.array([float(self.tt['t'][i]) for i in range(len(self.tt))])
        self.t = (self.ps_t / self.pd_td) * self.t_nondim
        self.surface_runoff_input = self.ps_m * self.ps_x**2 * np.array(
            [float(self.tt['E'][i]) for i in range(len(self.tt))]
        )
        self.q_out = self.ps_Q * (
            np.array([float(self.tt['Q_outQ'][i]) for i in range(len(self.tt))]) +
            np.array([float(self.tt['Q_outq'][i]) for i in range(len(self.tt))])
        )

        n_time = len(self.tt)
        self.pts_ni = self.oo['pts_ni'].item().astype(int) - 1
        self.pts_hs = np.column_stack([self.tt['pts_hs'][i] for i in range(n_time)]) * self.ps_h
        self.pts_hb = np.column_stack([self.tt['pts_hb'][i] for i in range(n_time)]) * self.ps_hb
        self.pts_u = np.column_stack([self.tt['pts_u'][i] for i in range(n_time)]) * self.ps_ub * self.pd_ty

        ni_gps = self.pp['ni_gps'].item().astype(int)
        gps_rows = []
        gps_matched = []
        pts_ni_source = self.oo['pts_ni'].item().astype(int)
        for ig, gps_index in enumerate(ni_gps):
            idx = np.where(pts_ni_source == gps_index)[0]
            if len(idx) > 0:
                gps_rows.append(idx[0])
                gps_matched.append(ig)
        self.gps_rows = np.array(gps_rows)
        self.gps_matched = np.array(gps_matched)

        ni_gps_0 = ni_gps - 1
        xx_flat = (self.ps_x / 1e3) * np.float64(self.gg['nx'].item()).flatten(order='F')
        yy_flat = (self.ps_x / 1e3) * np.float64(self.gg['ny'].item()).flatten(order='F')
        gps_x_km = xx_flat[ni_gps_0]
        gps_y_km = yy_flat[ni_gps_0]
        gps_cluster = np.zeros(len(ni_gps), dtype=int)
        cluster_centers = np.array([[0.0, 0.0], [20.0, -15.0], [40.0, -25.0]])
        for ig in range(len(ni_gps)):
            dists = np.sum((cluster_centers - [gps_x_km[ig], gps_y_km[ig]]) ** 2, axis=1)
            gps_cluster[ig] = np.argmin(dists) + 1
        dist_c2 = np.sqrt((gps_x_km - 17.0) ** 2 + (gps_y_km + 25.0) ** 2)
        gps_cluster[np.argmin(dist_c2)] = 4
        dist_c1 = np.sqrt((gps_x_km - 5.0) ** 2 + (gps_y_km + 10.0) ** 2)
        gps_cluster[np.argmin(dist_c1)] = 5

        self.gps_cluster = gps_cluster[self.gps_matched]
        self.gps_x_km = gps_x_km[self.gps_matched]
        self.gps_y_km = gps_y_km[self.gps_matched]

        self.hb_hs = self.pts_hb + self.pts_hs
        self.t_drainage = np.asarray(self.pp['t_drainage'].item()) * self.ps_t / self.pd_td
        self.t_duration = np.asarray(self.pp['t_duration'].item()) * self.ps_t / self.pd_td
        self.V_l = np.asarray(self.pp['V_l'].item(), dtype=float)

        lake_volume_m3 = self.ps_Q0 * self.ps_t * self.V_l
        lake_duration_days = np.asarray(self.t_duration, dtype=float)
        lake_center_days = np.asarray(self.t_drainage, dtype=float)
        valid_lakes = (
            np.isfinite(lake_volume_m3) &
            np.isfinite(lake_duration_days) &
            np.isfinite(lake_center_days) &
            (lake_duration_days > 0)
        )
        if np.any(valid_lakes):
            duration_seconds = lake_duration_days[valid_lakes] * 24.0 * 60.0 * 60.0
            time_offset = (
                self.t[np.newaxis, :] - lake_center_days[valid_lakes, np.newaxis]
            ) / lake_duration_days[valid_lakes, np.newaxis]
            gaussian_prefactor = lake_volume_m3[valid_lakes, np.newaxis] / (
                np.sqrt(2.0 * np.pi) * duration_seconds[:, np.newaxis]
            )
            self.lake_drainage_inflow = np.sum(
                gaussian_prefactor * np.exp(-0.5 * time_offset ** 2),
                axis=0,
            )
        else:
            self.lake_drainage_inflow = np.zeros_like(self.t)

        if self.right_column_mode == 'speed':
            dt_days = np.diff(self.t)
            if np.any(dt_days <= 0):
                raise ValueError('Expected t to be strictly increasing for dhdt_model computation.')
            self.t_mid = 0.5 * (self.t[:-1] + self.t[1:])
            self.right_model = np.diff(self.hb_hs, axis=1) / dt_days[np.newaxis, :] * 365.0
        else:
            t_ref = 151.0
            self.t_ref = t_ref
            z_model_ref = np.array([np.interp(t_ref, self.t, hb_hs_row) for hb_hs_row in self.hb_hs])
            self.right_model = self.hb_hs - z_model_ref[:, np.newaxis]

        self._loaded = True

    @staticmethod
    def _smooth_series(time_values, series_values, max_window):
        valid = np.isfinite(time_values) & np.isfinite(series_values)
        time_values = time_values[valid]
        series_values = series_values[valid]
        if time_values.size == 0:
            return None, None
        sort_index = np.argsort(time_values)
        time_sorted = time_values[sort_index]
        series_sorted = series_values[sort_index]
        window = max(1, min(max_window, len(series_sorted)))
        series_smooth = np.convolve(series_sorted, np.ones(window) / window, mode='same')
        return time_sorted, series_smooth

    @staticmethod
    def _station_vertical_displacement(t_station, w_station, t_ref):
        t_station = np.atleast_1d(t_station)
        w_station = np.atleast_1d(w_station)
        valid = np.isfinite(t_station) & np.isfinite(w_station)
        if np.sum(valid) <= 1:
            return None, None

        t_station = t_station[valid]
        w_station = w_station[valid]
        sort_index = np.argsort(t_station)
        t_station = t_station[sort_index]
        w_station = w_station[sort_index]

        t_station, unique_idx = np.unique(t_station, return_index=True)
        w_station = w_station[unique_idx]
        if t_station.size <= 1:
            return None, None

        z_station = cumulative_trapezoid(w_station / 365.0, t_station, initial=0.0)
        z_station = z_station - np.interp(t_ref, t_station, z_station)
        return t_station, z_station

    def _annotate_panel_a_speedup(self, ax, members, members_obs):
        target_day = 245.0
        series_time = None
        series_speed = None

        if len(members) > 0:
            series_time = self.t
            series_speed = np.mean(self.pts_u[members, :], axis=0)
        elif len(members_obs) > 0:
            t_all = np.concatenate([self.stations[im]['t_doy'] for im in members_obs])
            u_all = np.concatenate([self.stations[im]['u_s'] for im in members_obs])
            series_time, series_speed = self._smooth_series(t_all, u_all, max_window=6)

        if series_time is not None and series_speed is not None:
            valid = np.isfinite(series_time) & np.isfinite(series_speed)
            if np.any(valid):
                series_time = np.asarray(series_time)[valid]
                series_speed = np.asarray(series_speed)[valid]
                if series_time.size > 1 and series_time[0] <= target_day <= series_time[-1]:
                    speed_at_target = float(np.interp(target_day, series_time, series_speed))
                    speed_at_target = 150
                    # ax.annotate(
                    #     'runoff only',
                    #     xy=(target_day, speed_at_target),
                    #     xytext=(248.0, min(speed_at_target + 45.0, 335.0)),
                    #     fontsize=7,
                    #     ha='left',
                    #     va='bottom',
                    #     arrowprops=dict(arrowstyle='->', lw=0.8, color='k'),
                    #     bbox=dict(boxstyle='round,pad=0.2', fc='white', ec='none', alpha=0.85),
                    # )

        band_mask = (self.t_drainage >= self.tmin) & (self.t_drainage <= self.tmax)
        if np.any(band_mask):
            band_centers = self.t_drainage[band_mask]
            band_center = float(band_centers[np.argmin(np.abs(band_centers - target_day))])
            # ax.annotate(
            #     'HF+moulin inputs+runoff',
            #     xy=(196, 150.0),
            #     xytext=(band_center - 60.0, 195.0),
            #     fontsize=6,
            #     ha='left',
            #     va='bottom',
            #     arrowprops=dict(arrowstyle='->', lw=0.8, color='k'),
            #     bbox=dict(boxstyle='round,pad=0.2', fc='white', ec='none', alpha=0.0),
            # )

    def _plot_horizontal_speed(self, ax, members, members_obs, label, col, col_light, col_dark, panel_id=''):
        for i_d in range(len(self.t_drainage)):
            td_d = self.t_drainage[i_d]
            ts_d = self.t_duration[i_d]
            ax.axvspan(td_d - 0.5 * ts_d, td_d + 0.5 * ts_d, color=(0.7, 0.85, 1.0), alpha=0.4)

        if len(members_obs) > 0:
            all_t = np.concatenate([np.atleast_1d(self.stations[im]['t_doy']) for im in members_obs])
            if len(all_t) > 0:
                t_grid = np.linspace(np.nanmin(all_t), np.nanmax(all_t), 1000)
                u_interps = []
                for im in members_obs:
                    t_im = np.atleast_1d(self.stations[im]['t_doy'])
                    u_im = np.atleast_1d(self.stations[im]['u_s'])
                    valid = np.isfinite(t_im) & np.isfinite(u_im)
                    if np.sum(valid) > 1:
                        sort_index = np.argsort(t_im[valid])
                        u_interps.append(
                            np.interp(t_grid, t_im[valid][sort_index], u_im[valid][sort_index], left=np.nan, right=np.nan)
                        )
                if len(u_interps) > 0:
                    u_interps = np.array(u_interps)
                    valid_cols = np.any(np.isfinite(u_interps), axis=0)
                    if np.any(valid_cols):
                        with np.errstate(invalid='ignore'):
                            u_min = np.nanmin(u_interps[:, valid_cols], axis=0)
                            u_max = np.nanmax(u_interps[:, valid_cols], axis=0)
                        ax.fill_between(t_grid[valid_cols], u_min, u_max, color=col_light, alpha=0.5, label=r'obs range (Stevens et al., 2026)')

        if len(members_obs) > 0:
            t_all = np.concatenate([self.stations[im]['t_doy'] for im in members_obs])
            u_all = np.concatenate([self.stations[im]['u_s'] for im in members_obs])
            t_sorted, u_smooth = self._smooth_series(t_all, u_all, max_window=6)
            if t_sorted is not None:
                ax.plot(t_sorted, u_smooth, '-', color=col, lw=1.2, label=f'{label} obs mean')

        if len(members) > 0:
            ax.plot(self.t, np.mean(self.pts_u[members, :], axis=0), '-', color=col_dark, lw=1.2, label='Model')

        ax.set_ylabel('U [m/yr]', fontsize=7)
        ax.annotate(
            f'({panel_id}) {label} horizontal speed',
            xy=(0.02, 0.85),
            xycoords='axes fraction',
            fontsize=8,
            fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.2', fc='white', ec='none', alpha=0.0),
        )
        ax.set_xlim(self.tmin, self.tmax)
        ax.set_ylim(50, 350)
        ax.tick_params(axis='both', which='major', labelsize=7)
        ax.grid(True, which='both', alpha=0.3)
        if panel_id == 'a':
            self._annotate_panel_a_speedup(ax, members, members_obs)
        ax.legend(loc='upper right', fontsize=8, framealpha=0)

    def _plot_right_speed(self, ax, members, members_obs, label, col, col_light, col_dark, panel_id=''):
        for i_d in range(len(self.t_drainage)):
            td_d = self.t_drainage[i_d]
            ts_d = self.t_duration[i_d]
            ax.axvspan(td_d - 0.5 * ts_d, td_d + 0.5 * ts_d, color=(0.7, 0.85, 1.0), alpha=0.4)

        t_all_w = []
        w_all = []
        if len(members_obs) > 0:
            all_t_w = []
            for im in members_obs:
                ws = self.stations[im]['w_s']
                tt_w = self.stations[im]['t_doy']
                if isinstance(ws, np.ndarray) and ws.size > 0:
                    all_t_w.append(np.atleast_1d(tt_w))
                    t_all_w.append(tt_w)
                    w_all.append(ws)
            if len(all_t_w) > 0:
                all_t_w = np.concatenate(all_t_w)
                t_grid_w = np.linspace(np.nanmin(all_t_w), np.nanmax(all_t_w), 1000)
                w_interps = []
                for im in members_obs:
                    ws = self.stations[im]['w_s']
                    tt_w = self.stations[im]['t_doy']
                    if isinstance(ws, np.ndarray) and ws.size > 0:
                        tt_w = np.atleast_1d(tt_w)
                        ws = np.atleast_1d(ws)
                        valid = np.isfinite(tt_w) & np.isfinite(ws)
                        if np.sum(valid) > 1:
                            sort_index = np.argsort(tt_w[valid])
                            w_interps.append(
                                np.interp(t_grid_w, tt_w[valid][sort_index], ws[valid][sort_index], left=np.nan, right=np.nan)
                            )
                if len(w_interps) > 0:
                    w_interps = np.array(w_interps)
                    valid_cols = np.any(np.isfinite(w_interps), axis=0)
                    if np.any(valid_cols):
                        with np.errstate(invalid='ignore'):
                            w_min = np.nanmin(w_interps[:, valid_cols], axis=0)
                            w_max = np.nanmax(w_interps[:, valid_cols], axis=0)
                        ax.fill_between(t_grid_w[valid_cols], w_min, w_max, color=col_light, alpha=0.5, label=r'obs range')

        if len(t_all_w) > 0:
            t_all_w = np.concatenate(t_all_w)
            w_all = np.concatenate(w_all)
            t_sorted, w_smooth = self._smooth_series(t_all_w, w_all, max_window=max(1, round(len(w_all) / 200)))
            if t_sorted is not None:
                ax.plot(t_sorted, w_smooth, '-', color=col, lw=1.2, label=f'{label} obs mean')

        if len(members) > 0:
            ax.plot(self.t_mid, np.mean(self.right_model[members, :], axis=0), '-', color=col_dark, lw=1.1, label='Model')

        ax.set_ylabel('w [m/yr]', fontsize=7)
        ax.annotate(
            f'({panel_id}) {label} vertical speed',
            xy=(0.02, 0.85),
            xycoords='axes fraction',
            fontsize=8,
            fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.2', fc='white', ec='none', alpha=0.0),
        )
        ax.set_xlim(self.tmin, self.tmax)
        ax.set_ylim(-50, 100)
        ax.tick_params(axis='both', which='major', labelsize=7)
        ax.grid(True, which='both', alpha=0.3)
        ax.legend(loc='upper right', fontsize=8, framealpha=0)

    def _plot_right_displacement(self, ax, members, members_obs, label, col, col_light, col_dark, panel_id=''):
        for i_d in range(len(self.t_drainage)):
            td_d = self.t_drainage[i_d]
            ts_d = self.t_duration[i_d]
            ax.axvspan(td_d - 0.5 * ts_d, td_d + 0.5 * ts_d, color=(0.7, 0.85, 1.0), alpha=0.4)

        z_interps = []
        disp_time_support = []
        for im in members_obs:
            t_station, z_station = self._station_vertical_displacement(
                self.stations[im]['t_doy'],
                self.stations[im]['w_s'],
                self.t_ref,
            )
            if t_station is None:
                continue
            disp_time_support.append(t_station)

            interp_start = max(t_station[0], self.tmin)
            interp_end = min(t_station[-1], self.tmax)
            if interp_end <= interp_start:
                continue
            z_interps.append((t_station, z_station))

        if len(z_interps) > 0:
            all_t_disp = np.concatenate(disp_time_support)
            t_grid_z = np.linspace(np.nanmin(all_t_disp), np.nanmax(all_t_disp), 1000)
            z_grid = []
            for t_station, z_station in z_interps:
                z_grid.append(np.interp(t_grid_z, t_station, z_station, left=np.nan, right=np.nan))
            z_grid = np.array(z_grid)
            valid_cols = np.any(np.isfinite(z_grid), axis=0)
            if np.any(valid_cols):
                with np.errstate(invalid='ignore'):
                    z_min = np.nanmin(z_grid[:, valid_cols], axis=0)
                    z_max = np.nanmax(z_grid[:, valid_cols], axis=0)
                    z_mean = np.nanmean(z_grid[:, valid_cols], axis=0)
                ax.fill_between(t_grid_z[valid_cols], z_min, z_max, color=col_light, alpha=0.5, label=r'obs range')
                ax.plot(t_grid_z[valid_cols], z_mean, '-', color=col, lw=1.2, label=f'{label} obs mean')

        if len(members) > 0:
            ax.plot(self.t, np.mean(self.right_model[members, :], axis=0), '-', color=col_dark, lw=1.1, label='Model')

        ax.set_ylabel('$\\Delta z$ [m]', fontsize=7)
        ax.annotate(
            f'({panel_id}) {label} vertical displacement',
            xy=(0.02, 0.85),
            xycoords='axes fraction',
            fontsize=8,
            fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.2', fc='white', ec='none', alpha=0.0),
        )
        ax.set_xlim(self.tmin, self.tmax)
        ax.set_ylim(-1, 1)
        ax.tick_params(axis='both', which='major', labelsize=7)
        ax.grid(True, which='both', alpha=0.3)
        ax.legend(loc='upper right', fontsize=8, framealpha=0)

    def _plot_surface_runoff_input(self, ax, show_ylabel=True, show_legend=True):
        for i_d in range(len(self.t_drainage)):
            td_d = self.t_drainage[i_d]
            ts_d = self.t_duration[i_d]
            ax.axvspan(td_d - 0.5 * ts_d, td_d + 0.5 * ts_d, color=(0.7, 0.85, 1.0), alpha=0.4)

        runoff = np.asarray(self.surface_runoff_input, dtype=float)
        lake_inflow = np.asarray(self.lake_drainage_inflow, dtype=float)
        q_out = np.asarray(self.q_out, dtype=float)
        runoff_plot = np.where(np.isfinite(runoff), runoff, np.nan)
        lake_inflow_plot = np.where(np.isfinite(lake_inflow), lake_inflow, np.nan)
        q_out_plot = np.where(np.isfinite(q_out), q_out, np.nan)
        plot_mask = (self.t >= self.tmin) & (self.t <= self.tmax)
        panel_values = np.concatenate([
            runoff_plot[plot_mask & np.isfinite(runoff_plot)],
            lake_inflow_plot[plot_mask & np.isfinite(lake_inflow_plot)],
            q_out_plot[plot_mask & np.isfinite(q_out_plot)],
        ])
        positive_values = panel_values[panel_values >= 0]
        if positive_values.size > 0:
            y_min = 0.0
            y_max = float(np.nanmax(positive_values)) * 1.1
            y_max = 8000
        else:
            y_min, y_max = 0.0, 1.0

        ax.plot(self.t, runoff_plot, '-', color='k', lw=1.2, label='surface runoff input')
        ax.plot(self.t, lake_inflow_plot, '-', color='tab:blue', lw=1.2, label='lake drainage input')
        ax.plot(self.t, q_out_plot, '-', color='tab:green', lw=1.2, label='outflow rate')
        ax.annotate(
            'Inflow and outflow rates',
            xy=(0.02, 0.85),
            xycoords='axes fraction',
            fontsize=7,
            fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.2', fc='white', ec='none', alpha=0.0),
        )
        ax.set_xlim(self.tmin, self.tmax)
        ax.set_ylim(y_min, y_max)
        if show_ylabel:
            ax.set_ylabel(r'Q [m$^3$/s]', fontsize=7)
        ax.locator_params(axis='y', nbins=4)
        ax.tick_params(axis='both', which='major', labelsize=7)
        ax.grid(True, which='both', alpha=0.3)
        if show_legend:
            ax.legend(loc='upper right', fontsize=7, framealpha=0)

    def plot(self, save_path=None):
        if not self._loaded:
            self.load_data()

        plt.style.use('classic')
        plt.rcParams.update({'font.size': 9, 'figure.dpi': 150})
        n_clust = len(self.clusters_to_plot)
        n_rows = n_clust + 1
        fig, axes = plt.subplots(
            n_rows,
            2,
            figsize=(11, 1.2 + 1.8 * n_clust),
            sharex=True,
            gridspec_kw={'hspace': 0.1, 'height_ratios': [0.65] + [1.0] * n_clust},
        )
        if n_rows == 1:
            axes = np.atleast_2d(axes)

        self._plot_surface_runoff_input(axes[0, 0], show_ylabel=True, show_legend=True)
        self._plot_surface_runoff_input(axes[0, 1], show_ylabel=True, show_legend=True)

        for ic, kc in enumerate(self.clusters_to_plot):
            col = self.cluster_colors[kc - 1]
            col_light = tuple(min(1.0, 0.4 + 0.6 * c) for c in col)
            col_dark = tuple(c * 0.4 for c in col)
            members = self.gps_rows[self.gps_cluster == kc]
            members_obs = np.where(self.gps_cluster == kc)[0]
            label = self.cluster_labels[kc - 1]

            panel_left = chr(ord('a') + 2 * ic)
            panel_right = chr(ord('a') + 2 * ic + 1)

            ax_u = axes[ic + 1, 0]
            self._plot_horizontal_speed(ax_u, members, members_obs, label, col, col_light, col_dark, panel_id=panel_left)

            ax_right = axes[ic + 1, 1]
            if self.right_column_mode == 'speed':
                self._plot_right_speed(ax_right, members, members_obs, label, col, col_light, col_dark, panel_id=panel_right)
            else:
                self._plot_right_displacement(ax_right, members, members_obs, label, col, col_light, col_dark, panel_id=panel_right)

            if ic == n_clust - 1:
                ax_u.set_xlabel('day of year', fontsize=7)
                ax_right.set_xlabel('day of year', fontsize=7)

        with warnings.catch_warnings():
            warnings.filterwarnings(
                'ignore',
                message='This figure includes Axes that are not compatible with tight_layout, so results might be incorrect.',
                category=UserWarning,
            )
            fig.tight_layout(rect=(0, 0, 1, 0.97))

        if save_path is None:
            if self.right_column_mode == 'speed':
                save_path = self.repo_root / 'figures' / f'{self.casename}_gps_comparison.png'
            else:
                save_path = self.repo_root / 'figures' / f'{self.casename}_gps_displacement_comparison.png'
        else:
            save_path = Path(save_path)
            if not save_path.is_absolute():
                save_path = self.repo_root / save_path

        fig.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f'Saved to {save_path}')
        return fig