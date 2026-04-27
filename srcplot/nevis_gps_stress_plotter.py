from pathlib import Path
import shutil
import subprocess

import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import ScalarFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable


class NevisGpsStressPlotter:
    Z95 = 1.959963984540054
    GPS_CLUSTER_COLORS = [
        (1.0, 0.45, 0.20),
        (0.55, 1.0, 0.55),
        (0.40, 0.75, 1.0),
        (0.85, 0.45, 1.0),
        (1.0, 0.85, 0.30),
    ]
    GPS_CLUSTER_LABELS = ['S1', 'S2', 'S3', 'S4', 'S5']
    STRESS_GROUPS = [
        (
            r'Membrane stress [kPa]',
            [('tauxx_gps', r'$\tau_{xx}$'), ('tauyy_gps', r'$\tau_{yy}$'), ('tauxy_gps', r'$\tau_{xy}$')],
            1e3,
        ),
        (
            r'Principal stress [kPa]',
            [('sigma1_gps', r'$\sigma_1$'), ('sigma2_gps', r'$\sigma_2$')],
            1e3,
        ),
        (
            r'Depth-integrated stress [kPa]',
            [('Txx_gps', r'$T_{xx}$'), ('Tyy_gps', r'$T_{yy}$'), ('Txy_gps', r'$T_{xy}$')],
            1e3,
        ),
        (
            r'Basal & driving stress [kPa]',
            [('tau_b_gps', r'$\tau_b$'), ('taudx_gps', r'$\tau_{dx}$'), ('taudy_gps', r'$\tau_{dy}$')],
            1e3,
        ),
    ]
    COLORS = {
        'tauxx_gps': '#b22222',
        'tauyy_gps': '#1f77b4',
        'tauxy_gps': '#2e8b57',
        'sigma1_gps': '#7f0000',
        'sigma2_gps': '#004c6d',
        'Txx_gps': '#8c564b',
        'Tyy_gps': '#9467bd',
        'Txy_gps': '#e377c2',
        'tau_b_gps': '#111111',
        'taudx_gps': '#ff7f0e',
        'taudy_gps': '#17becf',
    }
    EVENT_COLORS = [
        '#D55E00',
        '#0072B2',
        '#009E73',
        '#CC79A7',
        '#E69F00',
        '#56B4E9',
        '#000000',
        '#F0E442',
    ]

    def __init__(
        self,
        casename,
        repo_root='.',
        matlab_path=None,
        stress_tmin_yr=0.0,
        stress_tmax_yr=1.0,
        auto_generate_stress=True,
    ):
        self.casename = casename
        self.repo_root = Path(repo_root).resolve()
        self.data_dir = self.repo_root / 'data'
        self.figure_dir = self.repo_root / 'figures'
        self.result_dir = self.repo_root / 'results' / casename
        self.mat_path = self.result_dir / f'{casename}.mat'
        self.stress_path = self.data_dir / f'gps_stresses_{casename}.mat'
        default_matlab = shutil.which('matlab')
        if default_matlab is None:
            fallback_matlab = Path('/Applications/MATLAB_R2025b.app/bin/matlab')
            if fallback_matlab.exists():
                default_matlab = str(fallback_matlab)
        self.matlab_path = matlab_path or default_matlab
        self.stress_tmin_yr = stress_tmin_yr
        self.stress_tmax_yr = stress_tmax_yr
        self.auto_generate_stress = auto_generate_stress
        self._stress = None
        self._raw = None

    @staticmethod
    def mat_array(value, dtype=float):
        while isinstance(value, np.ndarray) and value.dtype == object and value.size == 1:
            value = value.item()
        return np.atleast_1d(np.asarray(value, dtype=dtype)).squeeze()

    @staticmethod
    def _quote_matlab_string(value):
        return str(value).replace("'", "''")

    @staticmethod
    def _add_panel_annotation(ax, text, x=0.02, y=0.95, fontsize=9):
        ax.text(
            x,
            y,
            text,
            transform=ax.transAxes,
            ha='left',
            va='top',
            fontsize=fontsize,
            clip_on=False,
        )

    @staticmethod
    def _add_marker_label(ax, x, y, text, dx=4, dy=4, fontsize=6):
        ax.annotate(
            text,
            xy=(x, y),
            xytext=(dx, dy),
            textcoords='offset points',
            ha='left',
            va='bottom',
            fontsize=fontsize,
            color='k',
            bbox={'facecolor': 'white', 'edgecolor': 'none', 'alpha': 0.75, 'pad': 0.2},
            zorder=9,
            clip_on=False,
        )

    @staticmethod
    def _add_panel_colorbar(ax, mappable, label, label_side='right', labelpad=None):
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='3%', pad=0.08)
        colorbar = plt.colorbar(mappable, cax=cax)
        colorbar.set_label(label, fontsize=7)
        colorbar.ax.yaxis.set_label_position(label_side)
        if labelpad is not None:
            colorbar.ax.yaxis.labelpad = labelpad
        colorbar.ax.tick_params(labelsize=7)
        return colorbar

    @staticmethod
    def _add_sigma1_quiver(ax, xx, yy, sigma1, t1, xlim, ylim, stride=6,
                           scale_factor=0.003, min_mag=5.0, color='k',
                           linewidth=0.2, alpha=0.7):
        """Overlay headless quiver showing sigma1 direction and magnitude."""
        nI, nJ = sigma1.shape
        # Build regular grid from the mesh arrays
        x_vals = xx[:, 0] if xx.ndim == 2 else np.unique(xx)
        y_vals = yy[0, :] if yy.ndim == 2 else np.unique(yy)
        # Subsample within view limits with stride
        ix = np.where((x_vals >= xlim[0]) & (x_vals <= xlim[1]))[0][::stride]
        iy = np.where((y_vals >= ylim[0]) & (y_vals <= ylim[1]))[0][::stride]
        if ix.size == 0 or iy.size == 0:
            return
        xq = x_vals[ix]
        yq = y_vals[iy]
        Xq, Yq = np.meshgrid(xq, yq, indexing='ij')
        S = sigma1[np.ix_(ix, iy)]
        T = t1[np.ix_(ix, iy)]
        mag = np.abs(S)
        mask = np.isfinite(mag) & (mag >= min_mag)
        if not np.any(mask):
            return
        dx = mag * np.cos(T) * scale_factor
        dy = mag * np.sin(T) * scale_factor
        # Draw headless lines: plot both +/- directions (σ₁ is axial, not vectorial)
        ax.quiver(
            Xq[mask], Yq[mask], dx[mask], dy[mask],
            angles='xy', scale_units='xy', scale=1,
            headwidth=0, headlength=0, headaxislength=0,
            color=color, linewidth=linewidth, alpha=alpha, zorder=4,
        )
        ax.quiver(
            Xq[mask], Yq[mask], -dx[mask], -dy[mask],
            angles='xy', scale_units='xy', scale=1,
            headwidth=0, headlength=0, headaxislength=0,
            color=color, linewidth=linewidth, alpha=alpha, zorder=4,
        )

    def _run_matlab_batch(self, matlab_code):
        if self.matlab_path is None:
            raise RuntimeError('MATLAB executable not found. Cannot run MATLAB helper scripts automatically.')
        result = subprocess.run(
            [self.matlab_path, '-batch', matlab_code],
            cwd=self.repo_root,
            text=True,
            capture_output=True,
            check=False,
        )
        if result.stdout:
            print(result.stdout.strip())
        if result.returncode != 0:
            stderr = result.stderr.strip() if result.stderr else 'MATLAB exited with a non-zero status.'
            raise RuntimeError(stderr)
        return result

    def ensure_stress_data(self, force=False):
        if self.stress_path.exists() and not force:
            return self.stress_path
        if not self.mat_path.exists():
            raise FileNotFoundError(f'Missing case result file: {self.mat_path}')

        matlab_code = (
            f"cd('{self._quote_matlab_string(self.repo_root)}'); "
            f"nevis_gps_stress('{self._quote_matlab_string(self.casename)}', "
            f"{self.stress_tmin_yr:.12g}, {self.stress_tmax_yr:.12g});"
        )
        print(f'Generating GPS stress file with MATLAB for {self.casename} ...')
        try:
            self._run_matlab_batch(matlab_code)
        except RuntimeError as error:
            raise RuntimeError(f'Failed to generate GPS stress data for {self.casename}:\n{error}') from error
        if not self.stress_path.exists():
            raise FileNotFoundError(f'MATLAB finished without creating {self.stress_path}')
        self._stress = None
        return self.stress_path

    def _load_data(self, required_fields=(), force_recompute=False):
        if force_recompute or not self.stress_path.exists():
            self.ensure_stress_data(force=True)
        if self._raw is None:
            self._raw = sio.loadmat(self.mat_path, squeeze_me=True)
        if self._stress is None:
            self._stress = sio.loadmat(self.stress_path, squeeze_me=True)

        missing = [field for field in required_fields if field not in self._stress]
        if missing and self.auto_generate_stress and not force_recompute:
            self.ensure_stress_data(force=True)
            self._stress = sio.loadmat(self.stress_path, squeeze_me=True)
            missing = [field for field in required_fields if field not in self._stress]
        if missing:
            raise KeyError(
                'Missing fields in GPS stress file: ' + ', '.join(missing) +
                f'. Re-run nevis_gps_stress for {self.casename}.'
            )
        return self._stress, self._raw

    def _shared_context(self, required_fields=(), force_recompute=False):
        stress, raw = self._load_data(required_fields=required_fields, force_recompute=force_recompute)
        ps = raw['ps']
        pd_r = raw['pd']
        pp = raw['pp']
        return stress, raw, ps, pd_r, pp

    @staticmethod
    def _time_bounds_days(t_days):
        t_days = np.atleast_1d(np.asarray(t_days, dtype=float))
        finite = t_days[np.isfinite(t_days)]
        if finite.size == 0:
            raise ValueError('Stress file does not contain any finite time samples.')
        return float(np.min(finite)), float(np.max(finite))

    def _ensure_stress_window(self, required_fields, tmin_days, tmax_days, force_recompute=False):
        stress, raw = self._load_data(required_fields=required_fields, force_recompute=force_recompute)
        available_min, available_max = self._time_bounds_days(stress['t_days'])
        needs_extension = tmin_days < available_min or tmax_days > available_max

        if needs_extension and self.auto_generate_stress and not force_recompute:
            previous_tmin_yr = float(self.stress_tmin_yr)
            previous_tmax_yr = float(self.stress_tmax_yr)
            pad_days = 1.0
            self.stress_tmin_yr = min(previous_tmin_yr, (tmin_days - pad_days) / 365.0)
            self.stress_tmax_yr = max(previous_tmax_yr, (tmax_days + pad_days) / 365.0)
            print(
                'Extending GPS stress time window from '
                f'[{365.0 * previous_tmin_yr:.1f}, {365.0 * previous_tmax_yr:.1f}] d to '
                f'[{365.0 * self.stress_tmin_yr:.1f}, {365.0 * self.stress_tmax_yr:.1f}] d.'
            )
            stress, raw = self._load_data(required_fields=required_fields, force_recompute=True)

        return stress, raw

    def _empty_window_error_message(self, requested_min_days, requested_max_days, available_t_days):
        available_min, available_max = self._time_bounds_days(available_t_days)
        return (
            'No stress samples fall within the requested plotting window. '
            f'Requested [{requested_min_days:.1f}, {requested_max_days:.1f}] d, '
            f'available [{available_min:.1f}, {available_max:.1f}] d.'
        )

    def _lake_to_gps(self, pp, ps, pd_r, stress):
        gps_x_km = np.atleast_1d(stress['gps_x_km']).astype(float)
        gps_y_km = np.atleast_1d(stress['gps_y_km']).astype(float)
        if gps_x_km.size == 0:
            raise ValueError('No GPS stations were saved in the stress file.')
        if 'x_l' not in pp.dtype.names or 'y_l' not in pp.dtype.names:
            return np.array([], dtype=int)
        lake_x_km = self.mat_array(pp['x_l']) * float(ps['x']) / 1e3
        lake_y_km = self.mat_array(pp['y_l']) * float(ps['x']) / 1e3
        if lake_x_km.size == 0:
            return np.array([], dtype=int)
        return np.argmin(
            (lake_x_km[:, None] - gps_x_km[None, :]) ** 2 +
            (lake_y_km[:, None] - gps_y_km[None, :]) ** 2,
            axis=1,
        )

    @staticmethod
    def _selected_lake_label_map(event_ids, t_drainage):
        event_ids = np.unique(np.asarray(event_ids, dtype=int))
        if event_ids.size == 0:
            return {}
        sort_order = np.argsort(t_drainage[event_ids], kind='stable')
        sorted_event_ids = event_ids[sort_order]
        return {
            int(event_idx): f'L{label_index + 1}'
            for label_index, event_idx in enumerate(sorted_event_ids)
        }

    @staticmethod
    def _event_labels(event_ids, label_map, one_based=False):
        labels = []
        for event_id in event_ids:
            zero_based_event_id = int(event_id) - 1 if one_based else int(event_id)
            labels.append(label_map.get(zero_based_event_id, f'E{zero_based_event_id + 1}'))
        return labels

    @staticmethod
    def _sort_events_by_lake_label(event_ids, label_map):
        return sorted(
            [int(event_id) for event_id in event_ids],
            key=lambda event_id: int(label_map.get(event_id, f'E{event_id + 1}')[1:]),
        )

    def _format_drainage_summary(self, event_idx, t_drainage, t_duration, label=None):
        drainage_duration = float(t_duration[event_idx])
        shaded_duration = float(2 * self.Z95 * t_duration[event_idx])
        event_label = label or f'E{event_idx + 1}'
        return (
            f'{event_label} at {t_drainage[event_idx]:.1f} d, '
            f'drainage duration={drainage_duration:.2f} d, '
            f'blue shading={shaded_duration:.2f} d'
        )

    @staticmethod
    def _cluster_gps_stations(gps_x_km, gps_y_km):
        if gps_x_km.size == 0:
            return np.array([], dtype=int)
        cluster_centers = np.array([[0.0, 0.0], [20.0, -15.0], [40.0, -25.0]])
        gps_cluster = np.zeros(gps_x_km.size, dtype=int)
        for index in range(gps_x_km.size):
            dists = np.sum((cluster_centers - [gps_x_km[index], gps_y_km[index]]) ** 2, axis=1)
            gps_cluster[index] = int(np.argmin(dists)) + 1
        dist_c2 = np.sqrt((gps_x_km - 17.0) ** 2 + (gps_y_km + 25.0) ** 2)
        gps_cluster[int(np.argmin(dist_c2))] = 4
        dist_c1 = np.sqrt((gps_x_km - 5.0) ** 2 + (gps_y_km + 10.0) ** 2)
        gps_cluster[int(np.argmin(dist_c1))] = 5
        return gps_cluster

    def _gps_station_positions(self, raw):
        ps = raw['ps']
        pp = raw['pp']
        gg = raw['gg']
        ni_gps = self.mat_array(pp['ni_gps'], dtype=int) - 1
        if ni_gps.size == 0:
            return np.array([], dtype=float), np.array([], dtype=float)
        ps_x = float(ps['x'])
        xx_flat = (ps_x / 1e3) * np.float64(gg['nx'].item()).flatten(order='F')
        yy_flat = (ps_x / 1e3) * np.float64(gg['ny'].item()).flatten(order='F')
        return xx_flat[ni_gps], yy_flat[ni_gps]

    def _grid_context(self, raw):
        ps = raw['ps']
        gg = raw['gg']
        ps_x = float(ps['x'])
        xx = ps_x / 1e3 * np.float64(gg['nx'].item())
        yy = ps_x / 1e3 * np.float64(gg['ny'].item())
        nI = int(gg['nI'].item())
        nJ = int(gg['nJ'].item())
        return xx, yy, nI, nJ

    @staticmethod
    def _select_representative_drainage_groups(event_ids, t_drainage, n_groups=3):
        selected_times = t_drainage[event_ids]
        sort_order = np.argsort(selected_times)
        sorted_ids = event_ids[sort_order]
        sorted_times = selected_times[sort_order]
        if n_groups < 1:
            raise ValueError('n_groups must be at least 1.')
        if sorted_ids.size < n_groups:
            raise ValueError(f'Need at least {n_groups} drainage events to define {n_groups} drainage moments.')
        if n_groups == 1:
            groups = [np.arange(sorted_ids.size)]
        else:
            gap_rank = np.argsort(np.diff(sorted_times))
            split_points = np.sort(gap_rank[-(n_groups - 1):] + 1)
            groups = np.split(np.arange(sorted_ids.size), split_points)

        representative_event_ids = []
        moment_times = []
        group_event_ids = []
        for group in groups:
            group_ids = sorted_ids[group]
            group_times = sorted_times[group]
            rep_local = int(np.argmin(np.abs(group_times - np.median(group_times))))
            representative_event_ids.append(int(group_ids[rep_local]))
            moment_times.append(float(group_times[rep_local]))
            group_event_ids.append((group_ids + 1).tolist())
        return representative_event_ids, moment_times, group_event_ids

    def _export_sigma1_frames(self, moment_times, raw):
        xx, yy, nI, nJ = self._grid_context(raw)
        ps = raw['ps']
        pd_r = raw['pd']
        oo = raw['oo']
        ps_t = float(ps['t'])
        pd_td = float(pd_r['td'])
        tspan = self.mat_array(oo['t_span']) * ps_t / pd_td
        frame_numbers = np.array([np.argmin(np.abs(tspan - td)) + 1 for td in moment_times], dtype=int)
        filenames = [self.result_dir / f'{frame_number:04d}.mat' for frame_number in frame_numbers]

        export_request = self.repo_root / 'export_ice_fields.mat'
        sio.savemat(
            export_request,
            {
                'filenames': np.array([str(filename) for filename in filenames], dtype=object),
                'casename': self.casename,
            },
        )
        print(f'Requesting MATLAB export for frames: {frame_numbers}')
        try:
            matlab_code = f"cd('{self._quote_matlab_string(self.repo_root)}'); run('src/nevis_export_ice_fields')"
            self._run_matlab_batch(matlab_code)
        finally:
            if export_request.exists():
                export_request.unlink()

        frames = []
        for filename in filenames:
            raw_frame = sio.loadmat(filename, squeeze_me=True)
            ice = raw_frame['ice_fields']
            ice_fields = set(ice.dtype.names or [])
            nout = ice['nout'].item().astype(int) - 1
            if 'sigma1_pa' in ice_fields:
                sigma1_pa = self.mat_array(ice['sigma1_pa']).astype(float).copy()
            else:
                stress_scale_kpa = float(ice['stress_scale'])
                sigma1 = np.atleast_1d(np.asarray(ice['sigma1'].item(), dtype=float)).copy()
                sigma1_pa = 1e3 * stress_scale_kpa * sigma1

            if 'sigma1_elastic_pa' in ice_fields:
                sigma1_elastic_pa = self.mat_array(ice['sigma1_elastic_pa']).astype(float).copy()
            else:
                sigma1_elastic_pa = np.full_like(sigma1_pa, np.nan)

            if 'sigma1_total_pa' in ice_fields:
                sigma1_total_pa = self.mat_array(ice['sigma1_total_pa']).astype(float).copy()
            else:
                sigma1_total_pa = np.full_like(sigma1_pa, np.nan)

            sigma1_kpa = sigma1_pa / 1e3
            sigma1_elastic_kpa = sigma1_elastic_pa / 1e3
            sigma1_total_kpa = sigma1_total_pa / 1e3
            sigma1_kpa[nout] = np.nan
            sigma1_elastic_kpa[nout] = np.nan
            sigma1_total_kpa[nout] = np.nan

            # Principal stress direction angles (radians)
            t1 = self.mat_array(ice['t1']).astype(float).copy() if 't1' in ice_fields else np.full_like(sigma1_kpa, np.nan)
            t1_elastic = self.mat_array(ice['t1_elastic']).astype(float).copy() if 't1_elastic' in ice_fields else np.full_like(sigma1_kpa, np.nan)
            t1_total = self.mat_array(ice['t1_total']).astype(float).copy() if 't1_total' in ice_fields else np.full_like(sigma1_kpa, np.nan)
            t1[nout] = np.nan
            t1_elastic[nout] = np.nan
            t1_total[nout] = np.nan

            frames.append(
                {
                    't_days': float(ice['t_days']),
                    'sigma1': sigma1_kpa.reshape(nI, nJ, order='F'),
                    'sigma1_elastic': sigma1_elastic_kpa.reshape(nI, nJ, order='F'),
                    'sigma1_total': sigma1_total_kpa.reshape(nI, nJ, order='F'),
                    't1': t1.reshape(nI, nJ, order='F'),
                    't1_elastic': t1_elastic.reshape(nI, nJ, order='F'),
                    't1_total': t1_total.reshape(nI, nJ, order='F'),
                }
            )
        return xx, yy, frames

    def plot_selected_region_sigma1_three_moments(
        self,
        selected_events,
        pad_x_km=8.0,
        pad_y_km=16.0,
        sigma1_clim=(0.0, 200.0),
        save_path=None,
        force_recompute=False,
    ):
        required_fields = ['lake_x_km', 'lake_y_km']
        stress, raw, ps, pd_r, pp = self._shared_context(required_fields, force_recompute=force_recompute)

        lake_x_km = np.atleast_1d(stress['lake_x_km']).astype(float)
        lake_y_km = np.atleast_1d(stress['lake_y_km']).astype(float)
        t_drainage = self.mat_array(pp['t_drainage']) * float(ps['t']) / float(pd_r['td'])

        event_ids = np.asarray(selected_events, dtype=int) - 1
        if np.any(event_ids < 0) or np.any(event_ids >= len(t_drainage)):
            raise ValueError('selected_events contains an event id outside the available range.')
        lake_label_by_event = self._selected_lake_label_map(event_ids, t_drainage)

        representative_event_ids, moment_times, group_event_ids = self._select_representative_drainage_groups(
            event_ids,
            t_drainage,
            n_groups=3,
        )
        xx, yy, frames = self._export_sigma1_frames(moment_times, raw)
        gps_x_km, gps_y_km = self._gps_station_positions(raw)
        gps_cluster = self._cluster_gps_stations(gps_x_km, gps_y_km)

        selected_lake_x = lake_x_km[event_ids]
        selected_lake_y = lake_y_km[event_ids]
        xlim = (float(np.min(selected_lake_x) - pad_x_km), float(np.max(selected_lake_x) + pad_x_km))
        ylim = (float(np.min(selected_lake_y) - pad_y_km + 10), float(np.max(selected_lake_y) + pad_y_km - 10))
        min_map_width_ratio = 1.2
        x_span = xlim[1] - xlim[0]
        y_span = ylim[1] - ylim[0]
        target_x_span = max(x_span, min_map_width_ratio * y_span)
        if target_x_span > x_span:
            extra_x = 0.5 * (target_x_span - x_span)
            xlim = (xlim[0] - extra_x, xlim[1] + extra_x)
        map_box_aspect = (ylim[1] - ylim[0]) / max(xlim[1] - xlim[0], np.finfo(float).eps)
        gps_region_mask = (
            (gps_x_km >= xlim[0]) & (gps_x_km <= xlim[1]) &
            (gps_y_km >= ylim[0]) & (gps_y_km <= ylim[1])
        )

        print(
            'Representative lakes:',
            [lake_label_by_event[event_id] for event_id in representative_event_ids],
        )
        print('Representative drainage times [d]:', [round(value, 2) for value in moment_times])

        plt.style.use('classic')
        plt.rcParams.update({
            'font.size': 9,
            'figure.dpi': 180,
            'mathtext.fontset': 'dejavusans',
        })
        fig, axes = plt.subplots(1, len(frames), figsize=(3.8 * len(frames), 3.6), constrained_layout=True)
        axes = np.atleast_1d(axes)
        panel_labels = [chr(ord('a') + index) for index in range(len(frames))]
        pcm = None

        for index, (ax, frame, rep_event_id, event_group) in enumerate(
            zip(axes, frames, representative_event_ids, group_event_ids)
        ):
            pcm = ax.pcolormesh(
                xx,
                yy,
                frame['sigma1'],
                cmap='jet',
                shading='gouraud',
                vmin=float(sigma1_clim[0]),
                vmax=float(sigma1_clim[1]),
            )

            for cluster_index, (cluster_color, cluster_label) in enumerate(
                zip(self.GPS_CLUSTER_COLORS, self.GPS_CLUSTER_LABELS),
                start=1,
            ):
                mask = (gps_cluster == cluster_index) & gps_region_mask
                if not np.any(mask):
                    continue
                ax.scatter(
                    gps_x_km[mask],
                    gps_y_km[mask],
                    s=11,
                    marker='s',
                    color=cluster_color,
                    edgecolors='k',
                    linewidths=0.3,
                    label=cluster_label if index == 0 else None,
                    zorder=5,
                )

            ax.scatter(
                selected_lake_x,
                selected_lake_y,
                s=26,
                facecolors='none',
                edgecolors='k',
                linewidths=0.6,
                zorder=6,
                label='Selected lakes' if index == 0 else None,
            )
            group_zero_based = np.asarray(event_group, dtype=int) - 1
            ax.scatter(
                lake_x_km[group_zero_based],
                lake_y_km[group_zero_based],
                s=46,
                facecolors='none',
                edgecolors='w',
                linewidths=1.1,
                zorder=7,
            )
            ax.scatter(
                lake_x_km[rep_event_id],
                lake_y_km[rep_event_id],
                s=65,
                marker='*',
                color='w',
                edgecolors='k',
                linewidths=0.5,
                zorder=8,
                label='Representative event' if index == 0 else None,
            )
            self._add_marker_label(
                ax,
                lake_x_km[rep_event_id],
                lake_y_km[rep_event_id],
                lake_label_by_event[rep_event_id],
            )

            event_label = ', '.join(self._event_labels(event_group, lake_label_by_event, one_based=True))
            ax.set_title(f'({panel_labels[index]}) t = {frame["t_days"]:.1f} d | {event_label}', fontsize=9)
            ax.set_xlim(*xlim)
            ax.set_ylim(*ylim)
            ax.set_aspect('equal')
            ax.set_xlabel('x (km)')
            if index == 0:
                ax.set_ylabel('y (km)')
                ax.legend(loc='upper left', fontsize=6, frameon=True, ncol=1)
            ax.tick_params(labelsize=7)

        if pcm is None:
            raise ValueError('No stress frames were exported for the requested drainage moments.')
        colorbar = fig.colorbar(pcm, ax=axes, shrink=0.95, pad=0.02)
        colorbar.set_label(r'$\sigma_1$ (kPa)', fontsize=8)
        colorbar.ax.tick_params(labelsize=7)

        if save_path is None:
            save_path = self.figure_dir / f'{self.casename}_selected_region_sigma1_three_moments.png'
        else:
            save_path = Path(save_path)
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f'Saved to {save_path}')
        return fig

    def plot_selected_lake_sigma1_and_stress_snapshots(
        self,
        selected_events,
        margin_days=6.0,
        n_stress_panels=3,
        pad_x_km=8.0,
        pad_y_km=8.0,
        stress_component='total',
        sigma1_clim=None,
        stress_panel_plan=None,
        reference_time_days=190.0,
        timeseries_tmin_days=None,
        timeseries_tmax_days=None,
        save_path=None,
        force_recompute=False,
    ):
        required_fields = ['t_days', 'lake_x_km', 'lake_y_km', 'sigma1_lake', 'sigma1_elastic_lake']
        stress, raw = self._load_data(required_fields=required_fields, force_recompute=force_recompute)
        ps = raw['ps']
        pd_r = raw['pd']
        pp = raw['pp']
        t_drainage = self.mat_array(pp['t_drainage']) * float(ps['t']) / float(pd_r['td'])
        t_duration = self.mat_array(pp['t_duration']) * float(ps['t']) / float(pd_r['td'])

        event_ids = np.asarray(selected_events, dtype=int) - 1
        if np.any(event_ids < 0) or np.any(event_ids >= len(t_drainage)):
            raise ValueError('selected_events contains an event id outside the available range.')
        lake_label_by_event = self._selected_lake_label_map(event_ids, t_drainage)

        if sigma1_clim is None:
            viscous_sigma1_clim = (-200.0, 200.0)
            elastic_sigma1_clim = (-500.0, 500.0)
            viscous_delta_clim = (-50.0, 50.0)
            elastic_delta_clim = (-300.0, 300.0)
        else:
            viscous_sigma1_clim = tuple(float(value) for value in sigma1_clim)
            elastic_sigma1_clim = tuple(float(value) for value in sigma1_clim)
            viscous_delta_clim = viscous_sigma1_clim
            elastic_delta_clim = elastic_sigma1_clim

        event_start = t_drainage[event_ids] - 0.5 * t_duration[event_ids]
        event_end = t_drainage[event_ids] + 0.5 * t_duration[event_ids]
        tmin = float(np.min(event_start) - 0.5 * margin_days)
        tmax = float(np.max(event_end) + margin_days)
        timeseries_xmin = tmin if timeseries_tmin_days is None else float(timeseries_tmin_days)
        timeseries_xmax = tmax if timeseries_tmax_days is None else float(timeseries_tmax_days)
        stress_window_min = min(tmin, timeseries_xmin)
        stress_window_max = max(tmax, timeseries_xmax)
        stress, raw = self._ensure_stress_window(
            required_fields,
            stress_window_min,
            stress_window_max,
            force_recompute=force_recompute,
        )
        ps = raw['ps']
        pd_r = raw['pd']
        pp = raw['pp']

        t_days = np.atleast_1d(stress['t_days']).astype(float)
        lake_x_km = np.atleast_1d(stress['lake_x_km']).astype(float)
        lake_y_km = np.atleast_1d(stress['lake_y_km']).astype(float)
        sigma1_viscous = np.atleast_2d(stress['sigma1_lake']).astype(float) / 1e3
        sigma1_elastic = np.atleast_2d(stress['sigma1_elastic_lake']).astype(float) / 1e3
        if sigma1_viscous.shape[0] <= np.max(event_ids) or sigma1_elastic.shape[0] <= np.max(event_ids):
            raise ValueError('Lake stress array is shorter than selected_events. Recompute the stress file.')
        plot_mask = (t_days >= tmin) & (t_days <= tmax)
        if not np.any(plot_mask):
            raise ValueError(self._empty_window_error_message(tmin, tmax, t_days))
        timeseries_plot_mask = (t_days >= timeseries_xmin) & (t_days <= timeseries_xmax)
        if not np.any(timeseries_plot_mask):
            raise ValueError(self._empty_window_error_message(timeseries_xmin, timeseries_xmax, t_days))

        def compute_panel_limits(panel_values):
            finite_values = panel_values[np.isfinite(panel_values)]
            if finite_values.size == 0:
                raise ValueError('Selected events contain no finite stress values in the requested plotting window.')
            y_min = float(np.min(finite_values))
            y_max = float(np.max(finite_values))
            y_pad = 0.10 * max(y_max - y_min, 20.0)
            return (y_min - y_pad, y_max + y_pad)

        def build_stress_panels():
            if stress_panel_plan is None:
                representative_event_ids, moment_times, group_event_ids = self._select_representative_drainage_groups(
                    event_ids,
                    t_drainage,
                    n_groups=n_stress_panels,
                )
                return [
                    {
                        'time_days': float(moment_time),
                        'event_ids': list(event_group),
                        'highlight_event_ids': [int(event_id) for event_id in list(event_group)[:2]],
                        'title_text': f'{float(moment_time):.1f} d',
                    }
                    for rep_event_id, moment_time, event_group in zip(
                        representative_event_ids,
                        moment_times,
                        group_event_ids,
                    )
                ]

            selected_event_ids = set((event_ids + 1).tolist())
            panels = []
            for panel in stress_panel_plan:
                panel_event_ids = [int(value) for value in panel.get('event_ids', [])]
                if panel_event_ids:
                    missing_event_ids = sorted(set(panel_event_ids) - selected_event_ids)
                    if missing_event_ids:
                        raise ValueError(
                            'stress_panel_plan event_ids must be a subset of selected_events. '
                            f'Missing: {missing_event_ids}'
                        )
                    group_zero_based = np.asarray(panel_event_ids, dtype=int) - 1
                    if np.any(group_zero_based < 0) or np.any(group_zero_based >= len(t_drainage)):
                        raise ValueError('stress_panel_plan contains an event id outside the available range.')
                else:
                    group_zero_based = np.array([], dtype=int)

                if panel.get('time_days') is not None:
                    moment_time = float(panel['time_days'])
                elif group_zero_based.size:
                    group_times = t_drainage[group_zero_based]
                    representative_local_index = int(np.argmin(np.abs(group_times - np.median(group_times))))
                    moment_time = float(group_times[representative_local_index])
                else:
                    raise ValueError('Each stress_panel_plan entry must define time_days or event_ids.')

                if group_zero_based.size:
                    sort_order = np.argsort(np.abs(t_drainage[group_zero_based] - moment_time))
                    highlight_zero_based = group_zero_based[sort_order[:min(2, group_zero_based.size)]]
                else:
                    sort_order = np.argsort(np.abs(t_drainage[event_ids] - moment_time))
                    highlight_zero_based = event_ids[sort_order[:min(2, event_ids.size)]]

                title_text = f'{moment_time:.1f} d'

                panels.append(
                    {
                        'time_days': moment_time,
                        'event_ids': panel_event_ids,
                        'highlight_event_ids': (np.asarray(highlight_zero_based, dtype=int) + 1).tolist(),
                        'title_text': title_text,
                    }
                )
            return panels

        stress_panels = build_stress_panels()
        moment_times = [panel['time_days'] for panel in stress_panels]
        all_export_times = [reference_time_days] + moment_times
        xx, yy, all_frames = self._export_sigma1_frames(all_export_times, raw)
        reference_frame = all_frames[0]
        frames = all_frames[1:]
        gps_x_km, gps_y_km = self._gps_station_positions(raw)
        gps_cluster = self._cluster_gps_stations(gps_x_km, gps_y_km)

        selected_lake_x = lake_x_km[event_ids]
        selected_lake_y = lake_y_km[event_ids]
        # Zoom in 2x on the centre of the selected lakes.
        full_x_min = float(np.nanmin(np.asarray(xx, dtype=float)))
        full_x_max = float(np.nanmax(np.asarray(xx, dtype=float)))
        full_y_min = float(np.nanmin(np.asarray(yy, dtype=float)))
        full_y_max = float(np.nanmax(np.asarray(yy, dtype=float)))
        center_x = 0.75 * (float(np.min(selected_lake_x)) + float(np.max(selected_lake_x)))
        center_y = 0.75 * (float(np.min(selected_lake_y)) + float(np.max(selected_lake_y)))
        half_x = 0.375 * (full_x_max - full_x_min)
        half_y = 0.375 * (full_y_max - full_y_min)
        xlim = (
            max(center_x - half_x, full_x_min),
            min(center_x + half_x, full_x_max),
        )
        ylim = (
            max(center_y - half_y, full_y_min),
            min(center_y + half_y, full_y_max),
        )
        map_box_aspect = (ylim[1] - ylim[0]) / max(xlim[1] - xlim[0], np.finfo(float).eps)
        gps_region_mask = (
            (gps_x_km >= xlim[0]) & (gps_x_km <= xlim[1]) &
            (gps_y_km >= ylim[0]) & (gps_y_km <= ylim[1])
        )

        viscous_limits = compute_panel_limits(sigma1_viscous[event_ids][:, timeseries_plot_mask])
        elastic_limits = compute_panel_limits(sigma1_elastic[event_ids][:, timeseries_plot_mask])
        colors = [self.EVENT_COLORS[i % len(self.EVENT_COLORS)] for i in range(len(event_ids))]
        event_color_by_id = {
            int(event_idx): color for event_idx, color in zip(event_ids, colors)
        }
        lake_fill_color = '#6fadd6'
        lake_edge_color = '#245c86'
        highlighted_group_fill = '#FFD700'

        print('Reference time [d]:', reference_time_days)
        print('Requested stress panel times [d]:', [round(value, 2) for value in moment_times])
        print(
            'Stress panel groups:',
            [self._event_labels(panel['event_ids'], lake_label_by_event, one_based=True) for panel in stress_panels],
        )
        print('Exported stress panel frame times [d]:', [round(float(frame['t_days']), 2) for frame in frames])

        plt.style.use('classic')
        plt.rcParams.update({
            'font.family': 'DejaVu Sans',
            'mathtext.fontset': 'dejavusans',
            'font.size': 8,
            'figure.dpi': 200,
            'savefig.dpi': 600,
        })
        n_map_rows = 1 + len(stress_panels)  # reference row + stress-panel rows
        n_rows = 1 + n_map_rows               # timeseries + map rows
        ts_ratio = 0.45
        fig_height = 2.0 * (ts_ratio + n_map_rows)
        fig = plt.figure(figsize=(8, fig_height), dpi=200)
        gs = fig.add_gridspec(
            n_rows, 2,
            hspace=0, wspace=0.25,
            height_ratios=[ts_ratio] + [1.0] * n_map_rows,
        )
        axes = np.empty((n_rows, 2), dtype=object)
        for r in range(n_rows):
            for c in range(2):
                axes[r, c] = fig.add_subplot(gs[r, c])
        fig.subplots_adjust(hspace=0, wspace=0.25)
        # Manually shift timeseries row upward to create gap
        for c in range(2):
            pos = axes[0, c].get_position()
            axes[0, c].set_position([pos.x0, pos.y0 + 0.02, pos.width, pos.height])

        timeseries_axes = axes[0, :]
        _DASH_ON_OFF = (4, 4)
        _DASH_CYCLE = sum(_DASH_ON_OFF)
        _DASH_PHASE_SHIFT = 0.2 * _DASH_CYCLE
        _N_EVENTS = max(len(event_ids), 1)
        top_panel_specs = [
            (timeseries_axes[0], r'$\bf{(a)}$ Viscous $\sigma_1$', sigma1_viscous, viscous_limits),
            (timeseries_axes[1], r'$\bf{(b)}$ Elastic (blister) $\sigma_1$', sigma1_elastic, elastic_limits),
        ]

        for ax, panel_title, panel_data, panel_limits in top_panel_specs:
            for ie, (color, event_idx) in enumerate(zip(colors, event_ids)):
                start = float(t_drainage[event_idx] - self.Z95 * t_duration[event_idx])
                end = float(t_drainage[event_idx] + self.Z95 * t_duration[event_idx])
                center = float(t_drainage[event_idx])
                phase = (_DASH_PHASE_SHIFT + ie * _DASH_CYCLE / _N_EVENTS) % _DASH_CYCLE
                ax.axvspan(start, end, color='#d9edf7', alpha=0.55, lw=0, zorder=0)
                ax.axvline(
                    center,
                    color=color,
                    ls=(phase, _DASH_ON_OFF),
                    lw=1.4,
                    alpha=1.0,
                    zorder=3,
                    dash_capstyle='butt',
                    solid_capstyle='butt',
                )
                ax.plot(
                    t_days[timeseries_plot_mask],
                    panel_data[event_idx, timeseries_plot_mask],
                    color=color,
                    lw=1.0,
                    zorder=2,
                )
            self._add_panel_annotation(ax, panel_title, fontsize=7)
            ax.set_xlim(timeseries_xmin, timeseries_xmax)
            ax.set_ylim(*panel_limits)
            ax.set_xlabel(r'day of year', fontsize=5, labelpad=1)
            ax.set_ylabel(r'$\sigma_1$ [kPa]')
            ax.grid(True, alpha=0.22)
            ax.tick_params(labelsize=6)
            # Use plain day numbers instead of offset notation like +1.9e2
            fmt = ScalarFormatter(useOffset=False)
            fmt.set_scientific(False)
            ax.xaxis.set_major_formatter(fmt)

        legend_handles = []
        for event_idx in self._sort_events_by_lake_label(event_ids, lake_label_by_event):
            legend_handles.append(
                Line2D([0], [0], color=event_color_by_id[event_idx], lw=1.2, label=lake_label_by_event[event_idx])
            )
        legend_handles.append(
            Line2D(
                [0], [0], color='#333333', lw=1.0, ls=(0, _DASH_ON_OFF),
                dash_capstyle='butt', solid_capstyle='butt', label='drainage'
            )
        )
        main_legend = timeseries_axes[0].legend(
            handles=legend_handles,
            loc='lower right',
            bbox_to_anchor=(1, 0),
            ncol=2,
            frameon=True,
            framealpha=0.9,
            edgecolor='none',
            fontsize=4,
            columnspacing=0.8,
            labelspacing=0.3,
            handletextpad=0.4,
            handlelength=3.0,
        )
        timeseries_axes[0].add_artist(main_legend)
        timeseries_axes[0].legend(
            handles=[
                Patch(
                    facecolor='#d9edf7',
                    edgecolor='none',
                    alpha=0.55,
                    label='drainage\nwindow',
                )
            ],
            loc='upper right',
            bbox_to_anchor=(1, 1.0),
            ncol=1,
            frameon=True,
            framealpha=0.9,
            edgecolor='none',
            fontsize=5,
            handletextpad=0.4,
            handlelength=1.5,
        )

        # --- Reference stress row (row 1): t = reference_time_days ---
        angle_key_map = {'sigma1': 't1', 'sigma1_elastic': 't1_elastic', 'sigma1_total': 't1_total'}
        quiver_scale_map = {'sigma1': 0.012, 'sigma1_elastic': 0.003, 'sigma1_total': 0.003}
        for col_index, (field_name, field_label, clim) in enumerate([
            ('sigma1', 'Viscous', viscous_sigma1_clim),
            ('sigma1_elastic', 'Elastic', elastic_sigma1_clim),
        ]):
            ax = axes[1, col_index]
            stress_field = reference_frame[field_name]
            pcm = ax.pcolormesh(
                xx, yy, stress_field,
                cmap='viridis', shading='gouraud',
                vmin=float(clim[0]), vmax=float(clim[1]),
            )
            self._add_panel_colorbar(
                ax,
                pcm,
                r'$\sigma_1$ (kPa)',
                label_side='right' if col_index == 0 else 'right',
                labelpad=2 if col_index == 0 else None,
            )
            # if field_name != 'sigma1':
            #     self._add_sigma1_quiver(
            #         ax, xx, yy,
            #         reference_frame[field_name],
            #         reference_frame[angle_key_map[field_name]],
            #         xlim, ylim,
            #         scale_factor=quiver_scale_map[field_name],
            #         linewidth=0.15,
            #     )

            for lake_idx, event_idx in enumerate(event_ids):
                ax.scatter(
                    lake_x_km[event_idx], lake_y_km[event_idx],
                    s=24, facecolors=colors[lake_idx],
                    edgecolors='w', linewidths=0.6, zorder=6,
                )

            ref_panel_label = chr(ord('c') + col_index)
            ref_title = rf'$\bf{{({ref_panel_label})}}$ {field_label} {reference_time_days:.0f} d'
            self._add_panel_annotation(ax, ref_title, x=0.02, fontsize=7)
            ax.set_xlim(*xlim)
            ax.set_ylim(*ylim)
            ax.set_aspect('equal', adjustable='box')
            ax.set_xlabel('')
            ax.tick_params(labelbottom=False)
            if col_index == 0:
                ax.set_ylabel('y (km)')
            else:
                ax.tick_params(labelleft=False)
            ax.tick_params(labelsize=7)

        stress_legend_handles = [
            Line2D(
                [0], [0], linestyle='None', marker='o', markersize=5,
                markerfacecolor=event_color_by_id[event_idx], markeredgecolor='w', markeredgewidth=0.6,
                label=lake_label_by_event[event_idx],
            )
            for event_idx in self._sort_events_by_lake_label(event_ids, lake_label_by_event)
        ]
        stress_legend_handles.append(
            Line2D(
                [0], [0], linestyle='None', marker='*', markersize=10,
                markerfacecolor='white', markeredgecolor='k',
                label='Drainage',
            )
        )

        panel_labels = [chr(ord('e') + index) for index in range(2 * len(stress_panels))]
        for row_index, (frame, panel_spec) in enumerate(zip(frames, stress_panels), start=2):
            event_group = panel_spec['event_ids']
            highlight_event_ids = panel_spec['highlight_event_ids']
            for col_index, (field_name, field_label, clim) in enumerate(
                [
                    ('sigma1', 'Viscous', viscous_delta_clim),
                    ('sigma1_elastic', 'Elastic', elastic_delta_clim),
                ]
            ):
                ax = axes[row_index, col_index]
                stress_field = frame[field_name] - reference_frame[field_name]
                pcm = ax.pcolormesh(
                    xx,
                    yy,
                    stress_field,
                    cmap='bwr',
                    shading='gouraud',
                    vmin=float(clim[0]),
                    vmax=float(clim[1]),
                )
                self._add_panel_colorbar(
                    ax,
                    pcm,
                    r'$\Delta\sigma_1$ (kPa)',
                    label_side='right' if col_index == 0 else 'right',
                    labelpad=2 if col_index == 0 else None,
                )
                # if field_name != 'sigma1':
                #     self._add_sigma1_quiver(
                #         ax, xx, yy,
                #         frame[field_name],
                #         frame[angle_key_map[field_name]],
                #         xlim, ylim,
                #         scale_factor=quiver_scale_map[field_name],
                #         linewidth=0.15,
                #     )

                for lake_idx, event_idx_lake in enumerate(event_ids):
                    ax.scatter(
                        lake_x_km[event_idx_lake],
                        lake_y_km[event_idx_lake],
                        s=24,
                        facecolors=colors[lake_idx],
                        edgecolors='w',
                        linewidths=0.6,
                        zorder=6,
                    )

                if event_group:
                    group_zero_based = np.asarray(event_group, dtype=int) - 1
                    ax.scatter(
                        lake_x_km[group_zero_based],
                        lake_y_km[group_zero_based],
                        s=36,
                        facecolors=highlighted_group_fill,
                        edgecolors='k',
                        linewidths=0.5,
                        zorder=7,
                    )

                if highlight_event_ids:
                    highlight_zero_based = np.asarray(highlight_event_ids, dtype=int) - 1
                    ax.scatter(
                        lake_x_km[highlight_zero_based],
                        lake_y_km[highlight_zero_based],
                        s=90,
                        marker='*',
                        c='white',
                        edgecolors='k',
                        linewidths=0.5,
                        zorder=8,
                    )
                    for event_index in highlight_zero_based:
                        self._add_marker_label(
                            ax,
                            lake_x_km[event_index],
                            lake_y_km[event_index],
                            lake_label_by_event[int(event_index)],
                        )

                panel_label = panel_labels[2 * (row_index - 2) + col_index]
                panel_title = rf'$\bf{{({panel_label})}}$ $\Delta$ {field_label} {panel_spec["title_text"]}'
                # annotation_x = 0.28 if row_index == 1 and col_index == 0 else 0.02
                annotation_x = 0.02
                self._add_panel_annotation(ax, panel_title, x=annotation_x, fontsize=7)
                ax.set_xlim(*xlim)
                ax.set_ylim(*ylim)
                ax.set_aspect('equal', adjustable='box')
                if row_index == n_rows - 1:
                    ax.set_xlabel('x (km)')
                else:
                    ax.set_xlabel('')
                    ax.tick_params(labelbottom=False)
                if col_index == 0:
                    ax.set_ylabel('y (km)')
                else:
                    ax.tick_params(labelleft=False)
                ax.tick_params(labelsize=7)

        axes[1, 0].legend(
            handles=stress_legend_handles,
            loc='upper left',
            bbox_to_anchor=(0.0, 0.90),
            fontsize=5,
            frameon=True,
            ncol=1,
            handletextpad=0.2,
            markerscale=0.7,
        )

        if save_path is None:
            save_path = self.figure_dir / f'{self.casename}_selected_event_sigma1_and_stress_5x2.png'
        else:
            save_path = Path(save_path)
        fig.savefig(save_path, dpi=600, bbox_inches='tight', pad_inches=0.1)
        print(f'Saved to {save_path}')
        return fig

    def plot_gps_stress_timeseries(self, tmin_yr=0.45, tmax_yr=0.8, save_path=None, force_recompute=False):
        required_fields = [
            't_days', 'gps_x_km', 'gps_y_km', 'n_gps',
            'tauxx_gps', 'tauyy_gps', 'tauxy_gps',
            'sigma1_gps', 'sigma2_gps',
            'Txx_gps', 'Tyy_gps', 'Txy_gps',
            'tau_b_gps', 'taudx_gps', 'taudy_gps',
        ]
        stress, _, ps, pd_r, pp = self._shared_context(required_fields, force_recompute=force_recompute)

        t_days = np.atleast_1d(stress['t_days']).astype(float)
        gps_x_km = np.atleast_1d(stress['gps_x_km']).astype(float)
        gps_y_km = np.atleast_1d(stress['gps_y_km']).astype(float)
        n_gps = int(np.atleast_1d(stress['n_gps']).astype(int)[0])
        ps_t = float(ps['t'])
        pd_td = float(pd_r['td'])

        if tmin_yr is None or tmax_yr is None:
            plot_mask = np.ones_like(t_days, dtype=bool)
            tmin = float(t_days.min())
            tmax = float(t_days.max())
        else:
            tmin = 365.0 * tmin_yr
            tmax = 365.0 * tmax_yr
            plot_mask = (t_days >= tmin) & (t_days <= tmax)
        if not np.any(plot_mask):
            raise ValueError('No stress data fall within the requested plotting window.')

        t_drainage = self.mat_array(pp['t_drainage']) * ps_t / pd_td
        t_duration = self.mat_array(pp['t_duration']) * ps_t / pd_td
        lake_to_gps = self._lake_to_gps(pp, ps, pd_r, stress)

        plt.style.use('classic')
        plt.rcParams.update({
            'font.size': 8,
            'figure.dpi': 160,
            'mathtext.fontset': 'dejavusans',
            'axes.titlesize': 9,
            'axes.labelsize': 9,
        })
        fig, axes = plt.subplots(
            n_gps,
            len(self.STRESS_GROUPS),
            figsize=(8, max(2.4 * n_gps, 4.0)),
            sharex=True,
            constrained_layout=True,
        )
        if n_gps == 1:
            axes = axes[None, :]

        for ig in range(n_gps):
            event_idx = np.where(lake_to_gps == ig)[0] if lake_to_gps.size else np.array([], dtype=int)
            for jg, (title, entries, scale) in enumerate(self.STRESS_GROUPS):
                ax = axes[ig, jg]
                ax.set_title(title)

                for ie in event_idx:
                    half_width = self.Z95 * t_duration[ie]
                    ax.axvspan(
                        t_drainage[ie] - half_width,
                        t_drainage[ie] + half_width,
                        color='#d7ecff', alpha=0.7, lw=0,
                    )

                for field_name, label in entries:
                    series = np.atleast_2d(stress[field_name]).astype(float)[ig, plot_mask] / scale
                    ax.plot(t_days[plot_mask], series, lw=1.2, color=self.COLORS[field_name], label=label)

                ax.grid(True, alpha=0.25)
                ax.set_xlim(tmin, tmax)
                if ig == 0:
                    ax.legend(loc='best', fontsize=7, frameon=True, framealpha=0.85, edgecolor='none')
                if jg == 0:
                    ax.set_ylabel(
                        f'GPS {ig + 1}\n' + rf'$({gps_x_km[ig]:.1f},\;{gps_y_km[ig]:.1f})$ km',
                    )

            axes[ig, -1].annotate(
                f'{len(event_idx)} nearby events',
                xy=(1.0, 0.5),
                xytext=(6, 0),
                xycoords='axes fraction',
                textcoords='offset points',
                va='center',
                ha='left',
                fontsize=7,
                fontstyle='italic',
            )

        for ax in axes[-1, :]:
            ax.set_xlabel(r'day of year')

        if save_path is None:
            save_path = self.figure_dir / f'{self.casename}_gps_stress_timeseries.png'
        else:
            save_path = Path(save_path)
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f'Saved to {save_path}')

        for ig in range(n_gps):
            event_idx = np.where(lake_to_gps == ig)[0] if lake_to_gps.size else np.array([], dtype=int)
            event_summary = [self._format_drainage_summary(idx, t_drainage, t_duration) for idx in event_idx]
            print(f'GPS {ig + 1}: {len(event_idx)} assigned drainage events -> {event_summary}')
        return fig

    def plot_selected_lake_sigma1(self, selected_events, margin_days=6.0, save_path=None, force_recompute=False):
        required_fields = ['t_days', 'lake_x_km', 'lake_y_km', 'sigma1_lake', 'sigma1_elastic_lake']
        stress, raw = self._load_data(required_fields=required_fields, force_recompute=force_recompute)
        ps = raw['ps']
        pd_r = raw['pd']
        pp = raw['pp']
        t_drainage = self.mat_array(pp['t_drainage']) * float(ps['t']) / float(pd_r['td'])
        t_duration = self.mat_array(pp['t_duration']) * float(ps['t']) / float(pd_r['td'])

        event_ids = np.asarray(selected_events, dtype=int) - 1
        if np.any(event_ids < 0) or np.any(event_ids >= len(t_drainage)):
            raise ValueError('selected_events contains an event id outside the available range.')
        lake_label_by_event = self._selected_lake_label_map(event_ids, t_drainage)

        event_start = t_drainage[event_ids] - 0.5 * t_duration[event_ids]
        event_end = t_drainage[event_ids] + 0.5 * t_duration[event_ids]
        tmin = float(np.min(event_start) - 0.5*margin_days)
        tmax = float(np.max(event_end) + margin_days)
        stress, raw = self._ensure_stress_window(required_fields, tmin, tmax, force_recompute=force_recompute)

        t_days = np.atleast_1d(stress['t_days']).astype(float)
        lake_x_km = np.atleast_1d(stress['lake_x_km']).astype(float)
        lake_y_km = np.atleast_1d(stress['lake_y_km']).astype(float)
        sigma1_viscous = np.atleast_2d(stress['sigma1_lake']).astype(float) / 1e3
        sigma1_elastic = np.atleast_2d(stress['sigma1_elastic_lake']).astype(float) / 1e3
        if sigma1_viscous.shape[0] <= np.max(event_ids) or sigma1_elastic.shape[0] <= np.max(event_ids):
            raise ValueError('Lake stress array is shorter than selected_events. Recompute the stress file.')
        plot_mask = (t_days >= tmin) & (t_days <= tmax)
        if not np.any(plot_mask):
            raise ValueError(self._empty_window_error_message(tmin, tmax, t_days))

        def compute_panel_limits(panel_values):
            finite_values = panel_values[np.isfinite(panel_values)]
            if finite_values.size == 0:
                raise ValueError('Selected events contain no finite stress values in the requested plotting window.')
            y_min = float(np.min(finite_values))
            y_max = float(np.max(finite_values))
            y_pad = 0.10 * max(y_max - y_min, 20.0)
            return (y_min - y_pad, y_max + y_pad)

        viscous_limits = compute_panel_limits(sigma1_viscous[event_ids][:, plot_mask])
        elastic_limits = compute_panel_limits(sigma1_elastic[event_ids][:, plot_mask])

        plt.style.use('classic')
        plt.rcParams.update({
            'font.size': 9,
            'figure.dpi': 160,
            'mathtext.fontset': 'dejavusans',
            'axes.titlesize': 11,
            'axes.labelsize': 10,
            'legend.fontsize': 8,
        })
        fig, axes = plt.subplots(
            1,
            2,
            figsize=(12, 2.8),
            sharex=True,
            sharey=False,
            constrained_layout=True,
        )
        fig.subplots_adjust(wspace=0.18)
        colors = [self.EVENT_COLORS[i % len(self.EVENT_COLORS)] for i in range(len(event_ids))]
        event_color_by_id = {
            int(event_idx): color for event_idx, color in zip(event_ids, colors)
        }
        panel_info = [
            (r'$\bf{(a)}$ Viscous $\sigma_1$', sigma1_viscous, viscous_limits),
            (r'$\bf{(b)}$ Elastic (blister) $\sigma_1$', sigma1_elastic, elastic_limits),
        ]

        _DASH_ON_OFF = (4, 4)
        _DASH_CYCLE = sum(_DASH_ON_OFF)
        _DASH_PHASE_SHIFT = 0.2 * _DASH_CYCLE
        _N_EVENTS = max(len(event_ids), 1)

        for ax, (panel_title, panel_data, panel_limits) in zip(axes, panel_info):
            for ie, (color, event_idx) in enumerate(zip(colors, event_ids)):
                start = float(t_drainage[event_idx] - self.Z95 * t_duration[event_idx])
                end = float(t_drainage[event_idx] + self.Z95 * t_duration[event_idx])
                center = float(t_drainage[event_idx])
                ax.axvspan(start, end, color='#d9edf7', alpha=0.55, lw=0, zorder=0)
                phase = (_DASH_PHASE_SHIFT + ie * _DASH_CYCLE / _N_EVENTS) % _DASH_CYCLE
                ax.axvline(center, color=color,
                           ls=(phase, _DASH_ON_OFF), lw=2.2, alpha=1.0, zorder=3,
                           dash_capstyle='butt', solid_capstyle='butt')
                ax.plot(
                    t_days[plot_mask],
                    panel_data[event_idx, plot_mask],
                    color=color,
                    lw=1.8,
                    zorder=2,
                )
            ax.text(
                0.02, 0.95, panel_title,
                transform=ax.transAxes, fontsize=10,
                va='top', ha='left',
            )
            ax.set_xlim(tmin, tmax)
            ax.set_ylim(*panel_limits)
            ax.set_xlabel(r'day of year')
            ax.grid(True, alpha=0.22)

        axes[0].set_ylabel(r'$\sigma_1$ [kPa]')
        axes[1].set_ylabel(r'$\sigma_1$ [kPa]')

        legend_handles = []
        for event_idx in self._sort_events_by_lake_label(event_ids, lake_label_by_event):
            legend_handles.append(
                Line2D(
                    [0],
                    [0],
                    color=event_color_by_id[event_idx],
                    lw=2.0,
                    label=lake_label_by_event[event_idx],
                )
            )
        legend_handles.append(
            Line2D(
                [0],
                [0],
                color='#333333',
                lw=1.6,
                ls=(0, _DASH_ON_OFF),
                dash_capstyle='butt',
                solid_capstyle='butt',
                label='Drainage event',
            )
        )
        legend_handles.append(
            Patch(
                facecolor='#d9edf7',
                edgecolor='none',
                alpha=0.55,
                label=r'drainage',
            )
        )
        axes[0].legend(
            handles=legend_handles,
            loc='upper right',
            bbox_to_anchor=(1, 1),
            ncol=min(len(legend_handles), 2),
            frameon=True,
            framealpha=0.9,
            edgecolor='none',
            fontsize=8,
            # title=r'Selected drainage events',
            # title_fontsize=9,
        )
        axes[1].legend(
            handles=legend_handles,
            loc='upper right',
            bbox_to_anchor=(1, 1),
            ncol=min(len(legend_handles), 2),
            frameon=True,
            framealpha=0.9,
            edgecolor='none',
            fontsize=8,
            # title=r'Selected drainage events',
            # title_fontsize=9,
        )

        if save_path is None:
            save_path = self.figure_dir / f'{self.casename}_selected_event_sigma1_two_panel.png'
        else:
            save_path = Path(save_path)
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f'Saved to {save_path}')
        print(f'Blue shading spans +/- z95 * duration, so total shaded width = {2 * self.Z95:.6f} * duration.')

        for event_idx in event_ids:
            print(
                f'{self._format_drainage_summary(event_idx, t_drainage, t_duration, label=lake_label_by_event[event_idx])} -> '
                f'lake node {event_idx + 1} ({lake_x_km[event_idx]:.1f}, {lake_y_km[event_idx]:.1f}) km'
            )
        return fig
