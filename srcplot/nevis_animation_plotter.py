import os, glob, subprocess
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import TwoSlopeNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.animation import FuncAnimation, FFMpegWriter

class NevisIceAnimation:
    """Animation of blister thickness (left spatial) and ice velocity (right spatial) vs time.
    Python port of nevis_2d_animation_ice.m, class-based following NevisAnimation pattern."""

    def __init__(self, casename, tmin_yr=0.4, tmax_yr=0.8,
                 clusters_to_plot=None, frame_step=1, fps=5,
                 playback_days_per_second=None):
        self.casename = casename
        self.tmin_yr = tmin_yr
        self.tmax_yr = tmax_yr
        self.tmin = tmin_yr * 365
        self.tmax = tmax_yr * 365
        self.clusters_to_plot = clusters_to_plot or [1, 2, 3, 4, 5]
        self.FRAME_STEP = frame_step
        self.FPS = fps
        self.playback_days_per_second = playback_days_per_second
        self.sk = 10  # quiver skip
        self.cluster_colors = [
            (1.0, 0.45, 0.20),  # C1 bright orange
            (0.55, 1.0, 0.55),  # C2 light green
            (0.40, 0.75, 1.0),  # C3 light blue
            (0.85, 0.45, 1.0),  # C4 light purple
            (1.0, 0.85, 0.30),  # C5 bright yellow
        ]
        self.cluster_labels = ['S1', 'S2', 'S3', 'S4', 'S5']
        self._load_data()

    # ──────────────────────────────────────────────────────────────
    # Data loading
    # ──────────────────────────────────────────────────────────────
    def _load_data(self):
        mat_path = os.path.join('results', self.casename, self.casename + '.mat')
        raw = sio.loadmat(mat_path, squeeze_me=True)
        self.ps = raw['ps'];  self.pd = raw['pd'];  self.pp = raw['pp']
        self.oo = raw['oo'];  self.gg = raw['gg'];  self.aa = raw['aa']
        self.tt = raw['tt']
        self._convert_parameters()
        self._compute_time_series()
        self._compute_gps_clusters()
        self._load_gps_observations()
        self._find_frames()
        self._load_background_topography()
        self._export_ice_fields()

    def _load_background_topography(self):
        topo_path = os.path.join('data', 'nevis_regional', 'morlighem_for_nevis_140km.mat')
        raw_topo = sio.loadmat(topo_path, squeeze_me=True)
        topo = raw_topo['morlighem_for_nevis_140km']
        self.X_terrain = topo['X_m'].item()
        self.Y_terrain = topo['Y_m'].item()
        self.B_terrain = topo['B_m'].item()

    def _load_gps_observations(self):
        obs_path = os.path.join('data', 'station_timeseries_2022.mat')
        raw_obs = sio.loadmat(obs_path, squeeze_me=True, struct_as_record=False)
        stations = raw_obs['station_data']
        if np.ndim(stations) == 0:
            stations = np.array([stations], dtype=object)
        self.stations = np.atleast_1d(stations)

    def _convert_parameters(self):
        ps, pd, gg = self.ps, self.pd, self.gg
        self.ps_t   = float(ps['t']);    self.pd_td = float(pd['td'])
        self.pd_ty  = float(pd['ty'])
        self.ps_h   = float(ps['h']);    self.ps_hb = float(ps['hb'])
        self.ps_ub  = float(ps['u_b']); self.ps_phi = float(ps['phi'])
        self.ps_x   = float(ps['x']);    self.ps_S  = float(ps['S'])
        self.ps_hs  = float(ps['hs']);   self.ps_Q  = float(ps['Q'])
        self.ps_Q0  = float(ps['Q0']);   self.ps_m  = float(ps['m'])
        self.nI = int(gg['nI'].item());  self.nJ = int(gg['nJ'].item())
        self.xx = (self.ps_x / 1e3) * np.float64(gg['nx'].item())
        self.yy = (self.ps_x / 1e3) * np.float64(gg['ny'].item())

    def _compute_time_series(self):
        tt, n = self.tt, len(self.tt)
        ps_t, pd_td, pd_ty = self.ps_t, self.pd_td, self.pd_ty
        self.n_time = n
        self.t_ts = (ps_t / pd_td) * np.array([float(tt['t'][i]) for i in range(n)])
        self.Q_out_Q = self.ps_Q  * np.array([float(tt['Q_outQ'][i]) for i in range(n)])
        self.Q_out_q = self.ps_Q  * np.array([float(tt['Q_outq'][i]) for i in range(n)])
        self.Q_out_b = self.ps_Q0 * np.array([float(tt['Q_outb'][i]) for i in range(n)])
        self.E_ts    = self.ps_m * self.ps_x**2 * np.array([float(tt['E'][i]) for i in range(n)])
        self.N_ts    = self.ps_phi * np.array([float(tt['N'][i]) for i in range(n)])
        self.hs_ts   = self.ps_x**2 * self.ps_h * np.array([float(tt['hs'][i]) for i in range(n)])
        self.S_ts    = self.ps_x * self.ps_S * np.array([float(tt['S'][i]) for i in range(n)])
        self.U_ts    = self.ps_ub * np.array([float(tt['U'][i]) for i in range(n)])
        self.V_b     = self.ps_x**2 * self.ps_hb * np.array([float(tt['Vb'][i]) for i in range(n)])
        self.A_total = self.ps_x**2 * np.sum(
            np.float64(self.gg['Dx'].item()) * np.float64(self.gg['Dy'].item()))
        # GPS point time series
        pts_ni = self.oo['pts_ni'].item().astype(int) - 1
        self.pts_ni = pts_ni
        self.pts_hb = np.column_stack([tt['pts_hb'][i] for i in range(n)]) * self.ps_hb
        self.pts_hs = np.column_stack([tt['pts_hs'][i] for i in range(n)]) * self.ps_h
        self.pts_u  = np.column_stack([tt['pts_u'][i]  for i in range(n)]) * self.ps_ub * pd_ty
        self.t_mid = 0.5 * (self.t_ts[:-1] + self.t_ts[1:])
        dt_days = np.diff(self.t_ts)
        hb_hs = self.pts_hb + self.pts_hs
        self.dhdt_model = np.diff(hb_hs, axis=1) / dt_days[np.newaxis, :] * 365.0
        # drainage events
        self.t_drainage = np.asarray(self.pp['t_drainage'].item()) * ps_t / pd_td
        self.t_duration = np.asarray(self.pp['t_duration'].item()) * ps_t / pd_td
        # lake positions
        self.x_l = np.asarray(self.pp['x_l'].item()) * (self.ps_x / 1e3)
        self.y_l = np.asarray(self.pp['y_l'].item()) * (self.ps_x / 1e3)
        self.V_l = np.asarray(self.pp['V_l'].item())
        self.Q_lake_in = np.array([
            np.sum(self._lake_input(float(tt['t'][i]))) for i in range(n)
        ])

    def _compute_gps_clusters(self):
        ni_gps = self.pp['ni_gps'].item().astype(int)
        ni_gps_0 = ni_gps - 1
        xx_flat = (self.ps_x / 1e3) * np.float64(self.gg['nx'].item()).flatten(order='F')
        yy_flat = (self.ps_x / 1e3) * np.float64(self.gg['ny'].item()).flatten(order='F')
        gps_x = xx_flat[ni_gps_0]
        gps_y = yy_flat[ni_gps_0]
        # match GPS stations to pts_ni
        gps_rows, gps_matched = [], []
        for ig, g in enumerate(ni_gps):
            idx = np.where(self.oo['pts_ni'].item().astype(int) == g)[0]
            if len(idx) > 0:
                gps_rows.append(idx[0])
                gps_matched.append(ig)
        self.gps_rows = np.array(gps_rows)
        gps_matched = np.array(gps_matched)
        # k-means style cluster assignment
        cluster_centers = np.array([[0, 0], [20, -15], [40, -25]])
        gps_cluster = np.zeros(len(ni_gps), dtype=int)
        for ig in range(len(ni_gps)):
            dists = np.sum((cluster_centers - [gps_x[ig], gps_y[ig]])**2, axis=1)
            gps_cluster[ig] = np.argmin(dists) + 1
        # two manual exceptions
        gps_cluster[np.argmin(np.sqrt((gps_x - 17)**2 + (gps_y + 25)**2))] = 4
        gps_cluster[np.argmin(np.sqrt((gps_x - 5)**2 + (gps_y + 10)**2))] = 5
        # filter to matched only
        self.gps_cluster = gps_cluster[gps_matched]
        self.gps_x_km = gps_x[gps_matched]
        self.gps_y_km = gps_y[gps_matched]
        print(f'GPS: {len(ni_gps)} stations, {len(self.gps_rows)} matched, '
              f'clusters: {np.unique(self.gps_cluster)}')

    def _target_frame_dt_days(self):
        if self.playback_days_per_second is None:
            return None
        if self.playback_days_per_second <= 0:
            raise ValueError('playback_days_per_second must be positive')
        if self.FPS <= 0:
            raise ValueError('fps must be positive when adaptive playback is enabled')
        return self.playback_days_per_second / self.FPS

    def _resample_frame_indices(self, frame_indices, frame_times):
        target_dt = self._target_frame_dt_days()
        if target_dt is None or len(frame_indices) <= 2:
            return frame_indices, frame_times

        kept_indices = [int(frame_indices[0])]
        kept_times = [float(frame_times[0])]
        last_time = float(frame_times[0])

        for idx, time_day in zip(frame_indices[1:-1], frame_times[1:-1]):
            if float(time_day) - last_time >= target_dt:
                kept_indices.append(int(idx))
                kept_times.append(float(time_day))
                last_time = float(time_day)

        if kept_indices[-1] != int(frame_indices[-1]):
            kept_indices.append(int(frame_indices[-1]))
            kept_times.append(float(frame_times[-1]))

        return np.asarray(kept_indices, dtype=int), np.asarray(kept_times, dtype=float)

    def _find_frames(self):
        tspan_d = (self.ps_t / self.pd_td) * self.oo['t_span']
        frame_indices = np.where(
            (tspan_d >= self.tmin) & (tspan_d <= self.tmax))[0]
        if frame_indices.size == 0:
            raise ValueError(
                f'No frames found in [{self.tmin:.0f}, {self.tmax:.0f}] days for {self.casename}')

        if self.FRAME_STEP > 1:
            frame_indices = frame_indices[::self.FRAME_STEP]

        frame_times = np.asarray(tspan_d[frame_indices], dtype=float)
        n_before_resample = len(frame_indices)
        frame_indices, frame_times = self._resample_frame_indices(frame_indices, frame_times)

        self.frame_indices = frame_indices
        self.frame_times_days = frame_times
        filepath = os.path.join('results', self.casename, '')
        self.all_filenames = [
            os.path.join(filepath, f'{idx+1:04d}.mat') for idx in self.frame_indices]

        if self.playback_days_per_second is None:
            print(f'Found {len(self.frame_indices)} frames in '
                  f'[{self.tmin:.0f}, {self.tmax:.0f}] days')
        else:
            target_dt = self._target_frame_dt_days()
            print(
                f'Adaptive playback: kept {len(self.frame_indices)} of {n_before_resample} '
                f'frames with target {target_dt:.2f} days/frame '
                f'({self.playback_days_per_second:.2f} days/s at {self.FPS} fps)'
            )

    def _export_ice_fields(self):
        need_export = []
        for fn in self.all_filenames:
            try:
                keys = sio.whosmat(fn)
                names = [k[0] for k in keys]
                if 'ice_fields' not in names:
                    need_export.append(fn)
            except Exception:
                need_export.append(fn)
        if not need_export:
            print('All frames already have ice_fields exported.')
            return
        print(f'Exporting ice fields for {len(need_export)} frames via MATLAB...')
        sio.savemat('export_ice_fields.mat', {
            'filenames': np.array(need_export, dtype=object),
            'casename': self.casename,
        })
        matlab_path = "/Applications/MATLAB_R2025b.app/bin/matlab"
        subprocess.run([
            matlab_path, "-nodisplay", "-nosplash", "-nodesktop",
            "-batch", "run('src/nevis_export_ice_fields')"
        ], check=True)
        os.remove('export_ice_fields.mat')
        print('MATLAB export done.')

    # ──────────────────────────────────────────────────────────────
    # Per-frame loader
    # ──────────────────────────────────────────────────────────────
    def _load_frame(self, fn):
        raw_f = sio.loadmat(fn, squeeze_me=True)
        ice = raw_f['ice_fields']; vv = raw_f['vv']
        vel_scale = float(ice['vel_scale'])
        nout = ice['nout'].item().astype(int) - 1
        uxn = np.float64(ice['uxn'].item())
        vyn = np.float64(ice['vyn'].item())
        Un  = np.float64(ice['Un'].item())
        vux = vel_scale * uxn.reshape(self.nI, self.nJ, order='F')
        vuy = vel_scale * vyn.reshape(self.nI, self.nJ, order='F')
        vUU = vel_scale * Un; vUU[nout] = np.nan
        vUU = vUU.reshape(self.nI, self.nJ, order='F')
        hb = self.ps_hb * np.float64(vv['hb'].item()).reshape(
            self.nI, self.nJ, order='F')
        pb = self.ps_phi * np.float64(vv['pb'].item()).reshape(
            self.nI, self.nJ, order='F')
        phi = self.ps_phi * np.float64(vv['phi'].item()).reshape(
            self.nI, self.nJ, order='F')
        hs_flat = self.ps_h * np.float64(vv['hs'].item())
        hs_flat[nout] = np.nan
        hs = hs_flat.reshape(self.nI, self.nJ, order='F')
        return dict(vux=vux, vuy=vuy, vUU=vUU, hb=hb, pb=pb, hs=hs, phi=phi,
                    t_days=float(ice['t_days']),
                    t_nondim=float(vv['t'].item()))

    def _lake_input(self, t_nondim):
        t_dr = np.asarray(self.pp['t_drainage'].item())
        t_du = np.asarray(self.pp['t_duration'].item())
        V_l_ = np.asarray(self.pp['V_l'].item())
        Q0_  = float(self.ps['Q0'])
        return Q0_ * V_l_ / (np.sqrt(2 * np.pi) * t_du) * \
               np.exp(-0.5 * ((t_nondim - t_dr) / t_du)**2)

    # ──────────────────────────────────────────────────────────────
    # Figure setup
    # ──────────────────────────────────────────────────────────────
    def _setup_figure(self):
        plt.style.use('classic')
        plt.rcParams.update({
            'font.family': 'DejaVu Sans',
            'mathtext.fontset': 'dejavusans',
            'font.size': 14,
            'axes.titlesize': 16,
            'axes.labelsize': 13,
            'xtick.labelsize': 11,
            'ytick.labelsize': 11,
            'legend.fontsize': 10,
            'savefig.dpi': 200,
        })
        self.fig = plt.figure(figsize=(16, 8), dpi=100)
        # Outer grid: LEFT = spatial, RIGHT = time series
        gs_outer = self.fig.add_gridspec(1, 2, width_ratios=[1.4, 1], wspace=0.04)
        gs_left  = gs_outer[0].subgridspec(2, 1, hspace=0.10)
        gs_right = gs_outer[1].subgridspec(3, 1, hspace=0.35)
        # Left: spatial panels
        self.ax_hb  = self.fig.add_subplot(gs_left[0])
        self.ax_vel = self.fig.add_subplot(gs_left[1])
        # Right: time-series panels
        self.ax_a = self.fig.add_subplot(gs_right[0])
        self.ax_b = self.fig.add_subplot(gs_right[1])
        self.ax_c = self.fig.add_subplot(gs_right[2])

        self.ts_axes = [self.ax_a, self.ax_b, self.ax_c]
        # load reference frame
        self.ref = self._load_frame(self.all_filenames[0])
        # populate panels
        self._setup_blister_panel()
        self._setup_velocity_panel()
        self._setup_flux_panel()
        self._setup_speed_panel()
        self._setup_gps_panel()
        # Timestamp below spatial panels
        self.time_text = self.fig.text(
            0.15, 0.86, f't = {self.ref["t_days"]:.1f} d',
            fontsize=14, fontweight='bold', ha='center', va='center')
        # hide x labels for non-bottom right panels
        for ax_ in [self.ax_b]:
            ax_.tick_params(labelbottom=False)
        # Keep only persistent background artists when refreshing contours
        self.keep_collections = {self.bg_terrain, self.pc_hb, self.lake_sc, *self.gps_markers}
        self.vlines = [self.vline_a, self.vline_d, self.vline_e]
        self.fig.subplots_adjust(left=0.04, right=0.96, top=0.96, bottom=0.06)
        print(f'Initial frame: t = {self.ref["t_days"]:.1f} d')

    # ── helpers: drainage shading ──
    def _add_drainage_shading(self, ax, legend_label=None):
        for i_d in range(len(self.t_drainage)):
            label = legend_label if i_d == 0 else None
            ax.axvspan(
                self.t_drainage[i_d] - 0.5*self.t_duration[i_d],
                self.t_drainage[i_d] + 0.5*self.t_duration[i_d],
                color=(0.7, 0.85, 1.0), alpha=0.4, label=label
            )

    # ── (c) Fluxes ──
    def _setup_flux_panel(self):
        ax = self.ax_a; td = self.ref['t_days']
        self._add_drainage_shading(ax, legend_label='drainage period')
        q_out = self.Q_out_Q + self.Q_out_q
        ax.plot(self.t_ts, q_out, '-', color=(0.00, 0.45, 0.00), lw=2.4,
                label=r'$Q_\mathrm{out}$')
        ax.plot(self.t_ts, self.Q_lake_in, '-', color=(0.05, 0.35, 0.85), lw=2.4,
                label='lake drainage inflow')
        ax.plot(self.t_ts, self.E_ts, '-', color='k', lw=1.8,
                label=r'$Q_\mathrm{in}$')
        self.vline_a = ax.axvline(td, color='k', ls='--', lw=1)
        ax.set_xlim(self.tmin, self.tmax)
        ax.set_ylim(0, 8000)
        ax.set_xlabel(r'time $(d)$')
        ax.set_ylabel(r'$Q$ (m$^3$/s)')
        ax.legend(ncol=1, loc='upper right', fontsize=9)
        ax.text(0.02, 0.85, r'(c) flux', transform=ax.transAxes, fontsize=14, fontweight='bold')
        ax.grid(True, which='major', alpha=0.3)

    # ── (d) Ice speed + blister volume ──
    def _setup_speed_panel(self):
        ax = self.ax_b; td = self.ref['t_days']
        self._add_drainage_shading(ax)
        ax.plot(self.t_ts, self.U_ts * self.pd_ty, 'b-', lw=1.2,
                label=r'$U$ (m/yr)')
        self.vline_d = ax.axvline(td, color='k', ls='--', lw=1)
        ax2 = ax.twinx()
        ax2.plot(self.t_ts, self.V_b, 'r-', lw=1.2, label=r'$V_b$')
        ax.set_xlim(self.tmin, self.tmax)
        ax.set_xlabel(r'time $(d)$')
        ax.set_ylim(50, 150)
        ax.set_ylabel('U (m/yr)', color='b')
        ax.tick_params(axis='y', colors='b')
        ax2.set_ylabel(r'$V_b$ (m$^3$)', color='r')
        ax2.tick_params(axis='y', colors='r')
        ax.text(0.02, 0.85, r'(d) ice speed & blister volume', transform=ax.transAxes, fontsize=14, fontweight='bold')
        ax.grid(True, which='both', alpha=0.3)

    # ── (e) horizontal speed at GPS sites ──
    def _setup_gps_panel(self):
        ax = self.ax_c; td = self.ref['t_days']
        self._add_drainage_shading(ax)
        self.vline_e = ax.axvline(td, color='k', ls='--', lw=1)
        legend_handles = []
        legend_labels = []
        n_obs = len(self.stations)
        for kc in self.clusters_to_plot:
            mask = self.gps_cluster == kc
            if not np.any(mask):
                continue
            members = self.gps_rows[mask]
            obs_ids = np.where(mask)[0]
            obs_ids = obs_ids[obs_ids < n_obs]
            color = self.cluster_colors[kc-1]
            dark_color = tuple(np.clip(0.55 * np.array(color), 0, 1))
            model_u = np.nanmean(self.pts_u[members, :], axis=0)
            if np.any(np.isfinite(model_u)):
                h, = ax.plot(self.t_ts, model_u, color=dark_color, lw=2.6, label=f'{self.cluster_labels[kc-1]} $U$')
                legend_handles.append(h)
                legend_labels.append(self.cluster_labels[kc-1])
            # for obs_id in obs_ids:
            #     station = self.stations[obs_id]
            #     t_obs = np.asarray(getattr(station, 't_doy', np.array([])), dtype=float).ravel()
            #     u_obs = np.asarray(getattr(station, 'u_s', np.array([])), dtype=float).ravel()
            #     valid_u = np.isfinite(t_obs) & np.isfinite(u_obs)
            #     if np.any(valid_u):
            #         ax.scatter(t_obs[valid_u], u_obs[valid_u], s=4, color=color,
            #                    alpha=0.05, edgecolors='none')
        ax.set_xlim(self.tmin, self.tmax)
        ax.set_ylim(0, 300)
        ax.set_xlabel(r'time $(d)$')
        ax.set_ylabel(r'$U_{gps}$ (m/yr)')
        ax.text(0.02, 0.88, r'(e) GPS horizontal speed', transform=ax.transAxes,
                fontsize=14, fontweight='bold')
        ax.grid(True, which='both', alpha=0.3)
        if legend_handles:
            ax.legend(legend_handles, legend_labels, ncol=2, loc='upper right', fontsize=10)

    # ── Left top: blister thickness ──
    def _setup_blister_panel(self):
        ax = self.ax_hb; ref = self.ref
        terrain_norm = mcolors.PowerNorm(gamma=1, vmin=-1.0, vmax=1.0)
        self.bg_terrain = ax.pcolormesh(
            self.X_terrain/1e3, self.Y_terrain/1e3, self.B_terrain/1e3,
            norm=terrain_norm, cmap='terrain', alpha=0.55, zorder=0
        )
        n_cmap = 256
        bwr_arr = np.vstack([
            np.column_stack([np.linspace(0,1,n_cmap),
                             np.linspace(0,1,n_cmap), np.ones(n_cmap)]),
            np.column_stack([np.ones(n_cmap),
                             np.linspace(1,0,n_cmap), np.linspace(1,0,n_cmap)]),
        ])
        self.bwr_cmap = mcolors.ListedColormap(bwr_arr)
        self.pc_hb = ax.pcolormesh(self.xx, self.yy, ref['hb'] + ref['hs'],
                                   cmap=self.bwr_cmap, shading='auto',
                                   vmin=-1.0, vmax=1.0, zorder=1, alpha=0.85)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="3%", pad=0.08)
        cb = plt.colorbar(self.pc_hb, cax=cax)
        cb.set_label(r'Uplift (m)')
        self.ct_hb = ax.contour(
            self.xx, self.yy, ref['phi'], levels=10,
            colors='k', linewidths=0.5, zorder=2
        )
        # lake markers
        self.Q_peak = np.max(float(self.ps['Q0']) * self.V_l / (
            np.sqrt(2*np.pi) * np.asarray(self.pp['t_duration'].item())))
        lake_flux = self._lake_input(ref['t_nondim'])
        norm_sz = lake_flux / (self.Q_peak + 1e-30)
        sz = 15 + 395 * norm_sz**0.3
        self.lake_sc = ax.scatter(self.x_l, self.y_l, s=sz,
                                  c=[(0.2,0.5,0.9)]*len(self.x_l),
                                  edgecolors='k', zorder=5, alpha=1.0, linewidths=0.5, label='lake')
        # annotate each drainage event marker with its event id
        self.lake_labels = []
        for i, (x, y) in enumerate(zip(self.x_l, self.y_l), start=1):
            txt = ax.annotate(
                f'E{i}', xy=(x, y), xytext=(3, 3), textcoords='offset points',
                fontsize=5, ha='left', va='bottom', zorder=7, color='k',
                bbox=dict(boxstyle='round,pad=0.12', fc='white', ec='none', alpha=0.7)
            )
            self.lake_labels.append(txt)

        # GPS station markers
        self.gps_markers = []
        for kc in [1, 2, 3, 4, 5]:
            mask = self.gps_cluster == kc
            if not np.any(mask):
                continue
            marker = ax.scatter(self.gps_x_km[mask], self.gps_y_km[mask],
                                s=30, marker='s', zorder=6,
                                color=self.cluster_colors[kc-1], edgecolors='k',
                                linewidths=1.0, label=self.cluster_labels[kc-1])
            self.gps_markers.append(marker)
        ax.set_xlim(np.nanmin(self.xx), np.nanmax(self.xx))
        ax.set_ylim(np.nanmin(self.yy), np.nanmax(self.yy))
        ax.set_aspect('equal')
        ax.set_ylabel('y (km)')
        ax.tick_params(axis='both', labelsize=10)
        ax.legend(fontsize=10, loc='center left', ncol=1)
        ax.text(0.02, 0.95, r'(a) $h_b$ and $\phi$ contour', transform=ax.transAxes,
                ha='left', va='top', fontsize=14, fontweight='bold')

    # ── Left bottom: ice velocity ──
    def _setup_velocity_panel(self):
        ax = self.ax_vel; ref = self.ref; sk = self.sk
        self.pc_vel = ax.pcolormesh(self.xx, self.yy, ref['vUU'],
                                    cmap='viridis', shading='gouraud',
                                    vmin=0, vmax=200)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="3%", pad=0.08)
        cb = plt.colorbar(self.pc_vel, cax=cax)
        cb.set_label('U (m/yr)')
        self.quiv = ax.quiver(self.xx[::sk, ::sk], self.yy[::sk, ::sk],
                              ref['vux'][::sk, ::sk], ref['vuy'][::sk, ::sk],
                              scale=4000, width=0.003, color='k', alpha=0.6)
        ax.set_xlim(np.nanmin(self.xx), np.nanmax(self.xx))
        ax.set_ylim(np.nanmin(self.yy), np.nanmax(self.yy))
        ax.set_aspect('equal')
        ax.set_ylabel('y (km)')
        ax.set_xlabel('x (km)')
        ax.tick_params(axis='both', labelsize=10)
        ax.text(0.02, 0.95, r'(b) ice speed U (m/yr)', transform=ax.transAxes,
                ha='left', va='top', fontsize=14, fontweight='bold')

    # ──────────────────────────────────────────────────────────────
    # Animation
    # ──────────────────────────────────────────────────────────────
    def _update(self, frame_i):
        f = self._load_frame(self.all_filenames[frame_i])
        td = f['t_days']
        # time markers
        for vl in self.vlines:
            vl.set_xdata([td, td])
        self.time_text.set_text(f't = {td:.1f} d')
        # blister spatial (2D for gouraud, mpl >=3.10)
        self.pc_hb.set_array(f['hb'] + f['hs'])
        # refresh contours
        for c in self.ax_hb.collections[:]:
            if c not in self.keep_collections:
                c.remove()
        self.ct_hb = self.ax_hb.contour(self.xx, self.yy, f['phi'], levels=10,
                                        colors='k', linewidths=0.5, zorder=2)
        # lake marker sizes
        lf = self._lake_input(f['t_nondim'])
        self.lake_sc.set_sizes(15 + 395 * (lf / (self.Q_peak + 1e-30))**0.3)
        # velocity spatial (2D for gouraud)
        self.pc_vel.set_array(f['vUU'])
        sk = self.sk
        self.quiv.set_UVC(f['vux'][::sk, ::sk], f['vuy'][::sk, ::sk])
        if frame_i % 20 == 0:
            print(f'  frame {frame_i}/{len(self.all_filenames)}, t={td:.1f} d')
        return []

    def create_animation(self, output_filename=None):
        self._setup_figure()
        os.makedirs('./results/videos', exist_ok=True)
        if output_filename is None:
            output_filename = f'./results/videos/{self.casename}_python.mp4'
        anim = FuncAnimation(self.fig, self._update,
                             frames=len(self.all_filenames), blit=False)
        writer = FFMpegWriter(fps=self.FPS,
                              metadata=dict(title=self.casename),
                              bitrate=20000)
        anim.save(output_filename, writer=writer)
        print(f'Video saved: {output_filename}')
        plt.close(self.fig)
        return anim
