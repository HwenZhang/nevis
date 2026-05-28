import os
import numpy as np
import scipy.io as sio
import scipy.ndimage as ndi
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
from matplotlib.artist import Artist
from matplotlib.colors import Normalize, TwoSlopeNorm
from matplotlib.offsetbox import AnchoredOffsetbox, HPacker, TextArea
from mpl_toolkits.axes_grid1 import make_axes_locatable
from cycler import cycler
import subprocess

# Enable retina / 4K rendering in notebook
try:
    from IPython.display import set_matplotlib_formats
    set_matplotlib_formats('retina')
except Exception:
    pass

class NevisIceFlowPlotter:
    """
    Generate multi-panel snapshots of ice velocity and a configurable right-column field
    at user-specified time slices, using MATLAB for field computation
    and Python/matplotlib for plotting.

    Layout:
        Left column  — velocity   (row 1 = absolute, rows 2+ = % change)
        Right column — stress or uplift, controlled by ``right_column_mode``
    """

    def __init__(self, casename, results_dir='results',
                 matlab_path="/Applications/MATLAB_R2026a.app/bin/matlab",
                 right_column_mode='uplift',
                 show_quivers=True,
                 flux_quiver_reference=None):
        self.casename = casename
        self.results_dir = results_dir
        self.matlab_path = matlab_path
        self.right_column_mode = right_column_mode.lower()
        self.show_quivers = show_quivers
        if self.right_column_mode not in {'stress', 'uplift'}:
            raise ValueError("right_column_mode must be 'stress' or 'uplift'")
        # --- time slices (days) to plot ---
        self.tslice_d = np.array([100, 150, 180, 210])
        # --- plotting limits ---
        self.vel_clim = [0, 200]         # m/yr   (absolute)
        self.dvel_clim = [-100, 100]     # %       (relative change)
        self.sigma1_clim = [0, 200]      # kPa    (absolute)
        self.dsigma1_clim = [-100, 100]  # %       (relative change)
        self.uplift_clim = [0.0, 1.0]    # m       (absolute uplift hb+hs)
        self.quiver_skip = 8
        self.flux_quiver_skip = 6
        self.flux_quiver_max_length = 4.0
        self.flux_quiver_clip_percentile = 99.0
        self.flux_quiver_reference = flux_quiver_reference
        # containers
        self.data_loaded = False

    # ------------------------------------------------------------------
    # Data loading
    # ------------------------------------------------------------------
    def load_data(self):
        """Load base .mat, call MATLAB to export ice fields, read back."""
        mat_path = os.path.join(self.results_dir, self.casename,
                                self.casename + '.mat')
        raw = sio.loadmat(mat_path, squeeze_me=True)
        self.ps = raw['ps']; self.pd = raw['pd']; self.pp = raw['pp']
        self.oo = raw['oo']; self.gg = raw['gg']; self.aa = raw['aa']
        ps_t = np.float64(self.ps['t']); pd_td = np.float64(self.pd['td'])
        ps_x = np.float64(self.ps['x'])
        self.ps_phi = float(self.ps['phi'])
        self.xx = ps_x / 1e3 * np.float64(self.gg['nx'].item())
        self.yy = ps_x / 1e3 * np.float64(self.gg['ny'].item())
        self.nI = int(self.gg['nI'].item())
        self.nJ = int(self.gg['nJ'].item())
        self.ps_hb = float(self.ps['hb'])
        self.ps_h = float(self.ps['h'])
        # --- GPS station positions and clustering ---
        ni_gps = self.pp['ni_gps'].item().astype(int) - 1  # 0-based
        xx_flat = (ps_x / 1e3) * np.float64(self.gg['nx'].item()).flatten(order='F')
        yy_flat = (ps_x / 1e3) * np.float64(self.gg['ny'].item()).flatten(order='F')
        self.gps_x = xx_flat[ni_gps]
        self.gps_y = yy_flat[ni_gps]
        # cluster by proximity to 3 centers (km)
        cluster_centers = np.array([[0, 0], [20, -15], [40, -25]])
        self.gps_cluster = np.zeros(len(ni_gps), dtype=int)
        for ig in range(len(ni_gps)):
            dists = np.sum((cluster_centers - [self.gps_x[ig], self.gps_y[ig]])**2, axis=1)
            self.gps_cluster[ig] = np.argmin(dists) + 1
        # two exceptions
        dist_c2 = np.sqrt((self.gps_x - 17)**2 + (self.gps_y + 25)**2)
        self.gps_cluster[np.argmin(dist_c2)] = 4
        dist_c1 = np.sqrt((self.gps_x - 5)**2 + (self.gps_y + 10)**2)
        self.gps_cluster[np.argmin(dist_c1)] = 5
        self.cluster_colors = [
            (1.0, 0.45, 0.20),    # S1 bright orange
            (0.55, 1.0, 0.55),    # S2 light green
            (0.40, 0.75, 1.0),    # S3 light blue
            (0.85, 0.45, 1.0),    # S4 light purple
            (1.0, 0.85, 0.30),    # S5 bright yellow
        ]
        self.cluster_labels = ['S1', 'S2', 'S3', 'S4', 'S5']
        # find frame indices closest to requested times
        tspan = ps_t / pd_td * self.oo['t_span']
        self.nums = np.array([np.argmin(np.abs(tspan - td)) + 1
                              for td in self.tslice_d])
        filenames = [os.path.join(self.results_dir, self.casename,
                                  f"{num:04d}.mat") for num in self.nums]
        print(f"Requesting MATLAB export for frames: {self.nums}")
        # save request for MATLAB helper
        sio.savemat('export_ice_fields.mat', {
            'filenames': np.array(filenames, dtype=object),
            'casename': self.casename
        })
        subprocess.run([
            self.matlab_path, "-nodisplay", "-nosplash", "-nodesktop",
            "-batch", "run('src/nevis_export_ice_fields')"
        ], check=True)
        os.remove('export_ice_fields.mat')
        # read exported fields
        self.frames = []
        for fn in filenames:
            raw_f = sio.loadmat(fn, squeeze_me=True)
            ice = raw_f['ice_fields']
            vv = raw_f['vv']
            f = {}
            f['t_days'] = float(ice['t_days'])
            vel_scale = float(ice['vel_scale'])
            stress_scale = float(ice['stress_scale'])
            nout = ice['nout'].item().astype(int) - 1   # 0-based
            uxn = np.float64(ice['uxn'].item())
            vyn = np.float64(ice['vyn'].item())
            Un = np.float64(ice['Un'].item())
            # velocity in m/yr, reshaped (Fortran order)
            vux = vel_scale * uxn.reshape(self.nI, self.nJ, order='F')
            vuy = vel_scale * vyn.reshape(self.nI, self.nJ, order='F')
            vUU_flat = vel_scale * Un
            vUU_flat[nout] = np.nan
            f['vux'] = vux
            f['vuy'] = vuy
            f['vUU'] = vUU_flat.reshape(self.nI, self.nJ, order='F')
            # stresses in kPa
            sigma1 = np.float64(ice['sigma1'].item())
            sigma2 = np.float64(ice['sigma2'].item())
            t1 = np.float64(ice['t1'].item())
            sigma1_kPa = stress_scale * sigma1
            sigma1_kPa[nout] = np.nan
            sigma2_kPa = stress_scale * sigma2
            sigma2_kPa[nout] = np.nan
            f['sigma1'] = sigma1_kPa.reshape(self.nI, self.nJ, order='F')
            f['sigma2'] = sigma2_kPa.reshape(self.nI, self.nJ, order='F')
            f['t1'] = t1.reshape(self.nI, self.nJ, order='F')
            # uplift and hydrology fields reconstructed from vv / exported data
            hb = self.ps_hb * np.float64(vv['hb'].item())
            hs = self.ps_h * np.float64(vv['hs'].item())
            phi = self.ps_phi * np.float64(vv['phi'].item())
            qnet = np.float64(ice['qnet'].item())
            hb[nout] = np.nan
            hs[nout] = np.nan
            phi[nout] = np.nan
            qnet[nout] = np.nan
            f['hb'] = hb.reshape(self.nI, self.nJ, order='F')
            f['hs'] = hs.reshape(self.nI, self.nJ, order='F')
            f['phi'] = phi.reshape(self.nI, self.nJ, order='F')
            f['qnet'] = qnet.reshape(self.nI, self.nJ, order='F')
            f['uplift'] = f['hb'] + f['hs']
            self.frames.append(f)

        flux_reference = self._resolve_flux_quiver_reference()
        for frame in self.frames:
            frame['flux_quiver'] = self._compute_flux_quiver(
                frame['phi'], frame['qnet'], flux_reference
            )

        self.data_loaded = True
        print(f"Loaded {len(self.frames)} snapshots")
        print(f"Using shared flux quiver reference: {flux_reference:.3e}")

    # ------------------------------------------------------------------
    # Plotting helpers
    # ------------------------------------------------------------------
    @staticmethod
    def _add_colorbar(ax, mappable, label):
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="3%", pad=0.08)
        cb = plt.colorbar(mappable, cax=cax)
        cb.set_label(label, fontsize=8)
        cb.ax.tick_params(labelsize=7)
        return cb

    @staticmethod
    def _add_panel_annotation(ax, text, x=0.0225, y=0.96):
        children: list[Artist]
        if ',' in text:
            prefix, suffix = text.split(',', 1)
            children = [
                TextArea(prefix, textprops={'fontsize': 6, 'fontweight': 'bold'}),
                TextArea(',' + suffix, textprops={'fontsize': 6}),
            ]
        else:
            children = [TextArea(text, textprops={'fontsize': 6, 'fontweight': 'bold'})]

        packed = HPacker(children=children, align='baseline', pad=0, sep=0)
        anchored = AnchoredOffsetbox(
            loc='upper left',
            child=packed,
            pad=0,
            borderpad=0,
            frameon=False,
            bbox_to_anchor=(x, y),
            bbox_transform=ax.transAxes,
        )
        anchored.set_clip_on(False)
        ax.add_artist(anchored)

    @staticmethod
    def _bwr_cmap():
        from matplotlib.colors import LinearSegmentedColormap
        colors = [(0, 0, 1), (1, 1, 1), (1, 0, 0)]
        return LinearSegmentedColormap.from_list('bwr_custom', colors, N=512)

    @staticmethod
    def _quiver_segment(ax, x, y, u, v, color, qkw):
        if x.size == 0:
            return
        ax.quiver(
            x, y, u, v,
            color=color,
            facecolor=color,
            edgecolor=color,
            linewidth=0.5,
            **qkw,
        )

    @staticmethod
    def _fill_nan_nearest(field):
        mask = ~np.isfinite(field)
        if not np.any(mask):
            return field.copy()
        if np.all(mask):
            return np.zeros_like(field, dtype=float)
        filled = field.copy()
        _, indices = ndi.distance_transform_edt(mask, return_indices=True)
        filled[mask] = filled[tuple(indices[:, mask])]
        return filled

    def _resolve_flux_quiver_reference(self):
        if self.flux_quiver_reference is not None:
            qref = float(self.flux_quiver_reference)
            if np.isfinite(qref) and qref > 0:
                return qref

        positive_flux = []
        for frame in self.frames:
            qfield = frame['qnet']
            mask = np.isfinite(qfield) & (qfield > 0)
            if np.any(mask):
                positive_flux.append(qfield[mask])

        if not positive_flux:
            self.flux_quiver_reference = 1.0
            return self.flux_quiver_reference

        positive_flux = np.concatenate(positive_flux)
        qref = np.nanpercentile(positive_flux, self.flux_quiver_clip_percentile)
        if not np.isfinite(qref) or qref <= 0:
            qref = np.nanmax(positive_flux)
        if not np.isfinite(qref) or qref <= 0:
            qref = 1.0
        self.flux_quiver_reference = float(qref)
        return self.flux_quiver_reference

    def _compute_flux_quiver(self, phi_field, qmag_field, qref):
        phi_filled = self._fill_nan_nearest(phi_field)
        xvec = np.asarray(self.xx[:, 0], dtype=float)
        yvec = np.asarray(self.yy[0, :], dtype=float)
        dphidx, dphidy = np.gradient(phi_filled, xvec, yvec)
        direction_norm = np.hypot(dphidx, dphidy)
        direction_norm[direction_norm < np.finfo(float).eps] = 1.0

        dirx = -dphidx / direction_norm
        diry = -dphidy / direction_norm

        if not np.isfinite(qref) or qref <= 0:
            qref = 1.0

        sk = self.flux_quiver_skip
        xq = self.xx[::sk, ::sk]
        yq = self.yy[::sk, ::sk]
        dirx = dirx[::sk, ::sk]
        diry = diry[::sk, ::sk]
        qsub = qmag_field[::sk, ::sk]
        valid = np.isfinite(phi_field[::sk, ::sk]) & np.isfinite(qsub) & (qsub > 0)
        scaled_length = self.flux_quiver_max_length * np.clip(qsub / qref, 0.0, 1.0)

        return {
            'x': xq[valid],
            'y': yq[valid],
            'u': (dirx * scaled_length)[valid],
            'v': (diry * scaled_length)[valid],
        }

    @staticmethod
    def _plot_flux_quiver(ax, flux_quiver):
        if flux_quiver['x'].size == 0:
            return
        ax.quiver(
            flux_quiver['x'],
            flux_quiver['y'],
            flux_quiver['u'],
            flux_quiver['v'],
            color='w',
            angles='xy',
            scale_units='xy',
            scale=0.3,
            width=0.0032,
            headwidth=4.5,
            headlength=6.0,
            headaxislength=5.2,
            pivot='tail',
            alpha=0.9,
            zorder=4,
        )

    def _plot_gps_stations(self, ax):
        """Overlay GPS stations colored by cluster."""
        for cl in range(1, 6):
            mask = self.gps_cluster == cl
            if not np.any(mask):
                continue
            ax.scatter(
                self.gps_x[mask], self.gps_y[mask],
                s=10, marker='s', zorder=5,
                color=self.cluster_colors[cl - 1],
                edgecolors='k', linewidths=0.3,
                label=self.cluster_labels[cl - 1],
            )

    def _plot_velocity_abs(self, ax, frame, annotation_text):
        """Absolute speed + quiver arrows."""
        pc = ax.pcolormesh(self.xx, self.yy, frame['vUU'],
                           cmap='viridis', shading='gouraud',
                           vmin=self.vel_clim[0], vmax=self.vel_clim[1])
        self._add_colorbar(ax, pc, r'$U$ (m/yr)')
        if self.show_quivers:
            sk = 5
            ax.quiver(self.xx[::sk, ::sk], self.yy[::sk, ::sk],
                      frame['vux'][::sk, ::sk], frame['vuy'][::sk, ::sk],
                      scale=3000, width=0.003, color='k', alpha=0.6)
        self._plot_gps_stations(ax)
        ax.legend(loc='upper left', fontsize=5, frameon=True, bbox_to_anchor=(0.0, 0.90),
                  handletextpad=0.2, ncol=1, markerscale=0.8)
        ax.set_aspect('equal')
        self._add_panel_annotation(ax, annotation_text)
        ax.tick_params(labelsize=6)

    def _plot_velocity_pct(self, ax, frame, ref_frame, annotation_text):
        """Percentage change in speed relative to ref_frame."""
        with np.errstate(divide='ignore', invalid='ignore'):
            pct = 100.0 * (frame['vUU'] - ref_frame['vUU']) / ref_frame['vUU']
        norm_v = TwoSlopeNorm(vmin=self.dvel_clim[0], vcenter=0,
                              vmax=self.dvel_clim[1])
        pc = ax.pcolormesh(self.xx, self.yy, pct, cmap=self._bwr_cmap(),
                           norm=norm_v, shading='gouraud')
        self._add_colorbar(ax, pc, r'$\Delta U / U_0$ (%)')
        ax.set_aspect('equal')
        self._add_panel_annotation(ax, annotation_text)
        ax.tick_params(labelsize=6)

    def _plot_stress_abs(self, ax, frame, annotation_text):
        """Absolute sigma_1 + principal stress crosses colored by sign."""
        pc = ax.pcolormesh(self.xx, self.yy, frame['sigma1'],
                           cmap='jet', shading='gouraud',
                           vmin=self.sigma1_clim[0], vmax=self.sigma1_clim[1])
        self._add_colorbar(ax, pc, r'$\sigma_1$ (kPa)')
        if self.show_quivers:
            sk = self.quiver_skip
            xq = self.xx[::sk, ::sk]
            yq = self.yy[::sk, ::sk]
            s1 = frame['sigma1'][::sk, ::sk]
            s2 = frame['sigma2'][::sk, ::sk]
            theta = frame['t1'][::sk, ::sk]
            dx1 = np.cos(theta)
            dy1 = np.sin(theta)
            dx2 = np.cos(theta + np.pi / 2)
            dy2 = np.sin(theta + np.pi / 2)
            qkw = dict(
                scale=3000,
                width=0.002,
                headwidth=0,
                headlength=0,
                headaxislength=0,
            )
            finite1 = np.isfinite(s1)
            finite2 = np.isfinite(s2)
            pos1 = finite1 & (s1 >= 0)
            neg1 = finite1 & (s1 < 0)
            pos2 = finite2 & (s2 >= 0)
            neg2 = finite2 & (s2 < 0)
            self._quiver_segment(ax, xq[pos1], yq[pos1],  (s1 * dx1)[pos1],  (s1 * dy1)[pos1], 'red', qkw)
            self._quiver_segment(ax, xq[pos1], yq[pos1], -(s1 * dx1)[pos1], -(s1 * dy1)[pos1], 'red', qkw)
            self._quiver_segment(ax, xq[neg1], yq[neg1],  (s1 * dx1)[neg1],  (s1 * dy1)[neg1], 'blue', qkw)
            self._quiver_segment(ax, xq[neg1], yq[neg1], -(s1 * dx1)[neg1], -(s1 * dy1)[neg1], 'blue', qkw)
            self._quiver_segment(ax, xq[pos2], yq[pos2],  (s2 * dx2)[pos2],  (s2 * dy2)[pos2], 'red', qkw)
            self._quiver_segment(ax, xq[pos2], yq[pos2], -(s2 * dx2)[pos2], -(s2 * dy2)[pos2], 'red', qkw)
            self._quiver_segment(ax, xq[neg2], yq[neg2],  (s2 * dx2)[neg2],  (s2 * dy2)[neg2], 'blue', qkw)
            self._quiver_segment(ax, xq[neg2], yq[neg2], -(s2 * dx2)[neg2], -(s2 * dy2)[neg2], 'blue', qkw)
        self._plot_gps_stations(ax)
        ax.set_aspect('equal')
        self._add_panel_annotation(ax, annotation_text)
        ax.tick_params(labelsize=6)

    def _plot_stress_pct(self, ax, frame, ref_frame, annotation_text):
        """Percentage change in sigma_1 relative to ref_frame (no quiver)."""
        with np.errstate(divide='ignore', invalid='ignore'):
            pct = 100.0 * (frame['sigma1'] - ref_frame['sigma1']) / ref_frame['sigma1']
        norm_s = TwoSlopeNorm(vmin=self.dsigma1_clim[0], vcenter=0,
                              vmax=self.dsigma1_clim[1])
        pc = ax.pcolormesh(self.xx, self.yy, pct, cmap=self._bwr_cmap(),
                           norm=norm_s, shading='gouraud')
        self._add_colorbar(ax, pc, r'$\Delta\sigma_1 / \sigma_{1,0}$ (%)')
        ax.set_aspect('equal')
        self._add_panel_annotation(ax, annotation_text)
        ax.tick_params(labelsize=6)

    def _plot_uplift_abs(self, ax, frame, annotation_text):
        """Absolute uplift field hb + hs with water-flux quiver overlay."""
        # norm = Normalize(vmin=self.uplift_clim[0], vmax=self.uplift_clim[1])
        norm = Normalize(vmin=0, vmax=0.5)
        pc = ax.pcolormesh(self.xx, self.yy, frame['uplift'],
                           cmap='viridis', norm=norm, shading='gouraud')
        self._add_colorbar(ax, pc, r'$h_b + h_s$ (m)')
        if self.show_quivers:
            self._plot_flux_quiver(ax, frame['flux_quiver'])
        self._plot_gps_stations(ax)
        ax.set_aspect('equal')
        self._add_panel_annotation(ax, annotation_text)
        ax.tick_params(labelsize=6)

    def _plot_right_column(self, ax, frame, ref_frame, annotation_text, row_index):
        if self.right_column_mode == 'stress':
            if row_index == 0:
                self._plot_stress_abs(ax, frame, annotation_text)
            else:
                self._plot_stress_pct(ax, frame, ref_frame, annotation_text)
            return
        self._plot_uplift_abs(ax, frame, annotation_text)

    def _right_column_annotation(self, frame, ref_frame, panel_label, row_index):
        td = frame['t_days']
        if self.right_column_mode == 'stress':
            if row_index == 0:
                return f'({panel_label}) $\\sigma_1$,  $t = {td:.0f}$ d'
            return f'({panel_label}) $\\Delta\\sigma_1/\\sigma_{{1,0}}$,  $t = {td:.0f}$ vs {ref_frame["t_days"]:.0f} d'
        return f'({panel_label}) $h_b+h_s$,  $t = {td:.0f}$ d'

    # ------------------------------------------------------------------
    # Main figure
    # ------------------------------------------------------------------
    def generate_figure(self, save_path=None):
        """
        Layout: n_snap rows × 2 columns
          Left column  = velocity   (row 0: absolute, rows 1+: % change)
          Right column = stress or uplift, controlled by right_column_mode
          Stress quiver is used only for the first stress row.
        """
        if not self.data_loaded:
            self.load_data()
        plt.style.use('classic')
        plt.rcParams.update({
            'font.family': 'DejaVu Sans',
            'mathtext.fontset': 'dejavusans',
            'font.size': 8,
            'figure.dpi': 200,
            'savefig.dpi': 600,
        })
        n_snap = len(self.frames)
        fig, axes = plt.subplots(n_snap, 2, figsize=(8, 2.1 * n_snap),
                                 dpi=200)
        if n_snap == 1:
            axes = axes.reshape(1, 2)
        ref = self.frames[0]
        letters = [chr(97 + i) for i in range(2 * n_snap)]  # a, b, c, ...
        for i, frame in enumerate(self.frames):
            td = frame['t_days']
            li = letters[2 * i]              # left  panel label
            ri = letters[2 * i + 1]          # right panel label
            if i == 0:
                self._plot_velocity_abs(
                    axes[i, 0], frame,
                    f'({li}) $U$,  $t={td:.0f}$ d')
            else:
                self._plot_velocity_pct(
                    axes[i, 0], frame, ref,
                    f'({li}) $\\Delta U/U_0$,  $t={td:.0f}$ vs {ref["t_days"]:.0f} d')
            self._plot_right_column(
                axes[i, 1], frame, ref,
                self._right_column_annotation(frame, ref, ri, i),
                i)
            axes[i, 0].set_ylabel('$y$ (km)', fontsize=8)
            axes[i, 1].tick_params(labelleft=False)
            if i == n_snap - 1:
                axes[i, 0].set_xlabel('$x$ (km)', fontsize=8)
                axes[i, 1].set_xlabel('$x$ (km)', fontsize=8)
        fig.subplots_adjust(hspace=0.0, wspace=0.20)
        if save_path:
            fig.savefig(save_path, dpi=600, bbox_inches='tight',
                        pad_inches=0.1)
            print(f"Saved to {save_path}")
        return fig

# ======================================================================
# Example usage
# ======================================================================
if __name__ == "__main__":
    casename = 'n2d_regional_racmo_eps1e_02_kappa5e_11_mu5e0_partition3e_01_k01e_01_spinup'
    plotter = NevisIceFlowPlotter(casename, right_column_mode='uplift', show_quivers=True)
    plotter.tslice_d = np.array([100, 180, 190, 200])  # days
    plotter.vel_clim = [0, 200]
    plotter.sigma1_clim = [0, 200]
    plotter.dsigma1_clim = [-100, 100]
    plotter.dvel_clim = [-100, 100]
    plotter.uplift_clim = [0.0, 1.0]
    plotter.quiver_skip = 8
    fig = plotter.generate_figure(
        save_path=f'./figures/{casename}_ice_flow_snapshots.png'
    )
    plt.show()
