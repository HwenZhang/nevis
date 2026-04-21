import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import os
import scipy.io as sio
from mpl_toolkits.axes_grid1 import make_axes_locatable


def _enable_retina_display():
    try:
        from IPython.display import set_matplotlib_formats

        set_matplotlib_formats('retina')
    except Exception:
        pass

class VolumeAnalysis:
    def __init__(self, casenames, vkappa):
        self.casenames = casenames
        self.vkappa = vkappa
        self.avg_dVdt = []
        # set up a color scheme
        self.colors = plt.cm.viridis(np.linspace(0, 1, len(casenames)))

    def process_case(self, casename, vkappa_value, i, ax1):
        mat_path = os.path.join('results', casename, casename + '.mat')
        raw = sio.loadmat(mat_path, squeeze_me=True)

        tt = raw['tt']
        ps = raw['ps']

        ps_hb = np.float64(ps['hb'])
        ps_x = np.float64(ps['x'])
        ps_S = np.float64(ps['S'])
        ps_h = np.float64(ps['h'])
        ps_t = np.float64(ps['t'])

        V_b = ps_x**2 * ps_hb * np.float64(tt['Vb'])
        dVdt = np.gradient(V_b, ps_t*(24*60*60))  # velocity in m/s

        tt_Q_outb = np.float64(tt['Q_outb'])
        t = (ps_t / (24*60*60)) * np.float64(tt['t'])

        idx_start = np.where((V_b > 0.99 * np.max(V_b)) & (dVdt < 0))[0][0]
        # idx_start = np.where(t > 750.25)[0][0]
        idx_end = np.where(t > t[idx_start] + 30)[0][0]
        print(f"Processing: {casename}, Start time: {t[idx_start]:.2f} d, End time: {t[idx_end]:.2f} d")
        # find the time when Q_outb starts to increase significantly
        idx_out = np.where(tt_Q_outb > 1e-6)[0]
        if idx_out.size > 0:
            print(f"  Q_outb starts to increase at t = {t[idx_out[0]]:.2f} d")
        else:
            print(f"  Q_outb does not increase significantly in this case.")

        Vb0 = V_b[idx_start]
        # Vb0 = 1e7
        Vbend = V_b[idx_end]
        print(f"  Volume at start: {Vb0:.2f} m^3")
        print(f"  Volume at end: {Vbend:.2f} m^3")
        print(f"  Volume difference: {Vbend - Vb0:.2f} m^3 over {t[idx_end] - t[idx_start]:.2f} d")
        V_b[idx_start]= np.nan
        dVdt0 = (V_b[1]-V_b[0])/(t[1]-t[0])
        if casename.startswith('n2d_0m3s'):
            self.avg_dVdt.append((Vbend - Vb0) / (t[idx_end] - t[idx_start])-dVdt0)
        else:
            self.avg_dVdt.append((Vbend - Vb0) / (t[idx_end] - t[idx_start])-dVdt0)
        

        exp = int(np.log10(vkappa_value))
        if casename.startswith('n2d_0m3s'):
            ax1.plot(t-t[idx_start], V_b/Vb0 - 0*dVdt0*(t-t[idx_start])/Vb0, linewidth=1.5, linestyle='-', color=self.colors[i])
            # label=rf'$Q_m=0~m^3~s$'
        else:
            ax1.plot(t-t[idx_start], V_b/Vb0 - 0*dVdt0*(t-t[idx_start])/Vb0, linewidth=1.5, linestyle='--', color=self.colors[i])
            # label=rf'$Q_m=10~m^3~s$'


    def analyze(self, ax1):
        for i, casename in enumerate(self.casenames):
            self.process_case(casename, self.vkappa[i], i, ax1)

if __name__ == "__main__":
    casenames = ['n2d_0m3s_kappa6e_10_mu1e03_hbreg5e_3_V1e7_drainage',\
                  'n2d_0m3s_kappa4e_10_mu1e03_hbreg5e_3_V1e7_drainage',\
                  'n2d_0m3s_kappa2e_10_mu1e03_hbreg5e_3_V1e7_drainage',\
                 'n2d_0m3s_kappa1e_10_mu1e3_hbreg5e_3_V1e7_drainage',\
                 'n2d_0m3s_kappa1e_11_mu1e3_hbreg5e_3_V1e7_drainage',\
                 'n2d_0m3s_kappa1e_12_mu1e3_hbreg5e_3_V1e7_drainage',\
                 'n2d_0m3s_kappa1e_13_mu1e3_hbreg5e_3_V1e7_drainage',\
                 'n2d_0m3s_kappa1e_14_mu1e3_hbreg5e_3_V1e7_drainage',\
                 'n2d_0m3s_kappa1e_15_mu1e3_hbreg5e_3_V1e7_drainage']
    casenames = casenames[::-1]

    casenames2 = ['n2d_10m3s_kappa6e_10_mu1e03_hbreg5e_3_V1e7_drainage',\
                  'n2d_10m3s_kappa4e_10_mu1e03_hbreg5e_3_V1e7_drainage',\
                  'n2d_10m3s_kappa2e_10_mu1e03_hbreg5e_3_V1e7_drainage',\
                 'n2d_10m3s_kappa1e_10_mu1e03_hbreg5e_3_V1e7_drainage',\
                 'n2d_10m3s_kappa1e_11_mu1e03_hbreg5e_3_V1e7_drainage',\
                 'n2d_10m3s_kappa1e_12_mu1e03_hbreg5e_3_V1e7_drainage',\
                 'n2d_10m3s_kappa1e_13_mu1e03_hbreg5e_3_V1e7_drainage',\
                 'n2d_10m3s_kappa1e_14_mu1e03_hbreg5e_3_V1e7_drainage',\
                 'n2d_10m3s_kappa1e_15_mu1e03_hbreg5e_3_V1e7_drainage']
    casenames2 = casenames2[::-1]

    fig, ax = plt.subplots(1, 2, figsize=(7, 3), dpi=100)
    plt.style.use('classic')
    plt.rcParams["figure.dpi"] = 600
    plt.rcParams['font.family'] = 'DejaVu Sans'
    plt.rcParams['mathtext.fontset'] = 'dejavusans'
    plt.rcParams['font.size'] = 8
    plt.rcParams['axes.linewidth'] = 1.0
    plt.rcParams['lines.linewidth'] = 0.6
    plt.rcParams['grid.linewidth'] = 0.3
    plt.rcParams['axes.labelsize'] = 8
    plt.rcParams['xtick.labelsize'] = 8
    plt.rcParams['ytick.labelsize'] = 8
    _enable_retina_display()

    ax1 = plt.subplot(121)

    vkappa = np.array([1e-15, 1e-14, 1e-13, 1e-12, 1e-11, 1e-10, 2e-10, 4e-10, 6e-10])
    fontsize_panel = 9

    volume_analysis = VolumeAnalysis(casenames, vkappa)
    volume_analysis.analyze(ax1)
    avg_dVdt1 = volume_analysis.avg_dVdt
    volume_analysis = VolumeAnalysis(casenames2, vkappa)
    volume_analysis.analyze(ax1)
    avg_dVdt2 = volume_analysis.avg_dVdt
    ax1.plot([], [], color='k', linewidth=1.0, linestyle='-', label=rf'$Q_m=0~m^3~s$')
    ax1.plot([], [], color='k', linewidth=1.0, linestyle='--', label=rf'$Q_m=10~m^3~s$')

    ax1.set_xlim(0, 300.0)
    ax1.set_ylim(-0.2, 1.2)
    ax1.set_xlabel(r'$t~(d)$', fontsize=fontsize_panel)
    ax1.set_ylabel(r'normalised volume', fontsize=fontsize_panel)
    ax1.tick_params(axis='both', labelsize=fontsize_panel)
    ax1.text(0.03, 0.92, r'(a) volume change', fontsize=fontsize_panel + 2, transform=ax1.transAxes)
    ax1.legend(loc='lower left', fontsize=fontsize_panel - 2, frameon=True, handletextpad=0.2, ncol=2, labelspacing=0.2, numpoints=1)
    ax1.grid(alpha=0.5)

    divider = make_axes_locatable(ax1)
    cax1 = divider.append_axes("right", size="3%", pad=0.1)
    sm = plt.cm.ScalarMappable(cmap=plt.cm.viridis, norm=plt.Normalize(vmin=np.log10(vkappa.min()), vmax=np.log10(vkappa.max())))
    sm.set_array([])
    cbar = plt.colorbar(sm, cax=cax1)
    cbar.set_label(r'$\log_{10}(\kappa)$', fontsize=fontsize_panel)
    cbar.ax.tick_params(labelsize=fontsize_panel)
    cbar.set_ticks(np.arange(np.floor(np.log10(vkappa.min())), np.ceil(np.log10(vkappa.max())) + 1, 1))

    ax2 = plt.subplot(122)
    ax2.set_visible(False)
    pos_ax2 = ax2.get_position()
    ax2_width = pos_ax2.x1 - pos_ax2.x0
    ax2_height = pos_ax2.y1 - pos_ax2.y0
    ax2_pos = [pos_ax2.x0 + 0.30 * ax2_width, pos_ax2.y0, ax2_width, ax2_height]
    ax2 = fig.add_axes(ax2_pos)

    avg_dVdt1 = np.array(avg_dVdt1)
    avg_dVdt2 = np.array(avg_dVdt2)

    cmap = plt.cm.viridis
    norm = plt.Normalize(vmin=np.log10(vkappa.min()), vmax=np.log10(vkappa.max()))
    ax2.plot(vkappa, -avg_dVdt1, linestyle='-', linewidth=1.0,
             color='blue', label=r'$Q_m=0~\text{m}^3~\text{s}^{-1}$')
    ax2.plot(vkappa, -avg_dVdt2, linestyle='-', linewidth=1.0,
             color='red', label=r'$Q_m=10~\text{m}^3~\text{s}^{-1}$')
    for i in range(len(vkappa)):
        color_i = cmap(norm(np.log10(vkappa[i])))
        ax2.plot(vkappa[i], -avg_dVdt1[i], marker='o', linestyle='',
                 markersize=6, color=color_i)
        ax2.plot(vkappa[i], -avg_dVdt2[i], marker='o', linestyle='',
                 markersize=6, color=color_i)

    ax2.plot(vkappa, 1e15 * vkappa**(1), linestyle='--', linewidth=1.0, color='gray', label=r'linear in $\kappa$')

    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_xlim(1e-14, 1e-9)
    ax2.set_ylim(1e1, 1e6)
    ax2.set_xlabel(r'$\kappa$', fontsize=fontsize_panel)
    ax2.set_ylabel(r'rate of volume decrease ($m^3$/d)', fontsize=fontsize_panel)
    ax2.tick_params(axis='both', labelsize=fontsize_panel)
    ax2.text(0.03, 0.92, r'(c) volume rate vs $\kappa$', fontsize=fontsize_panel + 2, transform=ax2.transAxes)
    ax2.legend(loc='lower right', fontsize=fontsize_panel - 1, frameon=True, handletextpad=0.2, ncol=1, labelspacing=0.4, numpoints=1)
    ax2.grid(alpha=0.5)

    fig.savefig(
        './figures/Fig6_volume_poster.png',
        dpi=600,
        transparent=False,
        bbox_inches='tight',
        pad_inches=0.1
    )

    plt.show()