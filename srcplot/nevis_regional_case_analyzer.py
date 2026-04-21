import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import os
import scipy.io as sio


def _enable_retina_display():
    try:
        from IPython.display import set_matplotlib_formats

        set_matplotlib_formats('retina')
    except Exception:
        pass

class NevisRegionalCaseAnalyzer:
    """
    A class for loading and analyzing NEVIS .mat result files.

    This class can process one or more casenames, calculate the front propagation speed
    for each case, and store the results for later use.
    """
    def __init__(self, results_dir='results'):
        """
        Initialize the analyzer.

        Args:
            results_dir (str): Root directory containing casename result folders.
        """
        self.results_dir = results_dir
        self.results = {}

    def _load_and_process_case(self, casename):
        """
        Load and process data for a single case.
        """
        mat_path = os.path.join(self.results_dir, casename, f"{casename}.mat")
        try:
            raw = sio.loadmat(mat_path, squeeze_me=True)
        except FileNotFoundError:
            print(f"Warning: File not found {mat_path}, skipping this case.")
            return None

        # Extract data
        tt = raw['tt']
        ps = raw['ps']
        
        # Type conversion
        ps_x = np.float64(ps['x'])
        ps_t = np.float64(ps['t'])
        ps_hb = np.float64(ps['hb'])
        
        # Calculate derived quantities
        t = (ps_t / (24 * 60 * 60)) * np.float64(tt['t'])
        # Core analysis logic
        L_bx = np.abs(ps_x * np.float64(tt['Lbx']) - 40e3)
        L_by = np.abs(ps_x * np.float64(tt['Lby']) + 15e3)
        L_b = np.sqrt(L_bx**2 + L_by**2)
        # L_b = L_bx
        V_b = ps_x**2 * ps_hb * np.float64(tt['Vb'])
        Q_b_out = np.float64(tt['Q_outb'])
        return t, L_bx, L_by, L_b, V_b, Q_b_out

    def run_analysis(self, casenames):
        """
        Run analysis for a list of casenames.
        """
        vt = []
        Lbx = []
        Lby = []
        Lb = []
        Vb = []
        Q_out_b = []
        for casename in casenames:
            result = self._load_and_process_case(casename)
            if result is not None:
                t, L_bx, L_by, L_b, V_b, Q_b_out = result
                vt.append(t)
                Lbx.append(L_bx)
                Lby.append(L_by)
                Lb.append(L_b)
                Vb.append(V_b)
                Q_out_b.append(Q_b_out)

        return vt, Lbx, Lby, Lb, Vb, Q_out_b

if __name__ == "__main__":
    fontsize_panel = 9

    casenames = ['n2d_region_test_meanperms1_Hreg1000_kappa1e_10_mu1e1_spinup',\
                 'n2d_region_test_meanperms1_Hreg1000_kappa1e_11_mu1e1_spinup',\
                 'n2d_region_test_meanperms1_Hreg1000_kappa5e_12_mu1e1_spinup',\
                 'n2d_region_test_meanperms1_Hreg1000_kappa1e_99_mu1e1_spinup']

    casenames = ['n2d_region_melt_meanperms1_Hreg1000_kappa1e_10_mu1e1_td30_V1e8_drainage',\
                 'n2d_region_melt_meanperms1_Hreg1000_kappa1e_11_mu1e1_td30_V1e8_drainage',\
                 'n2d_region_melt_meanperms1_Hreg1000_kappa5e_12_mu1e1_td30_V1e8_drainage',\
                 'n2d_region_melt_meanperms1_Hreg1000_kappa1e_99_mu1e1_td30_V1e8_drainage']

    vkappa = [1e-10, 1e-11, 5e-12, 0]
    vmu = [10, 10, 10, 10]
    fig, ax = plt.subplots(1, 2, figsize=(6, 2.5), dpi=100)
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

    t_drainage = 30
    analyzer = NevisRegionalCaseAnalyzer()
    vt, Lbx, Lby, Lb, Vb, Q_out_b = analyzer.run_analysis(casenames)
    colors = ['b', 'r', 'g', 'k']
    t_Qout = []
    for i in range(len(casenames)):
        Q = Q_out_b[i]
        idx = np.where(Q > 1e-3)[0]

        if len(idx) == 0:
            t_Qout.append(vt[i][-1] - t_drainage)
        else:
            t_Qout.append(vt[i][idx[0]] - t_drainage)
        print(f"Case {casenames[i]}: t_Qout = {t_Qout[i]:.2f} days")

    ax = plt.subplot(121)
    for i in range(len(casenames)):
        Vb[i] = np.where((vt[i] < t_drainage), np.nan, Vb[i])
        ax.plot(vt[i] - t_drainage - 0.5, Vb[i] / 1.8e8, linestyle='-', linewidth=1.0, color=colors[i], label=r'$\kappa$=' + str(vkappa[i]))
        valid_indices = np.where(vt[i] > t_Qout[i] + t_drainage)[0]
        if len(valid_indices) > 0:
            last_valid_index = valid_indices[0]
            ax.plot(vt[i][last_valid_index] - t_drainage - 0.5, Vb[i][last_valid_index] / 1.8e8, 'o', color=colors[i], markersize=6)

    ax.set_xlim(0, 40)
    ax.set_ylim(0, 1.2)
    plt.legend(loc='upper right', fontsize=6, frameon=True, ncol=1,
              labelspacing=0.3, handletextpad=0.3, borderpad=0.2)
    ax.set_xlabel('$t~(d)$')
    ax.set_ylabel(r'$V_b/V_l$')
    ax.text(0.03, 0.95, '(a) volume', transform=ax.transAxes,
            ha='left', va='top', fontsize=fontsize_panel + 1, clip_on=False)
    ax.grid(alpha=0.5)

    ax = plt.subplot(122)
    for i in range(len(casenames)):
        window_size = 5
        Lb[i] = np.convolve(Lb[i], np.ones(window_size) / window_size, mode='same')
        Lb[i] = np.where((vt[i] < t_drainage + 0.5) | (vt[i] > t_Qout[i] + t_drainage), np.nan, Lb[i])
        ax.plot(vt[i] - t_drainage - 0.5, (Lb[i]) / 1e3, linestyle='-', linewidth=1.25, color=colors[i], label=r'$\kappa$=' + str(vkappa[i]))

        valid_indices = np.where(~np.isnan(Lb[i]))[0]
        if len(valid_indices) > 0:
            last_valid_index = valid_indices[-1]
            ax.plot(vt[i][last_valid_index] - t_drainage - 0.5, Lb[i][last_valid_index] / 1e3, 'o', color=colors[i], markersize=6)

    ax.set_xlim(0, 40)
    ax.set_ylim(0, 150)
    plt.legend(loc='upper right', fontsize=6, frameon=True, ncol=1,
              labelspacing=0.3, handletextpad=0.3, borderpad=0.2)
    ax.set_xlabel('$t~(d)$')
    ax.set_ylabel(r'$L_b~(km)$')
    ax.yaxis.set_label_coords(-0.10, 0.5)

    ax.text(0.03, 0.95, '(b) front', transform=ax.transAxes,
            ha='left', va='top', fontsize=fontsize_panel + 1, clip_on=False)
    ax.grid(alpha=0.5)

    plt.savefig('./figures/Fig8_regional_front_volume.png', dpi=600, transparent=False, bbox_inches='tight', pad_inches=0.1)