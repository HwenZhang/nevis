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


class NevisCaseAnalyzer:
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

        tt = raw['tt']
        ps = raw['ps']

        ps_x = np.float64(ps['x'])
        ps_t = np.float64(ps['t'])

        R_b = np.asarray(ps_x * np.float64(tt['Rb']), dtype=float)
        t = np.asarray((ps_t / (24 * 60 * 60)) * np.float64(tt['t']), dtype=float)

        R_b = np.where(R_b > 25e3, 0.0, R_b)
        if not np.any(R_b > 0):
            print(f"Warning: In {casename}, R_b has no positive values, unable to calculate speed.")
            return None

        R_max = np.max(R_b)
        start_indices = np.where((R_b > 0.2 * R_max) & (R_b < 0.3 * R_max))[0]
        end_indices = np.where((R_b > 0.8 * R_max) & (R_b < R_max))[0]

        if len(start_indices) == 0 or len(end_indices) == 0:
            print(f"Warning: Unable to determine start or end points in {casename}.")
            return None

        idx_start = start_indices[0]
        idx_end = end_indices[0]

        delta_R = R_b[idx_end] - R_b[idx_start]
        delta_t = t[idx_end] - t[idx_start]

        if delta_t == 0:
            print(f"Warning: Time difference is zero in {casename}, unable to calculate speed.")
            return None

        avg_speed = (delta_R / delta_t) / (24 * 3600)
        if avg_speed < 1e-3 or avg_speed > 1e2:
            avg_speed = np.nan
        print(f"Processing: {casename}, Average Speed: {avg_speed:.4f} m/s")
        return avg_speed

    def run_analysis(self, casenames):
        """
        Run analysis for a list of casenames.

        Args:
            casenames (list): List of casename strings to analyze.

        Returns:
            list: List of calculated average speeds.
        """
        avg_speeds = []
        for casename in casenames:
            print(f"Analyzing case: {casename}")
            speed = self._load_and_process_case(casename)
            if speed is not None:
                avg_speeds.append(speed)

        self.results['avg_speeds'] = avg_speeds
        return np.array(avg_speeds)


if __name__ == "__main__":
    fontsize_panel = 9

    casenames1 = ['n2d_0m3s_kappa1e_20_mu1e_3_hbreg5e_3_V1e6_drainage',\
                 'n2d_0m3s_kappa1e_20_mu1e_2_hbreg5e_3_V1e6_drainage',\
                 'n2d_0m3s_kappa1e_20_mu1e_1_hbreg5e_3_V1e6_drainage',\
                 'n2d_0m3s_kappa1e_20_mu1e0_hbreg5e_3_V1e6_drainage',\
                 'n2d_0m3s_kappa1e_20_mu1e1_hbreg5e_3_V1e6_drainage',\
                 'n2d_0m3s_kappa1e_20_mu1e2_hbreg5e_3_V1e6_drainage',\
                 'n2d_0m3s_kappa1e_20_mu1e3_hbreg5e_3_V1e6_drainage']

    casenames2 = ['n2d_0m3s_kappa1e_20_mu1e_3_hbreg5e_3_V1e7_drainage',\
                 'n2d_0m3s_kappa1e_20_mu1e_2_hbreg5e_3_V1e7_drainage',\
                 'n2d_0m3s_kappa1e_20_mu1e_1_hbreg5e_3_V1e7_drainage',\
                 'n2d_0m3s_kappa1e_20_mu1e0_hbreg5e_3_V1e7_drainage',\
                 'n2d_0m3s_kappa1e_20_mu1e1_hbreg5e_3_V1e7_drainage',\
                 'n2d_0m3s_kappa1e_20_mu1e2_hbreg5e_3_V1e7_drainage',\
                 'n2d_0m3s_kappa1e_20_mu1e3_hbreg5e_3_V1e7_drainage']

    casenames3 = ['n2d_0m3s_kappa1e_20_mu1e_3_hbreg5e_3_V1e8_drainage',\
                 'n2d_0m3s_kappa1e_20_mu1e_2_hbreg5e_3_V1e8_drainage',\
                 'n2d_0m3s_kappa1e_20_mu1e_1_hbreg5e_3_V1e8_drainage',\
                 'n2d_0m3s_kappa1e_20_mu1e0_hbreg5e_3_V1e8_drainage',\
                 'n2d_0m3s_kappa1e_20_mu1e1_hbreg5e_3_V1e8_drainage',\
                 'n2d_0m3s_kappa1e_20_mu1e2_hbreg5e_3_V1e8_drainage',\
                 'n2d_0m3s_kappa1e_20_mu1e3_hbreg5e_3_V1e8_drainage']

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

    vmu = np.array([1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3])

    ax = plt.subplot(121)
    analyzer = NevisCaseAnalyzer()
    avg_speed1 = analyzer.run_analysis(casenames1)
    avg_speed2 = analyzer.run_analysis(casenames2)
    avg_speed3 = analyzer.run_analysis(casenames3)

    mask = np.arange(len(vmu)) != 4
    ax.plot(vmu, avg_speed1, marker='o', markersize=5, linestyle='-', color='b', label=r'$V_l=10^6~m^3$')
    ax.plot(vmu[mask], avg_speed2[mask], marker='o', markersize=5, linestyle='-', color='r', label=r'$V_l=10^7~m^3$')
    ax.plot(vmu, avg_speed3, marker='o', markersize=5, linestyle='-', color='g', label=r'$V_l=10^8~m^3$')
    ax.plot(vmu[4], avg_speed2[4], marker='^', markersize=7, linestyle='None', color='red', label='wintertime reference')

    vmutheory = np.linspace(1e-4, 1e3, 10)
    ax.plot(vmutheory, 6e0 * vmutheory**(-1 / 11), linewidth=1.5, linestyle='--', color='gray')
    ax.plot(vmutheory, vmutheory**(-1), linewidth=1.5, linestyle=(0, (12, 3, 2, 3)), color='gray')
    ax.text(0.55, 0.82, r'$\mu^{-1/11}$(bending)', fontsize=fontsize_panel, color='black', transform=ax.transAxes, rotation=-180 / np.pi * 1 / 11)
    ax.text(0.67, 0.13, r'$\mu^{-1}$(gravity)', fontsize=fontsize_panel, color='black', transform=ax.transAxes, rotation=-0.8 * 180 / np.pi)

    ax.axhspan(0.01, 1.3, color='royalblue', alpha=0.2, label='observations')

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$\mu_{\text{eff}}~(Pa~s)$', fontsize=fontsize_panel)
    ax.set_ylabel('Average Speed (m/s)', fontsize=fontsize_panel)
    ax.tick_params(axis='both', labelsize=fontsize_panel - 1)
    ax.set_ylim(1e-4, 1e2)
    ax.set_xlim(1e-4, 1e3)
    ax.text(0.05, 0.9, r'(a) $\bar{v}\left(\mu\right)$', fontsize=fontsize_panel + 3, transform=ax.transAxes)
    ax.grid(alpha=0.5)
    ax.legend(loc='lower left', fontsize=fontsize_panel - 3, numpoints=1)

    ax = plt.subplot(122)
    vVl = np.array([1e6, 1e7, 1e8])

    casenames1 = ['n2d_0m3s_kappa1e_20_mu1e_1_hbreg5e_3_V1e6_drainage',\
                 'n2d_0m3s_kappa1e_20_mu1e_1_hbreg5e_3_V1e7_drainage',\
                 'n2d_0m3s_kappa1e_20_mu1e_1_hbreg5e_3_V1e8_drainage']

    casenames2 = ['n2d_0m3s_kappa1e_20_mu1e1_hbreg5e_3_V1e6_drainage',\
                 'n2d_0m3s_kappa1e_20_mu1e1_hbreg5e_3_V1e7_drainage',\
                 'n2d_0m3s_kappa1e_20_mu1e1_hbreg5e_3_V1e8_drainage']

    casenames3 = ['n2d_0m3s_kappa1e_20_mu1e_3_hbreg5e_3_V1e6_drainage',\
                 'n2d_0m3s_kappa1e_20_mu1e_3_hbreg5e_3_V1e7_drainage',\
                 'n2d_0m3s_kappa1e_20_mu1e_3_hbreg5e_3_V1e8_drainage']

    avg_speed1 = analyzer.run_analysis(casenames1)
    avg_speed2 = analyzer.run_analysis(casenames2)
    avg_speed3 = analyzer.run_analysis(casenames3)

    mask = np.arange(len(vVl)) != 1
    ax.plot(vVl[mask], avg_speed2[mask], marker='o', markersize=5, linestyle='-', color='b', label=r'$\mu=10^{1}~Pa~s$')
    ax.plot(vVl, avg_speed1, marker='o', markersize=5, linestyle='-', color='r', label=r'$\mu=10^{-1}~Pa~s$')
    ax.plot(vVl, avg_speed3, marker='o', markersize=5, linestyle='-', color='g', label=r'$\mu=10^{-3}~Pa~s$')
    ax.plot(vVl[1], avg_speed2[1], marker='^', markersize=7, linestyle='None', color='blue', label='wintertime reference')

    ax.plot(vVl, 1e-10 * vVl**(5 / 4), linewidth=1.5, linestyle=(0, (12, 3, 2, 3)), color='gray')
    ax.text(0.60, 0.40, r'$V_l^{5/4}$(gravity)', fontsize=fontsize_panel, color='black', transform=ax.transAxes, rotation=180 / np.pi * 1 / 2)
    ax.plot(vVl, 0.25e0 * vVl**(5 / 22), linewidth=1.5, linestyle='--', color='gray')
    ax.text(0.55, 0.85, r'$V_l^{5/22}$(bending)', fontsize=fontsize_panel, color='black', transform=ax.transAxes, rotation=180 / np.pi * 2 / 22)

    ax.axhspan(0.01, 1.3, color='royalblue', alpha=0.2, label='observations')

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$V_l~(m^3)$')
    ax.set_xlim(1e6, 1e8)
    ax.set_ylim(1e-4, 1e2)
    ax.text(0.05, 0.9, r'(b) $\bar{v}\left(V_l\right)$', fontsize=fontsize_panel + 3, transform=ax.transAxes)
    ax.grid(alpha=0.5)
    ax.legend(loc='lower right', fontsize=fontsize_panel - 3, numpoints=1)
    fig.savefig(
        './figures/Fig5_velocity.png',
        dpi=600,
        transparent=False,
        bbox_inches='tight',
        pad_inches=0.1
    )

    plt.show()