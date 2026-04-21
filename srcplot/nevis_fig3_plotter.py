import os
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import PowerNorm
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from cycler import cycler
import subprocess

class NevisFig3Plotter:
    """
    A class to generate Figure 3 (summertime) plot for NEVIS simulation data.
    This class packs all the plotting code from the first cell into a single class.
    """
    
    def __init__(self, casename, results_dir='results', matlab_path="/Applications/MATLAB_R2025a.app/bin/matlab"):
        """
        Initialize the plotter with casename and paths.
        
        Args:
            casename: Name of the simulation case
            results_dir: Directory containing results
            matlab_path: Path to MATLAB executable
        """
        self.casename = casename
        self.results_dir = results_dir
        self.matlab_path = matlab_path
        
        # Flag to control data import
        self.SKIP_DATA_IMPORT = False
        self.LOAD_NEW_DATA = 1  # Set to True to load new data, False to use existing data
        
        # -------------------- Parameters for plotting ---------------------------------------------
        okabe_ito = ["#E69F00","#56B4E9","#009E73","#F0E442",
                     "#0072B2","#D55E00","#CC79A7","#000000"]
        self.colors = okabe_ito
        self.colors = plt.get_cmap('tab10').colors  

        self.xa = 10 
        self.xb = 50
        self.ya = 5
        self.yb = 15
        self.plot_hb = False
        self.plot_qs = True
        self.h_min = 0  # min water thickness on the color bar (m)
        self.h_max = 0.1  # max water thickness on the color bar (m)
        self.q_min = 1e-3 # min water flux on the color bar (m^2/s)
        self.q_max = 5e-2 # max water flux on the color bar (m^2/s)
        # norm = PowerNorm(gamma=1, vmin=q_min, vmax=q_max)
        self.norm = LogNorm(vmin=self.q_min, vmax=self.q_max)

        self.tdrainge_d = 2.0 * 365 + 20
        self.tmin = self.tdrainge_d - 5
        self.tmax = self.tdrainge_d + 5
        self.tmin_d = self.tmin 
        self.tmax_d = self.tmax
        self.tslice = self.tdrainge_d + np.array([-1, 0, 1, 10])
        
        # Data containers
        self.data_loaded = False
        
    def load_data(self):
        """Load and process all data"""
        if (not self.SKIP_DATA_IMPORT) or (self.LOAD_NEW_DATA):
            print("Loading data...")
            # ------------------------- Import data -----------------------------------

            mat_path = os.path.join(self.results_dir, self.casename, self.casename + '.mat')
            raw = sio.loadmat(mat_path, squeeze_me=True)

            self.tt = raw['tt']
            self.ps = raw['ps']
            self.pp = raw['pp']
            self.oo = raw['oo']
            self.pd = raw['pd']
            self.aa = raw['aa']
            self.gg = raw['gg']

            self.ps_t = np.float64(self.ps['t'])
            self.pd_td = np.float64(self.pd['td'])
            self.ps_hb = np.float64(self.ps['hb'])
            self.ps_h = np.float64(self.ps['h'])
            self.ps_x = np.float64(self.ps['x'])
            self.ps_m = np.float64(self.ps['m'])
            self.ps_S = np.float64(self.ps['S'])
            self.ps_phi = np.float64(self.ps['phi'])
            self.ps_Q0 = np.float64(self.ps['Q0'])
            self.ps_Q = np.float64(self.ps['Q'])
            self.ps_qs = np.float64(self.ps['qs'])
            self.ps_qb = np.float64(self.ps['qb'])

            self.pp_c8 = np.float64(self.pp['c8'])

            self.tt_t = np.float64(self.tt['t'])
            self.tt_Qb_in = np.float64(self.tt['Qb_in'])
            self.tt_Qb_dec = np.float64(self.tt['Qb_dec'])
            self.tt_Q_outb = np.float64(self.tt['Q_outb'])
            self.tt_Q_outQ = np.float64(self.tt['Q_outQ'])
            self.tt_Q_outq = np.float64(self.tt['Q_outq'])
            self.tt_E = np.float64(self.tt['E'])

            # --- Compute time vectors -----------------------------------------------
            self.t = (self.ps_t / (24*60*60)) * self.tt_t        # days

            # --- Build derived quantities -------------------------------------------
            self.Q_b_in = np.float64(self.pd['Q_0']) * self.tt_Qb_in
            self.Q_b_dec = self.ps_h * self.ps_x**2 / self.ps_t * self.tt_Qb_dec

            self.Q_out_b = self.ps_Q0 * self.tt_Q_outb
            self.Q_out_b = self.ps_Q0 * self.tt_Q_outb
            self.Q_out_Q = self.ps_Q * self.tt_Q_outQ
            self.Q_out_q = self.ps_Q * self.tt_Q_outq
            self.E = self.ps_m * self.ps_x**2 * self.tt_E

            self.pts_ni = np.int64(self.oo['pts_ni'].item())
            self.phi_0 = self.ps_phi * np.float64(self.aa['phi_0'].item())
            self.pts_phi_lake = self.ps_phi * np.array([np.float64(sub_arr[0]) for sub_arr in self.tt['pts_phi']])
            self.pts_phi_moulin = self.ps_phi * np.array([np.float64(sub_arr[1]) for sub_arr in self.tt['pts_phi']])
            self.pts_hb_lake = self.ps_hb * np.array([np.float64(sub_arr[0]) for sub_arr in self.tt['pts_hb']])
            self.pts_hb_moulin = self.ps_hb * np.array([np.float64(sub_arr[1]) for sub_arr in self.tt['pts_hb']])
            self.pts_pb_lake = self.ps_phi * np.array([np.float64(sub_arr[0]) for sub_arr in self.tt['pts_pb']])
            self.pts_hc_lake = self.ps_h * np.array([np.float64(sub_arr[0]) for sub_arr in self.tt['pts_hc']])
            self.pts_hs_lake = self.ps_h * np.array([np.float64(sub_arr[0]) for sub_arr in self.tt['pts_hs']])

            self.pts_hb_gps1 = self.ps_hb * np.array([np.float64(sub_arr[2]) for sub_arr in self.tt['pts_hb']])
            self.pts_hb_gps2 = self.ps_hb * np.array([np.float64(sub_arr[3]) for sub_arr in self.tt['pts_hb']])
            self.pts_hb_gps3 = self.ps_hb * np.array([np.float64(sub_arr[4]) for sub_arr in self.tt['pts_hb']])
            self.pts_hb_gps4 = self.ps_hb * np.array([np.float64(sub_arr[5]) for sub_arr in self.tt['pts_hb']])

            self.pts_hc_gps1 = self.ps_h * np.array([np.float64(sub_arr[2]) for sub_arr in self.tt['pts_hc']])
            self.pts_hc_gps2 = self.ps_h * np.array([np.float64(sub_arr[3]) for sub_arr in self.tt['pts_hc']])
            self.pts_hc_gps3 = self.ps_h * np.array([np.float64(sub_arr[4]) for sub_arr in self.tt['pts_hc']])
            self.pts_hc_gps4 = self.ps_h * np.array([np.float64(sub_arr[5]) for sub_arr in self.tt['pts_hc']])

            self.pts_hs_gps1 = self.ps_h * np.array([np.float64(sub_arr[2]) for sub_arr in self.tt['pts_hs']])
            self.pts_hs_gps2 = self.ps_h * np.array([np.float64(sub_arr[3]) for sub_arr in self.tt['pts_hs']])
            self.pts_hs_gps3 = self.ps_h * np.array([np.float64(sub_arr[4]) for sub_arr in self.tt['pts_hs']])
            self.pts_hs_gps4 = self.ps_h * np.array([np.float64(sub_arr[5]) for sub_arr in self.tt['pts_hs']])

            self.N = self.ps_phi * np.float64(self.tt['N'])  # MPa
            self.pts_N_lake = np.float64(self.phi_0[self.pts_ni][0] - self.pts_phi_lake)
            self.pts_N_moulin = np.float64(self.phi_0[self.pts_ni][1] - self.pts_phi_moulin)

            self.V_b = self.ps_x**2 * self.ps_hb * np.float64(self.tt['Vb'])
            self.V_c = self.ps_x * self.ps_S * np.float64(self.tt['S'])
            self.V_s = self.ps_x**2 * self.ps_h * np.float64(self.tt['hs'])
            self.R_b = self.ps_x * np.float64(self.tt['Rb'])

            self.xx = self.ps_x/1e3 * np.float64(self.gg['nx'].item())
            self.yy = self.ps_x/1e3 * np.float64(self.gg['ny'].item())
            self.nI = int(self.gg['nI'].item())
            self.nJ = int(self.gg['nJ'].item())
            self.eI = int(self.gg['eI'].item())
            self.eJ = int(self.gg['eJ'].item())
            self.cI = int(self.gg['cI'].item())
            self.cJ = int(self.gg['cJ'].item())

            self.nmeanx = self.gg['nmeanx'].item()
            self.nmeany = self.gg['nmeany'].item()
            self.nmeans = self.gg['nmeans'].item()
            self.nmeanr = self.gg['nmeanr'].item()
            self.ein = self.gg['ein'].item()
            self.fin = self.gg['fin'].item()
            self.cin = self.gg['cin'].item()
            self.Dx = np.float64(self.gg['Dx'].item())
            self.Dy = np.float64(self.gg['Dy'].item())
            self.Ds = np.float64(self.gg['Ds'].item())
            self.Dr = np.float64(self.gg['Dr'].item())

            # --- get spatial pattern at different times -----------------------------
            tspan = self.ps_t/self.pd_td*self.oo['t_span']
            self.nums = np.array([np.argmin(np.abs(tspan - t_val))+1 for t_val in self.tslice])

            # run matlab script to extend the field for times
            filenames = [os.path.join(self.results_dir, self.casename, f"{num:04d}.mat") for num in self.nums]
            print(f"Processing files: {filenames}")
            # save the filenames to a text file
            # Save filenames and casename to a .mat file for MATLAB to load
            filenames_array = np.array(filenames, dtype=object)
            sio.savemat('extend_field.mat', {
                'filenames': filenames_array,
                'casename': self.casename
            })
            # Use MATLAB's -batch option to run the script and auto‐exit when done
            subprocess.run([
                self.matlab_path,
                "-nodisplay",
                "-nosplash",
                "-nodesktop",
                "-batch",
                "run('src/nevis_nodedischarge_python')"
            ], check=True)
            # Clear the temporary .mat file
            os.remove('extend_field.mat')

            self.vhb = []
            self.vphi1 = []
            self.vhs = []
            self.vhc = []
            self.vt = []
            self.vqs = []
            self.vqe = []
            self.vqQ = []
            self.vqb = []

            for i in range(len(self.nums)):
                con_path = os.path.join(self.results_dir, self.casename, f"{self.nums[i]:04d}.mat")
                raw_con = sio.loadmat(con_path, squeeze_me=True)
                vv = raw_con['vv']
                vv2 = raw_con['vv2']

                # geometry values
                phi1_val = np.float64(vv['phi'].item())
                hb_val = np.float64(vv['hb'].item())
                hs_val = np.float64(vv['hs'].item())
                pb_val = np.float64(vv['pb'].item())

                Sx_val = self.nmeanx.dot(np.float64(vv['Sx'].item()))
                Sy_val = self.nmeany.dot(np.float64(vv['Sy'].item()))
                Ss_val = self.nmeans.dot(np.float64(vv['Ss'].item()))
                Sr_val = self.nmeanr.dot(np.float64(vv['Sr'].item()))

                self.vphi1.append(np.array(np.reshape(phi1_val, (self.nI, self.nJ), order='F'), dtype=np.float64))
                self.vhb.append(np.array(np.reshape(hb_val, (self.nI, self.nJ), order='F'), dtype=np.float64))
                vpb = np.array(np.reshape(pb_val, (self.nI, self.nJ), order='F'), dtype=np.float64)
                self.vhs.append(np.array(np.reshape(hs_val, (self.nI, self.nJ), order='F'), dtype=np.float64))
                vSx = np.array(np.reshape(Sx_val, (self.nI, self.nJ), order='F'), dtype=np.float64)
                vSy = np.array(np.reshape(Sy_val, (self.nI, self.nJ), order='F'), dtype=np.float64)
                vSs = np.array(np.reshape(Ss_val, (self.nI, self.nJ), order='F'), dtype=np.float64)
                vSr = np.array(np.reshape(Sr_val, (self.nI, self.nJ), order='F'), dtype=np.float64)
                self.vhc.append(np.reshape(self.pp_c8*(Sx_val*self.Dx + Sy_val*self.Dy + Ss_val*self.Ds + Sr_val*self.Dr)/(self.Dx*self.Dy), (self.nI, self.nJ), order='F'))
                self.vt.append((self.ps_t / self.pd_td) * np.float64(vv['t'].item()))

                # flux values
                qs_val = np.float64(vv2['qs'].item())
                qe_val = np.float64(vv2['qe'].item())
                qQ_val = np.float64(vv2['qQ'].item())
                qb_val = np.float64(vv2['qb'].item())
                self.vqs.append(np.array(np.reshape(qs_val, (self.nI, self.nJ), order='F'), dtype=np.float64))
                self.vqe.append(np.array(np.reshape(qe_val, (self.nI, self.nJ), order='F'), dtype=np.float64))
                self.vqQ.append(np.array(np.reshape(qQ_val, (self.nI, self.nJ), order='F'), dtype=np.float64))
                self.vqb.append(np.array(np.reshape(qb_val, (self.nI, self.nJ), order='F'), dtype=np.float64))
            self.vt = np.array(self.vt, dtype=np.float64)-self.tdrainge_d
            
            self.data_loaded = True
            print("Data loading completed")

    def setup_plotting_style(self):
        """Set up all matplotlib plotting styles"""
        # ================= Plotting setup ======================================
        plt.style.use('classic')
        plt.rc('axes', prop_cycle=cycler('color', self.colors))
        plt.rcParams["figure.dpi"] = 600
        # plt.style.use('classic')
        plt.rcParams['font.family'] = 'DejaVu Sans'
        plt.rcParams['mathtext.fontset'] = 'dejavusans'
        plt.rcParams['font.size'] = 8

    def plot_spatial_panel(self, ax, panel_idx, panel_label, add_colorbar=False):
        """Plot a single spatial panel"""
        cmap = plt.cm.Blues
        
        if add_colorbar:
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="3%", pad=0.1)
        
        if self.plot_hb:
            pc = ax.pcolormesh(self.xx, self.yy,
                            self.ps_hb * self.vhb[panel_idx] + self.ps_h * self.vhs[panel_idx] + self.ps_h * self.vhc[panel_idx],
                            cmap=cmap, norm=self.norm, shading='gouraud')
            pc.set_edgecolor('face')
        elif self.plot_qs:
            pc = ax.pcolormesh(self.xx, self.yy,
                            self.ps_qs * self.vqs[panel_idx] + self.ps_qb * self.vqb[panel_idx] + self.ps_qs * self.vqQ[panel_idx],
                            cmap=cmap, norm=self.norm, shading='gouraud')
            pc.set_edgecolor('face')
        else:
            print("No valid plot type selected. Please set plot_hb or plot_qs to True.")
            return ax

        # Plot contour lines for the blister alone
        contour_lines = ax.contour(self.xx, self.yy, self.ps_hb * self.vhb[panel_idx], levels=[0.01], colors='blue', linewidths=1.0, linestyles='-')
        contour_lines_phi = ax.contour(self.xx, self.yy, self.ps_phi/1e6*self.vphi1[panel_idx], colors='k', linewidths=0.5, linestyles='-')
        
        # Add markers and labels only for first panel
        if panel_idx == 0:
            ax.plot(12.5, 10, 'wo', markersize=6, label='moulin')
            ax.plot(25.0, 10, 'bo', markersize=6, label='lake')
            ax.plot(30.0, 10, 'r^', markersize=6, label='GPS1')
            ax.plot(35.0, 10, 'g^', markersize=6, label='GPS2')
            ax.plot(40.0, 10, 'y^', markersize=6, label='GPS3')
            ax.legend(loc='lower left', fontsize=6, frameon=True, handletextpad=0.2, ncol=5, labelspacing=0.2, numpoints=1)
        else:
            ax.plot(25.0, 10, 'bo', markersize=6, label='lake')
        
        if add_colorbar:
            cbar = plt.colorbar(pc, cax=cax, orientation='vertical')
            cbar.set_label('thickness$\\;(m)$' if self.plot_hb else 'flux$\\;(m^2/s)$')
            cbar.ax.tick_params(labelsize=10)

        ax.text(0.025, 0.9, f'{panel_label} $t=$ ${self.vt[panel_idx]:.1f}$ d', transform=ax.transAxes,
                ha='left', va='top', fontsize=8, clip_on=False, color='black')
        ax.set_xlim(self.xa, self.xb)
        ax.set_ylim(self.ya, self.yb)
        ax.set_aspect('equal')
        ax.grid()

        manual_positions = [(15, 10), (25, 12.5), (35, 12.5), (40, 12.5), (42.5, 17.5)]
        labels = ax.clabel(contour_lines_phi, inline=False, fontsize=4, fmt='%1.1f MPa', manual=manual_positions)
        import matplotlib.patheffects as pe
        for txt in labels:
            txt.set_bbox({'facecolor':'none', 'edgecolor':'none', 'pad':0})
            txt.set_path_effects([pe.Stroke(linewidth=1.2, foreground='white'), pe.Normal()])
            txt.set_zorder(10)
        
        return ax

    def generate_figure(self, save_path=None):
        """Generate the complete Figure 3 plot"""
        if not self.data_loaded:
            self.load_data()
            
        self.setup_plotting_style()
        
        fig = plt.figure(figsize=(10, 4), dpi=100)
        gs = gridspec.GridSpec(3, 2, figure=fig, wspace=0.1, hspace=0.3)

        # ================= Panel (a): blister at time a =================
        ax1 = fig.add_subplot(gs[0, 0])
        self.plot_spatial_panel(ax1, 0, '(a)')
        ax1.set_ylabel('$y~(km)$')

        # ================= Panel (b): blister at time b ================
        ax2 = fig.add_subplot(gs[0, 1])
        self.plot_spatial_panel(ax2, 1, '(b)', add_colorbar=True)

        # ================= Panel (c): blister at time c ================
        ax3 = fig.add_subplot(gs[1, 0])
        self.plot_spatial_panel(ax3, 2, '(c)')
        ax3.set_xlabel('$x~(km)$')
        ax3.xaxis.set_label_coords(0.5, -0.2)
        ax3.set_ylabel('$y~(km)$')

        # ================= Panel (d): blister at time d ================
        ax4 = fig.add_subplot(gs[1, 1])
        self.plot_spatial_panel(ax4, 3, '(d)', add_colorbar=True)
        ax4.set_xlabel('$x~(km)$')
        ax4.xaxis.set_label_coords(0.5, -0.2)

        # ============== Panel (e): Flux (last row: spanning both columns) ==============
        # Panel (e) — combined h/S at the lake (bottom merged)
        ax5 = fig.add_subplot(gs[2, 0])
        ax5.set_visible(False)  # Hide original axes

        # Get positions of the axes above to align with
        pos_left = ax3.get_position()
        pos_right = ax4.get_position()
        pos_bottom = gs[2, :].get_position(fig)

        # Compute full available width
        full_width = pos_left.x1 - pos_left.x0

        # Shrink factor (e.g. 70% of full width)
        shrink_factor = 1
        ax5_width = full_width * shrink_factor

        # Center the shrunken panel horizontally
        ax5_x0 = pos_left.x0 + (full_width - ax5_width)/2
        ax5_pos = [ax5_x0, pos_bottom.y0, ax5_width, pos_bottom.height-0.02]

        # Colorbar axes to the right of the shrunken panel
        pad = 0.01
        colorbar_width = 0.015
        cax_pos = [ax5_x0 + ax5_width + pad, pos_bottom.y0, colorbar_width, pos_bottom.height]

        # Create new axes with the calculated positions
        ax5 = fig.add_axes(ax5_pos)

        ax5.plot(self.t-self.tdrainge_d, self.Q_b_in, '-', linewidth=1.0, label='lake influx')
        ax5.plot(self.t-self.tdrainge_d, self.Q_b_dec,'-', linewidth=1.0, label='leakage')
        ax5.plot(self.t-self.tdrainge_d, self.Q_out_b + self.Q_out_Q + self.Q_out_q - self.E, '-', linewidth=1.0, label=r'outflow variation')

        if '10m3s' in self.casename:
            ax5.plot(self.t-self.tdrainge_d, self.E, '-', linewidth=1.0, label='moulin influx')

        ax5.set_yscale('log')
        ax5.set_xlim(self.tmin_d-self.tdrainge_d, self.tmax_d-self.tdrainge_d+10)
        ax5.set_xlabel('$t~(d)$')
        ax5.xaxis.set_label_coords(0.5, -0.2)
        ax5.set_ylim(1e0, 1e3)
        ax5.set_xlim(-2, 10)
        ax5.set_ylabel(r'flux$~(m^3/s)$')
        ax5.text(0.025, 0.9, '(e) fluxes', transform=ax5.transAxes,
                ha='left', va='top', fontsize=8, clip_on=False)
        ax5.legend(loc='upper right', bbox_to_anchor=(1.0, 1.06), frameon=True, fontsize=7, ncol=2,
                  labelspacing=0.4, handletextpad=0.3, borderpad=0.2)
        ax5.grid(alpha=0.5)

        # ============== Panel (f): total volume and thickness ==============
        ax6 = fig.add_subplot(gs[2, 1])
        ax6.set_visible(False)  # Hide original axes

        # Get positions of the axes above to align with
        pos_right = ax4.get_position()
        pos_bottom = gs[2, 1].get_position(fig)

        # Compute full available width
        full_width = pos_right.x1 - pos_right.x0

        # Shrink factor (e.g. 70% of full width)
        shrink_factor = 1
        ax6_width = full_width * shrink_factor

        # Center the shrunken panel horizontally
        ax6_x0 = pos_right.x0 + (full_width - ax6_width)/2
        ax6_pos = [ax6_x0, pos_bottom.y0, ax6_width*0.95, pos_bottom.height-0.02]

        # Colorbar axes to the right of the shrunken panel
        pad = 0.01
        colorbar_width = 0.015
        cax_pos = [ax6_x0 + ax6_width + pad, pos_bottom.y0, colorbar_width, pos_bottom.height]

        # Create new axes with the calculated positions
        ax6 = fig.add_axes(ax6_pos)
        # ============== Plot thickness ==============
        ax6.plot(self.t-self.tdrainge_d, self.ps_hb*self.pts_hb_lake + self.ps_h*(self.pts_hc_gps1+self.pts_hs_gps1), 'b-', linewidth=1.0, label='lake', alpha=1.0)
        ax6.plot(self.t-self.tdrainge_d, self.ps_hb*self.pts_hb_gps2 + self.ps_h*(self.pts_hc_gps2+self.pts_hs_gps2), 'r-', linewidth=1.0, label='GPS1', alpha=1.0)
        ax6.plot(self.t-self.tdrainge_d, self.ps_hb*self.pts_hb_gps3 + self.ps_h*(self.pts_hc_gps3+self.pts_hs_gps3), 'g-', linewidth=1.0, label='GPS2', alpha=1.0)
        ax6.plot(self.t-self.tdrainge_d, self.ps_hb*self.pts_hb_gps4 + self.ps_h*(self.pts_hc_gps4+self.pts_hs_gps4), 'y-', linewidth=1.0, label='GPS3', alpha=1.0)

        ax6.set_xlim(self.tmin_d+3-self.tdrainge_d, self.tmax_d+10-self.tdrainge_d)
        ax6.set_xlabel('$t~(d)$')
        ax6.xaxis.set_label_coords(0.5, -0.2)
        ax6.set_ylim(0, 1.0)
        ax6.set_xlim(-2, 10)
        ax6.set_ylabel(r'thickness$~(m)$')
        ax6.text(0.025, 0.9, '(f)', transform=ax6.transAxes,
                ha='left', va='top', fontsize=8, clip_on=False)
        ax6.grid(alpha=1.0)

        # ============== Plot volume ==============
        ax6_twin = ax6.twinx()
        ax6_twin.plot(self.t-self.tdrainge_d, self.V_b/1e7, '--', color='dimgray', linewidth=1.0, label='$V_b/V_l$', alpha=0.8)
        ax6_twin.plot(self.t-self.tdrainge_d, self.R_b/5e4, '-', color='dimgray', linewidth=1.0, label='$R_b/L$', alpha=0.8)

        ax6_twin.set_xlabel('$t~(d)$')
        # ax6_twin.set_ylabel(r'$V_b/V_l$', fontsize=8)
        ax6_twin.set_ylim(0, 1.0)
        ax6_twin.tick_params(axis='y', labelsize=8)
        ax6_twin.grid(alpha=0.5)
        # make axis and label blue
        ax6_twin.spines['right'].set_color('dimgray')
        ax6_twin.yaxis.label.set_color('dimgray')
        ax6_twin.tick_params(axis='y', colors='dimgray')

        # combine legend entries from ax6 and its twin
        handles1, labels1 = ax6.get_legend_handles_labels()
        handles2, labels2 = ax6_twin.get_legend_handles_labels()
        plt.legend(handles1 + handles2, labels1 + labels2, loc='upper right',
            bbox_to_anchor=(1.0, 1.06), fontsize=7, ncol=3,
            labelspacing=0.1, handletextpad=0.3, borderpad=0.2, columnspacing=1.0,
            frameon=True, framealpha=0.5)

        if save_path:
            fig.savefig(save_path, dpi=600, transparent=False, bbox_inches='tight', pad_inches=0.1)
            print(f"Figure saved to {save_path}")
        
        return fig
    def generate_cross_section_figure(self, save_path=None):
        """Generate the cross-section figure (Figure 4)"""
        tslice_cross = self.tdrainge_d + np.array([-1, 0, 1, 2, 4])
        
        if not self.data_loaded or not np.array_equal(self.tslice, tslice_cross):
            self.tslice = tslice_cross
            self.data_loaded = False
            self.load_data()
            
        self.setup_plotting_style()
        
        # ================= Plotting setup ======================================
        fig = plt.figure(figsize=(4.5, 3), dpi=100)
        gs = gridspec.GridSpec(2, 2, width_ratios=[1, 0.05], wspace=0.05, hspace=0.3)

        ax1 = plt.subplot(gs[0, 0])
        ax2 = plt.subplot(gs[1, 0])
        cax = plt.subplot(gs[:, 1])

        pos1 = ax1.get_position()
        pos2 = ax2.get_position()
        poscax = cax.get_position()

        vertical_gap = -0.008
        # Set ax2's bottom to be just below ax1's top
        new_pos2_bottom = pos1.y0 - 0.8*pos2.height - vertical_gap
        ax2.set_position([pos2.x0, new_pos2_bottom, pos2.width, 0.8*pos2.height])
        cax.set_position([poscax.x0, pos2.y0+0.17, poscax.width, pos1.y0+ pos1.height - pos2.y0 - 0.2])

        plt.rcParams["figure.dpi"] = 100
        plt.rcParams['axes.linewidth'] = 0.5
        plt.rcParams['lines.linewidth'] = 0.6
        plt.rcParams['grid.linewidth'] = 0.3
        plt.rcParams['axes.labelsize'] = 7
        plt.rcParams['xtick.labelsize'] = 7
        plt.rcParams['ytick.labelsize'] = 7

        cmap = plt.cm.plasma
        norm = PowerNorm(gamma=1.0, vmin=np.min(self.vt), vmax=np.max(self.vt))
        sm = cm.ScalarMappable(norm=norm, cmap=cmap)
        sm.set_array([])

        # ================= Panel (a): blister profile (plan view) =================
        for i in range(len(self.nums)):
            color_i = cmap(norm(self.vt[i]))
            CS = ax1.contour(self.xx, self.yy, self.ps_hb * self.vhb[i], levels=[0.01], colors=[color_i])

        cbar = fig.colorbar(sm, cax=cax, orientation='vertical')
        cbar.set_label('time $(d)$')

        ax1.set_xlim(10, 50)
        ax1.set_ylim(5, 15)
        ax1.tick_params(axis='x', labelbottom=False)  # hide x ticks for top panel
        ax1.set_ylabel('$y~(km)$')
        ax1.set_aspect('equal', adjustable='box')
        ax1.text(0.025, 0.8, '(a) plan-view', transform=ax1.transAxes,
                 ha='left', va='bottom', fontsize=8, clip_on=False)
        ax1.grid(alpha=0.5)

        # ================= Panel (b): blister profile (cross section) =================
        p = []
        for i in range(len(self.nums)):
            color_i = cmap(norm(self.vt[i]))
            p.append(ax2.plot(self.xx[:,0], self.ps_hb * self.vhb[i][:,self.nJ//2], color=color_i))  # Plot each time step

        ax2.set_xlabel('$x~(km)$')
        ax2.xaxis.set_label_coords(0.5, -0.175)
        ax2.set_ylabel('$h_b~(m)$')
        ax2.set_xlim(10, 50)
        ax2.set_ylim(0, 0.75)
        ax2.text(0.025, 0.9, '(b) cross-section', transform=ax2.transAxes,
                 ha='left', va='top', fontsize=8, clip_on=False, color='black')
        ax2.grid(alpha=0.5)

        fig.tight_layout(rect=[0, 0, 0.95, 1]) # Adjust layout to make space for colorbar

        # if save_path:
        #     fig.savefig(save_path, dpi=600, transparent=False, bbox_inches='tight', pad_inches=0.1)
        #     print(f"Cross-section figure saved to {save_path}")

        return fig
    
# Example usage:
if __name__ == "__main__":

    plotter = NevisFig3Plotter('n2d_10m3s_kappa1e_10_mu1e01_hbreg5e_3_V1e7_drainage')
    fig1 = plotter.generate_figure(save_path='./figures/Fig3_summertime_poster.png')

    plotter = NevisFig3Plotter('n2d_0m3s_kappa1e_10_mu1e01_hbreg5e_3_V1e7_drainage')
    fig2 = plotter.generate_figure(save_path='./figures/Fig2_wintertime_poster.png')
    fig3 = plotter.generate_cross_section_figure(save_path='./figures/Fig4_cross_section_poster.png')
    # Show the plot
    plt.show()