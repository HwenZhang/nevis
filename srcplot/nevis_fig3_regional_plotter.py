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

        self.xa = -75
        self.xb = 60
        self.ya = -50
        self.yb = 25
        self.plot_hb = False
        self.plot_qs = True
        self.h_min = 0  # min water thickness on the color bar (m)
        self.h_max = 0.5  # max water thickness on the color bar (m)
        self.q_min = 1e-4  # min water flux on the color bar (m2/s)
        self.q_max = 1e-1  # max water flux on the color bar (
        # self.norm = PowerNorm(gamma=1, vmin=self.h_min, vmax=self.h_max)
        # self.norm = PowerNorm(gamma=1, vmin=self.q_min, vmax=self.q_max)
        self.norm = LogNorm(vmin=self.q_min, vmax=self.q_max)
        # self.norm = LogNorm(vmin=self.h_min, vmax=self.h_max)

        self.tdrainge_d = 30
        self.tmin = self.tdrainge_d - 5
        self.tmax = self.tdrainge_d + 40
        self.tmin_d = self.tmin 
        self.tmax_d = self.tmax
        self.tslice = self.tdrainge_d + np.array([0, 2, 20, 35])
        
        # Data containers
        self.data_loaded = False
        
    def load_terrain(self):
        """Load terrain data from file"""
        mat_path = os.path.join('./data/nevis_regional/morlighem_for_nevis_140km.mat')
        raw = sio.loadmat(mat_path, squeeze_me=True)
        morlighem_for_nevis_140km = raw['morlighem_for_nevis_140km']
        self.X_terrain = morlighem_for_nevis_140km['X_km'].item()
        self.Y_terrain = morlighem_for_nevis_140km['Y_km'].item()
        self.B_terrain = morlighem_for_nevis_140km['B_km'].item()
        self.S_terrain = morlighem_for_nevis_140km['S_km'].item()

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
            self.aa = raw['aa']

            self.ps_t = np.float64(self.ps['t'])
            self.pd_td = np.float64(self.pd['td'])
            self.ps_hb = np.float64(self.ps['hb'])
            self.ps_h = np.float64(self.ps['h'])
            self.ps_x = np.float64(self.ps['x'])
            self.ps_z = np.float64(self.ps['z'])
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
            self.pts_phi_lake = self.ps_phi * np.array([np.float64(sub_arr[-1]) for sub_arr in self.tt['pts_phi']])
            self.pts_phi_moulin = self.ps_phi * np.array([np.float64(sub_arr[-2]) for sub_arr in self.tt['pts_phi']])
            self.pts_hb_lake = self.ps_hb * np.array([np.float64(sub_arr[0]) for sub_arr in self.tt['pts_hb']])
            self.pts_hb_moulin = self.ps_hb * np.array([np.float64(sub_arr[1]) for sub_arr in self.tt['pts_hb']])
            self.pts_pb_lake = self.ps_phi * np.array([np.float64(sub_arr[0]) for sub_arr in self.tt['pts_pb']])
            self.pts_hc_lake = self.ps_h * np.array([np.float64(sub_arr[0]) for sub_arr in self.tt['pts_hc']])
            self.pts_hs_lake = self.ps_h * np.array([np.float64(sub_arr[0]) for sub_arr in self.tt['pts_hs']])

            self.pts_hb_gps1 = self.ps_hb * np.array([np.float64(sub_arr[-1]) for sub_arr in self.tt['pts_hb']])
            self.pts_hb_gps2 = self.ps_hb * np.array([np.float64(sub_arr[-2]) for sub_arr in self.tt['pts_hb']])
            self.pts_hb_gps3 = self.ps_hb * np.array([np.float64(sub_arr[-3]) for sub_arr in self.tt['pts_hb']])
            self.pts_hb_gps4 = self.ps_hb * np.array([np.float64(sub_arr[-4]) for sub_arr in self.tt['pts_hb']])
            self.pts_hb_gps5 = self.ps_hb * np.array([np.float64(sub_arr[-5]) for sub_arr in self.tt['pts_hb']])

            self.pts_hc_gps1 = self.ps_h * np.array([np.float64(sub_arr[-1]) for sub_arr in self.tt['pts_hc']])
            self.pts_hc_gps2 = self.ps_h * np.array([np.float64(sub_arr[-2]) for sub_arr in self.tt['pts_hc']])
            self.pts_hc_gps3 = self.ps_h * np.array([np.float64(sub_arr[-3]) for sub_arr in self.tt['pts_hc']])
            self.pts_hc_gps4 = self.ps_h * np.array([np.float64(sub_arr[-4]) for sub_arr in self.tt['pts_hc']])
            self.pts_hc_gps5 = self.ps_h * np.array([np.float64(sub_arr[-5]) for sub_arr in self.tt['pts_hc']])

            self.pts_hs_gps1 = self.ps_h * np.array([np.float64(sub_arr[-1]) for sub_arr in self.tt['pts_hs']])
            self.pts_hs_gps2 = self.ps_h * np.array([np.float64(sub_arr[-2]) for sub_arr in self.tt['pts_hs']])
            self.pts_hs_gps3 = self.ps_h * np.array([np.float64(sub_arr[-3]) for sub_arr in self.tt['pts_hs']])
            self.pts_hs_gps4 = self.ps_h * np.array([np.float64(sub_arr[-4]) for sub_arr in self.tt['pts_hs']])
            self.pts_hs_gps5 = self.ps_h * np.array([np.float64(sub_arr[-5]) for sub_arr in self.tt['pts_hs']])

            self.N = self.ps_phi * np.float64(self.tt['N'])  # MPa
            self.pts_N_lake = np.float64(self.phi_0[self.pts_ni][-1] - self.pts_phi_lake)
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

            # auxiliary fields
            H_val = np.float64(self.aa['H'].item())
            self.H = np.array(np.reshape(H_val, (self.nI, self.nJ), order='F'), dtype=np.float64)

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
        plt.rcParams['axes.linewidth'] = 1.2

    def plot_spatial_panel(self, ax, panel_idx, panel_label, add_colorbar=False):
        """Plot a single spatial panel"""
        cmap = plt.cm.Blues
        from matplotlib.colors import TwoSlopeNorm
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

        # Add the thickness map self.H*ps.z as the background
        norm = PowerNorm(gamma=1, vmin=-1.0, vmax=1.0)
        if panel_idx == 0:
            terrain = ax.pcolormesh(self.X_terrain/1e3, self.Y_terrain/1e3, self.B_terrain/1e3 + 0.0*(self.S_terrain-self.B_terrain)/1e3, norm=norm, cmap='terrain', alpha=1.0)
        else:
            terrain = ax.pcolormesh(self.X_terrain/1e3, self.Y_terrain/1e3, self.B_terrain/1e3 + 0.0*(self.S_terrain-self.B_terrain)/1e3, norm=norm, cmap='terrain', alpha=0.3)

        from mpl_toolkits.axes_grid1.inset_locator import inset_axes
        cbax = inset_axes(ax,
                 width="30%",   
                 height="7%",  
                 loc='lower left',
                 bbox_to_anchor=(0.65, 0.90, 1.0, 1.0),
                 bbox_transform=ax.transAxes,
                 borderpad=0)
        cbar_terr = plt.colorbar(terrain, cax=cbax, orientation='horizontal')
        cbar_terr.set_label('Bed elevation (km)')
        cbar_terr.set_ticks(np.array([-1.0, -0.5, 0.0, 0.5, 1.0]))
        cbar_terr.ax.tick_params(labelsize=6)

        # Plot contour lines for the blister alone
        # contour_lines = ax.contour(self.xx, self.yy, self.ps_hb * self.vhb[panel_idx], levels=[0.01], colors='blue', linewidths=1.0, linestyles='-')
        contour_lines = ax.contour(self.xx, self.yy, self.ps_phi/1e6*self.vphi1[panel_idx],
                       colors='k', linewidths=0.5, linestyles='-', label='Contour lines')
        
        ax.contour(
            self.xx, self.yy,
            self.ps_hb * self.vhb[panel_idx],
            levels=[0.1],
            colors='blue',
            linewidths=1.0,
            linestyles='-'
        )

        # Add markers and labels only for first panel
        if panel_idx == 0:
            ax.plot(40, -20, 'bo', markersize=6, label='lake')
            # ax.plot(40, -5, 'r^', markersize=6, label='GPS1')
            ax.plot(20, -15, 'r^', markersize=6, label='S1')
            ax.plot(0, -25, 'y^', markersize=6, label='S2')
            ax.plot(-40, -30, 'g^', markersize=6, label='S3')
        else:
            ax.plot(40, -20, 'bo', markersize=6)
            # ax.plot(40, -5, 'r^', markersize=6)
            ax.plot(20, -15, 'r^', markersize=6)
            ax.plot(0, -25, 'y^', markersize=6)
            ax.plot(-40, -30, 'g^', markersize=6)

        if add_colorbar:
            cbar = plt.colorbar(pc, cax=cax, orientation='vertical')
            cbar.set_label('$h\\;(m)$' if self.plot_hb else '$q\\;(m^2/s)$')
            cbar.ax.tick_params(labelsize=10)
            from matplotlib.ticker import ScalarFormatter
            # formatter = ScalarFormatter(useMathText=True)
            # formatter.set_scientific(True)
            # formatter.set_powerlimits((0, 0))  # 强制使用科学计数法
            # cbar.ax.yaxis.set_major_formatter(formatter)
            # cbar.ax.tick_params(labelsize=10)

        ax.text(0.025, 0.9, f'{panel_label} $t=$ ${self.vt[panel_idx]:.1f}$ d', transform=ax.transAxes,
                ha='left', va='top', fontsize=10, clip_on=False, color='black')
        ax.set_xlim(self.xa, self.xb)
        ax.set_ylim(self.ya, self.yb)
        ax.set_aspect('equal')
        ax.grid()
        if panel_idx == 0:
            ax.legend(loc='center left', fontsize=8, frameon=True, 
                    handletextpad=0.3, ncol=1, labelspacing=0.3, numpoints=1)
        # else:
            # ax.legend(loc='lower right', fontsize=8, frameon=True, 
            #         handletextpad=0.3, ncol=2, labelspacing=0.3, numpoints=1)
        return ax

    def generate_figure(self, save_path=None):
        """Generate the complete Figure 3 plot"""
        if not self.data_loaded:
            self.load_data()
            
        self.setup_plotting_style()
        
        self.load_terrain()
        fig = plt.figure(figsize=(10, 8), dpi=300)
        gs = gridspec.GridSpec(3, 2, figure=fig, wspace=0.05, hspace=0.2)

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
        ax3.xaxis.set_label_coords(0.5, -0.1)
        ax3.set_ylabel('$y~(km)$')

        # ================= Panel (d): blister at time d ================
        ax4 = fig.add_subplot(gs[1, 1])
        self.plot_spatial_panel(ax4, 3, '(d)', add_colorbar=True)
        ax4.set_xlabel('$x~(km)$')
        ax4.xaxis.set_label_coords(0.5, -0.1)

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
        ax5_pos = [ax5_x0, pos_bottom.y0+0.3*pos_bottom.height, ax5_width, 0.65*pos_bottom.height]

        # Create new axes with the calculated positions
        ax5 = fig.add_axes(ax5_pos)

        ax5.plot(self.t-self.tdrainge_d, self.Q_b_in, '-', linewidth=1.5, label='injection')
        # ax5.plot(self.t-self.tdrainge_d, self.Q_b_dec,'-', linewidth=1.0, label='leakage')
        ax5.plot(self.t-self.tdrainge_d, 0*self.Q_out_q+self.Q_out_Q, '-', linewidth=1.5, label='baseflow')
        ax5.plot(self.t-self.tdrainge_d, self.Q_out_b, '-', linewidth=1.5, label='blister outflow')

        drainage_start_ind = np.where(self.Q_b_in > 1e-6)[0][0]
        drainage_end_ind = np.where(self.Q_b_in > 1e-6)[0][-1]
        # check if empty
        if np.where(self.Q_out_b > 1e-1)[0].size > 0:
            reaching_boundary_ind = np.where(self.Q_out_b > 1e-1)[0][0]
        else:
            reaching_boundary_ind = 0

        ax5.axvspan(self.t[drainage_start_ind]-self.tdrainge_d, self.t[drainage_end_ind]-self.tdrainge_d, color='blue', alpha=0.2)
        ax5.axvspan(self.t[drainage_end_ind]-self.tdrainge_d, self.t[reaching_boundary_ind]-self.tdrainge_d, color='orange', alpha=0.2)
        # ax5.text(self.t[drainage_start_ind]-self.tdrainge_d-2, 2e0, 'drainage', color='blue', fontsize=8)
        ax5.text(self.t[drainage_end_ind]-self.tdrainge_d+10, 2e0, 'propagation', color='orange', fontsize=8)
        ax5.text(self.t[reaching_boundary_ind]-self.tdrainge_d+8, 2e0, 'outflow', color='green', fontsize=8)
        
        ax5.axvspan(self.t[reaching_boundary_ind]-self.tdrainge_d, self.tmax_d, color='green', alpha=0.2)

        ax5.set_yscale('log')
        ax5.set_xlim(self.tmin_d-self.tdrainge_d, self.tmax_d-self.tdrainge_d+10)
        ax5.set_xlabel('$t~(d)$')
        ax5.xaxis.set_label_coords(0.5, -0.2)
        ax5.set_ylim(1e0, 1e2 )
        ax5.set_xlim(-5, 60)
        ax5.set_ylabel(r'$Q~(m^3/s)$')
        ax5.text(0.025, 0.9, '(e) fluxes', transform=ax5.transAxes,
                ha='left', va='top', fontsize=10, clip_on=False)
        ax5.legend(loc='upper right', bbox_to_anchor=(1.015, 1.03), frameon=True, fontsize=7, ncol=1,
                  labelspacing=0.3, handletextpad=0.3, borderpad=0.2)
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
        shrink_factor = 0.95
        ax6_width = full_width * shrink_factor

        # Center the shrunken panel horizontally
        ax6_x0 = pos_right.x0 + 0.05*ax6_width
        ax6_pos = [ax6_x0, pos_bottom.y0+0.3*pos_bottom.height, ax6_width*0.95, 0.65*pos_bottom.height]

        # Colorbar axes to the right of the shrunken panel
        pad = 0.01
        colorbar_width = 0.015
        cax_pos = [ax6_x0 + ax6_width + pad, pos_bottom.y0, colorbar_width, pos_bottom.height]

        # Create new axes with the calculated positions
        ax6 = fig.add_axes(ax6_pos)
        # ============== Plot thickness ==============
        ax6.plot(self.t-self.tdrainge_d,\
                self.ps_hb*self.pts_hb_lake + 0*self.ps_h*(self.pts_hc_gps1+self.pts_hs_gps1),\
                'b-', linewidth=1.0, label='lake', alpha=1.0)
        # ax6.plot(self.t-self.tdrainge_d,\
        #         self.ps_hb*self.pts_hb_gps4 + 0*self.ps_h*(self.pts_hc_gps4+self.pts_hs_gps4),\
        #         'r-', linewidth=1.0, label='GPS1', alpha=1.0)
        ax6.plot(self.t-self.tdrainge_d,\
                self.ps_hb*self.pts_hb_gps3 + 0*self.ps_h*(self.pts_hc_gps3+self.pts_hs_gps3),\
                'r-', linewidth=1.0, label='S1', alpha=1.0)
        ax6.plot(self.t-self.tdrainge_d,\
                self.ps_hb*self.pts_hb_gps2 + 0*self.ps_h*(self.pts_hc_gps2+self.pts_hs_gps2),\
                'y-', linewidth=1.0, label='S2', alpha=1.0)
        ax6.plot(self.t-self.tdrainge_d,\
                self.ps_hb*self.pts_hb_gps1 + 0*self.ps_h*(self.pts_hc_gps1+self.pts_hs_gps1),\
                'g-', linewidth=1.0, label='S3', alpha=1.0)
        
        # ax6.plot(self.t-self.tdrainge_d,\
        #         self.ps_hb*self.pts_hb_gps5 + self.ps_h*(self.pts_hc_gps5+self.pts_hs_gps5),\
        #         'y-', linewidth=1.0, label='GPS5', alpha=1.0)

        ax6.set_xlabel('$t~(d)$')   
        ax6.xaxis.set_label_coords(0.5, -0.2)
        ax6.set_ylim(0, 2.0)
        ax6.set_xlim(-5, 40)
        ax6.set_ylabel(r'$h~(m)$')
        ax6.text(0.025, 0.9, '(f) time series', transform=ax6.transAxes,
                ha='left', va='top', fontsize=10, clip_on=False)
        ax6.grid(alpha=1.0)

        # ============== Plot volume ==============
        ax6_twin = ax6.twinx()
        ax6_twin.plot(self.t-self.tdrainge_d, self.V_b/(1.8e8), '-', color='grey', linewidth=1.0, label='$V_b/V_l$', alpha=0.8)
        # ax6_twin.plot(self.t-self.tdrainge_d, self.R_b/5e4, '-', color='dimgray', linewidth=1.0, label='$R_b/L$', alpha=0.8)

        ax6_twin.set_xlabel('$t~(d)$')
        ax6_twin.set_ylim(0, 2.0)
        ax6_twin.set_yticks([0, 0.5, 1.0])
        ax6_twin.tick_params(axis='y', labelsize=8)
        # ax6_twin.grid(alpha=0.5)
        # make axis and label blue
        ax6_twin.spines['right'].set_color('grey')
        ax6_twin.yaxis.label.set_color('grey')
        ax6_twin.tick_params(axis='y', colors='grey')

        # combine legend entries from ax6 and its twin
        handles1, labels1 = ax6.get_legend_handles_labels()
        handles2, labels2 = ax6_twin.get_legend_handles_labels()
        plt.legend(handles1 + handles2, labels1 + labels2, loc='upper right',
            bbox_to_anchor=(1.0, 0.95), fontsize=7, ncol=3,
            labelspacing=0.1, handletextpad=0.3, borderpad=0.2, columnspacing=1.0,
            frameon=True, framealpha=0.5)

        if save_path:
            fig.savefig(save_path, dpi=600, transparent=False, bbox_inches='tight', pad_inches=0.1)
            print(f"Figure saved to {save_path}")
        
        return fig
    
# Example usage:
if __name__ == "__main__":
    # casename = 'n2d_region_test_meanperms1_Hreg1000_kappa1e_99_mu2e1_spinup'
    casename = 'n2d_region_melt_meanperms1_Hreg1000_kappa1e_99_mu2e1_td30_V1e8_drainage'
    plotter = NevisFig3Plotter(casename)
    fig1 = plotter.generate_figure(save_path='./figures/Fig7_regional.png')#
 