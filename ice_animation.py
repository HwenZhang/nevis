## Animation: blister thickness (left) and ice velocity (right) vs time
## Python port of nevis_2d_animation_ice.m

import os, glob, subprocess
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import TwoSlopeNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.animation import FuncAnimation, FFMpegWriter

# ══════════════════════════════════════════════════════════════════════
# Settings
# ══════════════════════════════════════════════════════════════════════
casename = 'n2d_regional_racmo_eps1e_02_kappa5e_11_mu5e0_partition5e_01_k01e_01_spinup'
tmin_yr = 0.0
tmax_yr = tmin_yr + 1.0
tmin = tmin_yr * 365  # days
tmax = tmax_yr * 365
clusters_to_plot = [1, 2, 3, 4, 5]
FRAME_STEP = 2          # skip every N frames for speed (1 = every frame)
FPS = 5

# ══════════════════════════════════════════════════════════════════════
# Load base simulation data
# ══════════════════════════════════════════════════════════════════════
mat_path = os.path.join('results', casename, casename + '.mat')
raw = sio.loadmat(mat_path, squeeze_me=True)
ps = raw['ps']; pd_r = raw['pd']; pp = raw['pp']
oo = raw['oo']; gg = raw['gg']; aa = raw['aa']; tt = raw['tt']

ps_t   = float(ps['t']);   pd_td  = float(pd_r['td']); pd_ty = float(pd_r['ty'])
ps_h   = float(ps['h']);   ps_hb  = float(ps['hb'])
ps_ub  = float(ps['u_b']); ps_phi = float(ps['phi']); ps_x = float(ps['x'])
ps_S   = float(ps['S']);   ps_hs  = float(ps['hs'])
ps_Q   = float(ps['Q']);   ps_Q0  = float(ps['Q0'])
ps_m   = float(ps['m'])
nI = int(gg['nI'].item()); nJ = int(gg['nJ'].item())

xx = (ps_x / 1e3) * np.float64(gg['nx'].item())  # (nI, nJ) in km
yy = (ps_x / 1e3) * np.float64(gg['ny'].item())

# time series (days)
n_time = len(tt)
t_ts = (ps_t / pd_td) * np.array([float(tt['t'][i]) for i in range(n_time)])

# flux time series
Q_out_Q = ps_Q * np.array([float(tt['Q_outQ'][i]) for i in range(n_time)])
Q_out_q = ps_Q * np.array([float(tt['Q_outq'][i]) for i in range(n_time)])
Q_out_b = ps_Q0 * np.array([float(tt['Q_outb'][i]) for i in range(n_time)])
E_ts    = ps_m * ps_x**2 * np.array([float(tt['E'][i]) for i in range(n_time)])

# effective pressure, sheet/channel, ice speed time series
N_ts  = ps_phi * np.array([float(tt['N'][i]) for i in range(n_time)])
hs_ts = ps_x**2 * ps_h * np.array([float(tt['hs'][i]) for i in range(n_time)])
S_ts  = ps_x * ps_S * np.array([float(tt['S'][i]) for i in range(n_time)])
U_ts  = ps_ub * np.array([float(tt['U'][i]) for i in range(n_time)])
V_b   = ps_x**2 * ps_hb * np.array([float(tt['Vb'][i]) for i in range(n_time)])

A_total = ps_x**2 * np.sum(np.float64(gg['Dx'].item()) * np.float64(gg['Dy'].item()))

# GPS pts at pts_ni
pts_ni = oo['pts_ni'].item().astype(int) - 1
pts_hb = np.column_stack([tt['pts_hb'][i] for i in range(n_time)]) * ps_hb
pts_hs = np.column_stack([tt['pts_hs'][i] for i in range(n_time)]) * ps_h
pts_u  = np.column_stack([tt['pts_u'][i]  for i in range(n_time)]) * ps_ub * pd_ty  # m/yr

# GPS cluster assignment
ni_gps = pp['ni_gps'].item().astype(int)
ni_gps_0 = ni_gps - 1
xx_flat = (ps_x / 1e3) * np.float64(gg['nx'].item()).flatten(order='F')
yy_flat = (ps_x / 1e3) * np.float64(gg['ny'].item()).flatten(order='F')
gps_x_km = xx_flat[ni_gps_0]
gps_y_km = yy_flat[ni_gps_0]

gps_rows = []
gps_matched = []
for ig, g in enumerate(ni_gps):
    idx = np.where(oo['pts_ni'].item().astype(int) == g)[0]
    if len(idx) > 0:
        gps_rows.append(idx[0])
        gps_matched.append(ig)
gps_rows = np.array(gps_rows)
gps_matched = np.array(gps_matched)

cluster_centers = np.array([[0, 0], [20, -15], [40, -25]])
gps_cluster = np.zeros(len(ni_gps), dtype=int)
for ig in range(len(ni_gps)):
    dists = np.sum((cluster_centers - [gps_x_km[ig], gps_y_km[ig]])**2, axis=1)
    gps_cluster[ig] = np.argmin(dists) + 1
dist_c2 = np.sqrt((gps_x_km - 17)**2 + (gps_y_km + 25)**2)
gps_cluster[np.argmin(dist_c2)] = 4
dist_c1 = np.sqrt((gps_x_km - 5)**2 + (gps_y_km + 10)**2)
gps_cluster[np.argmin(dist_c1)] = 5
# Filter to only GPS stations found in pts_ni
gps_cluster = gps_cluster[gps_matched]
gps_x_km = gps_x_km[gps_matched]
gps_y_km = yy_flat[ni_gps_0][gps_matched]

print(f'GPS: {len(ni_gps)} stations, {len(gps_rows)} matched pts_ni, clusters: {np.unique(gps_cluster)}')

cluster_colors = [
    (1.0, 0.45, 0.20),   # C1 bright orange
    (0.55, 1.0, 0.55),   # C2 light green
    (0.40, 0.75, 1.0),   # C3 light blue
    (0.85, 0.45, 1.0),   # C4 light purple
    (1.0, 0.85, 0.30),   # C5 bright yellow
]
cluster_labels = ['C1', 'C2', 'C3', 'C4', 'C5']

# drainage event times
t_drainage = np.asarray(pp['t_drainage'].item()) * ps_t / pd_td
t_duration = np.asarray(pp['t_duration'].item()) * ps_t / pd_td

# lake positions (dimensional km)
x_l = np.asarray(pp['x_l'].item()) * (ps_x / 1e3)
y_l = np.asarray(pp['y_l'].item()) * (ps_x / 1e3)
V_l = np.asarray(pp['V_l'].item())

def lake_input(t_nondim, pp_obj, ps_obj):
    """Gaussian lake input flux (dimensional)."""
    t_dr = np.asarray(pp_obj['t_drainage'].item())
    t_du = np.asarray(pp_obj['t_duration'].item())
    V_l_ = np.asarray(pp_obj['V_l'].item())
    Q0_  = float(ps_obj['Q0'])
    return Q0_ * V_l_ / (np.sqrt(2 * np.pi) * t_du) * \
           np.exp(-0.5 * ((t_nondim - t_dr) / t_du)**2)

# ══════════════════════════════════════════════════════════════════════
# Find frames in time range
# ══════════════════════════════════════════════════════════════════════
tspan_d = (ps_t / pd_td) * oo['t_span']
frame_indices = np.where((tspan_d >= tmin) & (tspan_d <= tmax))[0]
frame_indices = frame_indices[::FRAME_STEP]
print(f'Found {len(frame_indices)} frames in [{tmin:.0f}, {tmax:.0f}] days')

filepath = os.path.join('results', casename, '')

# ══════════════════════════════════════════════════════════════════════
# Export ice fields via MATLAB (batch)
# ══════════════════════════════════════════════════════════════════════
# Collect all filenames for this range
all_filenames = [os.path.join(filepath, f'{idx+1:04d}.mat') for idx in frame_indices]

# Check which ones already have ice_fields exported
need_export = []
for fn in all_filenames:
    try:
        keys = sio.whosmat(fn)
        names = [k[0] for k in keys]
        if 'ice_fields' not in names:
            need_export.append(fn)
    except Exception:
        need_export.append(fn)

if need_export:
    print(f'Exporting ice fields for {len(need_export)} frames via MATLAB...')
    sio.savemat('export_ice_fields.mat', {
        'filenames': np.array(need_export, dtype=object),
        'casename': casename,
    })
    matlab_path = "/Applications/MATLAB_R2025b.app/bin/matlab"
    subprocess.run([
        matlab_path, "-nodisplay", "-nosplash", "-nodesktop",
        "-batch", "run('src/nevis_export_ice_fields')"
    ], check=True)
    os.remove('export_ice_fields.mat')
    print('MATLAB export done.')
else:
    print('All frames already have ice_fields exported.')

# ══════════════════════════════════════════════════════════════════════
# Read first frame to set up reference fields
# ══════════════════════════════════════════════════════════════════════
def load_frame(fn):
    """Read exported ice_fields + vv from a .mat file."""
    raw_f = sio.loadmat(fn, squeeze_me=True)
    ice = raw_f['ice_fields']
    vv = raw_f['vv']
    vel_scale = float(ice['vel_scale'])
    stress_scale = float(ice['stress_scale'])
    nout = ice['nout'].item().astype(int) - 1
    uxn = np.float64(ice['uxn'].item())
    vyn = np.float64(ice['vyn'].item())
    Un = np.float64(ice['Un'].item())
    vux = vel_scale * uxn.reshape(nI, nJ, order='F')
    vuy = vel_scale * vyn.reshape(nI, nJ, order='F')
    vUU = vel_scale * Un; vUU[nout] = np.nan
    vUU = vUU.reshape(nI, nJ, order='F')
    # blister
    hb = ps_hb * np.float64(vv['hb'].item()).reshape(nI, nJ, order='F')
    pb = ps_phi * np.float64(vv['pb'].item()).reshape(nI, nJ, order='F')
    t_days = float(ice['t_days'])
    t_nondim = float(vv['t'].item())
    return dict(vux=vux, vuy=vuy, vUU=vUU, hb=hb, pb=pb,
                t_days=t_days, t_nondim=t_nondim)

# ══════════════════════════════════════════════════════════════════════
# Build figure layout
# ══════════════════════════════════════════════════════════════════════
plt.style.use('classic')
plt.rcParams.update({
    'font.family': 'DejaVu Sans',
    'mathtext.fontset': 'dejavusans',
    'font.size': 8,
    'figure.dpi': 120,
})

# BWR colormap for blister
n_cmap = 256
bwr_arr = np.vstack([
    np.column_stack([np.linspace(0, 1, n_cmap), np.linspace(0, 1, n_cmap), np.ones(n_cmap)]),
    np.column_stack([np.ones(n_cmap), np.linspace(1, 0, n_cmap), np.linspace(1, 0, n_cmap)]),
])
bwr_cmap = mcolors.ListedColormap(bwr_arr)

ref = load_frame(all_filenames[0])

fig = plt.figure(figsize=(16, 9), dpi=120)
gs = fig.add_gridspec(7, 5, hspace=0.45, wspace=0.45)

# ── Left column: time series panels ──
ax_a = fig.add_subplot(gs[0, 0:2])   # (a) fluxes
ax_b = fig.add_subplot(gs[1, 0:2])   # (b) effective pressure
ax_c = fig.add_subplot(gs[2, 0:2])   # (c) sheet thickness
ax_d = fig.add_subplot(gs[3, 0:2])   # (d) ice speed + blister volume
ax_e = fig.add_subplot(gs[4, 0:2])   # (e) hb+hs at GPS clusters
ax_f = fig.add_subplot(gs[5, 0:2])   # (f) U at GPS clusters
ax_g = fig.add_subplot(gs[6, 0:2])   # (g) placeholder / blank

ts_axes = [ax_a, ax_b, ax_c, ax_d, ax_e, ax_f, ax_g]

# ── Right column: spatial panels ──
ax_hb = fig.add_subplot(gs[0:4, 2:5])   # blister hb
ax_vel = fig.add_subplot(gs[4:7, 2:5])  # ice velocity

# ── (a) Fluxes ──
for i_d in range(len(t_drainage)):
    ax_a.axvspan(t_drainage[i_d] - 0.5*t_duration[i_d],
                 t_drainage[i_d] + 0.5*t_duration[i_d],
                 color=(0.7, 0.85, 1.0), alpha=0.4)
ax_a.semilogy(t_ts, Q_out_Q + Q_out_q, '-', color=(0, 0.25, 0), lw=2, label=r'$Q_\mathrm{out}$')
ax_a.semilogy(t_ts, Q_out_b, 'b--', lw=1.2, label=r'$Q_\mathrm{out,b}$')
ax_a.semilogy(t_ts, E_ts, 'k-.', lw=1.2, label=r'$Q_\mathrm{in}$')
vline_a = ax_a.axvline(ref['t_days'], color='k', ls='--', lw=1)
ax_a.set_xlim(tmin, tmax); ax_a.set_ylim(1e-2, 1e4)
ax_a.set_ylabel(r'$Q$ (m$^3$/s)', fontsize=7)
ax_a.legend(fontsize=6, ncol=3, loc='upper right')
ax_a.text(0.02, 0.85, '(a) flux', transform=ax_a.transAxes, fontsize=7)
ax_a.grid(True, which='both', alpha=0.3)

# ── (b) Effective pressure ──
for i_d in range(len(t_drainage)):
    ax_b.axvspan(t_drainage[i_d] - 0.5*t_duration[i_d],
                 t_drainage[i_d] + 0.5*t_duration[i_d],
                 color=(0.7, 0.85, 1.0), alpha=0.4)
ax_b.plot(t_ts, N_ts / 1e6, 'b-', lw=1.2)
vline_b = ax_b.axvline(ref['t_days'], color='k', ls='--', lw=1)
ax_b.set_xlim(tmin, tmax); ax_b.set_ylim(0.5, 2.5)
ax_b.set_ylabel('N (MPa)', fontsize=7)
ax_b.text(0.02, 0.85, '(b) effective pressure', transform=ax_b.transAxes, fontsize=7)
ax_b.grid(True, which='both', alpha=0.3)

# ── (c) Sheet thickness ──
for i_d in range(len(t_drainage)):
    ax_c.axvspan(t_drainage[i_d] - 0.5*t_duration[i_d],
                 t_drainage[i_d] + 0.5*t_duration[i_d],
                 color=(0.7, 0.85, 1.0), alpha=0.4)
ax_c.plot(t_ts, hs_ts / A_total, 'b-', lw=1.2, label=r'$h_\mathrm{cav}$')
ax_c2 = ax_c.twinx()
ax_c2.plot(t_ts, S_ts / A_total, 'r-', lw=1.2, label=r'$S$')
vline_c = ax_c.axvline(ref['t_days'], color='k', ls='--', lw=1)
ax_c.set_xlim(tmin, tmax)
ax_c.set_ylabel('Avg h (m)', fontsize=7)
ax_c2.set_ylabel('S/A (m)', fontsize=7, color='r')
ax_c.text(0.02, 0.85, '(c) sheet & channel', transform=ax_c.transAxes, fontsize=7)
ax_c.grid(True, which='both', alpha=0.3)

# ── (d) Ice speed + blister volume ──
for i_d in range(len(t_drainage)):
    ax_d.axvspan(t_drainage[i_d] - 0.5*t_duration[i_d],
                 t_drainage[i_d] + 0.5*t_duration[i_d],
                 color=(0.7, 0.85, 1.0), alpha=0.4)
ax_d.plot(t_ts, U_ts * pd_ty, 'b-', lw=1.2, label=r'$U$ (m/yr)')
vline_d = ax_d.axvline(ref['t_days'], color='k', ls='--', lw=1)
ax_d2 = ax_d.twinx()
ax_d2.plot(t_ts, V_b, 'r-', lw=1.2, label=r'$V_b$')
ax_d.set_xlim(tmin, tmax); ax_d.set_ylim(50, 150)
ax_d.set_ylabel('U (m/yr)', fontsize=7)
ax_d2.set_ylabel(r'$V_b$ (m$^3$)', fontsize=7, color='r')
ax_d.text(0.02, 0.85, '(d) ice speed & blister volume', transform=ax_d.transAxes, fontsize=7)
ax_d.grid(True, which='both', alpha=0.3)

# ── (e) hb+hs at GPS clusters ──
for i_d in range(len(t_drainage)):
    ax_e.axvspan(t_drainage[i_d] - 0.5*t_duration[i_d],
                 t_drainage[i_d] + 0.5*t_duration[i_d],
                 color=(0.7, 0.85, 1.0), alpha=0.4)
for kc in clusters_to_plot:
    members = gps_rows[gps_cluster == kc]
    if len(members) == 0:
        continue
    ax_e.plot(t_ts, np.mean(pts_hb[members, :] + pts_hs[members, :], axis=0),
              '-', color=cluster_colors[kc-1], lw=1.2, label=cluster_labels[kc-1])
vline_e = ax_e.axvline(ref['t_days'], color='k', ls='--', lw=1)
ax_e.set_xlim(tmin, tmax)
ax_e.set_ylabel(r'$h_b + h_s$ (m)', fontsize=7)
ax_e.text(0.02, 0.85, r'(e) $h_b + h_s$ at GPS', transform=ax_e.transAxes, fontsize=7)
ax_e.legend(fontsize=5, ncol=5, loc='upper right')
ax_e.grid(True, which='both', alpha=0.3)

# ── (f) U at GPS clusters ──
for i_d in range(len(t_drainage)):
    ax_f.axvspan(t_drainage[i_d] - 0.5*t_duration[i_d],
                 t_drainage[i_d] + 0.5*t_duration[i_d],
                 color=(0.7, 0.85, 1.0), alpha=0.4)
for kc in clusters_to_plot:
    members = gps_rows[gps_cluster == kc]
    if len(members) == 0:
        continue
    ax_f.plot(t_ts, np.mean(pts_u[members, :], axis=0),
              '-', color=cluster_colors[kc-1], lw=1.2, label=cluster_labels[kc-1])
vline_f = ax_f.axvline(ref['t_days'], color='k', ls='--', lw=1)
ax_f.set_xlim(tmin, tmax); ax_f.set_ylim(0, 300)
ax_f.set_ylabel('U (m/yr)', fontsize=7)
ax_f.set_xlabel('t (d)', fontsize=7)
ax_f.text(0.02, 0.85, '(f) U at GPS', transform=ax_f.transAxes, fontsize=7)
ax_f.legend(fontsize=5, ncol=5, loc='upper right')
ax_f.grid(True, which='both', alpha=0.3)

# ── (g) blank / timestamp ──
ax_g.axis('off')

# hide x labels for non-bottom left panels
for ax_ in [ax_a, ax_b, ax_c, ax_d, ax_e]:
    ax_.tick_params(labelbottom=False)

# ── Right: blister thickness ──
pc_hb = ax_hb.pcolormesh(xx, yy, ref['hb'], cmap=bwr_cmap,
                          shading='gouraud', vmin=-1.0, vmax=1.0)
divider_hb = make_axes_locatable(ax_hb)
cax_hb = divider_hb.append_axes("right", size="3%", pad=0.08)
cb_hb = plt.colorbar(pc_hb, cax=cax_hb)
cb_hb.set_label(r'$h_b$ (m)', fontsize=7)
# pb contours
ct_hb = ax_hb.contour(xx, yy, ref['pb'], colors='k', linewidths=0.5)
# lake markers
Q_peak = np.max(float(ps['Q0']) * V_l / (np.sqrt(2*np.pi) * np.asarray(pp['t_duration'].item())))
lake_flux = lake_input(ref['t_nondim'], pp, ps)
norm_sz = lake_flux / (Q_peak + 1e-30)
sz = 5 + 295 * norm_sz**0.3
lake_sc = ax_hb.scatter(x_l, y_l, s=sz, c=[(0.2, 0.5, 0.9)]*len(x_l),
                         edgecolors='k', zorder=5, alpha=1.0)
# GPS stations
for kc in clusters_to_plot:
    mask = gps_cluster == kc
    if not np.any(mask):
        continue
    ax_hb.scatter(gps_x_km[mask], gps_y_km[mask], s=12, marker='s', zorder=6,
                  color=cluster_colors[kc-1], edgecolors='k', linewidths=0.3,
                  label=cluster_labels[kc-1])
ax_hb.set_aspect('equal')
ax_hb.set_ylabel('y (km)', fontsize=7)
ax_hb.set_title(r'$h_b$ and $p_b$ contour', fontsize=8)
ax_hb.tick_params(labelsize=6)
ax_hb.legend(fontsize=5, loc='upper left', ncol=1)

# ── Right: ice velocity ──
pc_vel = ax_vel.pcolormesh(xx, yy, ref['vUU'], cmap='viridis',
                            shading='gouraud', vmin=0, vmax=200)
divider_vel = make_axes_locatable(ax_vel)
cax_vel = divider_vel.append_axes("right", size="3%", pad=0.08)
cb_vel = plt.colorbar(pc_vel, cax=cax_vel)
cb_vel.set_label('U (m/yr)', fontsize=7)
sk = 5
quiv = ax_vel.quiver(xx[::sk, ::sk], yy[::sk, ::sk],
                     ref['vux'][::sk, ::sk], ref['vuy'][::sk, ::sk],
                     scale=3000, width=0.003, color='k', alpha=0.6)
ax_vel.set_aspect('equal')
ax_vel.set_ylabel('y (km)', fontsize=7)
ax_vel.set_xlabel('x (km)', fontsize=7)
ax_vel.set_title('Ice speed U (m/yr)', fontsize=8)
ax_vel.tick_params(labelsize=6)

# timestamp text
time_text = ax_g.text(0.5, 0.5, f't = {ref["t_days"]:.1f} d',
                      transform=ax_g.transAxes, fontsize=14,
                      fontweight='bold', ha='center', va='center')

vlines = [vline_a, vline_b, vline_c, vline_d, vline_e, vline_f]
# Snapshot collections that must survive the contour-removal loop in update()
keep_collections = set(ax_hb.collections)
print(f'Initial frame: t = {ref["t_days"]:.1f} d')

# ══════════════════════════════════════════════════════════════════════
# Animation update function
# ══════════════════════════════════════════════════════════════════════
def update(frame_i):
    fn = all_filenames[frame_i]
    f = load_frame(fn)
    td = f['t_days']
    # update time markers
    for vl in vlines:
        vl.set_xdata([td, td])
    time_text.set_text(f't = {td:.1f} d')
    # update spatial: blister (2D required for gouraud in mpl >=3.10)
    pc_hb.set_array(f['hb'])
    # update contours — remove old, draw new (preserve GPS + lake + pcolormesh)
    for c in ax_hb.collections[:]:
        if c not in keep_collections:
            c.remove()
    ax_hb.contour(xx, yy, f['pb'], colors='k', linewidths=0.5)
    # update lake marker size
    lake_flux_i = lake_input(f['t_nondim'], pp, ps)
    norm_sz_i = lake_flux_i / (Q_peak + 1e-30)
    lake_sc.set_sizes(5 + 295 * norm_sz_i**0.3)
    # update velocity (2D required for gouraud in mpl >=3.10)
    pc_vel.set_array(f['vUU'])
    quiv.set_UVC(f['vux'][::sk, ::sk], f['vuy'][::sk, ::sk])
    if frame_i % 20 == 0:
        print(f'  frame {frame_i}/{len(all_filenames)}, t={td:.1f} d')
    return []

# ══════════════════════════════════════════════════════════════════════
# Run animation and save
# ══════════════════════════════════════════════════════════════════════
os.makedirs('./results/videos', exist_ok=True)
anim = FuncAnimation(fig, update, frames=len(all_filenames), blit=False)
video_path = f'./results/videos/{casename}_python.mp4'
writer = FFMpegWriter(fps=FPS, metadata=dict(title=casename), bitrate=5000)
anim.save(video_path, writer=writer)

print(f'Video saved: {video_path}')
plt.close(fig)