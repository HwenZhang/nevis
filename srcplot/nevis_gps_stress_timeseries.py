from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio


Z95 = 1.959963984540054
STRESS_GROUPS = [
    (
        'Membrane stresses [kPa]',
        [('tauxx_gps', 'tau_xx'), ('tauyy_gps', 'tau_yy'), ('tauxy_gps', 'tau_xy')],
        1e3,
    ),
    (
        'Principal stresses [kPa]',
        [('sigma1_gps', 'sigma_1'), ('sigma2_gps', 'sigma_2')],
        1e3,
    ),
    (
        'Depth-integrated stresses [kPa]',
        [('Txx_gps', 'T_xx'), ('Tyy_gps', 'T_yy'), ('Txy_gps', 'T_xy')],
        1e3,
    ),
    (
        'Basal and driving [kPa]',
        [('tau_b_gps', 'tau_b'), ('taudx_gps', 'tau_dx'), ('taudy_gps', 'tau_dy')],
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


def mat_array(value, dtype=float):
    while isinstance(value, np.ndarray) and value.dtype == object and value.size == 1:
        value = value.item()
    return np.atleast_1d(np.asarray(value, dtype=dtype)).squeeze()


def plot_gps_stress_timeseries(
    casename,
    repo_root='.',
    tmin_yr=0.45,
    tmax_yr=0.8,
    save_path=None,
):
    repo_root = Path(repo_root).resolve()
    stress_path = repo_root / 'data' / f'gps_stresses_{casename}.mat'
    mat_path = repo_root / 'results' / casename / f'{casename}.mat'

    if not stress_path.exists():
        raise FileNotFoundError(
            f'Missing {stress_path}. Run analysis/nevis_gps_stresses.m first to generate the GPS stress timeseries file.'
        )

    stress = sio.loadmat(stress_path, squeeze_me=True)
    raw = sio.loadmat(mat_path, squeeze_me=True)
    ps = raw['ps']
    pd_r = raw['pd']
    pp = raw['pp']

    ps_t = float(ps['t'])
    pd_td = float(pd_r['td'])
    ps_x = float(ps['x'])

    t_days = np.atleast_1d(stress['t_days']).astype(float)
    gps_x_km = np.atleast_1d(stress['gps_x_km']).astype(float)
    gps_y_km = np.atleast_1d(stress['gps_y_km']).astype(float)
    n_gps = int(np.atleast_1d(stress['n_gps']).astype(int)[0])

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

    lake_x_km = mat_array(pp['x_l']) * ps_x / 1e3
    lake_y_km = mat_array(pp['y_l']) * ps_x / 1e3
    t_drainage = mat_array(pp['t_drainage']) * ps_t / pd_td
    t_duration = mat_array(pp['t_duration']) * ps_t / pd_td

    if not (len(lake_x_km) == len(lake_y_km) == len(t_drainage) == len(t_duration)):
        raise ValueError('Lake location and drainage arrays do not have matching lengths.')

    lake_to_gps = np.argmin(
        (lake_x_km[:, None] - gps_x_km[None, :]) ** 2 + (lake_y_km[:, None] - gps_y_km[None, :]) ** 2,
        axis=1,
    )

    plt.style.use('classic')
    plt.rcParams.update({'font.size': 8, 'figure.dpi': 160})
    fig, axes = plt.subplots(
        n_gps,
        len(STRESS_GROUPS),
        figsize=(14, max(2.4 * n_gps, 4.0)),
        sharex=True,
        constrained_layout=True,
    )

    if n_gps == 1:
        axes = axes[None, :]

    for ig in range(n_gps):
        event_idx = np.where(lake_to_gps == ig)[0]

        for jg, (title, entries, scale) in enumerate(STRESS_GROUPS):
            ax = axes[ig, jg]
            ax.set_title(title, fontsize=9)

            for ie in event_idx:
                half_width = Z95 * t_duration[ie]
                start = t_drainage[ie] - half_width
                end = t_drainage[ie] + half_width
                ax.axvspan(start, end, color='#d7ecff', alpha=0.7, lw=0)

            for field_name, label in entries:
                series = np.atleast_2d(stress[field_name])[ig, plot_mask] / scale
                ax.plot(t_days[plot_mask], series, lw=1.2, color=COLORS[field_name], label=label)

            ax.grid(True, alpha=0.25)
            ax.set_xlim(tmin, tmax)
            if ig == 0:
                ax.legend(loc='upper right', fontsize=7, frameon=False)
            if jg == 0:
                ax.set_ylabel(f'GPS {ig + 1}\n({gps_x_km[ig]:.1f}, {gps_y_km[ig]:.1f}) km')

        axes[ig, -1].text(
            1.02,
            0.5,
            f'{len(event_idx)} nearby events',
            transform=axes[ig, -1].transAxes,
            va='center',
            fontsize=8,
        )

    for ax in axes[-1, :]:
        ax.set_xlabel('t [d]')

    if save_path is None:
        save_path = repo_root / 'figures' / f'{casename}_gps_stress_timeseries.png'
    else:
        save_path = Path(save_path)
        if not save_path.is_absolute():
            save_path = repo_root / save_path
    fig.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f'Saved to {save_path}')

    for ig in range(n_gps):
        event_idx = np.where(lake_to_gps == ig)[0]
        event_summary = [
            f'{t_drainage[idx]:.1f} d, width95={2 * Z95 * t_duration[idx]:.2f} d'
            for idx in event_idx
        ]
        print(f'GPS {ig + 1}: {len(event_idx)} assigned drainage events -> {event_summary}')

    return fig
