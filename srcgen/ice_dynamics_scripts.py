#!/usr/bin/env python3
"""Generate regional ice-dynamics NEVIS case configs.

This replaces the old MATLAB script generator. It writes configs into ./cases;
run them with n2d_regional_template(casename), or set CASE_WORKFLOW=regional
when using nevis_run.sh.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path


CASE_LABEL = "test"

PARAMETER_COMBINATIONS = [
    # partition_ratio, eps_reg, kappa, mu, V, k0, surface_runoff, racmo_runoff, distributed_input
    (0.5, 0.01, 5e-11, 5e0, 1e1, 1e-1, 1, 1, 0),
    (0.4, 0.01, 5e-11, 5e0, 1e1, 1e-1, 1, 1, 0),
    (0.3, 0.01, 5e-11, 5e0, 1e1, 1e-1, 1, 1, 0),
    (0.2, 0.01, 5e-11, 5e0, 1e1, 1e-1, 1, 1, 0),
    (0.1, 0.01, 5e-11, 5e0, 1e1, 1e-1, 1, 1, 0),
    (0.3, 0.01, 5e-11, 1e1, 1e1, 1e-1, 1, 1, 0),
    (0.3, 0.01, 5e-11, 1e0, 1e1, 1e-1, 1, 1, 0),
    (0.3, 0.01, 1e-11, 1e1, 1e1, 1e-1, 1, 1, 0),
]


def matlab_number(value: float | int) -> str:
    return f"{value:.12g}"


def sci_name(value: float | int) -> str:
    text = f"{float(value):.0e}"
    mantissa, exponent = text.split("e")
    exponent = exponent.replace("+", "").replace("-", "_")
    return f"{mantissa}e{exponent}"


def runoff_name(racmo_runoff: int) -> str:
    return "racmo" if racmo_runoff else "analytic"


def spinup_case_name(
    partition_ratio: float,
    eps_reg: float,
    kappa: float,
    mu: float,
    k0: float,
    racmo_runoff: int,
) -> str:
    label = f"{CASE_LABEL}_" if CASE_LABEL else ""
    return (
        f"n2d_regional_{label}{runoff_name(racmo_runoff)}"
        f"_eps{sci_name(eps_reg)}_kappa{sci_name(kappa)}_mu{sci_name(mu)}"
        f"_partition{sci_name(partition_ratio)}_k0{sci_name(k0)}_spinup"
    )


def drainage_case_name(spinup_case: str, volume: float) -> str:
    volume_log10 = int(round(math.log10(volume))) if volume > 0 else 0
    return spinup_case.replace("_spinup", f"_V1e{volume_log10}_drainage_highelev")


def inversion_config_block(mode: str, file: str) -> str:
    if mode == "uniform":
        return """cfg.inversion.mode = 'uniform';
cfg.inversion.file = '';
"""
    if mode != "file":
        raise ValueError(f"Unsupported inversion mode: {mode}")
    return f"""cfg.inversion.mode = 'file';
cfg.inversion.file = '{file}';
"""


def initial_hydrology_config_block(
    mode: str,
    file: str,
    variable: str,
    k_factor: float,
    result_case: str | None = None,
    timestep_file: str = "0365.mat",
) -> str:
    if mode == "file":
        return f"""cfg.initial_hydrology.mode = 'file';
cfg.initial_hydrology.file = '{file}';
cfg.initial_hydrology.variable = '{variable}';
"""
    if mode == "k_factor":
        return f"""cfg.initial_hydrology.mode = 'k_factor';
cfg.initial_hydrology.k_factor = {matlab_number(k_factor)};
"""
    if mode == "result_timestep":
        if result_case is None:
            raise ValueError("result_timestep initial hydrology requires result_case")
        return f"""cfg.initial_hydrology.mode = 'result_timestep';
cfg.initial_hydrology.result_case = '{result_case}';
cfg.initial_hydrology.timestep_file = '{timestep_file}';
"""
    raise ValueError(f"Unsupported initial hydrology mode: {mode}")


def common_overrides(
    casename: str,
    partition_ratio: float,
    eps_reg: float,
    kappa: float,
    mu: float,
    k0: float,
    surface_runoff: int,
    racmo_runoff: int,
    distributed_input: int,
    inversion_mode: str,
    inversion_file: str,
) -> str:
    return f"""run(fullfile(fileparts(mfilename('fullpath')), 'nevis_regional.m'));

cfg.casename = '{casename}';
cfg.case_root = './cases';

{inversion_config_block(inversion_mode, inversion_file)}

cfg.oo.include_ice = 1;
cfg.oo.mean_perms = 1;
cfg.oo.modified_mean_perms = 0;
cfg.oo.display_residual = 0;
cfg.oo.visualize_vel = 0;
cfg.oo.N_coupling = 1;
cfg.oo.U_coupling = 1;
cfg.oo.boundary_method = 'stress_l_vel_tbl';
cfg.oo.mask_boundary_method = 'stress_free';
cfg.oo.plot_residual = 0;
cfg.oo.max_iter_new = 50;
cfg.oo.partition_ratio = {matlab_number(partition_ratio)};
cfg.oo.iterative_solver = 0;

cfg.pd.kappa_b = {matlab_number(kappa)};
cfg.pd.mu = {matlab_number(mu)};
cfg.pd.k_s = {matlab_number(k0)};
cfg.pd.melt = @(pd) (pd.G+(pd.u_b*pd.tau_b))/pd.rho_w/pd.L;
cfg.pd.meltinterior = @(pd) 0*((pd.G+((100/pd.ty)*pd.tau_b))/pd.rho_w/pd.L)*1e3;
cfg.pd.B_reg = @(pd) pd.Ye*(1000)^3/(12*(1-0.33)^2);

cfg.pp.eps_reg = {matlab_number(eps_reg)};
cfg.pp.Ub_reg = 1e-16;
cfg.pp.N_slide_reg = 1e-16;
cfg.pp.taud_reg = 1e-16;
cfg.pp.C2 = 0;

cfg.runoff.surface_runoff = logical({int(surface_runoff)});
cfg.runoff.RACMO_runoff = logical({int(racmo_runoff)});
cfg.runoff.distributed_input = logical({int(distributed_input)});
cfg.runoff.meltE_mm_per_day = 0;
cfg.runoff.meltE_ramp_days = 20;
cfg.runoff.smooth_days = 1;
cfg.runoff.subtract_lake_volume = true;
"""


def spinup_config(
    casename: str,
    partition_ratio: float,
    eps_reg: float,
    kappa: float,
    mu: float,
    k0: float,
    surface_runoff: int,
    racmo_runoff: int,
    distributed_input: int,
    inversion_mode: str,
    inversion_file: str,
    initial_hydrology_mode: str,
    initial_hydrology_file: str,
    initial_hydrology_variable: str,
    initial_hydrology_k_factor: float,
) -> str:
    volume_scale = 1 if racmo_runoff else 0
    return f"""% Auto-generated regional NEVIS ice-dynamics spinup case config.
% Generated by ./srcgen/ice_dynamics_scripts.py.
% Run with: n2d_regional_template('{casename}')

{common_overrides(casename, partition_ratio, eps_reg, kappa, mu, k0, surface_runoff, racmo_runoff, distributed_input, inversion_mode, inversion_file)}
{initial_hydrology_config_block(initial_hydrology_mode, initial_hydrology_file, initial_hydrology_variable, initial_hydrology_k_factor)}

cfg.lakes.mode = 'environs_lakes_catalogue';
cfg.lakes.file = 'mechanistic_drainage_catalogue_2022-main/catalogues/environs_lakes_2022B_250505_archive.mat';
cfg.lakes.variable = 'environs_lakes';
cfg.lakes.spinup_volume_scale = {volume_scale};

cfg.stations.mode = 'file';
cfg.stations.file = 'station_timeseries_2022.mat';
cfg.stations.variable = 'station_data';
cfg.stations.x_field = 'x_m';
cfg.stations.y_field = 'y_m';

cfg.run.t_span_days = [(1:1:180) (180.1:0.1:239.9) (240:1:365)];
cfg.run.save_timesteps = true;
cfg.run.save_pts_all = true;
"""


def drainage_config(
    casename: str,
    spinup_case: str,
    partition_ratio: float,
    eps_reg: float,
    kappa: float,
    mu: float,
    volume: float,
    k0: float,
    surface_runoff: int,
    racmo_runoff: int,
    distributed_input: int,
    inversion_mode: str,
    inversion_file: str,
    initial_hydrology_mode: str,
    initial_hydrology_file: str,
    initial_hydrology_variable: str,
    initial_hydrology_k_factor: float,
    initial_hydrology_timestep_file: str,
) -> str:
    if racmo_runoff:
        lake_block = """cfg.lakes.mode = 'environs_lakes_catalogue';
cfg.lakes.file = 'mechanistic_drainage_catalogue_2022-main/catalogues/environs_lakes_2022B_250505_archive.mat';
cfg.lakes.variable = 'environs_lakes';
cfg.lakes.spinup_volume_scale = 1;
"""
    else:
        lake_block = f"""cfg.lakes.mode = 'manual';
cfg.lakes.x_m = 40e3;
cfg.lakes.y_m = -15e3;
cfg.lakes.volume_m3 = {matlab_number(volume)};
cfg.lakes.drainage_after_start_days = 10;
cfg.lakes.duration_days = 0.25;
"""

    return f"""% Auto-generated regional NEVIS ice-dynamics drainage case config.
% Generated by ./srcgen/ice_dynamics_scripts.py.
% Run with: n2d_regional_template('{casename}')

{common_overrides(casename, partition_ratio, eps_reg, kappa, mu, k0, surface_runoff, racmo_runoff, distributed_input, inversion_mode, inversion_file)}
{initial_hydrology_config_block(initial_hydrology_mode, initial_hydrology_file, initial_hydrology_variable, initial_hydrology_k_factor, spinup_case, initial_hydrology_timestep_file)}

{lake_block}
cfg.stations.mode = 'file';
cfg.stations.file = 'station_timeseries_2022.mat';
cfg.stations.variable = 'station_data';
cfg.stations.x_field = 'x_m';
cfg.stations.y_field = 'y_m';

cfg.run.t_span_after_start_days = [(1:1:180) (180.1:0.1:239.9) (240:1:365)];
cfg.run.save_timesteps = true;
cfg.run.save_pts_all = true;
"""


def write_file(path: Path, text: str, overwrite: bool) -> bool:
    if path.exists() and not overwrite:
        return False
    path.write_text(text, encoding="utf-8")
    return True


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", default="./cases", help="Directory for generated case configs.")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing configs.")
    parser.add_argument("--spinup-only", action="store_true", help="Generate only spinup configs.")
    parser.add_argument("--drainage-only", action="store_true", help="Generate only drainage configs.")
    parser.add_argument(
        "--inversion-mode",
        choices=("file", "uniform"),
        default="file",
        help="Sliding coefficient source written to cfg.inversion.",
    )
    parser.add_argument(
        "--inversion-file",
        default="C_inversion_results.mat",
        help="Dataset-relative inversion file when --inversion-mode=file.",
    )
    parser.add_argument(
        "--spinup-initial-hydrology-mode",
        choices=("file", "k_factor"),
        default="file",
        help="Initial hydrology source for generated spinup configs.",
    )
    parser.add_argument(
        "--spinup-initial-hydrology-file",
        default="velocity_inverted.mat",
        help="Dataset-relative initial hydrology file for spinup file mode.",
    )
    parser.add_argument(
        "--spinup-initial-hydrology-variable",
        default="vv_hydro",
        help="Variable name in the spinup initial hydrology file.",
    )
    parser.add_argument(
        "--spinup-initial-hydrology-k-factor",
        type=float,
        default=0.9,
        help="k_factor value for spinup k_factor initial hydrology mode.",
    )
    parser.add_argument(
        "--drainage-initial-hydrology-mode",
        choices=("result_timestep", "file", "k_factor"),
        default="result_timestep",
        help="Initial hydrology source for generated drainage configs.",
    )
    parser.add_argument(
        "--drainage-initial-hydrology-file",
        default="velocity_inverted.mat",
        help="Dataset-relative initial hydrology file for drainage file mode.",
    )
    parser.add_argument(
        "--drainage-initial-hydrology-variable",
        default="vv_hydro",
        help="Variable name in the drainage initial hydrology file.",
    )
    parser.add_argument(
        "--drainage-initial-hydrology-k-factor",
        type=float,
        default=0.9,
        help="k_factor value for drainage k_factor initial hydrology mode.",
    )
    parser.add_argument(
        "--drainage-initial-hydrology-timestep-file",
        default="0365.mat",
        help="Result timestep file used for drainage result_timestep mode.",
    )
    args = parser.parse_args()

    if args.spinup_only and args.drainage_only:
        parser.error("--spinup-only and --drainage-only are mutually exclusive")

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    write_spinup = not args.drainage_only
    write_drainage = not args.spinup_only
    written_spinups = []
    written_drainage = []
    skipped = []

    unique_spinups: dict[str, tuple[float, float, float, float, float, int, int, int]] = {}
    for params in PARAMETER_COMBINATIONS:
        partition_ratio, eps_reg, kappa, mu, _volume, k0, surface_runoff, racmo_runoff, distributed_input = params
        name = spinup_case_name(partition_ratio, eps_reg, kappa, mu, k0, racmo_runoff)
        unique_spinups[name] = (partition_ratio, eps_reg, kappa, mu, k0, surface_runoff, racmo_runoff, distributed_input)

    if write_spinup:
        for name, params in sorted(unique_spinups.items()):
            path = out_dir / f"{name}.m"
            text = spinup_config(
                name,
                *params,
                args.inversion_mode,
                args.inversion_file,
                args.spinup_initial_hydrology_mode,
                args.spinup_initial_hydrology_file,
                args.spinup_initial_hydrology_variable,
                args.spinup_initial_hydrology_k_factor,
            )
            if write_file(path, text, args.overwrite):
                written_spinups.append(path.name)
            else:
                skipped.append(path.name)

    if write_drainage:
        for params in PARAMETER_COMBINATIONS:
            partition_ratio, eps_reg, kappa, mu, volume, k0, surface_runoff, racmo_runoff, distributed_input = params
            spinup = spinup_case_name(partition_ratio, eps_reg, kappa, mu, k0, racmo_runoff)
            name = drainage_case_name(spinup, volume)
            path = out_dir / f"{name}.m"
            text = drainage_config(
                name,
                spinup,
                partition_ratio,
                eps_reg,
                kappa,
                mu,
                volume,
                k0,
                surface_runoff,
                racmo_runoff,
                distributed_input,
                args.inversion_mode,
                args.inversion_file,
                args.drainage_initial_hydrology_mode,
                args.drainage_initial_hydrology_file,
                args.drainage_initial_hydrology_variable,
                args.drainage_initial_hydrology_k_factor,
                args.drainage_initial_hydrology_timestep_file,
            )
            if write_file(path, text, args.overwrite):
                written_drainage.append(path.name)
            else:
                skipped.append(path.name)

    print(f"Generated {len(written_spinups)} spinup configs and {len(written_drainage)} drainage configs in {out_dir}.")
    if skipped:
        print(f"Skipped {len(skipped)} existing configs; use --overwrite to replace them.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
