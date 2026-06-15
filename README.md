# NEVIS: A model for subglacial hydrology with blisters

The repository contains the MATLAB/Python code for simulating the subglacial hydrological system that incorporates the effects of supraglacial lake drainage and subglacial blisters. The model is described in [Zhang et al. (2025)](). It builds upon the subglacial hydrology model of [Hewitt (2013)](https://doi.org/10.1016/j.epsl.2013.04.022), extending it to include the formation and evolution of subglacial blisters caused by supraglacial lake drainage events.

The schematic of the model setup is shown below:
<p align="center">
  <img src="./figures/Fig1_schematic.png" alt="Model Schematic" width="500">
</p>

# Contents
```
.
├── README.md
├── test.m                     # Reference regional workflow: build -> invert -> forward -> drainage
├── cases/
│   ├── README.md              # Case/dataset contract and per-type data requirements
│   ├── templates/             # Starting points (not run as production)
│   │   ├── nevis_idealized.m              # Idealised template
│   │   ├── nevis_regional.m               # Forward regional template
│   │   ├── nevis_regional_inversion_spinup.m  # Hydrology-only spinup template
│   │   └── nevis_inversion.m              # Inversion template
│   └── <case_name>/           # One self-contained subdir per regional run
│       ├── <case_name>_inversion.m        # Inversion config (workflow-managed)
│       ├── <case_name>_inversion_spinup.m # Spinup case   (workflow-managed)
│       ├── <case_name>_forward.m          # Forward case  (workflow-managed)
│       └── <case_name>_drainage.m         # Drainage case (workflow-managed)
├── data/
│   ├── raw/                   # Raw inputs (BedMachine, MEaSUREs, RACMO, GNSS)
│   ├── pipelines/             # Low-level builders for gridded regional products
│   └── datasets/<name>/       # Self-contained data packages, one per case_name
├── src/                       # Core model source
│   └── workflow/              # Regional workflow drivers used by test.m
├── inversion/                 # Adjoint-based inversion entry points
├── srcgen/                    # Python parameter-sweep config generators
├── srcplot/                   # Python/MATLAB visualisation scripts
├── analysis/                  # MATLAB/notebook analysis utilities
├── results/<casename>/        # Per-case simulation outputs
├── logs/, .jobs/              # nevis_run.sh logs and job state
├── src/nevis_run_idealized.m  # Idealised runner
├── src/nevis_run_regional.m   # Regional runner
└── nevis_run.sh               # Batch runner (parallel jobs, spinup→drainage chain)
```

# Quick start: regional reference workflow

For regional simulations, [test.m](test.m) is the reference workflow that
matches the current case and data structure. It drives the whole pipeline for
a single `case_name`, end-to-end:

```
build dataset -> inversion (auto hydrology-only spinup) -> forward (spinup-like) -> drainage
   from raw       C_hat + velocity_inverted                 lakes off, year 2022    lakes on
```

Edit the **USER SETTINGS** block at the top of `test.m` and run it. Each
stage has an on/off toggle, so you can build once and re-run later stages
without rebuilding.

| Setting | Meaning |
|---|---|
| `case_name` | Names the dataset package **and** the per-case config directory `cases/<case_name>/`. |
| `bbox_km` | `[]` for the default region, or `[xmin, xmax, ymin, ymax]` in km, local NEVIS frame (M1 moulin at (0, 0)). Default is `[-80.10, 60.15, -50.07, 30.18]`. |
| `run_dataset_build` | Build `data/datasets/<case_name>/` from raw inputs. |
| `run_inversion` | Run the inversion; auto-runs the hydrology-only spinup if its outputs are missing. |
| `run_forward` | Run the forward regional model (lakes off — a spinup-like year). |
| `run_drainage` | Run the drainage case (lakes on), chained from the forward run's day-365 state. |
| `reinit_case` | `true` regenerates `cases/<case_name>/` from the current templates (overwrites local edits). |
| `sliding_law`, `spinup_days`, `spinup_timestep_file` | First-run-only seeds baked into the generated configs; edit the per-case file afterwards. |

`test.m` creates and uses two matched products for each `case_name`:

| Location | Produced by | Role |
|---|---|---|
| `data/datasets/<case_name>/` | `nevis_build_demo_dataset` and the inversion stage | Self-contained regional data package: geometry, velocity, runoff, stations, moulins, lakes, manifest, source states, and inversion outputs. |
| `cases/<case_name>/` | `nevis_prepare_demo_workflow` | Self-contained configs for the inversion spinup, inversion, forward run, and drainage run. |

Everything below explains what each stage does and where to tune it.

# How the workflow is organised

## Cases: templates vs per-case directories

A case is **resolved by name** through [src/nevis_resolve_case_config.m](src/nevis_resolve_case_config.m):
it looks for `cases/<casename>/<casename>.m`, then `cases/*/<casename>.m`,
then `cases/<casename>.m`.

- **Templates** live in [cases/templates/](cases/templates/) and are starting
  points, not production runs:
  [nevis_regional.m](cases/templates/nevis_regional.m),
  [nevis_regional_inversion_spinup.m](cases/templates/nevis_regional_inversion_spinup.m),
  [nevis_inversion.m](cases/templates/nevis_inversion.m),
  [nevis_idealized.m](cases/templates/nevis_idealized.m).
- **Per-case directories** `cases/<case_name>/` are materialized by the
  reference workflow on the first run for a new `case_name`. Each generated
  file inlines its template verbatim and appends the case-specific bindings —
  it is fully self-contained, with no `run(template)` chain. Subsequent runs
  leave these files untouched, so your edits survive; `reinit_case = true`
  regenerates them.

See [cases/README.md](cases/README.md) for the full case/dataset contract and
the per-type data requirements.

## Stage 1 — Build the dataset package

Each `case_name` gets its own **self-contained** dataset under
`data/datasets/<case_name>/`, built entirely from raw inputs. There is no
shared dataset to point at; a run never reads another case's data.

```matlab
addpath('data/pipelines');
report = build_nevis_dataset();                                  % default region
% report = build_nevis_dataset('data/datasets/my_region', [], [-100, 80, -60, 80]);
```

`build_nevis_dataset(output_dir, raw_root, bbox_km)` builds the four gridded
products end-to-end from raw inputs. With no arguments it reproduces the
original 140 km region under `data/datasets/nevis_2022_140km_rebuilt/`.
`bbox_km` is `[xmin, xmax, ymin, ymax]` in km, local NEVIS frame; default
`[-80.10, 60.15, -50.07, 30.18]`.

| Product file | Builder | Raw source |
|---|---|---|
| `geometry.mat` | `build_geometry.m` | BedMachine v5 NetCDF (auto-discovered under `data/raw/bedmachine/`) |
| `velocity.mat`  | `build_velocity.m`  | MEaSUREs 250 m GeoTIFFs |
| `racmo_runoff_2022.mat`      | `build_racmo_runoff_2022.m`      | RACMO2.3p2 1 km daily NetCDFs |
| `station_timeseries_2022.mat`   | `build_station_timeseries_2022.m`   | GNSS archives |

The reference workflow ([test.m](test.m) via
[src/workflow/nevis_build_demo_dataset.m](src/workflow/nevis_build_demo_dataset.m))
wraps that orchestrator and additionally makes the package complete for the
inversion and forward runs:

- `moulins_2022_coordinates.mat` — moulin coordinates projected into the local
  frame from the raw 2022 shapefile (`cfg.moulins.mode = 'coordinates'`, so
  catchments are rebuilt on whatever grid the bbox produces).
- `environs_lakes_2022B_250416.mat` — the lake catalogue, copied from raw into
  the package so runtime data stays isolated.
- `dataset_manifest.m` — the dataset contract, declaring every file, its
  variable, and expected fields.

See [data/DATA_PIPELINES.md](data/DATA_PIPELINES.md) for raw-input provenance,
key parameters (domain extent, resolution, grid sizes), and field-by-field
requirements.

## Stage 2 — Inversion

```matlab
[inv, vv_hydro, summary] = nevis_run_inversion('./cases/<case_name>/<case_name>_inversion.m');
```

The dispatcher reads `cfg.sliding_law` and calls:

- `weertman` → [inversion/nevis_inv_C.m](inversion/nevis_inv_C.m) (fixed-N adjoint inversion for `tau_b = C * u_b^m`).
- `regularized_coulomb` → [inversion/nevis_inv_C1.m](inversion/nevis_inv_C1.m) (iterative C–N coupled inversion).

For `cfg.initial_hydrology.mode = 'file'`, the runner checks
`cfg.source.state_file` and `cfg.initial_hydrology.file`. If either is missing,
it runs the configured hydrology-only spinup case (`cfg.spinup.casename`,
generated as `cases/<case_name>/<case_name>_inversion_spinup.m`), exports the
selected timestep (`cfg.spinup.timestep_file`), and then starts the inversion.
It never falls back to another case.

Outputs land under `data/datasets/<case_name>/inversions/<cfg.name>/`:
`C_inversion_results.mat` (the inverted field) and `velocity_inverted.mat`
(the reusable hydrology state). The forward case consumes these via
`cfg.inversion.file` and `cfg.initial_hydrology.file`.

The inverted field can optionally be split into pressure-sensitive and
pressure-insensitive sliding components for coupled forward runs by
[inversion/nevis_inv_partition.m](inversion/nevis_inv_partition.m). See
[inversion/INVERSION.md](inversion/INVERSION.md) for the full contract, and
[docs/nevis_inversion.tex](docs/nevis_inversion.tex) for the mathematics.

## Stage 3 — Forward and drainage

```matlab
[tt, vv, info, state] = nevis_run_regional('<case_name>_forward');    % lakes off
[tt, vv, info, state] = nevis_run_regional('<case_name>_drainage');   % lakes on
```

1. **Forward (spinup-like)** — `<case_name>_forward.m` runs the coupled
   hydrology/ice model over 2022 with lake drainage off
   (`cfg.lakes.spinup_volume_scale = 0`), starting from the inversion's
   hydrology state. It saves a day-365 timestep.
2. **Drainage** — `<case_name>_drainage.m` turns lakes on and chains from the
   forward run's day-365 state via
   `cfg.initial_hydrology.mode = 'result_timestep'`, a full-state restore of
   both hydrology and ice velocity.

## Idealised runs (analytical, no dataset)

Idealised runs need no input files — domain, bed, surface, moulins, lakes, and
forcing are all generated from `cfg`:

```matlab
[tt, vv, info, state] = nevis_run_idealized('<casename>');
```

Start from [cases/templates/nevis_idealized.m](cases/templates/nevis_idealized.m).

# Parameter sweeps

Sweeps spanning many configs can be auto-generated by the Python scripts in
`srcgen/`:

```bash
python3 srcgen/reference_scripts.py        --out-dir cases --overwrite   # idealised sweeps
python3 srcgen/regional_study_scripts.py   --out-dir cases --overwrite   # regional sweeps
python3 srcgen/convergence_test_scripts.py --out-dir cases --overwrite   # grid-resolution tests
python3 srcgen/ice_dynamics_scripts.py     --out-dir cases --overwrite \
    --inversion-file inversions/weertman_ns3_default/C_inversion_results.mat \
    --spinup-initial-hydrology-file inversions/weertman_ns3_default/velocity_inverted.mat
```

To run a sweep set, use the bash runner. It discovers `cases/n[12]d*.m`, runs
spinups first, and automatically queues each matching drainage case after its
spinup succeeds:

```bash
CASE_WORKFLOW=idealized MAX_PARALLEL=8 ./nevis_run.sh   # idealised sweep
CASE_WORKFLOW=regional  MAX_PARALLEL=4 ./nevis_run.sh   # regional sweep
SKIP_SPINUP=true CASE_WORKFLOW=regional ./nevis_run.sh  # drainage only (reuse spinups)
```

Configurable environment variables: `MAX_PARALLEL`, `MAX_RETRIES`, `JOB_TIMEOUT`, `POLL_INTERVAL`, `CASE_DIR`, `CASE_WORKFLOW`. Per-job logs go to `logs/<casename>.log` and job state to `.jobs/`.

# Results

Each run writes to `results/<casename>/`:
- `<casename>.mat` — final consolidated state (`pp`, `pd`, `ps`, `gg`, `aa`, `vv`, `vv2`, `oo`, `tt`, `info`, and `region`/`idealized`).
- Per-timestep state files (when `cfg.oo.save_timesteps = 1`).
- Time series at moulins, lakes, and stations (when `cfg.oo.save_pts_all = 1`).

# Postprocess and visualise

| Tool | Purpose |
|---|---|
| [srcplot/make_plot.ipynb](srcplot/make_plot.ipynb) | Static figures (e.g. Figures 2–7) |
| [srcplot/animation.ipynb](srcplot/animation.ipynb) | Time-evolution animations |
| [srcplot/schematic.ipynb](srcplot/schematic.ipynb) | Schematic figure (Figure 1) |
| `srcplot/nevis_*_plotter.py`, `nevis_*_analyzer.py` | Region-, GPS-, stress-, and ice-flow-specific plots |
| `srcplot/nevis_2d_animation*.m`, `nevis_1d_animation.m` | MATLAB animations (controllable time window/step) |
| `analysis/` | Mass balance, blister, power-law, convergence diagnostics |

# Ice-dynamics extension

The repository also includes a depth-averaged shallow-shelf (SSA) ice-velocity solver coupled to the hydrology/blister system. Core files:

```
src/nevis_velocity.m                 # SSA-style ice velocity solve
src/nevis_stresses.m                 # membrane, basal, and driving stress diagnostics
src/nevis_principal_stress.m         # principal stress magnitudes and directions
src/nevis_update_parameters_ice.m    # ice-dynamics parameter setup
src/nevis_label_ice.m                # ice velocity boundary labels
src/nevis_label_ice_test.m           # regional/test boundary labels
src/nevis_timesteps_ice.m            # coupled hydrology/blister/ice time stepping
src/nevis_export_ice_fields.m        # exporters used by animation scripts
```

Enable ice coupling in a regional case by setting `cfg.oo.include_ice = 1` (the reference workflow's forward/drainage cases already do). Coupled parameter sweeps are generated by `srcgen/ice_dynamics_scripts.py` (see [Parameter sweeps](#parameter-sweeps)). The pipeline is the same as the reference workflow: build dataset → run inversion → forward with ice coupling → drainage with ice coupling.

# License
This repository is licensed under the CC BY-NC 4.0 License. See [LICENSE](./LICENSE).

# References
Please cite the following paper when using this code.

Stevens, L. A. (2026). Ice-sheet hydro-fracture not advanced inland by lower-elevation lake drainages in Kalaallit Nunaat (v1.0.0). Zenodo. [doi.org/10.5281/zenodo.19387821](https://doi.org/10.5281/zenodo.19387821)

Zhang, H., Stevens, L. A., Hewitt, I. J., & Stuart, H. (in prep). A unified blister and subglacial hydrology framework for supraglacial lake drainage events. Journal of Geophysical Research: Earth Surface.

Hewitt, I. J. (2013, June). Seasonal changes in ice sheet motion due to melt water lubrication. Earth and Planetary Science Letters, 371–372, 16–25. [doi: 10.1016/j.epsl.2013.04.022](https://doi.org/10.1016/j.epsl.2013.04.022)

Stevens, L. A., Hewitt, I. J., Das, S. B., & Behn, M. D. (2018, September). Relationship Between Greenland Ice Sheet Surface Speed and Modeled Effective Pressure. Journal of Geophysical Research: Earth Surface, 123(9), 2258–2278. [doi: 10.1029/2017JF004581](https://doi.org/10.1029/2017JF004581)
