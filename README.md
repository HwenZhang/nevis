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
├── cases/                     # Lightweight case configs (cfg.casename = file name)
│   ├── nevis_idealized.m      # Idealised template
│   ├── nevis_regional.m       # Regional template
│   └── nevis_inversion.m      # Inversion template
├── data/
│   ├── raw/                   # Raw inputs (BedMachine, MEaSUREs, RACMO, GNSS)
│   ├── pipelines/             # Scripts that build the shared dataset
│   └── datasets/<name>/       # Curated, reusable regional data packages
├── src/                       # Core model source
├── inversion/                 # Adjoint-based inversion entry points
├── srcgen/                    # Python parameter-sweep config generators
├── srcplot/                   # Python/MATLAB visualisation scripts
├── analysis/                  # MATLAB/notebook analysis utilities
├── results/<casename>/        # Per-case simulation outputs
├── logs/, .jobs/              # nevis_run.sh logs and job state
├── n2d_idealized_template.m   # Idealised runner
├── n2d_regional_template.m    # Regional runner
└── nevis_run.sh               # Batch runner (parallel jobs, spinup→drainage chain)
```

# End-to-end workflow

The same five-step pattern applies to all run types. Only step 1 (data) and step 3 (which runner) differ between **idealised**, **inversion**, and **regional** runs.

```
  raw data  ──▶  build dataset  ──▶  write case config  ──▶  run model  ──▶  postprocess
   (1)            data/pipelines/      cases/<name>.m         step 3            srcplot/
                  data/datasets/<name>                                          analysis/
```

| Stage | Idealised | Inversion | Regional |
|---|---|---|---|
| 1. Data | None — geometry is analytical | Shared dataset (geometry, observed velocity, initial hydrology) | Shared dataset + inversion output |
| 2. Config | `cases/<name>.m` from `nevis_idealized.m` | `cases/<name>.m` from `nevis_inversion.m` | `cases/<name>.m` from `nevis_regional.m` |
| 3. Runner | `n2d_idealized_template('<name>')` | `nevis_run_inversion('./cases/<name>.m')` | `n2d_regional_template('<name>')` |
| 4. Output | `results/<name>/` | `data/datasets/<dataset>/inversions/<name>/` | `results/<name>/` |
| 5. Postprocess | `srcplot/`, `analysis/` | Used by regional runs as initial conditions | `srcplot/`, `analysis/` |

## Step 1 — Generate the data package (regional/inversion only)

Idealised runs need no input files; skip to step 2.

Regional and inversion runs read from a shared, reusable **dataset package** under `data/datasets/<dataset_name>/`. The shipped package is `nevis_2022_140km`. To rebuild it from raw data, run:

```matlab
addpath('data/pipelines');
report = build_nevis_dataset();                                  % default 140 km bbox
% report = build_nevis_dataset([], [], [-100, 80, -60, 80]);     % any other bbox (km)
```

`build_nevis_dataset(output_dir, raw_root, bbox_km)` rebuilds geometry, velocity, runoff, and stations end-to-end from raw inputs. With no arguments it reproduces the original 140 km region under `data/datasets/nevis_2022_140km_rebuilt/`. The `bbox_km` argument is `[xmin, xmax, ymin, ymax]` in km in the local NEVIS frame (M1 moulin at (0, 0)); default `[-79.95, 60.00, -49.92, 30.03]`.

| Product file | Builder | Raw source |
|---|---|---|
| `geometry.mat` | `build_geometry.m` | BedMachine v5 NetCDF (auto-discovered under `data/raw/bedmachine/`) |
| `velocity.mat`  | `build_velocity.m`  | MEaSUREs 250 m GeoTIFFs |
| `racmo_runoff_2022.mat`      | `build_racmo_runoff_2022.m`      | RACMO2.3p2 1 km daily NetCDFs |
| `station_timeseries_2022.mat`   | `build_station_timeseries_2022.m`   | GNSS archives |

The orchestrator also writes a `dataset_manifest.m` next to the data files (only if one does not already exist) so the rebuild can be bound to a regional case via `cfg.dataset.root` and consumed by the inversion. Auxiliary inputs the inversion needs but the pipeline does not generate (moulins, lakes, source model state, initial hydrology) must be supplied separately and declared in the manifest.

See [data/DATA_PIPELINES.md](data/DATA_PIPELINES.md) for raw-input provenance, key parameters (domain extent, resolution, grid sizes), and field-by-field requirements.

## Step 2 — Write a case config

Every run is driven by a single MATLAB config file at `cases/<casename>.m` in which `cfg.casename` must match the file name. Start from the template that matches the run type:

- **Idealised**: [cases/nevis_idealized.m](cases/nevis_idealized.m) — analytical geometry, moulin/lake placement, and forcing. See [cases/nevis_idealized_DATA_REQUIREMENTS.md](cases/nevis_idealized_DATA_REQUIREMENTS.md).
- **Inversion**: [cases/nevis_inversion.m](cases/nevis_inversion.m) — sliding law (`weertman` or `regularized_coulomb`), observed velocity, prior, and regularisation.
- **Regional**: [cases/nevis_regional.m](cases/nevis_regional.m) — binds to a `cfg.dataset` package and an inversion package; sets hydrology, ice-coupling, and forcing options. See [cases/nevis_regional_DATA_REQUIREMENTS.md](cases/nevis_regional_DATA_REQUIREMENTS.md).

Create a new case:

```bash
cp cases/nevis_regional.m cases/<new_casename>.m
# edit <new_casename>.m: set cfg.casename = '<new_casename>'
# adjust cfg.oo / cfg.pd / cfg.pp / cfg.run / cfg.lakes / cfg.runoff as needed
```

Parameter sweeps spanning many configs can be auto-generated by the Python scripts in `srcgen/`:

```bash
python3 srcgen/reference_scripts.py        --out-dir cases --overwrite   # idealised sweeps
python3 srcgen/regional_study_scripts.py   --out-dir cases --overwrite   # regional sweeps
python3 srcgen/convergence_test_scripts.py --out-dir cases --overwrite   # grid-resolution tests
python3 srcgen/ice_dynamics_scripts.py     --out-dir cases --overwrite \
    --inversion-file inversions/weertman_ns3_default/C_inversion_results.mat \
    --spinup-initial-hydrology-file inversions/weertman_ns3_default/velocity_inverted.mat
```

## Step 3 — Run the model

### Idealised

```matlab
[tt, vv, info, state] = n2d_idealized_template('<casename>');
```

Domain, bed, surface, moulins, lakes, and forcing are all generated analytically from `cfg`. No external data files are read.

### Inversion

```matlab
[inv, vv_hydro, summary] = nevis_run_inversion('./cases/<casename>.m');
```

The dispatcher reads `cfg.sliding_law` and calls:
- `weertman` → [inversion/nevis_inv_C.m](inversion/nevis_inv_C.m) (fixed-N adjoint inversion for `tau_b = C * u_b^m`).
- `regularized_coulomb` → [inversion/nevis_inv_C1.m](inversion/nevis_inv_C1.m) (iterative C–N coupled inversion).

Outputs (`C_inversion_results.mat` and `velocity_inverted.mat`) are written under `data/datasets/<dataset>/inversions/<cfg.name>/`. Regional cases later consume these files via `cfg.inversion.file` and `cfg.initial_hydrology.file`.

The inverted field can optionally be split into pressure-sensitive and pressure-insensitive sliding components for coupled forward runs by [inversion/nevis_inv_partition.m](inversion/nevis_inv_partition.m). For mathematical details see [docs/nevis_inversion.tex](docs/nevis_inversion.tex).

### Regional

```matlab
[tt, vv, info, state] = n2d_regional_template('<casename>');
```

Every regional case consists of two stages tied by file name convention:

1. **Spinup** — `<casename>_spinup.m` brings the system to a quasi-steady state with the lake-drainage volume turned off (`cfg.lakes.spinup_volume_scale = 0`).
2. **Drainage** — `<casename>_..._drainage*.m` reuses the spinup result as the initial condition and applies prescribed lake-drainage events.

### Batch runs

For parameter-sweep sets, use the bash runner. It discovers `cases/n[12]d*.m`, runs spinups first, and automatically queues each matching drainage case after its spinup succeeds:

```bash
# idealised sweep
CASE_WORKFLOW=idealized MAX_PARALLEL=8 ./nevis_run.sh

# regional sweep
CASE_WORKFLOW=regional  MAX_PARALLEL=4 ./nevis_run.sh

# skip spinup and run drainage configs only (uses pre-existing spinup results)
SKIP_SPINUP=true CASE_WORKFLOW=regional ./nevis_run.sh
```

Configurable environment variables: `MAX_PARALLEL`, `MAX_RETRIES`, `JOB_TIMEOUT`, `POLL_INTERVAL`, `CASE_DIR`, `CASE_WORKFLOW`. Per-job logs are written to `logs/<casename>.log` and job state to `.jobs/`.

## Step 4 — Results

Each run writes to `results/<casename>/`:
- `<casename>.mat` — final consolidated state (`pp`, `pd`, `ps`, `gg`, `aa`, `vv`, `vv2`, `oo`, `tt`, `info`, and `region`/`idealized`).
- Per-timestep state files (when `cfg.oo.save_timesteps = 1`).
- Time series at moulins, lakes, and stations (when `cfg.oo.save_pts_all = 1`).

## Step 5 — Postprocess and visualise

| Tool | Purpose |
|---|---|
| [srcplot/make_plot.ipynb](srcplot/make_plot.ipynb) | Static figures (e.g. Figures 2–7) |
| [srcplot/animation.ipynb](srcplot/animation.ipynb) | Time-evolution animations |
| [srcplot/schematic.ipynb](srcplot/schematic.ipynb) | Schematic figure (Figure 1) |
| `srcplot/nevis_*_plotter.py`, `nevis_*_analyzer.py` | Region-, GPS-, stress-, and ice-flow-specific plots |
| `srcplot/nevis_2d_animation*.m`, `nevis_1d_animation.m` | MATLAB animations |
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

Enable ice coupling in a regional case by setting `cfg.oo.include_ice = 1`. Coupled parameter sweeps are generated by `srcgen/ice_dynamics_scripts.py` (see Step 2). The pipeline is: build dataset → run inversion → spinup with ice coupling → drainage with ice coupling.

# License
This repository is licensed under the CC BY-NC 4.0 License. See [LICENSE](./LICENSE).

# References
Please cite the following paper when using this code.

Stevens, L. A. (2026). Ice-sheet hydro-fracture not advanced inland by lower-elevation lake drainages in Kalaallit Nunaat (v1.0.0). Zenodo. [doi.org/10.5281/zenodo.19387821](https://doi.org/10.5281/zenodo.19387821)

Zhang, H., Stevens, L. A., Hewitt, I. J., & Stuart, H. (in prep). A unified blister and subglacial hydrology framework for supraglacial lake drainage events. Journal of Geophysical Research: Earth Surface.

Hewitt, I. J. (2013, June). Seasonal changes in ice sheet motion due to melt water lubrication. Earth and Planetary Science Letters, 371–372, 16–25. [doi: 10.1016/j.epsl.2013.04.022](https://doi.org/10.1016/j.epsl.2013.04.022)

Stevens, L. A., Hewitt, I. J., Das, S. B., & Behn, M. D. (2018, September). Relationship Between Greenland Ice Sheet Surface Speed and Modeled Effective Pressure. Journal of Geophysical Research: Earth Surface, 123(9), 2258–2278. [doi: 10.1029/2017JF004581](https://doi.org/10.1029/2017JF004581)
