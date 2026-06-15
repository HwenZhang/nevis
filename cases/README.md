# NEVIS Cases

This directory holds the per-run case configs that drive the NEVIS
solvers and the templates they are derived from. It contains no heavy
data; runtime MAT/NetCDF/TIFF files live under [`./data/datasets/`](../data/datasets/).

## Layout

```
cases/
  README.md                              this file
  templates/                             starting points (do not run as production)
    nevis_inversion.m                    inversion config template (regional)
    nevis_regional.m                     forward regional case template
    nevis_regional_inversion_spinup.m    hydrology-only spinup case template
    nevis_idealized.m                    idealized analytical-domain template
  <case_name>/                           one subdir per regional run, e.g. nevis_2022_140km_rebuilt/
    <case_name>_inversion.m              fully inlined inversion config (workflow-managed)
    <case_name>_inversion_spinup.m       fully inlined spinup case      (workflow-managed)
    <case_name>_forward.m                fully inlined forward case     (workflow-managed)
    <case_name>_drainage.m               fully inlined drainage case    (workflow-managed)
    <case_name>_<variant>.m              additional hand-written variants
```

The workflow-managed files are materialized by `nevis_prepare_demo_workflow`
and respected on subsequent runs (only `reinit_case = true` regenerates them).
Any other `*.m` next to them is treated as a hand-written sibling — the
resolver picks it up via `cases/*/<casename>.m`, and the workflow never
touches it. See [Adding a hand-written case variant](#adding-a-hand-written-case-variant) below.

A case is resolved by name via [src/nevis_resolve_case_config.m](../src/nevis_resolve_case_config.m):
the resolver checks `cases/<casename>/<casename>.m`, then `cases/*/<casename>.m`,
then `cases/<casename>.m`. So `nevis_run_regional('<case_name>_forward')`
finds the per-case file in `cases/<case_name>/`, and
`nevis_run_regional('nevis_regional')` finds the template in
`cases/templates/`.

## Case types

There are two top-level case types: **regional** (real geometry + datasets) and
**idealized** (analytical domain). Each type has a template under
[`templates/`](templates/) that documents the contract. Run a template
directly to validate it against the bundled reference dataset;
production runs should use a per-case subdirectory.

### Regional

```matlab
[tt, vv, info, state] = nevis_run_regional('<casename>');
```

A regional case must declare:

- `cfg.casename` — must equal the requested name; the resolver guards this.
- `cfg.dataset` — binding to a heavy data package under `./data/datasets/`:

  ```matlab
  cfg.dataset.name = 'nevis_2022_140km';
  cfg.dataset.root = './data/datasets/nevis_2022_140km';
  cfg.dataset.manifest = 'dataset_manifest.m';
  ```

- `cfg.oo` — model and solver options.
- `cfg.pd` — dimensional physical parameters.
- `cfg.pp` — nondimensional overrides applied after `nevis_nondimension`.
- `cfg.inversion.file` — inverted basal field, relative to the dataset root.
- `cfg.initial_hydrology.file` — initial effective pressure state, relative to the dataset root.

The regional validator (`nevis_validate_region_config`) refuses to run if
the manifest's `dataset.name` does not match `cfg.dataset.name` or if any
declared file is missing from `cfg.dataset.root`.

### Idealized

```matlab
[tt, vv, info, state] = nevis_run_idealized('<casename>');
```

Idealized cases are fully self-contained — there are no heavy data
requirements. The config declares analytical geometry, initial
hydrology, moulin/lake placement, forcing, and run control. The runner
generates every gridded field from `cfg`.

Supported analytical modes:

- `cfg.geometry.mode = 'hewitt_2013_sqrt'` — built-in flat-bed Hewitt-style geometry.
- `cfg.geometry.mode = 'function'` — custom `bed_function` and `surface_function`.
- `cfg.initial_hydrology.mode = 'k_factor_cavity_equilibrium'` — built-in pressure/sheet initialization.
- `cfg.initial_hydrology.mode = 'function'` — custom function returning `[vv, aa]`.
- `cfg.moulins.mode = 'coordinates_fractional'` — x/y fractions of domain length/width.
- `cfg.moulins.mode = 'coordinates_m'` — explicit x/y positions in meters.
- `cfg.moulins.mode = 'function'` — custom position function.

Each moulin can carry its own input: set `cfg.moulins.input_m3_per_s` to
a scalar or per-moulin vector, or set `cfg.moulins.input_function` to a
function of model time.

## Creating a new case

The recommended path is the demo workflow ([test.m](../test.m)). It
materializes a per-case subdirectory containing fully self-contained
configs — no `run(template)` chain at runtime:

```matlab
% in test.m:
case_name = '<new_case_name>';
bbox_km = [];          % or a custom [xmin, xmax, ymin, ymax] in km
test
```

This generates `cases/<new_case_name>/` with four files (forward,
drainage, spinup, inversion), each carrying its template content inlined
plus the case-specific bindings appended. Edit those files to tune solver,
physics, and run options. Set `reinit_case = true` in `test.m` to
regenerate from the current templates (overwrites local edits in that
dir).

To create a parameter-sweep variant of an existing regional case
manually:

1. Copy `cases/<case_name>/<case_name>_forward.m` to
   `cases/<new_case_name>/<new_case_name>_forward.m`.
2. Update `cfg.casename` to `<new_case_name>_forward`.
3. Keep `cfg.dataset` pointing at the same dataset package if
   topography, forcing, inversion, initial hydrology, moulins, lakes,
   and stations are unchanged.

For a new idealized case, copy `cases/templates/nevis_idealized.m` to
`cases/<new_case_name>.m`, set `cfg.casename = '<new_case_name>'`, edit
the analytical config, and run
`nevis_run_idealized('<new_case_name>')`.

### Adding a hand-written case variant

Within an existing per-case directory you can drop additional
`*.m` files that the workflow does not own. The resolver finds them by
casename via `cases/*/<casename>.m`, so they run with
`nevis_run_regional('<variant_name>')` just like the workflow-managed
files. Use this for variants that share the same dataset / inversion but
need a different forcing or chaining strategy.

The reference per-case dir ships with one such variant:

| File | Casename | Purpose |
|---|---|---|
| `nevis_2022_140km_rebuilt_drainage.m` | `nevis_2022_140km_rebuilt_drainage` | Year-2022 forward run with lake drainage events ON, chained from the spinup-like `*_forward` run's day-365 timestep. |

Pattern to follow:

1. Copy `<case_name>_forward.m` to `<case_name>_<variant>.m` in the same dir.
2. Replace the `% ====== workflow overrides ======` block at the bottom with your variant's overrides:
   - `cfg.casename = '<case_name>_<variant>';`
   - Keep `cfg.dataset.*` and `cfg.inversion.*` pointing at the same
     dataset/inverted-basal field as the workflow-managed forward case.
   - Either keep `cfg.initial_hydrology.mode = 'file'` (start from the
     inversion's hydrology) or switch to `'result_timestep'` with
     `result_case`/`timestep_file` to chain from a saved timestep of a
     previous run. `result_timestep` is a full-state restore: both
     hydrology (`phi`, `hs`) and ice velocity (`u`, `v`) are inherited
     from the saved `vv`, so the chained run continues seamlessly from a
     coupled (`include_ice = 1`) predecessor.
   - Override anything else you want different (e.g.
     `cfg.lakes.spinup_volume_scale`, `cfg.run.t_span_days`).
3. Run with `nevis_run_regional('<case_name>_<variant>')`.

This file is not managed by the workflow; `reinit_case = true` does not
overwrite it.

## Inversion

```matlab
[inv, vv_hydro, summary] = nevis_run_inversion( ...
    './cases/<case_name>/<case_name>_inversion.m');
```

Set `cfg.sliding_law` in the inversion config to `'weertman'`
(→ [inversion/nevis_inv_C.m](../inversion/nevis_inv_C.m)) or
`'regularized_coulomb'`
(→ [inversion/nevis_inv_C1.m](../inversion/nevis_inv_C1.m)).

See [inversion/INVERSION.md](../inversion/INVERSION.md) for the
end-to-end contract (required inputs, output products, common
pitfalls).

Ice-dynamics parameter sweeps generated by
[`srcgen/ice_dynamics_scripts.py`](../srcgen/ice_dynamics_scripts.py)
can target an explicit inversion package and initial hydrology source:

```bash
python3 ./srcgen/ice_dynamics_scripts.py --overwrite \
  --inversion-file inversions/weertman_ns3_default/C_inversion_results.mat \
  --spinup-initial-hydrology-file inversions/weertman_ns3_default/velocity_inverted.mat
```

## Data requirements: regional case dataset

Regional cases read all heavy data through `cfg.dataset.root/dataset_manifest.m`.
The manifest declares each file, the variable inside it, and any
expected fields. The runtime ([nevis_run_regional.m](../src/nevis_run_regional.m)
→ `nevis_setup_region_case`) refuses to run if anything declared is
missing, mismatched, or wrong size. Each case dataset under
`./data/datasets/<case_name>/` should contain its own complete set of
files; runtime configs must not read auxiliary data from another case
dataset.

| Config section | File | Required variable / fields | Purpose |
|---|---|---|---|
| `geometry` | `geometry.mat` | variable `geometry` with `X_m`, `Y_m`, `B_m`, `S_m` | Model grid, bed, surface, ice-mask, hydrology/blister/boundary labels |
| `velocity` | `velocity.mat` | variable `velocity` with `u_obs_dim`, `v_obs_dim` | Observed surface velocity (m/yr) mapped to model edges |
| `inversion` | `C_inversion_results.mat` | `C_hat`, `u_obs_noisy`, `v_obs_noisy`, `N_current`, `aa`, `pp`, `gg`, `oo`, `ps` | Inverted basal friction field for `nevis_inv_partition` |
| `initial_hydrology` | `velocity_inverted.mat` | variable `vv_hydro` with `N`; optionally `u`, `v` | Initial effective pressure and velocity |
| `moulins` | `moulins_2022_coordinates.mat` | variable `moulins` with `x_m`, `y_m` | Moulin coordinates; catchments are rebuilt on the current grid |
| `lakes` | `environs_lakes_2022B_250416.mat` | variable `environs_lakes` | Case-local lake positions, drainage type, timing, volume |
| `runoff` | `racmo_runoff_2022.mat` | variable `racmo_runoff_2022` | Daily surface runoff, `ntime × gg.nIJ`, mm w.e./day |
| `stations` | `station_timeseries_2022.mat` | variable `station_data` with `x_m`, `y_m` | GPS/station points for time-series output |
| manifest | `dataset_manifest.m` | variable `dataset` | Declares the files, variables, and expected fields above |

### How to obtain or rebuild each file

The active pipelines under [`data/pipelines/`](../data/pipelines/)
produce the four data products (`geometry`, `velocity`, runoff,
stations) from raw sources. See
[data/DATA_PIPELINES.md](../data/DATA_PIPELINES.md) for full builder
contracts; this section captures only the per-file expectations.

- **`geometry.mat`** — Source: BedMachine/Morlighem bed and surface.
  Project/crop to a regular x/y grid in meters; store `X_m`, `Y_m`,
  `B_m`, `S_m`. `cfg.geometry.skip` controls model resolution and must
  match downstream gridded products.
- **`velocity.mat`** — Source: MEaSUREs surface velocity mosaic (or
  equivalent observed product). Interpolate onto the model node grid;
  store `u_obs_dim`, `v_obs_dim` in m/yr. NaNs are allowed where
  observations are unavailable — the importer fills with
  nearest-neighbor before mapping to edges.
- **`C_inversion_results.mat`** — Output of `inversion/nevis_inv_C.m`
  (or `nevis_inv_C1.m`) on the same geometry, mask, and velocity. Save
  all variables required by `nevis_inv_partition`. Alternative
  inversion packages may live at
  `inversions/<name>/C_inversion_results.mat` and be selected
  explicitly with `cfg.inversion.file`.
- **`velocity_inverted.mat`** — Hydrology/velocity state saved after
  inversion or a coupled forward initialization. Required:
  `vv_hydro.N` (one value per model node); optional: `vv_hydro.u`
  (size `gg.eIJ`), `vv_hydro.v` (size `gg.fIJ`). Alternative initial
  states pair with their inversion package under
  `inversions/<name>/velocity_inverted.mat` and are selected via
  `cfg.initial_hydrology.file`.
- **`moulins_2022_coordinates.mat`** — Source: raw moulin shapefile
  under `data/raw/gnss_2022/NatalieTurner_Moulins_2022/`. Store
  `moulins.x_m` and `moulins.y_m` in meters in the same local NEVIS
  coordinate frame as the geometry. Set `cfg.moulins.mode =
  'coordinates'`; `nevis_import_region_moulins` rebuilds `ni_m` and
  `sum_m` on the current grid, so this works when bbox changes.
- **Lake catalogue** — Source: mechanistic drainage catalogue for the
  target year/region, copied from raw data into the current case dataset.
  Store `environs_lakes`, or configure a generic table via
  `cfg.lakes.mode = 'table'`. Coordinates in meters share the regional
  origin used by geometry/runoff/moulins/stations. Timing fields must be
  day-of-year indices consistent with
  `cfg.run.t_span_days`.
- **`racmo_runoff_2022.mat`** — Source: RACMO (or equivalent surface
  runoff product) preprocessed to the model node grid. One row per
  time sample, one column per node; mm water equivalent per day unless
  `cfg.runoff.units` and the importer are updated consistently. The
  raw RACMO NetCDFs are useful provenance but are not read at runtime.
- **`station_timeseries_2022.mat`** — Source: GNSS/GPS processing
  workflow. Store `station_data`; each station must have `x_m`, `y_m`
  in the same regional coordinate system as the model grid.

## Data requirements: idealized

Idealized cases have **no** required topography, velocity, inversion,
runoff, moulin, lake, or station data files. Everything is generated
analytically from `cfg`. The only on-disk requirement is the case
config itself (e.g. `cases/<casename>.m`).

## Files moved out of this directory

Stale flat-layout files generated by earlier versions of the workflow,
and historical runtime products, have been archived under
[`data/redundant/`](../data/redundant/) so this directory stays
isolated and minimal. The notable ones:

- `data/redundant/cases/` — old generated `<case_name>_forward.m` /
  `<case_name>_inversion_spinup.m` stubs and unused smoke-test
  fixtures.
- `data/redundant/nevis_regional/` — historical regional MAT products
  no longer read by `nevis_run_regional`, including older runoff
  packages, raw RACMO archives, the legacy moulin file, and prior
  inversion outputs. These are kept for provenance, preprocessing
  history, and older analysis scripts.
