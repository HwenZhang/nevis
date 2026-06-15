# NEVIS Inversion

This note explains what the inversion does, how to run it, what it
expects on disk, and how to point it at a custom dataset built with
`build_nevis_dataset`.

## What it does

Given a model geometry, a domain mask, and an observed surface velocity
field, the inversion recovers a spatially varying **basal sliding
parameter** `C_hat` (and, in the regularized-Coulomb mode, also updates
the effective pressure `N`). The result is the bridge between
observation and the forward hydrology model: it makes the modelled basal
shear stress consistent with the observed surface flow.

Two sliding laws are supported:

| `cfg.sliding_law` | Solver | What gets updated |
|---|---|---|
| `weertman` | [nevis_inv_C.m](nevis_inv_C.m) | `C_hat` only; `N` is fixed by the initial hydrology source |
| `regularized_coulomb` | [nevis_inv_C1.m](nevis_inv_C1.m) | Outer-loop updates of `C_hat` and `N` together |

## Entry point

```matlab
addpath('inversion');
% Normal path: the per-case config materialized by the demo workflow.
[inv, vv_hydro, summary] = nevis_run_inversion('./cases/<case_name>/<case_name>_inversion.m');
% The template under cases/templates/ can also be run directly to validate
% against the reference dataset it is bound to.
```

[nevis_run_inversion.m](nevis_run_inversion.m) loads the config file,
checks `cfg.sliding_law`, ensures file-mode source/hydrology inputs exist,
and dispatches to the matching solver. If `cfg.initial_hydrology.mode =
'file'` and either `cfg.source.state_file` or
`cfg.initial_hydrology.file` is missing, it runs the configured
hydrology-only spinup case and exports the selected timestep first.
Returns:

- `inv` — struct with `C_hat`, observed/noisy velocity, current `N`, and
  the solver state needed for downstream cases.
- `vv_hydro` — a reusable hydrology state (with at least `N`, plus `u`/`v`
  when available). Forward cases pick this up via
  `cfg.initial_hydrology.file`.
- `summary` — diagnostics: cost-function history, schedule used,
  sliding-law tag.

## Required inputs

The inversion is driven entirely by a single config file (`.m`, `.mat`, or
`.json`). The template is [cases/templates/nevis_inversion.m](../cases/templates/nevis_inversion.m);
the relevant sections are:

1. **Dataset binding** — `cfg.dataset.name`, `cfg.dataset.root`,
   `cfg.dataset.manifest`. The dataset directory must contain the four
   pipeline products (`geometry.mat`, `velocity.mat`,
   `racmo_runoff_2022.mat`, `station_timeseries_2022.mat`) plus the
   auxiliary files declared in its `dataset_manifest.m`
   (moulins, lakes catalogue, inversion outputs, initial hydrology). The
   regional validator (`nevis_validate_region_config`) refuses to run if
   the manifest's `dataset.name` does not equal `cfg.dataset.name`.

2. **Source model state** — the file used to initialize grid, parameters,
   and geometry:

   ```matlab
   cfg.source.state_file = 'source_states/<regional_spinup_case>_source.mat';
   ```

   In file mode, `nevis_run_inversion` creates it automatically from
   `cfg.spinup.casename` and `cfg.spinup.timestep_file` if it is missing.
   The inversion uses **this file's `gg`** as the model grid, not the
   dataset's geometry. The dataset supplies observed velocity and forcing.

3. **Observed velocity binding** — `cfg.velocity.file`,
   `cfg.velocity.variable`, `cfg.velocity.u_field`,
   `cfg.velocity.v_field`. Defaults to
   `velocity.mat` → variable `velocity` with fields
   `u_obs_dim`, `v_obs_dim` (m/yr).

4. **Initial hydrology**, via `cfg.initial_hydrology.mode`:

   | mode | required fields | what it does |
   |---|---|---|
   | `file` | `file`, `variable` | loads `vv.N` or `vv.phi` from an existing MAT file |
   | `k_factor` | `k_factor` (scalar 0…1) | `phi = phi_a + k_factor * (phi_0 - phi_a)`; no external file |
   | `constant_N` | `N` (scalar) | uniform dimensionless effective pressure; no external file |

   `k_factor` and `constant_N` are useful bootstraps for the spinup run,
   but they are usually too crude as direct inversion initial conditions.

5. **Output paths**, relative to `cfg.dataset.root`:

   - `cfg.output.inversion_file` — destination for `C_inversion_results.mat`
     (consumed by forward regional cases via `cfg.inversion.file`).
   - `cfg.output.initial_hydrology_file` — destination for
     `velocity_inverted.mat`. The forward case picks this up via
     `cfg.initial_hydrology.file` and bootstraps from `vv_hydro.N`.
   - Optional intermediate files when
     `cfg.output.save_intermediate = true`.

6. **Solver / regularization / optimizer schedules** — see the template
   for the full list. No hidden defaults; every option used must be
   declared.

## Outputs

After a successful run the dataset directory contains, under
`inversions/<cfg.name>/`:

- `C_inversion_results.mat` — the inverted sliding field plus everything
  required to reproduce the result (`C_hat`, `u_obs_noisy`,
  `v_obs_noisy`, `N_current`, `aa`, `pp`, `gg`, `oo`, `ps`).
- `velocity_inverted.mat` — the reusable `vv_hydro` state. Forward
  regional cases consume this as their initial hydrology.

A forward regional run bound to the same dataset then references these
two files via `cfg.inversion.file` and `cfg.initial_hydrology.file`
inside its case config — see
[cases/templates/nevis_regional.m](../cases/templates/nevis_regional.m).

## Recommended Spinup-Then-Inversion Workflow

For each domain/bbox, configure a hydrology-only spinup, then let the
inversion entry point ensure the needed files:

1. Build the dataset package for the domain.
2. Create or derive a regional spinup config from
   [cases/templates/nevis_regional.m](../cases/templates/nevis_regional.m). The default
   example is [cases/templates/nevis_regional_inversion_spinup.m](../cases/templates/nevis_regional_inversion_spinup.m):

   ```matlab
   cfg.oo.include_ice = 0;
   cfg.inversion.mode = 'uniform';
   cfg.initial_hydrology.mode = 'k_factor';
   cfg.lakes.spinup_volume_scale = 0;
   cfg.run.t_span_days = 1:365;
   ```

3. Configure the spinup source and the inversion file inputs:

   ```matlab
   cfg.spinup.casename = 'nevis_regional_inversion_spinup';
   cfg.spinup.timestep_file = '0365.mat';
   cfg.source.state_file = 'source_states/nevis_regional_inversion_spinup_source.mat';

   cfg.initial_hydrology.mode = 'file';
   cfg.initial_hydrology.file = 'initial_hydrology/nevis_regional_inversion_spinup_0365.mat';
   cfg.initial_hydrology.variable = 'vv';
   ```

4. Run inversion:

   ```matlab
   [inv, vv_hydro, summary] = nevis_run_inversion('./cases/templates/nevis_inversion.m');
   ```

If the configured source/hydrology files are missing, `nevis_run_inversion`
checks that `cfg.spinup.casename` is hydrology-only (`cfg.oo.include_ice =
0`), runs it if the spinup output/timestep is absent, exports the files,
and then starts inversion. There is no fallback to another case.

## Running inversion on a custom-bbox dataset

The supported path is the demo workflow ([test.m](../test.m)). Set
`case_name` and `bbox_km`, keep `run_inversion = true`, and run it. The
workflow:

1. Builds a self-contained dataset under `data/datasets/<case_name>/`,
   including `moulins_2022_coordinates.mat` (projected from the raw 2022
   shapefile onto the new domain), the lake catalogue, and a full
   `dataset_manifest.m`. See [data/DATA_PIPELINES.md](../data/DATA_PIPELINES.md).
2. Materializes `cases/<case_name>/<case_name>_inversion.m` and
   `cases/<case_name>/<case_name>_inversion_spinup.m`, both bound to the new
   dataset, with the source-state and initial-hydrology paths pre-wired.
3. Calls `nevis_run_inversion` on that config. Because
   `cfg.initial_hydrology.mode = 'file'` and the source/hydrology files do not
   exist yet, the runner first runs the hydrology-only spinup case on the new
   grid, exports the selected timestep, and then starts the inversion.

So a custom bbox needs no manual file copying — moulins, lakes, manifest,
spinup, and source state are all produced on the new grid automatically.

If you drive `nevis_run_inversion` outside `test.m`, reproduce the same
bindings by hand in a per-case `.m` config:

```matlab
cfg.dataset.name = '<case_name>';
cfg.dataset.root = './data/datasets/<case_name>';
cfg.spinup.casename = '<case_name>_inversion_spinup';   % hydrology-only, include_ice = 0
cfg.spinup.timestep_file = '0365.mat';
cfg.source.state_file = 'source_states/<case_name>_inversion_spinup_source.mat';
cfg.initial_hydrology.mode = 'file';
cfg.initial_hydrology.file = 'initial_hydrology/<case_name>_inversion_spinup_0365.mat';
cfg.initial_hydrology.variable = 'vv';
cfg.output.inversion_file = 'inversions/<cfg.name>/C_inversion_results.mat';
cfg.output.initial_hydrology_file = 'inversions/<cfg.name>/velocity_inverted.mat';
```

then `nevis_run_inversion('./cases/<case_name>/<case_name>_inversion.m')`. The
runner generates the source state from the spinup if it is missing; there is no
fallback to another case.

## Common pitfalls

- **`DatasetNameMismatch`** — `cfg.dataset.name` and the manifest's
  `dataset.name` must agree exactly. The auto-generated manifest names
  the dataset after the output-directory basename; rename either side to
  match.
- **Grid mismatch** — the inversion uses `gg` from
  `cfg.source.state_file`, so `pp.ni_m` from a different grid will not
  index the right nodes. Re-run the spinup or pick a source state that
  matches the current dataset's grid.
- **Wrong velocity units** — the importer expects m/yr (`cfg.velocity.units
  = 'm_per_year'`). If your dataset is in m/s, change the units field
  and ensure the importer supports it.
- **Initial hydrology absent** — with `mode = 'file'`,
  `nevis_run_inversion` generates `cfg.initial_hydrology.file` from
  `cfg.spinup.casename` if it is missing. `k_factor` or `constant_N`
  should be used for the spinup bootstrap, not as the default final
  inversion input.
