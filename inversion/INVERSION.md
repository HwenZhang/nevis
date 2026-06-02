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
[inv, vv_hydro, summary] = nevis_run_inversion('./cases/nevis_inversion.m');
```

[nevis_run_inversion.m](nevis_run_inversion.m) loads the config file,
reads `cfg.sliding_law`, and dispatches to the matching solver. Returns:

- `inv` — struct with `C_hat`, observed/noisy velocity, current `N`, and
  the solver state needed for downstream cases.
- `vv_hydro` — a reusable hydrology state (with at least `N`, plus `u`/`v`
  when available). Forward cases pick this up via
  `cfg.initial_hydrology.file`.
- `summary` — diagnostics: cost-function history, schedule used,
  sliding-law tag.

## Required inputs

The inversion is driven entirely by a single config file (`.m`, `.mat`, or
`.json`). The template is [cases/nevis_inversion.m](../cases/nevis_inversion.m);
the relevant sections are:

1. **Dataset binding** — `cfg.dataset.name`, `cfg.dataset.root`,
   `cfg.dataset.manifest`. The dataset directory must contain the four
   pipeline products (`geometry.mat`, `velocity.mat`,
   `racmo_runoff_2022.mat`, `station_timeseries_2022.mat`) plus the
   auxiliary files declared in its `dataset_manifest.m`
   (moulins, lakes catalogue, inversion outputs, initial hydrology). The
   regional validator (`nevis_validate_region_config`) refuses to run if
   the manifest's `dataset.name` does not equal `cfg.dataset.name`.

2. **Source model state** — `cfg.source.state_file`. A `.mat` file
   containing `pp`, `pd`, `ps`, `gg`, `aa`, and `oo` from a prior
   regional run. The inversion uses **this file's `gg`** as the model
   grid, not the dataset's geometry. The dataset only supplies observed
   velocity and forcing; the grid topology is fixed by the source state.

3. **Observed velocity binding** — `cfg.velocity.file`,
   `cfg.velocity.variable`, `cfg.velocity.u_field`,
   `cfg.velocity.v_field`. Defaults to
   `velocity.mat` → variable `velocity` with fields
   `u_obs_dim`, `v_obs_dim` (m/yr).

4. **Initial hydrology**, one of three modes via
   `cfg.initial_hydrology.mode`:

   | mode | required fields | what it does |
   |---|---|---|
   | `k_factor` | `k_factor` (scalar 0…1) | `phi = phi_a + k_factor * (phi_0 - phi_a)`; no external file |
   | `constant_N` | `N` (scalar) | uniform dimensionless effective pressure; no external file |
   | `file` | `file`, `variable` | loads `vv.N` (or `phi`) from disk; use a prior `velocity_inverted.mat` to warm-start |

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
[cases/nevis_regional.m](../cases/nevis_regional.m).

## Running inversion on a custom-bbox dataset

`build_nevis_dataset(out_dir, [], bbox_km)` produces the four core
products plus an auto-generated `dataset_manifest.m`, but **does not**
generate:

- a source model state on the new grid (`cfg.source.state_file`)
- moulin / lake catalogue files matching the new domain
- prior inversion / initial-hydrology products

These have to be brought in before the inversion can run. Concretely:

1. Run `build_nevis_dataset` to produce `<out_dir>/` with `geometry.mat`,
   `velocity.mat`, `racmo_runoff_2022.mat`,
   `station_timeseries_2022.mat`, `dataset_manifest.m`.
2. Open `<out_dir>/dataset_manifest.m` and append the auxiliary sections
   modelled on
   [data/datasets/nevis_2022_140km/dataset_manifest.m](../data/datasets/nevis_2022_140km/dataset_manifest.m):
   `dataset.moulins`, `dataset.lakes`, `dataset.inversion`,
   `dataset.initial_hydrology`.
3. Place the auxiliary files at the paths declared in step 2 (copy or
   symlink the moulins/lakes files from the reference dataset; produce a
   source state aligned with the new grid).
4. Derive a config from the template at runtime (this is what `test.m`
   does after a custom-bbox build):

   ```matlab
   clear cfg
   run('./cases/nevis_inversion.m');
   [~, dataset_name] = fileparts(output_dir);
   cfg.dataset.name = dataset_name;
   cfg.dataset.root = output_dir;
   cfg.name = sprintf('%s_%s', cfg.name, dataset_name);
   cfg.output.inversion_file = fullfile('inversions', cfg.name, 'C_inversion_results.mat');
   cfg.output.initial_hydrology_file = fullfile('inversions', cfg.name, 'velocity_inverted.mat');
   cfg.source.state_file = '<path to source state on the new grid>';
   inv_config = fullfile(output_dir, 'nevis_inversion.mat');
   save(inv_config, 'cfg');
   ```

   Alternatively, copy `cases/nevis_inversion.m` to a new `.m` config and
   edit the same fields by hand.
5. Run:

   ```matlab
   nevis_run_inversion(inv_config);
   ```

Until step 3's source state is available, the inversion will error out
when `nevis_inv_C` loads `cfg.source.state_file`. Reuse the reference
dataset for end-to-end testing first; only switch to a custom bbox once
the rebuild and the source state are aligned.

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
- **Initial hydrology absent** — `mode = 'file'` requires the file to
  exist before the run; `mode = 'k_factor'` or `'constant_N'` is the
  safe choice for the first pass.
