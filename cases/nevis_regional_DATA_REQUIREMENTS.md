# Data Requirements for `nevis_regional`

This directory is a strict case package for the case named `nevis_regional`.
The workflow function and `./cases/nevis_regional.m` are bound by
`cfg.casename`; a run with another `casename` must use its own
`./cases/<casename>.m` config file.

No missing file, missing variable, or size mismatch should fall back to data in
another directory.

Heavy reusable regional data live in the shared dataset package declared by
`cfg.dataset`:

```matlab
cfg.dataset.name = 'nevis_2022_140km';
cfg.dataset.root = './data/datasets/nevis_2022_140km';
cfg.dataset.manifest = 'dataset_manifest.m';
```

All case-specific model parameters live in `nevis_regional.m`:

- `cfg.oo` for model and solver options.
- `cfg.pd` for dimensional physical parameters.
- `cfg.pp` for nondimensional parameter overrides applied after
  `nevis_nondimension`.

The standardized workflow function should not set case-specific `oo`, `pd`, or
`pp` values directly.

## Case vs Dataset

Parameter sweep cases should copy only this config file and edit the copied
`./cases/<casename>.m`. They should point to the same `cfg.dataset` if they use the
same topography, forcing, inversion, initial state, moulins, lakes, and
stations.

The case directory should contain only:

| File | Purpose |
|---|---|
| `<casename>.m` | Declares case identity, dataset binding, parameters, and run options. |
| `DATA_REQUIREMENTS.md` | Documents the case/dataset contract. |

The dataset directory should contain the manifest plus the heavy files below.
The regional validator checks that each configured file is declared in the
manifest and exists inside `cfg.dataset.root`.

The original runtime MAT files for this template were moved to
`./data/datasets/nevis_2022_140km/` so new
parameter-sweep cases can reuse them without duplicating storage.

## Shared Dataset Runtime Files

These files are required by the current shared dataset
`./data/datasets/nevis_2022_140km/dataset_manifest.m` and are read by
`n2d_regional_template(casename)` through `nevis_setup_region_case`.

| Config section | File | Required variable/fields | Purpose |
|---|---|---|---|
| `geometry` | `geometry.mat` | variable `geometry` with `X_m`, `Y_m`, `B_m`, `S_m` | Builds the model grid, bed, surface, ice thickness mask, hydrology labels, blister labels, and ice boundary labels. |
| `velocity` | `velocity.mat` | variable `velocity` with `u_obs_dim`, `v_obs_dim` | Provides observed surface velocity, in m/yr, mapped onto model edges. |
| `inversion` | `C_inversion_results.mat` | `C_hat`, `u_obs_noisy`, `v_obs_noisy`, `N_current`, `aa`, `pp`, `gg`, `oo`, `ps` | Provides the inverted basal friction/sliding field used by `nevis_inv_partition`. |
| `initial_hydrology` | `velocity_inverted.mat` | variable `vv_hydro` with `N`; optionally `u`, `v` | Provides initial effective pressure and, when present, initial velocity fields. |
| `moulins` | `nevis_170207a.mat` | variable `pp` with `ni_m`, `sum_m` | Provides moulin node indices and catchment matrix. Must match the model grid node count. |
| `lakes` | `mechanistic_drainage_catalogue_2022-main/catalogues/environs_lakes_2022B_250505_archive.mat` | variable `environs_lakes` | Provides lake positions, drainage type, timing, and volume. |
| `runoff` | `racmo_runoff_2022.mat` | variable `racmo_runoff_2022` | Provides daily surface runoff as `ntime x gg.nIJ`, in mm water equivalent per day. |
| `stations` | `station_timeseries_2022.mat` | variable `station_data` with `x_m`, `y_m` | Provides GPS/station points to save model time series. |
| manifest | `dataset_manifest.m` | variable `dataset` | Declares this dataset's files, variables, and expected fields. |

## How To Obtain or Rebuild Each File

### `geometry.mat`

Source: Morlighem/BedMachine-style bed and surface products already subset to
the target region.

Required processing:

1. Project or crop the source data onto a regular x/y grid in meters.
2. Store a struct with fields `X_m`, `Y_m`, `B_m`, and `S_m`.
3. Set `cfg.geometry.skip` so that the downsampled grid has the intended model
   resolution and matches downstream gridded products.

### `velocity.mat`

Source: MEaSUREs surface velocity mosaic or another observed velocity product.

Required processing:

1. Interpolate velocity components onto the same downsampled node grid used by
   the model.
2. Store the result as `u_obs_dim` and `v_obs_dim` in m/yr.
3. Keep NaNs only where observations are unavailable; the importer fills NaNs
   with nearest-neighbor values before mapping to edges.

### `C_inversion_results.mat`

Source: output from the inversion workflow, normally `inversion/nevis_inv_C.m`.

Required processing:

1. Run the inversion on the same geometry, mask, and velocity observations as
   the forward case.
2. Save the variables required by `nevis_inv_partition`: `C_hat`,
   `u_obs_noisy`, `v_obs_noisy`, `N_current`, `aa`, `pp`, `gg`, `oo`, and `ps`.
3. Keep this file inside the shared dataset directory, not inside each
   parameter-sweep case directory. Alternative inversion packages may be stored
   below `inversions/<name>/C_inversion_results.mat` and selected explicitly in
   the case config with `cfg.inversion.file`.

### `velocity_inverted.mat`

Source: velocity/hydrology state saved after inversion or a coupled forward
initialization.

Required processing:

1. Save variable `vv_hydro`.
2. `vv_hydro.N` must have one value per model node.
3. Optional `vv_hydro.u` and `vv_hydro.v` must have sizes `gg.eIJ` and
   `gg.fIJ`.
4. Alternative initial states paired with inversion packages may be stored
   below `inversions/<name>/velocity_inverted.mat` and selected explicitly in
   the case config with `cfg.initial_hydrology.file`.

### `nevis_170207a.mat`

Source: moulin coordinates processed through `nevis_moulins`, or a saved
regional spinup state containing moulin mappings.

Required processing:

1. Use the same grid as the model case.
2. Store `pp.ni_m` and `pp.sum_m`.
3. `size(pp.sum_m, 2)` must equal `gg.nIJ`.

For a new region, prefer storing raw moulin coordinates in a new file and set
`cfg.moulins.mode = 'coordinates'`, with matching `x_field` and `y_field`.

### Lake Catalogue File

Source: mechanistic drainage catalogue for the target year/region.

Required processing:

1. Store variable `environs_lakes`, or configure a generic table via
   `cfg.lakes.mode = 'table'`.
2. Coordinates must be in meters relative to the same regional origin used by
   geometry, runoff, moulins, and stations.
3. Timing fields must be day-of-year indices consistent with the configured
   model `t_span_days`.

### `racmo_runoff_2022.mat`

Source: RACMO or another surface runoff product, preprocessed onto the model
node grid.

Required processing:

1. Convert/crop/project the raw climate product to the model region.
2. Interpolate or aggregate to the same model node grid.
3. Store a matrix with one row per day/time sample and one column per model
   node.
4. Use units of mm water equivalent per day unless `cfg.runoff.units` and the
   importer are updated consistently.

The raw RACMO NetCDF files are useful provenance/intermediate files, but they
are not read by the current runtime script.

### `station_timeseries_2022.mat`

Source: GNSS/GPS station processing workflow.

Required processing:

1. Store variable `station_data`.
2. Each station must have `x_m` and `y_m` in the same regional coordinate
   system as the model grid.

## Files Moved to `../redundant`

The following files used to live in this case directory but are not referenced
by the current `nevis_regional.m` and are not read by `n2d_regional_template`.
They were moved under `./data/redundant/nevis_regional/` to keep this runtime
case package isolated and minimal. They may still be useful for provenance,
preprocessing, older experiments, or other scripts.

Small metadata or local-system files:

- `.DS_Store`

Old model state or historical regional products:

- `0365.mat`
- `mosaicV1.mat`
- `moulin_joughin.mat`
- `CA2011.mat`
- `CA2012.mat`
- `CA2013.mat`
- `CA2014.mat`

Historical or alternate processed runoff products:

- `runoff_2008_nevis140.mat`
- `runoff_2009_nevis140.mat`
- `runoff_2010_nevis140.mat`
- `runoff_2022_nevis.mat`
- `runoff_2022_nevis_noSK_300m.mat`
- `runoff_2023_nevis140.mat`
- `runoff_2023_nevis_noSK_300m.mat`
- `runoff_for_nevis_noSK_2022_300m.m`
- `runoff_for_nevis_noSK_2023_300m.m`

Raw RACMO archive/intermediate files not read at runtime:

- `2022.zip`
- `2023.zip`
- `2022/runoff.2022_AMJ.BN_RACMO2.3p2_ERA5_3h_FGRN055.1km.DD.nc`
- `2022/runoff.2022_JAS.BN_RACMO2.3p2_ERA5_3h_FGRN055.1km.DD.nc`
- `2022/runoff.2022_JFM.BN_RACMO2.3p2_ERA5_3h_FGRN055.1km.DD.nc`
- `2022/runoff.2022_OND.BN_RACMO2.3p2_ERA5_3h_FGRN055.1km.DD.nc`
- `2022/t2m.2022_AMJ.BN_RACMO2.3p2_ERA5_3h_FGRN055.1km.DD.nc`
- `2022/t2m.2022_JAS.BN_RACMO2.3p2_ERA5_3h_FGRN055.1km.DD.nc`
- `2022/t2m.2022_JFM.BN_RACMO2.3p2_ERA5_3h_FGRN055.1km.DD.nc`
- `2022/t2m.2022_OND.BN_RACMO2.3p2_ERA5_3h_FGRN055.1km.DD.nc`
- `2023/runoff.2023_AMJ.BN_RACMO2.3p2_ERA5_3h_FGRN055.1km.DD.nc`
- `2023/runoff.2023_JAS.BN_RACMO2.3p2_ERA5_3h_FGRN055.1km.DD.nc`
- `2023/runoff.2023_JFM.BN_RACMO2.3p2_ERA5_3h_FGRN055.1km.DD.nc`
- `2023/runoff.2023_OND.BN_RACMO2.3p2_ERA5_3h_FGRN055.1km.DD.nc`
- `2023/t2m.2023_AMJ.BN_RACMO2.3p2_ERA5_3h_FGRN055.1km.DD.nc`
- `2023/t2m.2023_JAS.BN_RACMO2.3p2_ERA5_3h_FGRN055.1km.DD.nc`
- `2023/t2m.2023_JFM.BN_RACMO2.3p2_ERA5_3h_FGRN055.1km.DD.nc`
- `2023/t2m.2023_OND.BN_RACMO2.3p2_ERA5_3h_FGRN055.1km.DD.nc`

The `./cases` directory should contain lightweight case config files and
documentation. Heavy files referenced by `<casename>.m` should live in the
declared shared dataset package.
