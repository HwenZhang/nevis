# Data Requirements for `nevis_idealized`

This directory is a strict analytical data package for the case named
`nevis_idealized`.

The only required case file is:

| File | Purpose |
|---|---|
| `nevis_idealized.m` | Defines `cfg`, including all analytical geometry, initial hydrology, moulin/lake placement, forcing, physical parameters, solver options, and run control. |

There are no required topography, velocity, inversion, runoff, moulin, lake, or
station data files for this idealized case. The workflow generates those pieces
analytically from `cfg`.

The config supports these analytical modes:

- `cfg.geometry.mode = 'hewitt_2013_sqrt'` for the built-in flat-bed Hewitt-style geometry.
- `cfg.geometry.mode = 'function'` for custom `bed_function` and `surface_function`.
- `cfg.initial_hydrology.mode = 'k_factor_cavity_equilibrium'` for the built-in pressure/sheet initialization.
- `cfg.initial_hydrology.mode = 'function'` for a custom function that returns `[vv, aa]`.
- `cfg.moulins.mode = 'coordinates_fractional'` for x/y fractions of domain length/width.
- `cfg.moulins.mode = 'coordinates_m'` for explicit x/y positions in meters.
- `cfg.moulins.mode = 'function'` for a custom position function.

Each moulin can have its own direct input. Set
`cfg.moulins.input_m3_per_s` to a scalar or a vector with one entry per active
moulin, or set `cfg.moulins.input_function` to return scalar/vector m^3/s as a
function of model time.

To create another idealized case:

1. Copy `./cases/nevis_idealized.m` to `./cases/<new_casename>.m`.
2. Edit the copied config file.
3. Set `cfg.casename = '<new_casename>'`.
4. Keep `cfg.case_root = './cases'`.
5. Edit only the config values and analytical functions in `<new_casename>.m`.
6. Run `n2d_idealized_template('<new_casename>')`.

The standardized workflow does not read case-specific parameters from the
runner function. Missing required config fields are errors.
