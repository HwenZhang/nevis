# NEVIS Tutorial: Build a Case and Vary Moulin / Lake / Runoff Inputs

This is a standalone, hands-on guide to creating a regional NEVIS run like the
one in [`test.m`](test.m) and then changing its **moulin**, **lake drainage**,
and **surface runoff** forcing by editing config files only — no core-model
edits, no manual data shuffling.

It assumes the raw inputs are already in place under `data/raw/` (BedMachine,
MEaSUREs, RACMO, GNSS, moulin shapefile, lake catalogue). If you need to
(re)build the gridded dataset products from raw sources, see
[`data/DATA_PIPELINES.md`](data/DATA_PIPELINES.md). The contract for a case
package is in [`cases/README.md`](cases/README.md).

---

## 1. The mental model

A run has two halves:

| Half | What it is | Where it lives |
|---|---|---|
| **Dataset package** | Heavy, reusable gridded data (geometry, velocity, runoff, stations) + moulin coordinates + lake catalogue + manifest | `data/datasets/<case_name>/` |
| **Case configs** | Lightweight `cfg.*` files that bind to a dataset and set physics / numerics / forcing knobs | `cases/<case_name>/` |

`test.m` is the orchestrator. For one `case_name` it:

1. Builds the dataset package (`run_dataset_build`).
2. Runs the inversion to obtain basal slipperiness in the sliding law. The inversion auto-runs a hydrology-only spinup (`run_inversion`).
3. Runs the forward model — lakes off, spinup-like year, with the inverted
   basal field and a hydrology-only spinup timestep as the hydrology state
   (`run_forward`).
4. Runs the drainage case — lakes on, chained from the forward run's day-365
   hydrology state (`run_drainage`).

The **first** time you run a new `case_name`, `nevis_prepare_demo_workflow`
materializes four **self-contained** config files into `cases/<case_name>/`:

```
cases/<case_name>/
  <case_name>_inversion.m          inversion config
  <case_name>_inversion_spinup.m   hydrology-only spinup case
  <case_name>_forward.m            forward regional case (lakes OFF)
  <case_name>_drainage.m           drainage case (lakes ON)
```

Each file inlines the relevant template
([`cases/templates/`](cases/templates/)) and appends the workflow's
dataset/path bindings. After that, **the workflow leaves them alone** — so to
change the physics or forcing for this case, you edit the per-case file. The
templates are only re-read on first materialization or when
`reinit_case = true`.

> **Key rule:** edit `cases/<case_name>/<case_name>_*.m`, not the templates,
> to change an existing case. Set `reinit_case = true` in `test.m` only when
> you want to throw away local edits and regenerate from the current templates.

---

## 2. Create a new case (the `test.m` path)

Open [`test.m`](test.m) and edit only the **USER SETTINGS** block:

```matlab
%% USER SETTINGS
case_name = 'my_run_2022';      % dataset + config prefix for this run
bbox_km   = [];                 % [] = default 140 km region;
                                %  or [xmin, xmax, ymin, ymax] in km

run_dataset_build = true;       % build geometry/velocity/runoff/stations + moulins/lakes
run_inversion     = true;       % inversion (auto-runs spinup if needed)
run_forward       = true;       % forward regional run (lakes off)
run_drainage      = true;       % drainage case (lakes on); needs forward 0365.mat

reinit_case       = false;      % true: regenerate cases/<case_name>/ from templates

% First-run-only seeds (baked into the per-case files):
sliding_law          = 'weertman';   % or 'regularized_coulomb'
spinup_days          = 1:365;        % hydrology-only spinup output days
spinup_timestep_file = '0365.mat';   % spinup timestep that initializes inversion
```

Then run it:

```matlab
test
```

The first run creates `data/datasets/my_run_2022/` (geometry, velocity,
`racmo_runoff_2022.mat`, stations, `moulins_2022_coordinates.mat`,
`environs_lakes_2022B_250416.mat`, and `dataset_manifest.m`) and
`cases/my_run_2022/` with the four config files.

**Re-running a stage only.** Once the dataset and configs exist, you can flip
individual `run_*` flags to `false` to skip stages. For example, after a full
run, iterate on the drainage stage alone with:

```matlab
run_dataset_build = false;
run_inversion     = false;
run_forward       = false;
run_drainage      = true;
```

`case_name` is validated as a MATLAB identifier (letters, digits,
underscores), so keep it simple: `my_run_2022`, not `my-run/2022`.

---

## 3. Where the forcing knobs live

Moulins, lakes, and runoff are configured in the **forward** and **drainage**
per-case files:

- `cases/<case_name>/<case_name>_forward.m` — lakes effectively off
  (`spinup_volume_scale = 0` from the inlined regional template; the forward
  override does not turn the catalogue on). It uses the inversion result for
  `cfg.inversion.file`, restores hydrology from the hydrology-only spinup
  timestep, and initializes coupled ice velocity from the inversion hydrology
  output.
- `cases/<case_name>/<case_name>_drainage.m` — lakes on
  (`cfg.lakes.mode = 'environs_lakes_catalogue'`, `spinup_volume_scale = 1`).
  It chains hydrology from the forward run's day-365 timestep and initializes
  coupled ice velocity from the inversion hydrology output.

Both files end in a workflow override block (`% ====== workflow overrides ======`
or `% ====== workflow drainage overrides ======`). That block is appended
**after** the inlined template, so **whatever it sets wins**. Two ways to
change a value:

1. **Edit the override block at the bottom** (recommended — it is the last
   word and easy to find), or
2. Edit the corresponding line in the inlined template section above, *as long
   as the override block does not also set it.*

The full set of `cfg.moulins`, `cfg.lakes`, and `cfg.runoff` fields — with
inline documentation of every option — is in the template
[`cases/templates/nevis_regional.m`](cases/templates/nevis_regional.m). The
sections below summarize the ones you will actually turn.

---

## 4. Set up the model parameters: `cfg.oo`, `cfg.pd`, `cfg.pp`

Physics and numerics live in three structs, each applied at a different stage
of model setup. They are documented field-by-field in
[`cases/templates/nevis_regional.m`](cases/templates/nevis_regional.m); this
section explains *how* they fit together and which knobs you actually turn.

| Struct | What it holds | When it is applied |
|---|---|---|
| `cfg.oo` | Model & solver **options** (mostly switches and method names) | After `nevis_defaults`, before case setup |
| `cfg.pd` | **Dimensional** physical parameters (SI units) | Scalars first, then any function handles |
| `cfg.pp` | **Non-dimensional** overrides | After `nevis_nondimension` (and after ice-parameter updates) |

The validator
([`src/nevis_validate_region_config.m`](src/nevis_validate_region_config.m))
requires the full set of fields listed in the template to be **present and
non-empty** in each struct. So tune values in place — **don't delete fields**.
The safest place to override is the workflow override block at the bottom of
the per-case file (it is the last word), but for `oo`/`pd`/`pp` the override
block usually leaves them alone, so editing the inlined template lines in the
body is fine too.

### 4a. `cfg.oo` — options and solver switches

Numeric `0`/`1` are used for legacy logical switches. The ones you most often
change:

```matlab
cfg.oo.include_ice        = 1;   % 1 couples hydrology to the ice-velocity solve; 0 hydrology-only
cfg.oo.evaluate_variables = 1;   % 1 expands derived variables after the final step (nevis_backbone)
cfg.oo.adjust_boundaries  = 1;   % 1 lets the timestepper update pressure boundary nodes
cfg.oo.N_coupling         = 1;   % 1 couples effective pressure into sliding/ice dynamics
cfg.oo.U_coupling         = 1;   % 1 couples velocity into hydrology opening/melt terms
cfg.oo.max_iter_new       = 50;  % max Newton iterations per nonlinear solve
cfg.oo.iterative_solver   = 0;   % 0 direct linear solver; 1 iterative path where supported
cfg.oo.partition_ratio    = 1.0; % fraction of inverted C assigned to C1 (0..1); 1 => C2=0
```

Ice-flow boundary treatment is chosen by name (see the template for the full
list of supported strings):

```matlab
cfg.oo.boundary_method      = 'stress_l_vel_tbl';  % e.g. 'default','per_tblr','vel_tblr','stress_tblr',...
cfg.oo.mask_boundary_method = 'stress_free';        % or 'velocity'
```

Diagnostics (keep off for production speed): `cfg.oo.display_residual`,
`cfg.oo.plot_residual`, `cfg.oo.visualize_vel`.

> A **hydrology-only** run (e.g. a spinup) sets `cfg.oo.include_ice = 0`. The
> coupled forward/drainage runs keep it `1`. The spinup case file
> (`<case_name>_inversion_spinup.m`) is where the hydrology-only `oo` lives.

### 4b. `cfg.pd` — dimensional physical parameters (SI)

These are physical quantities in SI units. Fields may be **scalars** or
**function handles** of `pd`; handles are evaluated *after* the scalar fields
are set, so a derived parameter can depend on earlier ones:

```matlab
cfg.pd.u_b     = @(pd) 100/pd.ty;   % sliding speed [m/s] (100 m/yr -> m/s)
cfg.pd.tau_b   = 60e3;              % characteristic basal shear stress [Pa]
cfg.pd.k_s     = 0.1;               % sheet flux coefficient [1/Pa/s]
cfg.pd.melt    = @(pd) (pd.G + pd.u_b*pd.tau_b)/pd.rho_w/pd.L;  % basal melt rate [m/s]: (G + u_b*tau_b)/(rho_w*L)
cfg.pd.kappa_b = 5e-11;             % blister leakage relaxation coefficient (dimensionless; used when relaxation_term = 1)
cfg.pd.alpha_b = 0;                 % relaxation rate of the blister [s^-1] (used when relaxation_term = 0); 0 disables it
cfg.pd.mu      = 5.0;               % water viscosity [Pa s] (blister leakage)
```

Physical meanings and units above are taken from
[`src/nevis_defaults.m`](src/nevis_defaults.m) and
[`src/nevis_nondimension.m`](src/nevis_nondimension.m) — the authoritative
source for every parameter. Inline comments in
[`cases/templates/nevis_regional.m`](cases/templates/nevis_regional.m) match
those.

> **Units matter.** `cfg.pd` is the only struct in dimensional SI. To set a
> rate in "per year," divide by `pd.ty` (seconds per year) inside a handle, as
> with `u_b` above. Melt/runoff conversions elsewhere use `pd.td` (seconds per
> day).

### 4c. `cfg.pp` — non-dimensional overrides

Applied after non-dimensionalization. Mostly regularization floors plus the
pressure-insensitive sliding coefficient:

```matlab
cfg.pp.eps_reg      = 0.01;   % ice-flow regularization / viscosity floor
cfg.pp.Ub_reg       = 1e-16;  % sliding-speed regularization floor
cfg.pp.N_slide_reg  = 1e-16;  % effective-pressure sliding regularization floor
cfg.pp.taud_reg     = 1e-16;  % driving-stress regularization floor
cfg.pp.C2           = 0;      % non-dim pressure-insensitive sliding coefficient (usually set by partitioning)
```

### 4d. A typical parameter-sweep recipe

To sweep one parameter while holding the dataset, inversion, and initial state
fixed, make a hand-written sibling case (so `reinit_case` can't clobber it) and
override just that field — see the worked example in §10. For example, a basal
shear-stress sweep is three sibling files differing only in:

```matlab
cfg.casename = 'my_run_2022_tau40';   cfg.pd.tau_b = 40e3;
cfg.casename = 'my_run_2022_tau60';   cfg.pd.tau_b = 60e3;
cfg.casename = 'my_run_2022_tau80';   cfg.pd.tau_b = 80e3;
```

> The inversion config (`<case_name>_inversion.m`) has its **own** parameter
> sections (`cfg.solver`, `cfg.regularization`, `cfg.prior`, `cfg.opts_inv`,
> `cfg.forward_hydrology`); see [`inversion/INVERSION.md`](inversion/INVERSION.md).
> The `oo`/`pd`/`pp` here are for the forward/drainage and spinup cases.

---

## 5. Vary the moulin input

Moulins are point sinks that collect catchment runoff. In a regional run you
vary them by **which coordinates** are supplied; catchments are rebuilt on the
current grid (`cfg.moulins.mode = 'coordinates'`), so this is robust to a
changed `bbox_km`.

### 5a. The local NEVIS coordinate frame (M1 origin)

Every coordinate field used below — `cfg.moulins.x_m/y_m`,
`cfg.lakes.x_m/y_m`, the runoff source points, `cfg.stations.x_m/y_m`, and the
`X_m/Y_m` arrays inside `geometry.mat` — is expressed in **metres in the local
NEVIS frame**. That frame is a translated Polar Stereographic projection whose
**origin is the M1 moulin** (Stevens et al. 2015), located at approximately
**(68.72° N, 49.53° W)** on the Greenland Ice Sheet. The frame is built with
the standard Greenland Polar Stereographic parameters:

| Parameter | Value | Meaning |
|---|---|---|
| origin (lat, lon) | `(68.72, -49.53)` | M1 moulin: `(x, y) = (0, 0)` |
| `radius` (`a`) | `6378137.0` m | WGS84 semi-major axis |
| `eccen` (`e`) | `0.08181919` | WGS84 eccentricity |
| `lat_true` | `70°` | Standard parallel |
| `lon_posy` | `-45°` | Meridian aligned with +Y |

So when a config says `x_m = -12000`, that point is 12 km west of M1 in the
Polar Stereographic plane; `y_m = 5000` is 5 km north of M1 along the +Y
meridian (`-45°`). The `bbox_km` in `test.m` (`[xmin, xmax, ymin, ymax]`) is
in the same frame.

To convert between latitude/longitude and the NEVIS frame, project with these
parameters and then translate by `(x_M1, y_M1)`:

```matlab
% lat/lon -> NEVIS (x_m, y_m)
[x_M1, y_M1] = polarstereo_fwd(68.72, -49.53, 6378137.0, 0.08181919, 70, -45);
[x_ps, y_ps] = polarstereo_fwd(lat,   lon,    6378137.0, 0.08181919, 70, -45);
x_m = x_ps - x_M1;   y_m = y_ps - y_M1;

% NEVIS (x_m, y_m) -> lat/lon
[lat, lon] = polarstereo_inv(x_m + x_M1, y_m + y_M1, 6378137.0, 0.08181919, 70, -45);
```

All gridded products in a dataset package (`geometry.mat`, `velocity.mat`,
`racmo_runoff_2022.mat`, `station_timeseries_2022.mat`) must share this origin
and these projection parameters, or the importer will place them at the wrong
location on the model grid. See
[`docs/GPS_DATA_LOCATIONS_AND_EXTRACTION.md`](docs/GPS_DATA_LOCATIONS_AND_EXTRACTION.md)
§§1, 8 for the full conversion chain used by the bundled MEaSUREs/BedMachine
pipelines.

### 5b. Use a different set of moulins

Point the case at a different coordinate file. Put the file in the case's
dataset directory and declare it in the manifest (see §9), then in
`<case_name>_forward.m` / `<case_name>_drainage.m`:

```matlab
cfg.moulins.file     = 'moulins_alt_set.mat';   % under data/datasets/<case_name>/
cfg.moulins.variable = 'moulins';
cfg.moulins.mode     = 'coordinates';
cfg.moulins.x_field  = 'x_m';   % easting in metres, local NEVIS frame
cfg.moulins.y_field  = 'y_m';   % northing in metres
```

The `.mat` must hold a struct (variable `moulins`) with `x_m` and `y_m`
vectors in metres in the same local NEVIS frame as `geometry.mat` — see §5a
(M1 moulin at the origin). To create a custom set in MATLAB:

```matlab
moulins.x_m = [ -12000;  3000;  21000 ];   % metres
moulins.y_m = [   5000; -8000;  14000 ];
save('data/datasets/my_run_2022/moulins_alt_set.mat', 'moulins');
```

### 5c. Add or remove a few moulins

Load the generated `moulins_2022_coordinates.mat`, edit the `x_m`/`y_m`
vectors, and save under a new name (keep the original for provenance), then
point the config at the new file as above.

### 5d. Domain-edge handling

Two knobs control what happens to moulins near/outside the domain:

```matlab
cfg.moulins.keep_all       = false;  % true: try to retain all configured moulins
cfg.moulins.move_to_domain = true;   % true: snap out-of-domain moulins to nearest valid node
```

> **Note:** in a regional run, moulin *discharge* is the routed surface runoff
> for each moulin's catchment — it is not set per-moulin in the config. To
> change how much water reaches the moulins, change the **runoff** forcing
> (§7). (Per-moulin prescribed input via `cfg.moulins.input_m3_per_s` /
> `input_function` is an *idealized*-case feature; see
> [`cases/README.md`](cases/README.md).)

---

## 6. Vary the lake drainage input

Lakes are the drainage events turned on in the **drainage** case. The importer
([`src/nevis_import_region_lakes.m`](src/nevis_import_region_lakes.m)) supports
three `cfg.lakes.mode` values.

### 6a. Catalogue mode (default)

```matlab
cfg.lakes.file                = 'environs_lakes_2022B_250416.mat';
cfg.lakes.variable            = 'environs_lakes';
cfg.lakes.mode                = 'environs_lakes_catalogue';
cfg.lakes.spinup_volume_scale = 1;    % 0 removes lake water; 1 = catalogue volumes
cfg.lakes.time_offset_days    = 0;    % shift ALL drainage times by N days
```

Quick, no-new-data variations on the bundled catalogue:

- **Turn lakes off** (make the drainage case behave like the forward case):
  set `cfg.lakes.spinup_volume_scale = 0`.
- **Scale all lake volumes**: `spinup_volume_scale = 0.5` halves every lake's
  injected volume.
- **Shift the drainage season**: `cfg.lakes.time_offset_days = -10` moves every
  event 10 days earlier. (Make sure the shifted times still fall inside
  `cfg.run.t_span_days`.)

The catalogue importer keeps only `drainage_type_num` 1 (fast/hydrofracture,
forced to a 2-hour duration) and 2 (slow), and reads volume/timing per lake
from the catalogue. To change *individual* lakes, build your own catalogue or
use table/manual mode below.

### 6b. Manual mode — hand-specified lakes (no data file)

Best for prescribing a small number of drainage events directly in the config.
Replace the lake override block with:

```matlab
cfg.lakes.mode          = 'manual';
cfg.lakes.x_m           = [ -5000;  12000 ];   % metres, local NEVIS frame
cfg.lakes.y_m           = [  8000;  -3000 ];
cfg.lakes.volume_m3     = [ 2.0e6;  5.0e6 ];   % scalar => same volume for all
cfg.lakes.drainage_day  = [ 195;    201    ];  % day-of-year (absolute)
cfg.lakes.duration_days = [ 2/24;   1.0    ];  % 2/24 = a 2-hour fast event
```

Notes:
- Use `cfg.lakes.drainage_day` for **absolute** day-of-year timing (matches
  `cfg.run.t_span_days`). Alternatively use
  `cfg.lakes.drainage_after_start_days` to time events relative to run start.
- `volume_m3` and `duration_days` may be scalars (broadcast to all lakes) or
  per-lake vectors.
- Make the run resolve the drainage window finely. The default drainage
  `t_span_days` already refines around mid-season:
  `[1:180, 180.01:0.01:220, 221:365]`. If your events fall outside that
  window, widen the fine band so the short drainage pulse is resolved.

### 6c. Table mode — generic catalogue from a file

For a different year/region catalogue stored as a struct or table:

```matlab
cfg.lakes.mode        = 'table';
cfg.lakes.file        = 'lakes_my_region.mat';
cfg.lakes.variable    = 'lakes_table';        % the struct/table variable name
cfg.lakes.x_field     = 'x_m';
cfg.lakes.y_field     = 'y_m';
cfg.lakes.volume_field = 'volume_m3';
cfg.lakes.start_field  = 'drainage_time_start'; % day-of-year
cfg.lakes.end_field    = 'drainage_time_end';   % duration = end - start
```

Coordinates in metres must share the regional origin used by
geometry/runoff/moulins; timing fields must be day-of-year indices consistent
with `cfg.run.t_span_days`. The file must be declared in the manifest (§9).

---

## 7. Vary the surface runoff input

Runoff is the distributed surface meltwater forcing, configured under
`cfg.runoff` and consumed by
[`src/nevis_import_region_runoff.m`](src/nevis_import_region_runoff.m).

### 7a. Swap the runoff product

```matlab
cfg.runoff.file     = 'racmo_runoff_2022.mat';   % under data/datasets/<case_name>/
cfg.runoff.variable = 'racmo_runoff_2022';
cfg.runoff.units    = 'mm_we_per_day';
cfg.runoff.mode     = 'file';
```

The variable must be a `ntime × gg.nIJ` matrix — **one row per time sample, one
column per model node** — in mm water-equivalent per day. The column count must
equal the number of model nodes for the current grid, or the importer errors
with a size mismatch. To use a different year or a regridded product, build that
matrix on the same grid (see [`data/DATA_PIPELINES.md`](data/DATA_PIPELINES.md)
for the RACMO builder), save it, declare it in the manifest, and point
`cfg.runoff.file` at it.

### 7b. Smoothing

```matlab
cfg.runoff.smooth_days  = 5;      % moving-average window (days); 0/1 disables
cfg.runoff.use_smoothed = false;  % true: run on smoothed runoff;
                                  % false: keep raw for the run, store smoothed for diagnostics
```

Smoothing is volume-conserving (it rescales each node's series back to its raw
total). Set `use_smoothed = true` to actually drive the model with the smoothed
field.

### 7c. Scaling / extra melt / routing toggles

```matlab
cfg.runoff.subtract_lake_volume = true;   % scale runoff down by total lake/runoff volume
                                          %  to avoid double-counting lake water
cfg.runoff.distributed_input    = true;   % apply runoff as distributed surface input
cfg.runoff.surface_runoff       = true;   % master enable for surface runoff forcing
cfg.runoff.meltE_mm_per_day     = 0;      % uniform extra melt (mm w.e./day) on top of runoff
cfg.runoff.meltE_ramp_days      = 20;     % e-folding ramp for the uniform extra melt
```

Common variations:

- **Turn surface runoff off entirely** (e.g. to isolate lake forcing):
  `cfg.runoff.surface_runoff = false`.
- **Add a uniform melt background**: `cfg.runoff.meltE_mm_per_day = 5` with a
  `meltE_ramp_days` ramp.
- **Don't discount runoff for lakes**: `cfg.runoff.subtract_lake_volume =
  false` (otherwise total runoff is reduced by the lake/runoff volume ratio so
  lake water isn't counted twice).

> The legacy `cfg.runoff.RACMO_runoff` switch stays `false`: runoff is already
> preprocessed into a node matrix, so the RACMO-specific pathway is not used.

---

## 8. Preparing the RACMO runoff input (format)

The runtime never reads raw RACMO. It reads a single preprocessed matrix that
[`data/pipelines/build_racmo_runoff_2022.m`](data/pipelines/build_racmo_runoff_2022.m)
produces. Understanding both ends — the **runtime format** and the **raw
sources** — lets you supply your own runoff product.

### 8a. The runtime format (what the model consumes)

`data/datasets/<case_name>/racmo_runoff_2022.mat` holds **one variable**:

| Property | Value |
|---|---|
| Variable name | `racmo_runoff_2022` (must match `cfg.runoff.variable`) |
| Shape | `ntime × nNodes` — **rows = time samples (days), columns = model nodes** |
| Node count | must equal `gg.nIJ` for the current grid (else the importer errors) |
| Node order | column `k` = the `k`-th entry of `reshape(X_m(1:skip:end,1:skip:end), [], 1)` and matching `Y_m` |
| Units | mm water-equivalent per day (`cfg.runoff.units = 'mm_we_per_day'`) |
| Time base | row `d` = day-of-year `d`; the importer indexes it by model time in days |

The importer
([`src/nevis_import_region_runoff.m`](src/nevis_import_region_runoff.m))
converts to the model's internal melt scale, optionally smooths
(`smooth_days`), optionally discounts lake volume (`subtract_lake_volume`), and
exposes it as `pp.runoff_function(t)`. **Column order and node count are tied to
the grid**, so a runoff matrix built for one `bbox_km`/`skip` cannot be reused
on a different grid — rebuild it.

### 8b. Build it from the bundled sources

The supported path is the dataset builder, which picks the right RACMO source
for your bbox automatically:

```matlab
addpath('data/pipelines');
raw_root = 'data/raw';
bbox_km = [];  % or [xmin, xmax, ymin, ymax] in km

% Build geometry first (it defines the grid), then runoff onto that grid:
build_geometry('data/datasets/my_run_2022', raw_root, bbox_km);     % writes geometry.mat
out = build_racmo_runoff_2022('data/datasets/my_run_2022', raw_root); % writes racmo_runoff_2022.mat

% ...or do all four gridded products + manifest at once:
% report = build_nevis_dataset('data/datasets/my_run_2022', raw_root, bbox_km);
```

`build_racmo_runoff_2022` requires `geometry.mat` in the output dir (it reads
`X_m`, `Y_m`, `skip` from it), then:

1. Loads a runoff source as scattered points `(Xv, Yv)` in metres (local NEVIS
   frame, M1 at the origin) with a `days × points` value matrix.
2. Builds the query grid `Xq = X_m(1:skip:end, 1:skip:end)` (same as every
   other gridded product).
3. For each day, interpolates with `scatteredInterpolant(..., 'natural','none')`
   and sets out-of-domain NaNs to `0`.
4. Saves the `days × nodes` matrix as `racmo_runoff_2022` (`-v7.3`).

Source priority (highest first):

1. **Regional 300 m cache** — `data/raw/racmo/2022/runoff_2022_nevis_noSK_300m.mat`,
   a struct `runoff_2022_nevis` with fields `X_m`, `Y_m` (metres) and `runoff`
   (`days × points`). Used when the geometry bbox fits within a 5 km tolerance.
2. **Raw pan-Greenland NetCDFs** — four 2022 seasonal files
   `runoff.2022_{JFM,AMJ,JAS,OND}.BN_RACMO2.3p2_ERA5_3h_FGRN055.1km.DD.nc`. The
   builder reads variable `runoffcorr` plus `x`/`y`, subtracts the M1 origin to
   move into the local frame, and crops to the bbox + 2 km pad.

### 8c. Supply your own runoff product (any source)

You have two clean options.

**Option 1 — produce the final matrix directly.** Build a `ntime × gg.nIJ`
matrix yourself (mm w.e./day) in the model node order described in §8a, save
it, declare it in the manifest (§9), and point the config at it:

```matlab
% my_runoff is days x nNodes, mm w.e./day, columns in model-node order
racmo_runoff_2022 = my_runoff;          % keep the variable name, or set cfg.runoff.variable
save('data/datasets/my_run_2022/runoff_myproduct.mat', 'racmo_runoff_2022', '-v7.3');
```

```matlab
cfg.runoff.file     = 'runoff_myproduct.mat';
cfg.runoff.variable = 'racmo_runoff_2022';
cfg.runoff.units    = 'mm_we_per_day';
cfg.runoff.mode     = 'file';
```

The easiest way to get the node order right is to interpolate onto
`X_m(1:skip:end, 1:skip:end)` / `Y_m(1:skip:end, 1:skip:end)` from
`geometry.mat` and `reshape` column-major — exactly what the builder does in
§8b.

**Option 2 — feed the builder.** If your product is daily fields on a regular
or scattered grid in metres (local NEVIS frame), wrap it in the cache struct
format (`X_m`, `Y_m`, `runoff` with `runoff` being `days × points`) and save it
as the cache file the builder loads, or adapt `load_runoff_raw` for your raw
format. Then run `build_racmo_runoff_2022` so the regridding/zero-fill matches
the reference pipeline.

> **Units / time length.** Keep mm w.e./day unless you also change the importer
> (the `'mm_we_per_day'` conversion is hard-wired). The matrix should have at
> least as many day-rows as the largest day in `cfg.run.t_span_days`; the
> bundled product is 365 rows for a full 2022 year. A different year is just a
> different matrix under a different filename, declared in the manifest.

---

## 9. When a new data file needs a manifest entry

The validator
([`src/nevis_validate_region_config.m`](src/nevis_validate_region_config.m))
checks that configured data files exist under `cfg.dataset.root` and, for
file-reading modes, are declared in `dataset_manifest.m`. Some modes skip the
manifest check (`'manual'`, `'coordinates'`, `'uniform'`, `'result_timestep'`,
etc.), but skipping the manifest check does **not** always mean skipping file
I/O: `cfg.moulins.mode = 'coordinates'` still reads `cfg.moulins.file`.

So you need a manifest entry **only** when you point a `file`-reading section
at a new file. That covers:

- A new **runoff** file (`cfg.runoff` is always `mode = 'file'`).
- A new **lake** file in `'environs_lakes_catalogue'` or `'table'` mode.
- A new **moulin** file — note `'coordinates'` mode still requires the file to
  exist, but skips the manifest match check. Matching the manifest keeps the
  package self-describing; declaring it is good practice.

To declare it, open `data/datasets/<case_name>/dataset_manifest.m` and set the
relevant section's `file` (and `files`/`directory` if you keep several), e.g.:

```matlab
dataset.runoff.file = 'racmo_runoff_2023.mat';
dataset.lakes.file  = 'lakes_my_region.mat';
```

The config `file` must match the manifest's `file` (or appear in its `files`
list, or sit under its declared `directory`), or the validator raises
`DatasetManifestMismatch`.

`'manual'`-mode lakes need **no** file and **no** manifest entry — that is the
lowest-friction way to experiment with lake events. Regional moulins are file
based in both supported modes (`'coordinates'` and `'node_indices'`), so a
custom moulin set still needs a `.mat` file even though manifest matching is
optional for `'coordinates'`.

---

## 10. Worked example: a runoff + lake-timing sensitivity case

Goal: same geometry/velocity/inversion as `my_run_2022`, but (a) smoothed
runoff with a uniform melt background, and (b) two hand-specified hydrofracture
lakes 10 days earlier than usual.

1. Run `test.m` once with `case_name = 'my_run_2022'` to build the dataset and
   the four configs (or reuse an existing case).

2. Make a hand-written drainage variant so the workflow won't touch it. Copy
   the managed drainage file:

   ```
   cp cases/my_run_2022/my_run_2022_drainage.m \
      cases/my_run_2022/my_run_2022_sensitivity.m
   ```

3. In `my_run_2022_sensitivity.m`, edit the override block at the bottom:

   ```matlab
   cfg.casename = 'my_run_2022_sensitivity';   % MUST match the file/casename

   % --- runoff: smoothed + uniform melt background ---
   cfg.runoff.smooth_days      = 7;
   cfg.runoff.use_smoothed     = true;
   cfg.runoff.meltE_mm_per_day = 4;
   cfg.runoff.meltE_ramp_days  = 20;

   % --- lakes: two manual hydrofracture events, 10 days early ---
   cfg.lakes.mode          = 'manual';
   cfg.lakes.x_m           = [ -5000; 12000 ];
   cfg.lakes.y_m           = [  8000; -3000 ];
   cfg.lakes.volume_m3     = [ 3.0e6; 6.0e6 ];
   cfg.lakes.drainage_day  = [ 185;   191   ];   % 10 days earlier than the catalogue
   cfg.lakes.duration_days = 2/24;               % fast events, broadcast to both

   % make sure the run resolves the early drainage window
   cfg.run.t_span_days = [1:170, 170.01:0.01:210, 211:365];
   ```

   Leave `cfg.dataset.*`, `cfg.inversion.*`, `cfg.initial_hydrology.*`, and
   `cfg.initial_ice_velocity.*` pointing at the same files the managed drainage
   case used (they were copied with the file) so it chains from the same
   forward/inversion state. In the copied drainage case, hydrology is chained
   from the forward run's saved timestep while initial ice velocity is loaded
   from the inversion hydrology output.

4. Run just this variant directly (no need to re-run `test.m`):

   ```matlab
   addpath('src'); addpath('src/workflow'); addpath('inversion');
   [tt, vv, info, state] = nevis_run_regional('my_run_2022_sensitivity');
   ```

Because lakes are `'manual'` and runoff reuses the existing dataset file, no
manifest edits are needed. Results land under `results/my_run_2022_sensitivity/`.

---

## 11. Quick reference

| I want to… | Edit (in `<case_name>_drainage.m` / `_forward.m`) | New data file? | Manifest entry? |
|---|---|---|---|
| Couple/decouple ice dynamics | `cfg.oo.include_ice` (1/0) | no | no |
| Change Newton iterations / solver | `cfg.oo.max_iter_new`, `cfg.oo.iterative_solver` | no | no |
| Change basal sliding speed / shear stress | `cfg.pd.u_b` (handle), `cfg.pd.tau_b` (Pa) | no | no |
| Change sheet conductivity | `cfg.pd.k_s` | no | no |
| Change a non-dim regularization floor | `cfg.pp.eps_reg`, `cfg.pp.*_reg` | no | no |
| Use a different moulin set | `cfg.moulins.file` (+ `x_field`/`y_field`) | yes | recommended |
| Add/remove a few moulins | edit `moulins` struct, save new file, repoint `cfg.moulins.file` | yes | recommended |
| Turn lakes off | `cfg.lakes.spinup_volume_scale = 0` | no | no |
| Scale all lake volumes | `cfg.lakes.spinup_volume_scale = <s>` | no | no |
| Shift all drainage times | `cfg.lakes.time_offset_days = <n>` | no | no |
| Hand-pick a few lakes | `cfg.lakes.mode = 'manual'` (+ x/y/vol/day/duration) | no | no |
| Use another lake catalogue | `cfg.lakes.mode = 'table'` (+ field map) or `'environs_lakes_catalogue'` | yes | yes |
| Swap runoff product/year | `cfg.runoff.file` (`ntime × nIJ`, mm w.e./day; see §8) | yes | yes |
| Build runoff matrix for a new grid | `build_racmo_runoff_2022` after `build_geometry` (§8b) | yes | yes |
| Smooth runoff | `cfg.runoff.smooth_days`, `use_smoothed` | no | no |
| Add uniform melt | `cfg.runoff.meltE_mm_per_day`, `meltE_ramp_days` | no | no |
| Turn runoff off | `cfg.runoff.surface_runoff = false` | no | no |
| Regenerate configs from templates | `reinit_case = true` in `test.m` | no | no |

---

## 12. Gotchas

- **Casename must match.** A per-case file's `cfg.casename` must equal the name
  passed to `nevis_run_regional`, and (for managed files) the folder name.
  The resolver guards this.
- **Don't edit templates to change one case.** Templates only seed *new* cases
  (or `reinit_case = true`). Edit the per-case file.
- **Don't delete `cfg.oo`/`cfg.pd`/`cfg.pp` fields.** The validator requires
  every field listed in the template to be present and non-empty — change
  values in place, don't remove lines.
- **`cfg.pd` is dimensional SI; `cfg.pp` is non-dimensional.** Mixing them up
  (e.g. setting a "per year" rate without dividing by `pd.ty` inside a handle)
  silently gives wrong physics.
- **`reinit_case = true` overwrites local edits** in `cases/<case_name>/`.
  Keep experiments in a hand-written sibling (e.g. `*_sensitivity.m`) the
  workflow never owns.
- **Timing must fall inside `cfg.run.t_span_days`** and be finely resolved
  where drainage pulses occur, or short fast events get smeared or missed.
- **Runoff column count must equal the model node count** for the current grid;
  changing `bbox_km` or `cfg.geometry.skip` changes that count, so rebuild the
  runoff matrix (rerun the dataset build) rather than reusing one from another
  grid.
- **Lake/runoff double counting:** with `subtract_lake_volume = true`, the
  runoff is discounted by the total lake/runoff volume ratio. If you swap in a
  very different lake set or runoff product, re-check that the discount is what
  you intend.
