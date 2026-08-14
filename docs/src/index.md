# GEMB.jl

```@meta
CurrentModule = GEMB
```

The **Glacier Energy and Mass Balance (GEMB)** model is a column model of firn processes for cryosphere research. GEMB.jl is a Julia implementation of the GEMB model, providing high-performance simulation of snow, firn, and ice evolution driven by surface climate forcing.

GEMB models grain growth, albedo, radiative transfer, thermodynamics, accumulation, melt, layer management, and densification within a vertical snow/firn/ice column.

## Installation

GEMB.jl can be installed from the Julia package manager:

```julia
using Pkg
Pkg.add("GEMB")
```

Or in the Pkg REPL (press `]`):

```
pkg> add GEMB
```

## Quick Start

Using GEMB requires four basic steps:

1. **Define Climate Forcing** -- Use [`initialize_forcing`](@ref) to create forcing from time series data, `simulate_climate_forcing` (from GEMB_ClimateForcing.jl) to generate synthetic test data, or use [GEMB_ClimateForcing.jl](https://github.com/alex-s-gardner/GEMB_ClimateForcing.jl) to download ERA5-Land data.
2. **Define Model Parameters** -- Use [`ModelParameters`](@ref) to set model configuration (densification model, albedo method, grid geometry, etc.).
3. **Initialize a Column** -- Use [`initialize_profile`](@ref) to create an initial profile of temperature, density, grid spacing, and other column properties.
4. **Run GEMB** -- Pass the profile, climate forcing, and model parameters to the [`gemb`](@ref) function.

```julia
using GEMB
using GEMB_ClimateForcing

# Initialize model parameters
mp = initialize_parameters(output_frequency=:daily)

# Generate synthetic climate forcing (3-hour time step) — returns DimStack
ds = simulate_climate_forcing("test_1", 3)
cf = GEMB.initialize_forcing(ds)   # convert to ClimateForcing

# Initialize the firn column profile
profile, mp = initialize_profile(mp, cf)

# Run GEMB
output = gemb(profile, cf, mp)
```

The output is a `DimStack` (from [DimensionalData.jl](https://github.com/rafaqz/DimensionalData.jl)) containing time series of surface fluxes and vertical profiles at the specified output frequency.

## Using Real Climate Data

For production runs with ERA5-Land reanalysis data, use the [GEMB_ClimateForcing.jl](https://github.com/alex-s-gardner/GEMB_ClimateForcing.jl) package which automatically downloads and formats climate data:

```julia
# Install GEMB_ClimateForcing (first time only)
using Pkg
Pkg.add(url="https://github.com/alex-s-gardner/GEMB_ClimateForcing.jl")

using GEMB
using GEMB_ClimateForcing

# Download ERA5-Land data for Summit Station, Greenland
forcing_data = climate_forcing(:era5land, 72.58, -38.48;
                                time_range=(DateTime(2020,1,1), DateTime(2020,12,31)),
                                token=ENV["CDS_API_KEY"])

# Convert to GEMB ClimateForcing (core initialize_forcing method)
cf = GEMB.initialize_forcing(forcing_data)

# Run GEMB
mp = initialize_parameters(output_frequency=:daily)
profile, mp = initialize_profile(mp, cf)
output = gemb(profile, cf, mp)
```

GEMB_ClimateForcing produces a `DimStack`, which the core `initialize_forcing(::DimStack)` method converts to GEMB's `ClimateForcing` type. `GEMB_ClimateForcing` is an optional companion package — one producer of a conforming DimStack.

## Spinup

For research applications, the column should be spun up to a quasi-steady state before running simulations with transient forcing. Use [`gemb_spinup`](@ref) to repeat the forcing over multiple cycles:

```julia
using GEMB

using GEMB_ClimateForcing

mp = initialize_parameters(output_frequency=:last)

# Generate synthetic climate forcing (returns DimStack)
ds = simulate_climate_forcing("test_1", 3)
cf = GEMB.initialize_forcing(ds)

# Create a single-year climatological average for spinup
cf_clim = forcing_climatology(cf)

# Initialize the column and spin up
profile, mp = initialize_profile(mp, cf_clim)
spun_up_profile = gemb_spinup(profile, cf_clim, mp; max_iterations=5,
                              convergence_delta_density=0.01)

# The spun-up profile carries provenance: which climatology years were averaged
# and how the spinup converged (`metadata` is re-exported from DimensionalData).
metadata(spun_up_profile)   # (spinup_cycles=…, spinup_converged=…, climatology_n_years=…, …)

# Now run with transient forcing
mp_run = initialize_parameters(output_frequency=:daily)
output = gemb(spun_up_profile, cf, mp_run)

# That provenance propagates onto the output. A profile run without spinup
# instead records `spinup_performed => false`.
metadata(output)
```

## Model Architecture

After climate forcing, model parameters, and the initial state of the column are defined, the `gemb` function calls `gemb_core` for each time step of the climate forcing. At each time step, `gemb_core` calls a series of physics functions that update the column grain size, albedo, shortwave radiation, temperature, accumulation, meltwater, and density. The `manage_layer_thickness` function merges and splits layers to keep thicknesses within their configured bounds and to hold the layer count fixed; total column depth is pinned separately at the end of each time step.

### Physics Modules

| Module | Description |
|--------|-------------|
| `calculate_grain_size` | Evolution of effective grain radius, dendricity, and sphericity |
| `calculate_albedo` | Snow, firn, and ice albedo from grain radius, density, cloud amount |
| `calculate_shortwave_radiation` | Vertical distribution of absorbed shortwave radiation |
| `calculate_temperature` | Temperature profile from energy absorption and thermal diffusion |
| `calculate_accumulation` | Precipitation and deposition added to the column |
| `calculate_melt` | Meltwater production, pore water content, grid adjustment |
| `calculate_density` | Snow/firn densification |
| `manage_layer_thickness` | Layer splitting and merging to maintain grid constraints |

## Output Variables

The output `DimStack` contains monolevel (1D time series) and profile (2D depth-time) variables.

Every layer carries CF-style metadata — `units`, `long_name`, `cell_methods`, and a
`standard_name` where a CF standard name applies exactly — and the stack carries
`Conventions = "CF-1.11"` alongside the run's provenance. The `cell_methods` column below
is that attribute: it records how each value was reduced over the output interval, which
is the difference between a per-interval total, a per-interval average, and a snapshot.
The table lives in [`GEMB_CF_ATTRIBUTES`](@ref); read one layer's attributes with
[`cf_attributes`](@ref).

### Monolevel Outputs (dimensions: `Ti`)

| Variable | Units | `cell_methods` | Description |
|----------|-------|----------------|-------------|
| `melt` | kg m⁻² | `time: sum` | Melt mass produced, including bare-ice melt |
| `runoff` | kg m⁻² | `time: sum` | Meltwater leaving the column |
| `refreeze` | kg m⁻² | `time: sum` | Meltwater refrozen within the column |
| `evaporation_condensation` | kg m⁻² | `time: sum` | Net vapour mass exchange: condensation/deposition (+, mass gain) or evaporation/sublimation (−) |
| `shortwave_net` | W m⁻² | `time: mean` | Net shortwave radiation, positive toward the surface |
| `longwave_net` | W m⁻² | `time: mean` | Net longwave radiation, positive toward the surface |
| `heat_flux_sensible` | W m⁻² | `time: mean` | Sensible heat flux, positive toward the surface |
| `heat_flux_latent` | W m⁻² | `time: mean` | Latent heat flux, positive toward the surface |
| `albedo_surface` | 1 | `time: mean` | Broadband surface albedo |
| `densification_from_compaction` | m | `time: sum` | Thickness lost to dry compaction |
| `densification_from_melt` | m | `time: sum` | Thickness lost to melt and wet compaction |
| `thickness_cumulative` | m | `time: mean` | Cumulative thickness change since the start of the run |
| `firn_air_content` | m | `time: mean` | Total air height in the firn column |
| `valid_profile_length` | 1 | `time: point` | Number of active vertical levels |
| `temperature_air` | K | `time: mean` | Near-surface air temperature (forcing summary) |
| `precipitation` | kg m⁻² | `time: sum` | Total precipitation (forcing summary) |
| `rain` | kg m⁻² | `time: sum` | Liquid fraction of precipitation; snowfall is `precipitation - rain` |

### Profile Outputs (dimensions: `Z x Ti`)

All profile fields are instantaneous snapshots at the output time (`cell_methods = "time: point"`),
top-justified with the surface at row 1. `Z` is a 1-based cell index, not a depth — recover
cell-centre heights with `dz2z(output[:dz])`.

| Variable | Units | Description |
|----------|-------|-------------|
| `temperature` | K | Layer temperature |
| `dz` | m | Layer thickness |
| `density` | kg m⁻³ | Layer bulk density |
| `water` | kg m⁻² | Layer liquid (pore) water content |
| `grain_radius` | mm | Effective grain radius |
| `grain_dendricity` | 1 | Grain dendricity (0--1) |
| `grain_sphericity` | 1 | Grain sphericity (0--1) |
| `albedo` | 1 | Layer direct-beam broadband albedo |
| `albedo_diffuse` | 1 | Layer diffuse-radiation albedo |

## Examples

Example scripts are provided in the `examples/` directory:

- **`synthetic_example.jl`**: Complete workflow using synthetic climate forcing (spinup + run)
- **`era5_example.jl`**: Workflow using ERA5 reanalysis data (with data download instructions)

## Citation

Please cite any use of GEMB as:

> Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass Balance (GEMB): a model of firn processes for cryosphere research, Geosci. Model Dev., 16, 2277--2302, [https://doi.org/10.5194/gmd-16-2277-2023](https://doi.org/10.5194/gmd-16-2277-2023), 2023.

## Author Information

The Glacier Energy and Mass Balance (GEMB) model was created by Alex Gardner, with contributions from Nicole-Jeanne Schlegel and Chad Greene. The Julia implementation (GEMB.jl) is available at [https://github.com/alex-s-gardner/GEMB.jl](https://github.com/alex-s-gardner/GEMB.jl).
