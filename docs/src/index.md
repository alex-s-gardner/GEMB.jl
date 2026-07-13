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

1. **Define Climate Forcing** -- Use [`initialize_forcing`](@ref) to create forcing from time series data, [`simulate_climate_forcing`](@ref) to generate synthetic test data, or use [GEMB_ClimateForcing.jl](https://github.com/alex-s-gardner/GEMB_ClimateForcing.jl) to download ERA5/MERRA-2 data.
2. **Define Model Parameters** -- Use [`ModelParameters`](@ref) to set model configuration (densification model, albedo method, grid geometry, etc.).
3. **Initialize a Column** -- Use [`initialize_profile`](@ref) to create an initial profile of temperature, density, grid spacing, and other column properties.
4. **Run GEMB** -- Pass the profile, climate forcing, and model parameters to the [`gemb`](@ref) function.

```julia
using GEMB

# Initialize model parameters
mp = ModelParameters(output_frequency="daily")

# Generate synthetic climate forcing (3-hour time step)
cf = simulate_climate_forcing("test_1", 3)

# Initialize the firn column profile
profile = initialize_profile(mp, cf)

# Run GEMB
output = gemb(profile, cf, mp)
```

The output is a `DimStack` (from [DimensionalData.jl](https://github.com/rafaqz/DimensionalData.jl)) containing time series of surface fluxes and vertical profiles at the specified output frequency.

## Using Real Climate Data

For production runs with ERA5, ERA5-Land, or MERRA-2 reanalysis data, use the [GEMB_ClimateForcing.jl](https://github.com/alex-s-gardner/GEMB_ClimateForcing.jl) package which automatically downloads and formats climate data:

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

# Convert to GEMB ClimateForcing (automatic via package extension)
cf = GEMB.ClimateForcing(forcing_data)

# Run GEMB
mp = ModelParameters(output_frequency="daily")
profile = initialize_profile(mp, cf)
output = gemb(profile, cf, mp)
```

GEMB.jl includes a package extension that automatically loads when both `GEMB` and `GEMB_ClimateForcing` are imported, providing seamless conversion from downloaded climate data to GEMB's `ClimateForcing` type.

## Spinup

For research applications, the column should be spun up to a quasi-steady state before running simulations with transient forcing. Use [`gemb_spinup`](@ref) to repeat the forcing over multiple cycles:

```julia
using GEMB

mp = ModelParameters(output_frequency="last")

# Generate or load climate forcing
cf = simulate_climate_forcing("test_1", 3)

# Initialize the column
profile = initialize_profile(mp, cf)

# Spin up over 5 cycles
spun_up_profile = gemb_spinup(profile, cf, mp, 5)

# Now run with transient forcing
mp_run = ModelParameters(output_frequency="daily")
output = gemb(spun_up_profile, cf, mp_run)
```

## Model Architecture

After climate forcing, model parameters, and the initial state of the column are defined, the `gemb` function calls `gemb_core` for each time step of the climate forcing. At each time step, `gemb_core` calls a series of physics functions that update the column grain size, albedo, shortwave radiation, temperature, accumulation, meltwater, and density. The `manage_layers` function adjusts the depth and number of vertical layers to ensure layer thicknesses remain within configured bounds.

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
| `manage_layers` | Layer splitting and merging to maintain grid constraints |

## Output Variables

The output `DimStack` contains monolevel (1D time series) and profile (2D depth-time) variables:

### Monolevel Outputs (dimensions: `Ti`)

| Variable | Units | Description |
|----------|-------|-------------|
| `melt` | kg m⁻² | Melt mass |
| `runoff` | kg m⁻² | Runoff mass |
| `refreeze` | kg m⁻² | Refrozen mass |
| `evaporation_condensation` | kg m⁻² | Evaporation (+) or condensation (-) |
| `shortwave_net` | W m⁻² | Net shortwave radiation |
| `longwave_net` | W m⁻² | Net longwave radiation |
| `heat_flux_sensible` | W m⁻² | Sensible heat flux |
| `heat_flux_latent` | W m⁻² | Latent heat flux |
| `albedo_surface` | fraction | Surface albedo |
| `densification_from_compaction` | m | Compaction from densification |
| `densification_from_melt` | m | Compaction from melt |
| `thickness_cumulative` | m | Cumulative thickness change |
| `firn_air_content` | m | Total air height in the firn column |
| `valid_profile_length` | integer | Number of active vertical levels |

### Profile Outputs (dimensions: `Z x Ti`)

| Variable | Units | Description |
|----------|-------|-------------|
| `temperature` | K | Column temperature |
| `dz` | m | Layer thickness |
| `density` | kg m⁻³ | Column density |
| `water` | kg m⁻² | Pore water content |
| `grain_radius` | mm | Effective grain radius |
| `grain_dendricity` | fraction | Grain dendricity (0--1) |
| `grain_sphericity` | fraction | Grain sphericity (0--1) |

## Examples

Example scripts are provided in the `examples/` directory:

- **`synthetic_example.jl`**: Complete workflow using synthetic climate forcing (spinup + run)
- **`era5_example.jl`**: Workflow using ERA5 reanalysis data (with data download instructions)

## Citation

Please cite any use of GEMB as:

> Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass Balance (GEMB): a model of firn processes for cryosphere research, Geosci. Model Dev., 16, 2277--2302, [https://doi.org/10.5194/gmd-16-2277-2023](https://doi.org/10.5194/gmd-16-2277-2023), 2023.

## Author Information

The Glacier Energy and Mass Balance (GEMB) model was created by Alex Gardner, with contributions from Nicole-Jeanne Schlegel and Chad Greene. The Julia implementation (GEMB.jl) is available at [https://github.com/alex-s-gardner/GEMB.jl](https://github.com/alex-s-gardner/GEMB.jl).
