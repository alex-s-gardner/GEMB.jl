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
2. **Define Model Parameters** -- Use [`initialize_parameters`](@ref) to set and validate model configuration (densification model, albedo method, grid geometry, etc.). It builds and checks a [`ModelParameters`](@ref).
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
profile = initialize_profile(mp, cf)

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
profile = initialize_profile(mp, cf)
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
profile = initialize_profile(mp, cf_clim)
spun_up_profile = gemb_spinup(profile, cf_clim, mp; max_iterations=5,
                              convergence_delta_density=0.01)

# Convergence is judged on the column-mean density. Two tests are available, and
# when both are given both must hold: `convergence_delta_density` bounds the change
# between consecutive cycles, while `convergence_drift_density` bounds the trend —
# the least-squares slope of column-mean density against cycle over the trailing
# `drift_window` cycles [kg m⁻³ per cycle]. The trend test is the stricter claim: a
# column creeping steadily at just under the delta tolerance passes the step test
# while still densifying.
spun_up_profile = gemb_spinup(profile, cf_clim, mp; max_iterations=200,
                              convergence_delta_density=0.01,
                              convergence_drift_density=0.005)

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

### Meltwater percolation scheme

GEMB percolates meltwater with a tipping-bucket scheme: water moves cell by cell from the
surface down, refreezing where cold content allows, retained up to the irreducible water
content ([`irreducible_saturation`](@ref)), and routed to runoff at a contiguous run of cells
at or above `impermeable_density` thicker than `impermeable_thickness`. There is no
preferential-flow (heterogeneous, "piping") domain and no Richards-equation matrix flow, and
none is planned.

That is a deliberate choice, not a missing feature. RetMIP (Vandecrux et al., 2020)
intercompared nine firn models at four Greenland sites and found that the three models with
explicit deep or preferential percolation (CFM-Cr, CFM-KM, UppsalaUniDeepPerc) performed worse
than the bucket schemes at three of the four sites — they infiltrated water too deeply and
carried a warm bias in firn temperature at the dry-snow site (Summit) and both percolation
sites (Dye-2 mean error +3.6 to +6.2 °C, KAN\_U +1.8 to +4.7 °C). At Dye-2 in 2016 the CFM
models percolated to 10 m against 2.5 m observed by upward-looking radar, and built
multi-metre near-surface ice slabs where none are observed. Their advantage was confined to
the firn-aquifer site, where only the deep-percolation schemes recharged the aquifer at all.
RetMIP's conclusion (their Sect. 5.2) is that until the physics of preferential flow in firn
is better constrained by field and laboratory observation, the more complex schemes do not
necessarily give better results than simple bucket schemes.

The levers RetMIP does identify for bucket schemes are *when* water is blocked and *how fast*
it then leaves, not how deep it goes. Both are exposed:

**The impermeability criterion** — `impermeable_density` and `impermeable_thickness`. See
[`calculate_melt`](@ref) and the physics notes in the README for the range the
participating models spanned and how it maps onto their skill at the ice-slab site.

**The runoff timescale** — `runoff_method`, with `surface_slope` as the driving hydraulic
gradient. Under the `:instantaneous` default blocked water leaves the
column within the timestep and no cell can ever hold more than its irreducible water. The
other two let it pond into the pore space above the barrier and drain laterally over a finite
timescale:

| `runoff_method` | Runoff law | Reference |
|---|---|---|
| `:instantaneous` | All blocked water leaves within the timestep | The default; what both comparison models and all three RetMIP bucket lineages use |
| `:ZuoOerlemans` | `drain = excess · min(1, Δt/τ)`, `τ = c₁ + c₂·exp(−c₃·S)` | Zuo and Oerlemans (1996) eqs. 21–22, coefficients via Langen et al. (2017) |
| `:Darcy` | `drain = min(excess, ρ_w·Δt·K_sat·K_rel·S)` | Calonne et al. (2012) eq. 6 with van Genuchten (1980) / Yamaguchi et al. (2012) relative permeability |

This is where RetMIP's evidence actually points. The two models with the lowest firn-temperature
error at the ice-slab site KAN\_U — DMIHH (−1.6 °C) and GEUS (+0.6 °C), against a spread
reaching +4.7 °C — were both bucket schemes that *delay* runoff rather than models that
percolate deeper; DMIHH uses the `:ZuoOerlemans` timescale and GEUS a Darcy flux to a virtual
downslope neighbour. At the other end, DTU runs water off immediately and produced runoff
unrealistic enough that RetMIP excluded it from their multi-model mean.

Delayed runoff is also what makes saturated firn representable at all: with the hard
irreducible clamp of `:instantaneous`, a saturated cell cannot exist, and RetMIP (their
Sect. 5.4) note that models so constrained "are incapable of modeling actual aquifers" —
the firn-aquifer site was dropped from their retention evaluation for exactly this reason.
Aquifers form bottom-up under either delayed method, and the `aquifer_thickness` and
`aquifer_depth` outputs report the resulting water table. Both are opt-in: they are inert at
the `:instantaneous` default.

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
| `albedo_broadband` | 1 | `time: mean` | Broadband surface albedo |
| `densification_from_compaction` | m | `time: sum` | Thickness lost to dry compaction |
| `densification_from_melt` | m | `time: sum` | Thickness lost to melt and wet compaction |
| `strain_thinning` | m | `time: sum` | Thickness lost to horizontal ice-dynamic strain (positive = thinning); zero unless `horizontal_strain_rate` is set |
| `thickness_cumulative` | m | `time: mean` | Cumulative thickness change since the start of the run |
| `firn_air_content` | m | `time: mean` | Total air height in the firn column |
| `firn_air_content_10m` | m | `time: mean` | Air height in the top 10 m, the depth-limited form the firn-core literature reports |
| `firn_air_content_20m` | m | `time: mean` | Air height in the top 20 m |
| `percolation_depth` | m | `time: maximum` | Deepest the wetting front reached during the interval; comparable to upward-looking-radar estimates |
| `ice_slab_thickness` | m | `time: point` | Total thickness of cells at or above `impermeable_density`, whether or not they block flow |
| `ice_slab_depth` | m | `time: point` | Depth to the top of the shallowest flow-blocking ice slab; `NaN` when none qualifies |
| `aquifer_thickness` | m | `time: point` | Total thickness of cells holding standing (super-irreducible) water; 0 unless `runoff_method` delays runoff |
| `aquifer_depth` | m | `time: point` | Depth to the water table, i.e. the shallowest cell above irreducible saturation; `NaN` when the column holds no standing water |
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
| `age` | d | Mass-weighted mean age of the cell's mass, from column initialization |

## Examples

Example scripts are provided in the `examples/` directory:

- **`synthetic_example.jl`**: Complete workflow using synthetic climate forcing (spinup + run)
- **`era5_example.jl`**: Workflow using ERA5 reanalysis data (with data download instructions)

## Citation

Please cite any use of GEMB as:

> Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass Balance (GEMB): a model of firn processes for cryosphere research, Geosci. Model Dev., 16, 2277--2302, [https://doi.org/10.5194/gmd-16-2277-2023](https://doi.org/10.5194/gmd-16-2277-2023), 2023.

## Author Information

The Glacier Energy and Mass Balance (GEMB) model was created by Alex Gardner, with contributions from Nicole-Jeanne Schlegel and Chad Greene. The Julia implementation (GEMB.jl) is available at [https://github.com/alex-s-gardner/GEMB.jl](https://github.com/alex-s-gardner/GEMB.jl).
