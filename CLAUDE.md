# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

GEMB.jl is a Julia implementation of the Glacier Energy and Mass Balance model. This is a physics-based snow/firn/ice model that simulates surface mass and energy balance processes including:
- Snow grain metamorphism
- Albedo evolution
- Shortwave radiation absorption with depth
- Temperature profile evolution
- Melt and refreezing
- Densification (compaction and wet compaction)
- Layer management

The model is a translation from MATLAB to Julia, maintaining fidelity to the original implementation while leveraging Julia's performance and type system.

## Quick Start

Basic workflow for running GEMB:

```julia
using GEMB
using GEMB_ClimateForcing  # synthetic/real forcing lives in the companion package

# 1. Initialize model parameters
params = initialize_parameters()

# 2. Create climate forcing (DimStack) and convert to a ClimateForcing
ds = simulate_climate_forcing("test_1", 3)  # 3-hourly synthetic data
forcing = GEMB.initialize_forcing(ds)

# 3. Initialize the vertical profile (returns possibly depth-adjusted params)
profile, params = initialize_profile(params, forcing)

# 4. Run the model
output = gemb(profile, forcing, params)

# 5. Extract surface temperature time series
T_surface = surface_timeseries(output.temperature)
```

Note: `simulate_climate_forcing`, the `fit_*` climate functions, and the humidity
conversions were moved out of GEMB.jl into the companion `GEMB_ClimateForcing.jl`
(commit a669134). GEMB.jl keeps only the physics model; a `DimStack` is the neutral
interface, and `initialize_forcing(ds)` (a core method, always available) converts it.
`forcing_climatology` is a core GEMB.jl method that operates on a `ClimateForcing`
(build one with `initialize_forcing` first).

## DimStack Forcing Interface and the GEMB_ClimateForcing Companion

A `DimStack` is GEMB's neutral, producer-agnostic forcing interface. The
conversion `initialize_forcing(::DimStack)` is a **core** method (a method of
`initialize_forcing` in `src/initialize_forcing.jl`, always available — no
extension required): it validates the required layers, a `DateTime`-indexed `Ti`
dimension, and the required metadata, then forwards to the vector method of
`initialize_forcing`.

Any source can produce a conforming `DimStack`. The
[GEMB_ClimateForcing.jl](https://github.com/alex-s-gardner/GEMB_ClimateForcing.jl)
companion is one such producer, loading real climate data (ERA5-Land) from CDS.
It is an optional dependency — install it only when you need real/synthetic
forcing (it is also a test dependency; see below).

Installation (GEMB_ClimateForcing is not yet in the General registry):
```julia
using Pkg
Pkg.add(url="https://github.com/alex-s-gardner/GEMB_ClimateForcing.jl")
```

Usage with the companion producer:
```julia
using GEMB
using GEMB_ClimateForcing

# Download ERA5-Land data for a location
forcing_data = climate_forcing(:era5land, 67.0, -50.0;
                                time_range=(DateTime(2020,1,1), DateTime(2020,12,31)),
                                token=ENV["CDS_API_KEY"])

# Convert DimStack → ClimateForcing (core initialize_forcing method)
cf = GEMB.initialize_forcing(forcing_data)

# Use with GEMB
mp = initialize_parameters()
profile, mp = initialize_profile(mp, cf)
output = gemb(profile, cf, mp)
```

## Development Commands

### Testing

**Prerequisite:** The test target depends on the unregistered companion package
`GEMB_ClimateForcing.jl`, resolved via `[sources]` in `test/Project.toml`
pointing at `../../GEMB_ClimateForcing.jl` (a sibling of this repo). `[sources]`
is only honored on **Julia ≥ 1.11**; on 1.10 `Pkg.test()` fails with
`expected package GEMB_ClimateForcing to be registered`. Check the companion out
as a sibling directory, and keep the CI matrix off the `'1'` alias (which can
resolve to a 1.10.x release).

```bash
# Run all tests
julia --project=. -e 'using Pkg; Pkg.test()'

# Run a specific test file (from project root)
julia --project=. -e 'using GEMB, Test; include("test/test_thermal_conductivity.jl")'

# Or interactively from Julia REPL
julia --project=.
using GEMB, Test
include("test/test_thermal_conductivity.jl")
```

### Documentation
```bash
# Build documentation locally
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

### Package Management
```bash
# Install dependencies
julia --project=. -e 'using Pkg; Pkg.instantiate()'

# Add a new dependency
julia --project=. -e 'using Pkg; Pkg.add("PackageName")'
```

## Architecture

### Core Structure

GEMB follows a modular physics-based architecture:

1. **Types** (`types.jl`): Defines three key structs:
   - `ModelParameters`: All model configuration. Method/option fields use `Symbol` types (e.g., `:Arthern`, `:GardnerSharp`, `:daily`)
   - `ClimateForcing`: Time-series meteorological forcing, a `DimensionalData.AbstractDimStack` subtype. The 13 forcing variables are stack layers (`DimArray`s); `time_step::Int` (seconds) and the scalar means/observation heights live in the stack `metadata` (a `NamedTuple`). Both are reachable by name via property access (`cf.temperature_air`, `cf.time_step`). Time-varying model parameters are stored as Fill arrays. Supports the DimStack API and time indexing: `cf[Ti=At(t)]` returns a `ClimateForcing` (a sub-stack), not a `ClimateForcingStep`
   - `ClimateForcingStep`: Single timestep forcing values (plain struct of scalars for the physics loop)

2. **Initialization** (`initialize_*.jl`):
   - `initialize_parameters()`: Creates `ModelParameters` with defaults or user overrides
   - `initialize_forcing()`: Loads/creates `ClimateForcing` time series
   - `initialize_profile()`: Sets up initial vertical profile of temperature, density, etc.

3. **Physics Modules** (`calculate_*.jl`): Individual physics processes called in sequence:
   - `calculate_grain_size()`: Snow grain metamorphism
   - `calculate_albedo()`: Surface albedo evolution
   - `calculate_shortwave_radiation()`: SW absorption profile
   - `calculate_temperature()`: Energy balance and temperature evolution
   - `calculate_accumulation()`: Snow/rain addition
   - `calculate_melt()`: Melt, runoff, and refreezing
   - `calculate_density()`: Non-melt densification
   - `manage_layers()`: Grid management (merging/splitting)

4. **Integration**:
   - `gemb_core(state, cfs, mp, verbose)`: Single timestep integration calling all physics modules. Accepts a `state` NamedTuple and returns `(state, flux)` where `state` carries forward and `flux` contains budget terms for output
   - `gemb(profile, climate_forcing, mp)` in `gemb_driver.jl`: Main driver function that loops over time, accumulates output

5. **Utilities**:
   - `gemb_spinup()` in `spinup.jl`: Cycles forcing to reach equilibrium (for multi-millennial spinups). Returns the final equilibrated profile and all spinup output.
   - `gemb_profile()`, `gemb_interp()` in `profile_extract.jl`: Extract/interpolate profiles at specific times/depths
   - `surface_timeseries()`: Extract surface values from column arrays
   - `dz2z()`: Convert grid spacing to depth coordinates
   - `forcing_climatology()` in `forcing_climatology.jl`: Averages complete years of a
     `ClimateForcing` (dropping leap day 366 and partial years) into a one-year climatological
     cycle, typically for `gemb_spinup`. Consumes and returns a `ClimateForcing`.

   Climate-forcing generation (`simulate_climate_forcing()`, the `fit_*` fitting functions,
   and humidity conversions like `dewpoint_to_vapor_pressure()`) now live in the companion
   `GEMB_ClimateForcing.jl`, not in this package. They return/consume a `DimStack`; use
   `GEMB.initialize_forcing(ds)` to convert into the `ClimateForcing` type GEMB consumes.

### Data Flow

```
Initialize Profile → Time Loop [
    forcing_step = climate_forcing[Ti=At(t)]
    state, flux = gemb_core(state, forcing_step, mp, verbose)
    Accumulate flux → Store to output at intervals
] → Return DimStack
```

### Key Design Principles

- **MATLAB Fidelity**: Physics match original MATLAB implementation. Comments reference MATLAB line numbers.
- **DimensionalData.jl**: All input/output arrays use `DimArray` with explicit dimensions (`Ti` for time, `Z` for vertical). Indexing uses keyword syntax: `output[:temperature][Z=1:10, Ti=At(t)]`
- **State as NamedTuple**: The column state (temperature, dz, density, etc.) is passed between timesteps as a plain NamedTuple of vectors — no DimStack overhead in the hot loop
- **Symbols for Options**: Model parameters that select methods use `Symbol` (e.g., `albedo_method=:GardnerSharp`, `output_frequency=:daily`)
- **ClimateForcing Indexing**: `ClimateForcing` is an `AbstractDimStack`, so `climate_forcing[Ti(a .. b)]` / `climate_forcing[Ti=At(t)]` slice it like any DimStack and return a `ClimateForcing`. `ClimateForcingStep` (the per-timestep scalar struct) is built separately in the `gemb` hot loop by integer-indexing the unwrapped forcing vectors. Time-varying model parameters (black carbon, cloud properties) are stored as `FillArrays.Fill`-backed DimArrays, ready to become truly time-varying
- **Immutable Parameters**: `ModelParameters` is immutable; create new instance for modifications
- **Energy/Mass Conservation**: When `verbose=true`, `gemb_core()` validates conservation laws each timestep
- **Output Padding**: Profile outputs include padding (`output_padding` parameter) to accommodate column growth without reallocation

### Test Structure

Tests validate against MATLAB reference data stored in `test/` directory (MAT files). Each physics module has a corresponding test file that:
1. Loads MATLAB reference inputs
2. Runs Julia function
3. Compares outputs with reference (typically to ~1e-12 relative tolerance)

The test suite in `runtests.jl` runs all module tests in order of dependency.

### Constants and Physical Models

Physical constants are defined in `constants.jl` (e.g., `C_ICE`, `LF`, `CtoK`) and match MATLAB conventions.

GEMB simulates:
- **Surface Energy Balance**: Radiative (shortwave/longwave) and turbulent (sensible/latent) fluxes using Monin-Obukhov similarity theory
- **Subsurface Thermodynamics**: Heat equation with phase change, meltwater percolation, refreezing
- **Firn Densification**: Compaction schemes (empirical/semi-empirical, e.g., Arthern, Herron-Langway)
- **Hydrology**: "Bucket" scheme for meltwater retention, percolation, refreezing (irreducible water content)
- **Dynamic Albedo**: Long-term memory albedo accounting for grain growth and specific surface area
- **Grid Management**: Lagrangian-style vertical grid that merges/splits layers dynamically

## Testing Against MATLAB

When modifying physics functions, always run the corresponding test to ensure MATLAB equivalence is maintained. Tests use the MAT.jl package to load reference data stored in `test/` directory (MAT files). Validation typically requires ~1e-12 relative tolerance to match MATLAB. If you add new physics, generate reference data using `test/generate_reference_data.m` (MATLAB script).

### Cross-Validation Examples
- `examples/synthetic_example.jl` - Simple synthetic forcing test
- `examples/era5_example.jl` - ERA5 reanalysis data example  
- `examples/compare_outputs.jl` - MATLAB/Julia cross-validation

## DimensionalData.jl Usage

- Always use `Ti` dimension for time indices
- Use `Z` dimension for vertical (depth) indices  
- Use keyword indexing syntax: `array[Ti=At(t)]`, `array[Z=1:n, Ti=At(t)]`
- When extracting values from DimArray for computation, use `Vector{Float64}(...)` or `parent()` to get plain Arrays
- Results are packaged back into DimStack for return
- DimArray indexing with `At()` has zero overhead — dimension info is compiled away

## Performance Considerations

- First run includes JIT compilation time; subsequent runs are 2-5x faster than MATLAB
- GEMB prioritizes computational efficiency for multi-millennial spinups required by deep firn columns
- Memory efficient with minimal allocations after warmup
- The hot loop (`gemb_core`) uses plain NamedTuple of vectors (not DimStack) to avoid overhead
- Use `@time` or BenchmarkTools.jl to measure performance after compilation
