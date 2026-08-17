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

GEMB.jl is the reference implementation of GEMB. It began as a translation of an earlier
MATLAB model, but as of v2.0.0 it is developed independently: the two are maintained
separately and are no longer expected to agree numerically. Changes to the physics are
judged on their own merits and against the published literature, not against the MATLAB
output (see [Changing the physics](#changing-the-physics)).

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

# 3. Initialize the vertical profile (its grid is sized to the forcing climate)
profile = initialize_profile(params, forcing)

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
profile = initialize_profile(mp, cf)
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

### Benchmarking

`bench/opt_bench.jl` is the canonical performance + numerical-fingerprint driver
(used by the `julia_optimize` skill). It runs the representative hot path — a
75-year climatological spinup (`gemb_spinup` on `forcing_climatology(cf)`)
followed by a multi-decade 3-hourly transient run — and emits a numerical
snapshot (per-field sum/min/max/count) so optimizations can be checked for
bit-level equivalence against a baseline. It needs the `GEMB_ClimateForcing.jl`
companion checked out as a sibling (same as the test target).

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
   - `manage_layer_thickness()`: Merge/split cells into their thickness bands; restore the fixed cell count
   - `apply_horizontal_strain!()` (in `grid_ops.jl`): Ice-dynamic layer thinning at constant density, `dz *= exp(-strain_rate*dt)`. Runs after `calculate_density` and before `trim_bottom!`; a no-op at the `horizontal_strain_rate = 0.0` default

4. **Integration**:
   - `gemb_core(state, cfs, mp, verbose)`: Single timestep integration calling all physics modules. Accepts a `state` NamedTuple and returns `(state, flux)` where `state` carries forward and `flux` contains budget terms for output
   - `gemb(profile, climate_forcing, mp)` in `gemb_driver.jl`: Main driver function that loops over time, accumulates output

5. **Plotting (extension)**:
   - `gemb_plot_output()` is exported but its methods live in `ext/GEMBMakieExt.jl`, a package extension weak-dep'd on `Makie` (triggered by `CairoMakie`/`GLMakie`/`WGLMakie`). Calling it without a Makie backend loaded errors with a load hint. `src/plotting.jl` holds the core (backend-agnostic) plotting stub/exports.

6. **Utilities**:
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

- **Cited physics**: Every scheme traces to a published source, cited in the docstring of the function that implements it. Some comments still reference line numbers in the MATLAB model the code was originally translated from; these are historical provenance, not a specification to conform to.
- **DimensionalData.jl**: All input/output arrays use `DimArray` with explicit dimensions (`Ti` for time, `Z` for vertical). Indexing uses keyword syntax: `output[:temperature][Z=1:10, Ti=At(t)]`
- **State as NamedTuple**: The column state (temperature, dz, density, etc.) is passed between timesteps as a plain NamedTuple of vectors — no DimStack overhead in the hot loop
- **Symbols for Options**: Most model parameters that select methods use `Symbol`, read in an `if`/`elseif` chain (e.g., `albedo_method=:GardnerSharp`, `output_frequency=:daily`). This is the older convention and still covers the majority of the option fields.
- **Type dispatch for algorithm options**: Where the option selects a whole *algorithm* rather than a coefficient set, it is a singleton type under an abstract supertype, dispatched on rather than branched on. `thermal_solver::S` (`ExplicitThermal()` / `ImplicitThermal() <: AbstractThermalSolver`) is the worked example: `ModelParameters{S<:AbstractThermalSolver}` is parameterized on it, `_thermal_solve!(::ExplicitThermal, ...)` and `_thermal_solve!(::ImplicitThermal, ...)` replace the branch, selection resolves at compile time, and an invalid value is unconstructible rather than caught in `validate_parameters`. A third scheme is a new subtype plus one method, with no central branch to edit. Prefer this form for new algorithm-level options.
- **Description and scratch are separate**: `ModelParameters` is a pure description and holds no mutable state, so one `mp` can be read by any number of threads stepping columns concurrently. This is a design requirement, not an incidental property — never put a buffer in `mp` or in a solver singleton. Reusable scratch lives in a per-run workspace threaded through the call chain: `ThermalWorkspace` (`src/types.jl`) is the worked example, created by `gemb`/`gemb_spinup` (or defaulted per call), carrying one buffer set per scheme, selected by `_solver_workspace(solver, ws)` dispatch and grown by the grow-only `_resize_workspace!`. Two rules make the reuse sound and are pinned by tests in `test/test_calculate_temperature.jl`: buffers are indexed by the column length (never `length(buffer)`), and every entry is written before it is read. A new scheme with scratch adds a buffer struct plus a `_solver_workspace` method.
- **ClimateForcing Indexing**: `ClimateForcing` is an `AbstractDimStack`, so `climate_forcing[Ti(a .. b)]` / `climate_forcing[Ti=At(t)]` slice it like any DimStack and return a `ClimateForcing`. `ClimateForcingStep` (the per-timestep scalar struct) is built separately in the `gemb` hot loop by integer-indexing the unwrapped forcing vectors. Time-varying model parameters (black carbon, cloud properties) are stored as `FillArrays.Fill`-backed DimArrays, ready to become truly time-varying
- **Immutable Parameters**: `ModelParameters` is immutable; create new instance for modifications
- **Energy/Mass Conservation**: When `verbose=true`, `gemb_core()` validates conservation laws each timestep
- **Fixed-length column**: The column holds a constant cell count and a constant total depth (`sum(profile[:dz])`, chosen at initialization within the `column_depth_max` ceiling) for the whole run, enforced each timestep by the two controllers in `grid_ops.jl`. Profile outputs are sized exactly to the column and top-justified (surface at row 1), so there is no padding and no NaN scanning

### Test Structure

Each physics module has a corresponding `test/test_*.jl` file. Validation comes from four
sources, in rough order of strength:

1. **Closed-form checks** — where a scheme has an analytic solution, the test computes it
   independently and compares to ~1e-12 (e.g. the Arthern grain-growth integral in
   `test_calculate_grain_size.jl`, checked in SI against the module's mm working units).
2. **Published values** — coefficients and reference points taken from the paper the
   scheme is cited to.
3. **Independent implementations** — agreement with the Community Firn Model and IMAU-FDM
   where the algebra is the same; see `docs/src/cfm_comparison.md` and
   `docs/src/imau_fdm_comparison.md`.
4. **Invariants** — conservation of mass and energy (`gemb_core` with `verbose=true`),
   monotonicity, bounds, and the physical limits of each scheme.

On top of those, `test/test_synthetic_regression.jl` pins three whole-model fields as a
self-consistency fingerprint. It is a *change detector*, not a physics validation: an
intended change to the defaults or to a scheme is expected to move it, and re-pinning is
part of making such a change. Every re-pin is recorded in that file with what moved it and
by how much, so the history stays auditable.

The test suite in `runtests.jl` runs all module tests in order of dependency.

### Constants and Physical Models

Physical constants are defined in `constants.jl` (e.g., `LF`, `CtoK`, `GRAVITY`), each with
its source in a comment. Numerical tolerances live there too (commit `a5247ea`) rather than
as literals at their use sites. New constants belong in that file, not inline.

GEMB simulates:
- **Surface Energy Balance**: Radiative (shortwave/longwave) and turbulent (sensible/latent) fluxes using Monin-Obukhov similarity theory
- **Subsurface Thermodynamics**: Heat equation with phase change, meltwater percolation, refreezing
- **Firn Densification**: Compaction schemes (empirical/semi-empirical, e.g., Arthern, Herron-Langway)
- **Hydrology**: "Bucket" scheme for meltwater retention, percolation, refreezing (irreducible water content)
- **Dynamic Albedo**: Long-term memory albedo accounting for grain growth and specific surface area
- **Grid Management**: Lagrangian-style vertical grid that merges/splits layers dynamically

## Changing the physics

When modifying a physics function, run its `test/test_*.jl` file and then the full suite.
Two things follow from GEMB.jl being independent of the MATLAB model:

- **The bar is the literature, not the old output.** Justify a change by citing the source
  it comes from, and put that citation in the docstring. Where the change touches a scheme
  the Community Firn Model or IMAU-FDM also implements, check it against theirs and record
  the comparison in the relevant `docs/src/*_comparison.md` page.
- **New behaviour is opt-in unless a default change is the point.** Prefer adding a new
  `Symbol` option value, which leaves every existing run bit-identical. Changing a default
  is a deliberate, separately-argued act — the four default changes in v2.0.0 are the
  worked example.

If a change moves model output, add a one-line entry to `docs/src/physics_notes.md` stating
what changed and what it affects, concisely and factually. Do not editorialize. If the change also alters the synthetic regression fingerprint beyond its
tolerances, re-pin `test/test_synthetic_regression.jl` and append the reason and the
measured deltas to the history at the top of that file.

If the change corrects a defect that also exists in the
[MATLAB model](https://github.com/alex-s-gardner/GEMB), opening an issue there is a
courtesy to its users, but the two codebases are no longer kept in step.

### Examples
- `examples/synthetic_example.jl` - Simple synthetic forcing test
- `examples/era5_example.jl` - ERA5 reanalysis data example

## DimensionalData.jl Usage

- Always use `Ti` dimension for time indices
- Use `Z` dimension for vertical (depth) indices  
- Use keyword indexing syntax: `array[Ti=At(t)]`, `array[Z=1:n, Ti=At(t)]`
- When extracting values from DimArray for computation, use `Vector{Float64}(...)` or `parent()` to get plain Arrays
- Results are packaged back into DimStack for return
- DimArray indexing with `At()` has zero overhead — dimension info is compiled away

## Performance Considerations

- First run includes JIT compilation time; measure only after warmup
- GEMB prioritizes computational efficiency for multi-millennial spinups required by deep firn columns
- Memory efficient with minimal allocations after warmup
- The hot loop (`gemb_core`) uses plain NamedTuple of vectors (not DimStack) to avoid overhead
- Use `@time` or BenchmarkTools.jl to measure performance after compilation
