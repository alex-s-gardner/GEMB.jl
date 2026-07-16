# GEMB.jl: Glacier Energy and Mass Balance Model

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://alex-s-gardner.github.io/GEMB.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://alex-s-gardner.github.io/GEMB.jl/dev/)
[![Build Status](https://github.com/alex-s-gardner/GEMB.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/alex-s-gardner/GEMB.jl/actions/workflows/CI.yml?query=branch%3Amain)

## Overview

GEMB.jl is a Julia implementation of the Glacier Energy and Mass Balance (GEMB, the "B" is silent) model - a comprehensive one-dimensional physical model designed to simulate the surface energy balance and vertical firn evolution of glaciers and ice sheets. It couples atmospheric forcing with subsurface thermodynamics and densification physics to resolve the evolution of temperature, density, water content, and grain properties over time.

GEMB is a column model (no horizontal communication) of intermediate complexity, prioritizing computational efficiency to accommodate the multi-millennial spin-ups required for initializing deep firn columns. It is used for interpreting satellite altimetry data, firn studies, surface mass balance inversion from satellite data, ice core studies, uncertainty quantification and model exploration in cryosphere research. A complete description of GEMB can be found in [*Gardner et al*., 2023](https://doi.org/10.5194/gmd-16-2277-2023).

This Julia implementation maintains fidelity to the original MATLAB version while leveraging Julia's performance, type system, and modern ecosystem including DimensionalData.jl for labeled arrays.

## Key Capabilities

GEMB simulates a wide range of physical processes critical to glacier health:

* **Surface Energy Balance (SEB):** Resolves radiative fluxes (shortwave/longwave) and turbulent heat fluxes (sensible/latent) using Monin-Obukhov similarity theory.
* **Subsurface Thermodynamics:** Solves the heat equation with phase change, meltwater percolation, and refreezing.
* **Firn Densification:** Simulates the compaction of snow into firn and ice using empirical or semi-empirical schemes.
* **Hydrology:** Tracks meltwater retention (irreducible water content), percolation, and refreezing using a "bucket" scheme.
* **Dynamic Albedo:** Models albedo evolution with long-term memory, accounting for grain growth and specific surface area.
* **Grid Management:** Utilizes a dynamic Lagrangian-style vertical grid that evolves with accumulation and ablation, automatically merging and splitting layers to maintain numerical stability.

## Installation

```julia
using Pkg
Pkg.add("GEMB")
```

Or for the development version:

```julia
using Pkg
Pkg.add(url="https://github.com/alex-s-gardner/GEMB.jl")
```

## Quick Start

```julia
using GEMB

# 1. Initialize model parameters
params = initialize_parameters()

# 2. Create or load climate forcing data
forcing = simulate_climate_forcing("test_1", 3)  # 3-hourly synthetic data

# 3. Initialize the vertical profile
profile = initialize_profile(params, forcing)

# 4. Run the model
output = gemb(profile, forcing, params)

# 5. Extract surface temperature time series
T_surface = surface_timeseries(output.temperature)
```

## Basic Workflow

Using GEMB requires four basic steps:

1. **Define Climate Forcing:** Use `initialize_forcing()` to create climate forcing from observed time series, or use `simulate_climate_forcing()` to create synthetic data for testing.

2. **Define Model Parameters:** Use `initialize_parameters()` to set model parameters such as which densification model is used, output frequency, and physics options.

3. **Initialize a Column:** Use `initialize_profile()` to create an initial profile of temperature, density, grid spacing, and other column properties.

4. **Run GEMB:** Pass the profile, climate forcing, and model parameters to the `gemb()` function.

## Examples

The `examples/` directory contains working examples:

* `examples/synthetic_example.jl` - Simple example using synthetic climate forcing
* `examples/era5_example.jl` - Example using ERA5 reanalysis data
* `examples/compare_outputs.jl` - Cross-validation between MATLAB and Julia implementations

## Key Functions

### Initialization
- `initialize_parameters()` - Create model parameters with defaults or overrides
- `initialize_forcing()` - Load/create climate forcing time series
- `initialize_profile()` - Set up initial vertical profile

### Main Model
- `gemb()` - Main driver function that runs the model
- `gemb_core()` - Single timestep integration
- `gemb_spinup()` - Cycle forcing to reach equilibrium

### Utilities
- `gemb_profile()` - Extract vertical profiles at specific times
- `gemb_interp()` - Interpolate profiles to specific depths
- `forcing_climatology()` - Create climatology from time series
- `surface_timeseries()` - Extract surface values from output arrays
- `dz2z()` - Convert grid spacing to center coordinates

### Climate Data Utilities
- `simulate_climate_forcing()` - Generate synthetic forcing data
- `vapor_pressure_to_relative_humidity()` - Convert vapor pressure to RH
- `relative_humidity_to_vapor_pressure()` - Convert RH to vapor pressure
- `dewpoint_to_vapor_pressure()` - Convert dewpoint to vapor pressure

### Climate Fitting Functions
- `fit_air_temperature()` - Fit temperature model to data
- `fit_precipitation()` - Fit precipitation model to data
- `fit_longwave_irradiance_delta()` - Fit longwave radiation delta
- `fit_seasonal_daily_noise()` - Fit seasonal noise patterns

## Documentation

Full documentation is available at:
- **Stable:** https://alex-s-gardner.github.io/GEMB.jl/stable/
- **Development:** https://alex-s-gardner.github.io/GEMB.jl/dev/

See also the comprehensive MATLAB documentation at the [original GEMB repository](https://github.com/alex-s-gardner/GEMB).

## Output Structure

GEMB returns a `DimStack` containing labeled arrays with explicit dimensions:
- `Ti` - Time dimension
- `Z` - Vertical (depth) dimension

Output variables include:
- `temperature` - Subsurface temperature (K)
- `density` - Snow/firn/ice density (kg/m³)
- `grain_radius` - Grain size (m)
- `water_content` - Liquid water content (-)
- `dz` - Grid spacing (m)
- Surface fluxes and mass balance components

## Testing

GEMB.jl includes a comprehensive test suite that validates against the original MATLAB implementation:

```julia
using Pkg
Pkg.test("GEMB")
```

Tests validate:
- Individual physics modules against MATLAB reference data (typically ~1e-12 relative tolerance)
- Full model integration
- Spinup functionality
- Profile extraction and interpolation

## Performance

Julia's JIT compilation and type inference, combined with allocation and compute
optimizations to the hot physics loop, provide substantial performance benefits over the
reference MATLAB implementation.

**Benchmark workload** (`examples/synthetic_example.jl` / `examples/GEMB_example_synthetic.m`):
a 75-year climatological spinup followed by a 32-year transient run — **~107 model-years at
3-hourly resolution** — on a single firn column. Total wall-clock time for the full
workflow (spinup + transient run), measured on an Apple M2 Max as the minimum of 3 runs
after warmup:

| Implementation | Total runtime | Speedup vs MATLAB |
|---|---|---|
| MATLAB (R2024b) | 98.8 s | 1× |
| Julia (current) | **7.0 s** | **~14×** |

Both were re-benchmarked on the same machine with the same protocol
(`examples/GEMB_example_synthetic.m` timed physics only, post-warmup, min of 3 runs).

The speedup came in three stages. First, a hot-loop **allocation** rewrite eliminated ~40
temporary `Vector`/`BitVector` objects per timestep across the grain-size, density,
albedo, shortwave, and temperature physics functions, replacing mask-broadcasting and
gather/scatter patterns with scalar `@inbounds for` loops and caller-owned scratch buffers.
Second, a **compute** pass on the thermal solver (`calculate_temperature`, ~55% of each
timestep, driven by a ~36-iteration sub-stepping loop): the diffusion stencil was made
branch-free by peeling the boundary cells out of the loop, the shortwave-penetration update
was fused into the stencil reads (removing one full-column pass per sub-step), and the
sub-step-invariant turbulent-flux terms (bulk coefficient, Exner pressure factor, roughness
logs) were hoisted out of the loop.

Third, a **type-stability and residual-allocation** pass. The driver's time loop was
extracted into a function barrier: `ClimateForcing`'s fields are typed `::DimArray` (a
`UnionAll`, not a concrete type), so indexing them per timestep inferred to `Any` and
dispatched at runtime — passing the forcing series into the inner loop as concrete arrays
lets the compiler specialize, eliminating the per-step dynamic dispatch (isolated driver
call: −11% runtime, −41% allocations). The remaining broadcast temporaries in
`calculate_melt` (pore-water refreeze, water-squeeze, melt/refreeze passes) and the
split-detection `BitVector` in `manage_layers` were also converted to scalar loops. This
cut total hot-path allocations by ~30% (27.1M → 19.1M objects).

Every optimization stage is numerically bit-identical to the prior implementation (max
relative difference 2.9e-16 across all output fields; full MATLAB-validation test suite
green).

The first call to any function includes JIT compilation overhead; subsequent calls use
compiled native code.

## Differences from MATLAB Version

GEMB.jl maintains high fidelity to the MATLAB implementation while embracing Julia idioms:

- **DimensionalData.jl:** All arrays use explicit dimension labels (`Ti`, `Z`) instead of implicit ordering
- **Immutable Parameters:** `ModelParameters` is immutable; create new instance for modifications
- **Multiple Dispatch:** Physics functions use multiple dispatch for type-based specialization
- **Broadcasting:** Uses Julia's `.` broadcasting syntax instead of implicit array operations
- **Module System:** Organized as a proper Julia package with explicit exports

## Prerequisites

GEMB.jl requires:
- Julia ≥ 1.9
- DimensionalData.jl
- Statistics (standard library)
- Dates (standard library)

For MATLAB cross-validation tests:
- MATLAB.jl (optional, for running MATLAB comparison tests)
- MAT.jl (for reading MATLAB .mat files)

## Citation

If you use GEMB.jl in your research, please cite:

Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass Balance (GEMB): a model of firn processes for cryosphere research, Geosci. Model Dev., 16, 2277–2302, [https://doi.org/10.5194/gmd-16-2277-2023](https://doi.org/10.5194/gmd-16-2277-2023), 2023.

If you use GEMB model outputs, please cite:

Schlegel, N.-J., & Gardner, A. (2025). Output from the Glacier Energy and Mass Balance (GEMB v1.0) forced with 3-hourly ERA5 fields and gridded to 10km, Greenland and Antarctica 1979-2024 (1.4) [Data set]. Zenodo. [https://doi.org/10.5281/zenodo.14714746](https://doi.org/10.5281/zenodo.14714746)

## Related Repositories

- [GEMB (MATLAB)](https://github.com/alex-s-gardner/GEMB) - Original MATLAB implementation
- [GEMB.jl Documentation](https://alex-s-gardner.github.io/GEMB.jl/)

## Contributing

Contributions are welcome! Please:
1. Maintain fidelity to the MATLAB implementation for physics functions
2. Add tests that validate against MATLAB reference data
3. Follow Julia style guidelines
4. Update documentation for new features

## License

See LICENSE file for details.

## Authors

GEMB was created by Alex Gardner, with contributions from Nicole-Jeanne Schlegel and Chad Greene. The Julia translation was developed by Alex Gardner.

For questions or issues, please open an issue on GitHub.
