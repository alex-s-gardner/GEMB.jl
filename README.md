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

GEMB.jl is not yet registered in the Julia General registry. Install directly from GitHub:

```julia
using Pkg
Pkg.add(url="https://github.com/alex-s-gardner/GEMB.jl")
```

## Quick Start

```julia
using GEMB
using GEMB_ClimateForcing  # synthetic/real forcing lives in the companion package

# 1. Initialize model parameters
params = initialize_parameters()

# 2. Create climate forcing (a DimStack) and convert to a ClimateForcing
ds = simulate_climate_forcing("test_1", 3)  # 3-hourly synthetic data
forcing = initialize_forcing(ds)

# 3. Initialize the vertical profile
profile = initialize_profile(params, forcing)

# 4. Run the model
output = gemb(profile, forcing, params)

# 5. Extract surface temperature time series
T_surface = surface_timeseries(output.temperature)
```

`simulate_climate_forcing` (and real-data loaders) live in the companion
[GEMB_ClimateForcing.jl](https://github.com/alex-s-gardner/GEMB_ClimateForcing.jl),
which produces a `DimStack`. GEMB's core `initialize_forcing(::DimStack)` method
converts it into the `ClimateForcing` the model consumes.

## Basic Workflow

Using GEMB requires four basic steps:

1. **Define Climate Forcing:** Use `initialize_forcing()` to create a `ClimateForcing` — either from observed time-series vectors, or from a `DimStack` produced by the companion `GEMB_ClimateForcing.jl` (e.g. `simulate_climate_forcing()` for synthetic test data, or the ERA5-Land loader for real data).

2. **Define Model Parameters:** Use `initialize_parameters()` to set model parameters such as which densification model is used, output frequency, and physics options.

3. **Initialize a Column:** Use `initialize_profile()` to create an initial profile of temperature, density, grid spacing, and other column properties.

4. **Run GEMB:** Pass the profile, climate forcing, and model parameters to the `gemb()` function.

## Examples

The `examples/` directory contains working examples:

* `examples/synthetic_example.jl` - Simple example using synthetic climate forcing
* `examples/era5_example.jl` - Example using ERA5 reanalysis data

## Key Functions

### Initialization
- `initialize_parameters()` - Create model parameters with defaults or overrides
- `initialize_forcing()` - Create a `ClimateForcing`, either from time-series vectors or from a forcing `DimStack`
- `initialize_profile()` - Set up initial vertical profile

### Main Model
- `gemb()` - Main driver function that runs the model
- `gemb_core()` - Single timestep integration
- `gemb_spinup()` - Cycle forcing to reach equilibrium

### Utilities
- `gemb_profile()` - Extract vertical profiles at specific times
- `gemb_interp()` - Interpolate profiles to specific depths
- `forcing_climatology()` - Average a `ClimateForcing` into a one-year climatological cycle (for spinup)
- `surface_timeseries()` - Extract surface values from output arrays
- `dz2z()` - Convert grid spacing to center coordinates

### Climate Data (companion package)

Climate-forcing generation and loading live in the companion
[GEMB_ClimateForcing.jl](https://github.com/alex-s-gardner/GEMB_ClimateForcing.jl),
which produces a `DimStack` for `initialize_forcing`. These are **not** exported
by GEMB.jl:

- `simulate_climate_forcing()` - Generate synthetic forcing data
- `climate_forcing()` - Download ERA5-Land reanalysis data
- `dewpoint_to_vapor_pressure()`, `vapor_pressure_to_relative_humidity()`, `relative_humidity_to_vapor_pressure()` - Humidity conversions
- `fit_air_temperature()`, `fit_precipitation()`, `fit_longwave_irradiance_delta()`, `fit_seasonal_daily_noise()` - Climate fitting functions

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

The Julia hot loop has been tuned for low allocations and type stability while remaining
numerically faithful to the MATLAB reference: outputs are bit-identical to the original
vectorized implementation (max relative difference ~3e-16 across all output fields), and the
full MATLAB-validation test suite passes.

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
- Julia ≥ 1.11
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
