# GEMB.jl: Glacier Energy and Mass Balance Model

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://alex-s-gardner.github.io/GEMB.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://alex-s-gardner.github.io/GEMB.jl/dev/)
[![Build Status](https://github.com/alex-s-gardner/GEMB.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/alex-s-gardner/GEMB.jl/actions/workflows/CI.yml?query=branch%3Amain)

## Overview

GEMB.jl is a Julia implementation of the Glacier Energy and Mass Balance (GEMB, the "B" is silent) model - a comprehensive one-dimensional physical model designed to simulate the surface energy balance and vertical firn evolution of glaciers and ice sheets. It couples atmospheric forcing with subsurface thermodynamics and densification physics to resolve the evolution of temperature, density, water content, and grain properties over time.

GEMB is a column model (no horizontal communication) of intermediate complexity, prioritizing computational efficiency to accommodate the multi-millennial spin-ups required for initializing deep firn columns. It is used for interpreting satellite altimetry data, firn studies, surface mass balance inversion from satellite data, ice core studies, uncertainty quantification and model exploration in cryosphere research. A complete description of GEMB can be found in [*Gardner et al*., 2023](https://doi.org/10.5194/gmd-16-2277-2023).

This Julia implementation began as a translation of the original MATLAB version and leverages Julia's performance, type system, and modern ecosystem including DimensionalData.jl for labeled arrays. It now deviates from the MATLAB version where the physics warranted it; those deviations are listed under [Differences from MATLAB Version](#differences-from-matlab-version).

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

# 3. Initialize the vertical profile (returns possibly depth-adjusted params)
profile, params = initialize_profile(params, forcing)

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
| Julia (current) | **5.7 s** | **~17×** |

Both were re-benchmarked on the same machine with the same protocol
(`examples/GEMB_example_synthetic.m` timed physics only, post-warmup, min of 3 runs).

The Julia hot loop has been tuned for low allocations and type stability. The
MATLAB-validation test suite passes (~1e-12 relative tolerance per module) for every
physics function except where a deviation is documented under
[Physics deviations](#physics-deviations).

Whole-column output is **not** bit-identical to the MATLAB reference, and cannot be: the
fixed-length vertical grid (see below) merges and splits cells at depth to hold the cell
count constant, which changes the deep discretization. Integrated surface quantities agree
closely with the pre-refactor Julia results — melt −0.17%, runoff +0.26%, refreeze −0.50%,
net shortwave −0.31%, surface albedo +0.01%, and latent heat flux +0.79% over the 107
model-year benchmark — but a cell-by-cell comparison against MATLAB will differ at depth by
more than floating-point rounding.

The first call to any function includes JIT compilation overhead; subsequent calls use
compiled native code.

## Differences from MATLAB Version

GEMB.jl embraces Julia idioms and deviates from the MATLAB implementation where justified:

- **DimensionalData.jl:** All arrays use explicit dimension labels (`Ti`, `Z`) instead of implicit ordering
- **Immutable Parameters:** `ModelParameters` is immutable; create new instance for modifications
- **Multiple Dispatch:** Physics functions use multiple dispatch for type-based specialization
- **Broadcasting:** Uses Julia's `.` broadcasting syntax instead of implicit array operations
- **Module System:** Organized as a proper Julia package with explicit exports
- **Fixed-length vertical grid:** the column holds a constant cell count and a constant total
  depth for the whole run, rather than growing and shrinking as cells are added and removed
  (see below)

### Physics deviations

Changes that alter model output relative to the MATLAB reference. Upstream issue numbers
refer to [the MATLAB repository](https://github.com/alex-s-gardner/GEMB/issues).

- **`firn_air_content` is metres of air**, `Σ dz (1 − ρ/ρ_ice)`. MATLAB normalizes the same
  mass deficit by 1000 rather than `density_ice`, which yields metres of water equivalent —
  9.0% lower for `density_ice = 910`. Affects the `firn_air_content` output only
  (upstream [#198](https://github.com/alex-s-gardner/GEMB/issues/198)).
- **`:ArthernB` densification is corrected and selectable.** MATLAB omits gravity and the
  per-second-to-per-year conversion from Arthern et al. (2010) eq. B1 — a combined factor of
  ~3.1e8, leaving the scheme inert — and applies the overlying cell's density to the whole
  overlying depth rather than integrating `Σ(ρⱼdzⱼ + waterⱼ)g`. Corrected here and validated
  against the Community Firn Model's `Arthern2010T` to 1e-8. MATLAB marks the scheme
  "DO NOT USE"; it is now a permitted `densification_method`
  (upstream [#200](https://github.com/alex-s-gardner/GEMB/issues/200)). The steady-state
  initial guess falls back to `:Arthern`, which the spinup then relaxes. `:LiZwally` and
  `:Helsen` remain gated.
- **Added `densification_method = :GSFC2020`** — Medley et al. (2022) GSFC-FDM v1.2.1
  eq. 18, the recalibrated successor to `:Arthern`: the same `dρ/dt = c(ρᵢ−ρ)` form with the
  mean accumulation raised to a fitted exponent (α₀ = 0.91, α₁ = 0.644) and a per-stage
  activation energy (59500 and 56870 J mol⁻¹ against Arthern's single 60000). Not present in
  MATLAB. Cross-checked against the Community Firn Model's `GSFC2020` to the 0.102% `g`
  offset. Because α ≠ 1 the accumulation units are load-bearing rather than absorbed into a
  prefactor: `precipitation_mean` is kg m⁻² yr⁻¹, which is the convention the exponents were
  fitted against.
- **Added `densification_method = :Simonsen2013`** — Simonsen et al. (2013), `:Arthern`'s form
  and activation energies retuned for Greenland: a constant factor 0.8 below 550 kg m⁻³, and
  1.25·γ above with `γ = 61.7/√A · exp(−3800/RT̄)`, where γ scales the second stage only. Not
  present in MATLAB. The form follows Lundin et al. (2017) FirnMICE eqs. A36–A37; that paper
  publishes no numeric tuning scalars, so the values 0.8 and 1.25 are the Community Firn
  Model's, against which the implementation is cross-checked to the 0.102% `g` offset. As for
  `:GSFC2020`, the accumulation units are load-bearing (`γ ∝ A^−1/2`): `precipitation_mean` is
  kg m⁻² yr⁻¹.
- **Added `densification_method = :Crocus` and `:CrocusPure`** — viscous settling of Vionnet
  et al. (2012) eqs. 5–9, `dD/D = −σ/η·dt` with `η = f1·f2·η₀(ρ/cη)·exp(aη(T_fus−T) + bη·ρ)`.
  Not present in MATLAB, and the only scheme in which liquid water affects densification:
  `f1 = (1 + 60·W_liq/(ρ_w·D))⁻¹` reduces viscosity by up to ~61× in wet cells, `f2 =
  min(4, exp(min(0.4 mm, gs−0.2 mm)/0.1 mm))` raises it for coarse angular grains. Overburden
  is the same `Σ(ρdz + water)·g` integral `:ArthernB` uses, taken at the cell midpoint (the
  paper's half-own-weight rule for the surface layer, applied uniformly). Integrated as
  `D·exp(−σ/η·dt)`, eq. 5's exact solution at constant σ and η, rather than through the
  `c(ρᵢ−ρ)` relaxation the other schemes share. Cross-checked against the Community Firn
  Model's `Crocus`, which differs in hardcoding `f2 = 4` and adopting van Kampenhout et al.
  (2017)'s retuned `cη = 358`; the paper's `f2` and `cη = 250` are used. Two departures from
  the published law: `f2` is floored at 1 (eq. 9 is unbounded below and would soften GEMB's
  fresh snow 7×, a regime the law is not defined on — Crocus applies eq. 9 only to
  non-dendritic snow with `gs ≥ 0.3 mm`), and `:Crocus` hands cells at or above 450 kg m⁻³
  to `:GSFC2020`. That handover exists because eq. 7 is fitted to a 1–2 m alpine snowpack and
  `exp(bη·ρ)` saturates in firn: the unblended law gives 0.7–0.8× `:Arthern`'s compaction
  rate in the top few metres but only 0.02–0.09× below 20 m. 450 kg m⁻³ is the threshold the
  Community Firn Model uses for the same purpose; the handover is a branch, not a weighted
  blend, since neither paper prescribes a blending function. `:CrocusPure` applies eq. 5 at
  every density.
- **`water_irreducible_saturation` now applies during percolation.** MATLAB declares a local
  `0.07` that overrides the documented parameter at every percolation retention site; only
  the pre-percolation squeeze read the parameter. Changes output wherever
  `water_irreducible_saturation != 0.07` (upstream
  [#199](https://github.com/alex-s-gardner/GEMB/issues/199)).
- **Added `water_irreducible_method = :ColeouLesaffre`** — density-dependent irreducible
  saturation, `S_wi = wmi/(1−wmi)·ρᵢρ/(ρ_w(ρᵢ−ρ))` with `wmi = 0.057(ρᵢ−ρ)/ρ + 0.017`
  (Coléou and Lesaffre 1998 eq. 3 via Langen et al. 2017 eq. 4), matching the Community Firn
  Model. Retention is zero at and above `DENSITY_PORE_CLOSEOFF`, as in CFM. Default
  `:constant` is unchanged and bit-identical.
- **Added `heat_capacity_method = :CuffeyPaterson`** — temperature-dependent specific heat,
  `c_p(T) = 152.5 + 7.122 T` (Cuffey and Paterson 2010 eq. 9.1). MATLAB uses the constant
  2102 J kg⁻¹ K⁻¹, its value at the melting point, which overstates the cold content of firn
  at 240 K by 12.9% and at 210 K by 27.5%. Default `:constant` reproduces the MATLAB value.
- **Internal energy is the enthalpy integral `∫c_p dT`, not `M·T·c_p`.** The latter is valid
  only for constant `c_p`; under `:CuffeyPaterson` it overstates enthalpy by `(b/2)T²`, which
  is 0.79·`LF` at the melting point. Cell mixing, the thermal solver, and every energy budget
  now work in enthalpy. For the default `:constant` the two forms are algebraically identical,
  so the change is arithmetic reordering (~1e-13 per operation); the 75-cycle spinup is not
  converged and amplifies this to ~0.3% in melt-related whole-run totals, so the synthetic
  regression pins were re-centered. Absolute column enthalpy is not comparable between the two
  methods (574 061 vs 307 385 J kg⁻¹ at 273.15 K), so reported thermal energy shifts when
  switching. The enthalpy solver costs ~17% more wall-clock than the constant-`c_p` form it
  replaced (4.9 s to 5.7 s on the 107-model-year benchmark); allocations fell.

### Fixed-length vertical grid

MATLAB GEMB lets the column length float: cells are pushed and popped as snow accumulates,
melts out, splits, and merges. GEMB.jl pins both the cell count and the total column depth via
two independent controllers — count is restored by exactly conservative merge/split operations
at depth, and depth by a continuous signed adjustment to the bottom cell, the model's only
basal mass and energy flux. Mass and energy are conserved throughout (checked every timestep
under `verbose=true`; the whole-run budget closes to ~1e-11 kg m⁻²).

What this changes for consumers of the output:

- **Profile arrays are sized exactly to the column** and are **top-justified** — the surface
  is always row 1. They were previously overpadded (default 1000 extra rows) and
  bottom-justified, so the surface row moved between timesteps and validity had to be
  recovered by scanning for `NaN`. Neither is needed now; `valid_profile_length` is retained
  and constant.
- **`output_padding` and `column_zmin` are removed**, and `column_zmax` is renamed
  **`column_depth`** — padding is obsolete, and depth is now pinned exactly rather than
  banded within `[column_zmin, column_zmax]`.
- **`thickness_cumulative` carries a real basal flux**, signed: mass leaves through the base
  under accumulation and enters under ablation. Read it as basal flux, not as glacier
  thickness change — the column is an Eulerian window on the firn, not a prognostic
  ice-thickness model. For thickness change, use the surface terms
  (`precipitation - runoff + evaporation_condensation`).

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
1. Deviate from the MATLAB implementation only where physically justified, and document the
   deviation under [Physics deviations](#physics-deviations)
2. Add tests that validate against MATLAB reference data
3. Follow Julia style guidelines
4. Update documentation for new features

## License

See LICENSE file for details.

## Authors

GEMB was created by Alex Gardner, with contributions from Nicole-Jeanne Schlegel and Chad Greene. The Julia translation was developed by Alex Gardner.

For questions or issues, please open an issue on GitHub.
