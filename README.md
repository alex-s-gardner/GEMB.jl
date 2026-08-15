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

# 3. Initialize the vertical profile (its grid is sized to the forcing climate)
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

- **Grain metamorphism runs whenever any scheme reads grain size**, not only when
  `albedo_method` is `:GardnerSharp` or `:BrunLefebre`. That gate is a performance
  optimization — grains coarsen regardless of the albedo scheme — but it does not cover the
  other consumers of `grain_radius`: `:ArthernB`/`:Crocus`/`:CrocusPure` densification and the
  `:grain_radius_threshold`/`:grain_radius_w_threshold` emissivity methods. Under
  MATLAB's gate those configurations run with `grain_radius` frozen at its initial profile,
  silently. Since `:ArthernB` compaction goes as `1/r²` and grains only coarsen, the frozen
  radius is always too small and compaction always too fast: over 32 years of synthetic
  forcing, final mean column density is 671.5 rather than 877.7 kg m⁻³. `grain_size_required`
  enumerates the consumers in one place. Defaults are unaffected — grain growth already ran
  under `:GardnerSharp` (upstream [#202](https://github.com/alex-s-gardner/GEMB/issues/202)).
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
  the published law: `f2` is floored at 1 (eq. 9 is bounded above but not below, so it
  *softens* fine-grained snow — 0.368, a 2.7× speedup, at GEMB's 0.05 mm fresh-snow radius.
  Crocus never reaches that regime because eq. 9 applies only to non-dendritic snow, whose
  `gs` the paper puts at 0.3–0.4 mm where eq. 9 gives 2.7–4; GEMB carries one grain radius
  for both regimes, so the floor keeps `f2` a stiffening correction across the domain eq. 9
  actually covers and inert below it), and `:Crocus` hands cells at or above 450 kg m⁻³
  to `:GSFC2020`. That handover exists because eq. 7 is fitted to a 1–2 m alpine snowpack and
  `exp(bη·ρ)` saturates in firn: the unblended law gives 0.7–0.8× `:Arthern`'s compaction
  rate in the top few metres but only 0.02–0.09× below 20 m. 450 kg m⁻³ is the threshold the
  Community Firn Model uses for the same purpose; the handover is a branch, not a weighted
  blend, since neither paper prescribes a blending function. `:CrocusPure` applies eq. 5 at
  every density.
- **Added `densification_method = :Barnola1991`** — Herron & Langway (1980) stage 1 below
  550 kg m⁻³, Barnola et al. (1991) pressure sintering above:
  `dρ/dt = ρ·A₀·exp(−Q/RT)·f(ρ)·σ³` (Tellus 43B, 83–90, eq. 2) with `A₀ = 2.54e4 MPa⁻³ s⁻¹`
  converted to Pa⁻³ s⁻¹, `Q = 60 kJ mol⁻¹`.
  Not present in MATLAB. The only scheme here that treats the firn–ice transition
  mechanistically rather than extrapolating a snow/firn fit: `f(ρ)` is a polynomial fitted to
  Pimienta and Duval (1987) below 800 kg m⁻³ and the analytic isolated-spherical-pore form
  `(3/16)(1−ρ/ρᵢ)/(1−(1−ρ/ρᵢ)^⅓)³` above, which vanishes with porosity, so the law
  self-limits as ρ → ρᵢ instead of being driven to zero by a `(ρᵢ−ρ)` factor. That
  self-limiting holds only while the closed-pore branch is reached at all: the 800 kg m⁻³
  handover is absolute while the polynomial carries no `ρᵢ`, so `:Barnola1991` requires
  `density_ice >= 900` (enforced in `validate_parameters`), against `[800, 950]` elsewhere.
  Overburden is
  the same `Σ(ρdz + water)·g` integral `:ArthernB` uses, at the cell midpoint as for
  `:Crocus`. Integrated with forward Euler on `dρ/dt`. Stage 1 is bit-for-bit
  `:HerronLangway`, sharing one kernel; the steady-state initial guess falls back to
  `:HerronLangway`, exact below 550 and approximate above. Cross-checked against the
  Community Firn Model's `Barnola1991`, which agrees on both branches and all four
  polynomial coefficients. `n = 3` throughout, as in the paper, whose `A₀` was fitted against
  it over 0.55–0.8 g cm⁻³; the paper's `n = 1` remark scopes itself to bulk ice below
  close-off (Pimienta and Duval 1987 torsion-tested 0.85 g cm⁻³ ice and report no
  densification rate), so it is not an unimplemented firn branch. Two departures, both
  documented in `calculate_density.jl`: (i) the paper's ΔP is overburden *minus bubble
  pressure*; only the overburden is applied, so effective stress is overstated above
  830 kg m⁻³ where pores close. This paper gives no expression for the bubble term; Goujon et
  al. (2003) eqs. A11–A12 do, so it belongs in a separate scheme rather than as a patch here.
  (ii) only the closed-pore branch scales with `density_ice`, and the branches meet in value
  at ρᵢ = 919.96 and in slope at 920.06 against the paper's stated C¹ matching, so the fit
  assumes ρᵢ ≈ 920 and at GEMB's default 910 the rate steps down 14% crossing 800 kg m⁻³ (a
  discontinuity in `dρ/dt`, not in ρ). The 550 kg m⁻³ handover is also a rate discontinuity, by
  construction in the paper: it joins a depth-independent accumulation-driven rate to a σ³ one,
  so on a 250 K, 300 kg m⁻² yr⁻¹ column sintering is ~1.3e4× slower than Herron–Langway just
  above 550 at 1 m depth and overtakes it only near 13 m. Calibrated over −14 to −57 °C and
  2.2–65 g cm⁻² yr⁻¹; no liquid-water term.
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
- **Added `impermeable_density` and `impermeable_thickness`**, exposing the bucket scheme's
  density-based impermeability criterion, which MATLAB hardcodes as `830` kg m⁻³ and `0.1` m in
  `calculate_melt`. Water is routed to runoff at a contiguous run of cells at or above
  `impermeable_density` thicker than `impermeable_thickness`. RetMIP (Vandecrux et al. 2020)
  recommends bucket schemes adopt such a criterion and shows model spread at ice-slab sites is
  dominated by this pair: the density threshold spans 810 kg m⁻³ (DMIHH, after Gregory et al.
  2014, which gave the lowest firn-temperature RMSE at KAN_U) to 917 (DTU, whose runoff was
  unrealistic enough to be excluded from the paper's multi-model mean). The `:ColeouLesaffre`
  gate stays on the constant `DENSITY_PORE_CLOSEOFF` — where capillary retention ceases for
  want of connected pore space is a different question from where a lens stops conducting flow
  — so lowering the flow criterion does not silently change retention. Defaults reproduce
  MATLAB bit-identically.
- **Added percolation-depth, ice-slab, and depth-limited firn-air-content diagnostics**
  (`percolation_depth`, `ice_slab_thickness`, `ice_slab_depth`, `firn_air_content_10m`,
  `firn_air_content_20m`), the quantities RetMIP evaluates models against — upward-looking
  radar wetting-front depth (Heilig et al. 2018), firn-core slab observations, and published
  FAC (Vandecrux et al. 2019), which is reported to a fixed depth rather than whole-column.
  `percolation_depth` is an interval maximum, matching what radar measures; the slab terms are
  instantaneous, so `ice_slab_depth = NaN` ("no slab") cannot poison an interval mean, and are
  scanned after the grid controllers run, so recomputing slab depth from the output `dz` and
  `density` reproduces them. All are read-only and take no part in the mass or energy budget,
  so output is unchanged.
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
- **Spinup convergence is judged on the mass-weighted mean density of the whole column**, not
  on a cubic-spline depth average over a separate `convergence_depth` (removed, along with its
  half-the-column default). The column depth is now fixed for the run, so the full column is
  the same domain at every cycle and the mean is exact without interpolation. Shifts the cycle
  at which a spinup exits — the synthetic example converges at 48 cycles rather than 63,
  changing whole-run totals by ~0.3% — but not the physics of any cycle. MATLAB has no
  equivalent criterion.
- **Added `convergence_drift_density`**, a trend-based spinup exit test: the least-squares slope
  of column-mean density against cycle over the trailing `drift_window` cycles (default 10),
  in kg m⁻³ per cycle. Where `convergence_delta_density` bounds the step between consecutive
  cycles — which a column creeping steadily at just under the tolerance passes while still
  densifying — this bounds the trend. When both are given both must hold. Not present in
  MATLAB; off by default, so it changes no existing run.
- **No initialized cell starts above the melt point.** `initialize_profile` clamps the
  temperature it fills to 273.15 K, with a warning. The climate-derived path was already
  clamped in `_steady_state_temperature`; this extends the ceiling to the two
  MATLAB-fidelity flags (`constant_temperature`, and both flags together), which filled
  `temperature_air_mean` verbatim as MATLAB's `model_initialize_profile` does. Affects only
  sites whose mean annual air temperature is above freezing, where the old column began as
  ice above its melting point holding no water — enthalpy the column has no state to carry.

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
  **`column_depth_max`** — padding is obsolete, and the column depth is now pinned exactly
  (to the depth `initialize_profile` derives from the climate, capped by `column_depth_max`)
  rather than banded within `[column_zmin, column_zmax]`.
- **`albedo_method = :Bougamont2005` is removed, and albedo is no longer carried state.** It
  was the only prognostic (time-decay) albedo scheme; the four remaining methods (`:None`,
  `:GardnerSharp`, `:BrunLefebre`, `:GreuellKonzelmann`) are all diagnostic functions of the
  current column, so `calculate_albedo` returns two scalars instead of overwriting per-cell
  vectors, and the parameters `albedo_wet_snow_t0`, `albedo_dry_snow_t0`, and `albedo_K` are
  gone. The per-layer `albedo` and `albedo_diffuse` profile outputs are removed, and the
  surface time series `albedo_surface` is renamed **`albedo_broadband`**; it is now recorded
  from the value used in that timestep's shortwave balance rather than from the post-
  accumulation column, so it differs slightly from the old series on snowfall steps. Output
  under the four remaining methods is otherwise unchanged.
- **`thickness_cumulative` carries a real basal flux**, signed: mass leaves through the base
  under accumulation and enters under ablation. Read it as basal flux, not as glacier
  thickness change — the column is an Eulerian window on the firn, not a prognostic
  ice-thickness model. For thickness change, use the surface terms
  (`precipitation - runoff + evaporation_condensation`).
- **Added `horizontal_strain_rate`** [yr⁻¹, default `0.0`], the trace of the horizontal
  strain-rate tensor `ε̇_xx + ε̇_yy`. By incompressibility it thins (positive, divergence) or
  thickens (negative, convergence) every cell at constant density by `exp(-D·dt)`; the mass
  it exports leaves laterally and is reported as the new `strain_thinning` output and in the
  mass budget, separately from the basal flux. MATLAB has no ice-dynamic term; at the default
  the two agree exactly.

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
