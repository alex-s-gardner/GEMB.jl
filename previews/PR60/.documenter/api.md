---
---

# API Reference {#API-Reference}



## Types {#Types}
<details class='jldocstring custom-block' open>
<summary><a id='GEMB.ModelParameters' href='#GEMB.ModelParameters'><span class="jlbinding">GEMB.ModelParameters</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
ModelParameters
```


All GEMB model configuration parameters with validation. Construct with keyword arguments; unspecified fields use defaults.

Defaults follow the recommendations of the two firn-model intercomparisons (Vandecrux et al., 2020, RetMIP; Lundin et al., 2017, FirnMICE) and, where those are silent, the shipped defaults of the Community Firn Model and IMAU-FDM. See the parameter table in `README.md`.

Parameterized on the thermal solver type so that `thermal_solver` is concretely typed and the scheme is resolved at compile time. `::ModelParameters` in a signature still matches every instance.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.ClimateForcing' href='#GEMB.ClimateForcing'><span class="jlbinding">GEMB.ClimateForcing</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
ClimateForcing <: DimensionalData.AbstractDimStack
```


Time-series surface meteorological forcing for GEMB, implemented as an `AbstractDimStack`. It therefore supports the full DimensionalData stack API — `keys`, `length`, `dims`, `map`, `layers`, and DimStack-style indexing such as date-range subsetting `cf[Ti(a .. b)]` and single-time selection `cf[Ti(At(t))]` — all returning a `ClimateForcing` (a sub-stack), never a scalar step.

**Layers (13 time-series `DimArray`s sharing a common `Ti` dimension)**

Required forcing: `temperature_air`, `pressure_air`, `precipitation`, `wind_speed`, `shortwave_downward`, `longwave_downward`, `vapor_pressure`. Time-varying model parameters (typically `Fill`-backed): `black_carbon_snow`, `black_carbon_ice`, `cloud_optical_thickness`, `solar_zenith_angle`, `shortwave_downward_diffuse`, `cloud_fraction`.

**Metadata (scalars, carried in the stack `metadata` as a `NamedTuple`)**

`time_step::Int` [s], plus `temperature_air_mean`, `wind_speed_mean`, `precipitation_mean`, `temperature_observation_height`, `wind_observation_height`, source provenance (`dataset`, `latitude`, `longitude`, `elevation_offset`) and climatology provenance (`climatology_window_start`, `climatology_window_stop`, `climatology_n_years`, `climatology_steps_per_year`; set by [`forcing_climatology`](/api#GEMB.forcing_climatology)).

Both layers and scalar metadata are reachable by name via property access (`cf.temperature_air` returns the layer `DimArray`; `cf.time_step` returns the concrete scalar), so existing `cf.<field>` code continues to work unchanged. Construct via the 19-argument positional constructor (13 `DimArray`s then the 6 scalars, in the order listed above) or, more commonly, via [`initialize_forcing`](/api#GEMB.initialize_forcing).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.ClimateForcingStep' href='#GEMB.ClimateForcingStep'><span class="jlbinding">GEMB.ClimateForcingStep</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
ClimateForcingStep
```


Single time-step forcing values extracted from ClimateForcing. Plain struct of scalars for efficient use in the physics loop.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.ClimateSummary' href='#GEMB.ClimateSummary'><span class="jlbinding">GEMB.ClimateSummary</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
ClimateSummary
```


Climatological summary of a [`ClimateForcing`](/api#GEMB.ClimateForcing), reduced to the handful of scalars needed to guess a steady-state firn column. Built by [`initialize_climate_summary`](/api#GEMB.initialize_climate_summary).

Every field is derived from a small, bounded number of passes over the forcing series (one for the sums and the harmonic, plus one per iteration of the albedo/melt fixed point — typically 1-6 in total), so building one is cheap relative to any model integration.

**Mass balance [kg m-2 yr-1]**
- `accumulation`: snowfall only — precipitation falling below `mp.rain_temperature_threshold`. Rain is excluded because it does not add mass to the column as snow.
  
- `rainfall`: the complement, precipitation falling as rain.
  
- `melt`: annual melt from a surface-energy-balance estimate (see [`initialize_climate_summary`](/api#GEMB.initialize_climate_summary)).
  
- `refreeze`: the part of `melt + rainfall` the column can refreeze, capped by one annual layer's cold content.
  
- `balance`: net annual mass balance, `accumulation + refreeze - melt`. **Its sign is the only regime discriminator**: positive buries snow and grows firn, non-positive exposes ice. Rainfall is not subtracted — see the note in [`initialize_climate_summary`](/api#GEMB.initialize_climate_summary).
  
- `accumulation_effective`: `accumulation + refreeze`, the mass flux that drives densification.
  

**Temperature [K]**
- `temperature_air_mean`: mean annual air temperature.
  
- `temperature_air_effective`: the forcing's Arrhenius-weighted mean, carried through from `cf.temperature_air_effective`. Used _only_ where a rate factor is evaluated (`DensificationCoeffs`, under `mp.mean_temperature_method = :arrhenius`), never where a temperature is a temperature — the marched thermal profile and the annual harmonic are built about `temperature_air_mean`, which is the column's actual mean state.
  
- `temperature_amplitude`, `temperature_phase`: amplitude [K] and phase [rad] of the fitted annual cycle, `T(t) = T̄ + A·cos(ωt − φ)`.
  
- `latent_warming`: `LF·refreeze/(c_p·accumulation_effective)`, the mean warming from refreezing latent heat.
  

**Other means**
- `wind_speed_mean` [m s-1], `pressure_air_mean` [Pa], and the observation heights, carried through for the physics the marcher reuses.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.AbstractThermalSolver' href='#GEMB.AbstractThermalSolver'><span class="jlbinding">GEMB.AbstractThermalSolver</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
AbstractThermalSolver
```


Supertype for the schemes that advance the subsurface temperature profile over one forcing timestep. The concrete subtype held by `ModelParameters.thermal_solver` selects the scheme by dispatch, so a new scheme is a new subtype plus one `GEMB._thermal_solve!` method — there is no central branch to extend.

The interface a subtype must implement is

```julia
_thermal_solve!(solver, temperature, dz, density, K, shortwave_flux,
                water_surface, grain_radius, sfc, cfs, mp, verbose)
    -> (longwave_upward, heat_flux_sensible, heat_flux_latent,
        heat_flux_basal, evaporation_condensation)
```


updating `temperature` in place and returning forcing-step averages (the first four in W m-2, the last in kg m-2 accumulated over the step). Two invariants bind every implementation: the bottom cell is a Dirichlet reservoir whose temperature must be returned bit-unchanged, and the returned fluxes must be the ones actually applied to the column, since `gemb_core` closes its energy and mass budgets against them.

Subtypes are singletons today. They are structs rather than symbols so that a scheme can later carry its own preallocated workspace as fields without touching any dispatch site.

See [`ExplicitThermal`](/api#GEMB.ExplicitThermal) and [`ImplicitThermal`](/api#GEMB.ImplicitThermal).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.ExplicitThermal' href='#GEMB.ExplicitThermal'><span class="jlbinding">GEMB.ExplicitThermal</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
ExplicitThermal()
```


Explicit finite-volume thermal scheme, sub-stepped to the Von Neumann stability limit (see `GEMB._max_safe_dt`). The default. Conserves energy to the last bit — diffusion is applied as one flux per face with opposite signs to the two adjacent cells — at the cost of a sub-step count set by the single stiffest cell in the column.

**References**
- Patankar, S. V. (1980). _Numerical Heat Transfer and Fluid Flow_, Ch. 3-4.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.ImplicitThermal' href='#GEMB.ImplicitThermal'><span class="jlbinding">GEMB.ImplicitThermal</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
ImplicitThermal()
```


Backward-Euler thermal scheme on a tridiagonal system, solved by the Thomas algorithm. Unconditionally stable, so it needs no stability sub-stepping and is insensitive to thin layers; `mp.dt_divisors` and `GEMB._max_safe_dt` are unused on this path.

**References**
- Patankar, S. V. (1980). _Numerical Heat Transfer and Fluid Flow_, Ch. 4.
  
- Thomas, L. H. (1949). _Elliptic Problems in Linear Difference Equations over a Network_.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>


## Initialization {#Initialization}
<details class='jldocstring custom-block' open>
<summary><a id='GEMB.initialize_parameters' href='#GEMB.initialize_parameters'><span class="jlbinding">GEMB.initialize_parameters</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
initialize_parameters(; kwargs...)
```


Create and validate a `ModelParameters` struct. Matches MATLAB's `model_initialize_parameters.m`.

All parameters have defaults matching the MATLAB version. Validation checks are performed after construction.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.initialize_forcing' href='#GEMB.initialize_forcing'><span class="jlbinding">GEMB.initialize_forcing</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
initialize_forcing(time, temperature_air, pressure_air, precipitation,
                   wind_speed, shortwave_downward, longwave_downward,
                   vapor_pressure; kwargs...)
```


Create a `ClimateForcing` struct from time-series vectors. Matches MATLAB's `model_initialize_forcing.m`.

**Arguments**
- `time::Vector{DateTime}`: Time stamps
  
- `temperature_air::Vector{Float64}`: Air temperature [K]
  
- `pressure_air::Vector{Float64}`: Air pressure [Pa]
  
- `precipitation::Vector{Float64}`: Precipitation per timestep [kg m-2]
  
- `wind_speed::Vector{Float64}`: Wind speed [m s-1]
  
- `shortwave_downward::Vector{Float64}`: Incoming shortwave [W m-2]
  
- `longwave_downward::Vector{Float64}`: Incoming longwave [W m-2]
  
- `vapor_pressure::Vector{Float64}`: Vapor pressure [Pa]
  

**Keyword Arguments**
- `temperature_air_mean::Float64`: Climatological mean temperature [K]
  
- `wind_speed_mean::Float64`: Climatological mean wind speed [m s-1]
  
- `precipitation_mean::Float64`: Climatological mean precipitation [kg m-2 yr-1]
  
- `temperature_observation_height::Float64`: Height of temperature observation [m]
  
- `wind_observation_height::Float64`: Height of wind observation [m]
  
- `accumulation_mean::Float64`: Climatological mean _snowfall_ [kg m-2 yr-1], i.e. `precipitation_mean` less the rain fraction. Derived from the record by partitioning precipitation on `rain_temperature_threshold` when not supplied. Read by `calculate_density` under `densification_accumulation = :accumulation`.
  
- `temperature_air_effective::Float64`: Arrhenius-weighted mean air temperature [K], `Eg/(R·log(mean(exp(Eg/(R·T)))))`. Derived from the record when not supplied. Read by `calculate_density` under `mean_temperature_method = :arrhenius`.
  
- `rain_temperature_threshold::Float64`: Rain/snow partition temperature [K] used only to derive `accumulation_mean`; keep it equal to the `ModelParameters` field of the same name.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>



```julia
initialize_forcing(stack::DimStack) -> ClimateForcing
```


Create a `ClimateForcing` struct from a forcing `DimStack`.

A `DimStack` is GEMB's neutral, producer-agnostic forcing interface: any source (e.g. the companion `GEMB_ClimateForcing.jl` package, or a hand-built stack) can supply forcing as long as it carries the required layers, a `DateTime`-indexed `Ti` dimension, and the required metadata. This method extracts the forcing vectors and metadata from `stack` and forwards them to the vector method of [`initialize_forcing`](/api#GEMB.initialize_forcing), which applies the same validation and normalization.

**Required Variables in DimStack**
- `temperature_air`: Air temperature (K)
  
- `pressure_air`: Surface pressure (Pa)
  
- `vapor_pressure`: Vapor pressure (Pa)
  
- `wind_speed`: Wind speed (m/s)
  
- `precipitation`: Precipitation rate (kg/m²/timestep)
  
- `shortwave_downward`: Downward shortwave radiation (W/m²)
  
- `longwave_downward`: Downward longwave radiation (W/m²)
  

**Required Metadata**
- `temperature_air_mean`: Mean air temperature (K)
  
- `wind_speed_mean`: Mean wind speed (m/s)
  
- `precipitation_mean`: Annual precipitation (kg/m²/year)
  
- `temperature_observation_height`: Height of temperature observations (m)
  
- `wind_observation_height`: Height of wind observations (m)
  

**Optional Metadata (defaults used if not present)**
- `accumulation_mean` (kg/m²/year; derived from the record by rain/snow partitioning)
  
- `temperature_air_effective` (K; derived from the record as the Arrhenius-weighted mean)
  
- `black_carbon_snow` (default: 0.0)
  
- `black_carbon_ice` (default: 0.0)
  
- `cloud_optical_thickness` (default: 0.0)
  
- `solar_zenith_angle` (degrees, default: 0.0)
  
- `shortwave_downward_diffuse` (W/m², default: 0.0)
  
- `cloud_fraction` (default: 0.1)
  

**Examples**

```julia
using GEMB
using GEMB_ClimateForcing  # one producer of a conforming DimStack

ds = simulate_climate_forcing("test_1", 3)
cf = initialize_forcing(ds)

mp = initialize_parameters()
profile = initialize_profile(mp, cf)
output = gemb(profile, cf, mp)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.forcing_climatology' href='#GEMB.forcing_climatology'><span class="jlbinding">GEMB.forcing_climatology</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
forcing_climatology(cf::ClimateForcing)
forcing_climatology(cf::ClimateForcing, datetime_range::Tuple{DateTime,DateTime})
forcing_climatology(cf::ClimateForcing, selector)
```


Compute climatological average forcing from a [`ClimateForcing`](/api#GEMB.ClimateForcing).

Creates a single-year average forcing by:
2. Optionally subsetting to a time window (`datetime_range` or `selector`)
  
3. Eliminating leap days (day 366 of the year)
  
4. Excluding partial years
  
5. Averaging every time-series field across the complete years
  

Returns a new `ClimateForcing` holding one year of climatological forcing. The time step and scalar metadata (`temperature_air_mean`, `wind_speed_mean`, `precipitation_mean`, `temperature_observation_height`, `wind_observation_height`) are carried forward unchanged.

The window can be given as a `(start, stop)` `DateTime` tuple, or — since `ClimateForcing` is an `AbstractDimStack` — as any DimensionalData `Ti` selector, e.g. a closed interval `DateTime(1950,1,1) .. DateTime(1980,12,31)`, `At`, `Near`, or an index range. Both forms subset via `cf[Ti(...)]` before averaging.

Typically used to build a repeating forcing cycle for [`gemb_spinup`](/api#GEMB.gemb_spinup).

Matches MATLAB's `forcing_climatology.m`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.initialize_profile' href='#GEMB.initialize_profile'><span class="jlbinding">GEMB.initialize_profile</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
initialize_profile(mp::ModelParameters, cf::ClimateForcing;
                   constant_density=false, constant_temperature=false)
    -> profile::DimStack
```


Initialize a GEMB firn column profile as a DimStack, as a steady-state guess derived entirely from the climate forcing.

The grid is sized to the depth this climate needs (see [`_derive_column_depth`](/internals_grid#GEMB._derive_column_depth-Tuple{ModelParameters,%20ClimateSummary})) rather than to the configured `mp.column_depth_max`, which acts as a ceiling. Nothing else has to be told about that: the derived depth lives in the returned `dz`, and both [`gemb`](/api#GEMB.gemb) and [`gemb_spinup`](/api#GEMB.gemb_spinup) take the column depth they hold fixed from `sum(profile[:dz])`, not from `mp`.

**The scheme**

[`initialize_climate_summary`](/api#GEMB.initialize_climate_summary) reduces the forcing to a handful of scalars in one pass: snowfall and rainfall (partitioned on `mp.rain_temperature_threshold`, so rain is not counted as accumulating mass), a surface-energy-balance melt estimate, the cold-content-capped refreeze, the fitted annual temperature harmonic, and the net annual mass balance `b = snowfall + refreeze − melt`.

[`steady_state_profile`](/api#GEMB.steady_state_profile) then marches one parcel of snow forward in age with `b` as its burial rate, recording every state variable as it is buried: density (relaxing toward `mp.density_ice` under the run's own `mp.densification_method`), temperature (a damped annual wave about the latent-heat-warmed mean), grain size (evolved by [`calculate_grain_size`](/internals_physics#GEMB.calculate_grain_size-Tuple{Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20ClimateForcingStep,%20ModelParameters}) itself) and irreducible water. No albedo is stored: [`calculate_albedo`](/internals_physics#GEMB.calculate_albedo-Tuple{Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Float64,%20ClimateForcingStep,%20ModelParameters}) diagnoses it from the column at the start of every timestep, including the first.

**There is no regime threshold.** The sign of `b` is the only discriminator, and it needs none: where `b ≤ 0` nothing is ever buried, the march terminates immediately, and the column is the ice it exposes. Sites near `b = 0` get intermediate profiles rather than a cliff between a deep firn column and a block of solid ice.

**Keyword arguments**

Both flags exist **only to reproduce the MATLAB `model_initialize_profile` initialization** for the fidelity and regression tests. They are not physically meaningful choices for a real run — the climate-derived guess is strictly better — and either one builds the grid on the configured `mp.column_depth_max`, bypassing the derived depth.
- `constant_density`: fill density with `mp.density_ice`, with the matching bare-ice grain state (`grain_dendricity = 0`, `grain_sphericity = 0`, `grain_radius = 2.5` mm) and no pore water.
  
- `constant_temperature`: fill temperature with `cf.temperature_air_mean`, clamped to the melt point as every other path is.
  

Setting **both** takes a separate early-return path ([`_uniform_ice_profile`](/internals_grid#GEMB._uniform_ice_profile-Tuple{ModelParameters,%20Real})) that reproduces the MATLAB column exactly: pure ice, uniform mean-annual temperature, no march and no climate summary.

Every path here returns temperature at or below the melt point (273.15 K): the climate-derived march clamps in [`_steady_state_temperature`](/internals_grid#GEMB._steady_state_temperature-Tuple{Real,%20Real,%20ClimateSummary,%20ModelParameters}), and the fidelity paths clamp `cf.temperature_air_mean` (with a warning) before filling.

Returns a DimStack with Z dimension containing:
- dz, temperature, density, water, grain_radius, grain_dendricity, grain_sphericity, age
  

`age` [d] follows `mp.initialize_age`. Under `:steady_state` (the default) it is the residence time the march itself integrated — the accumulated burial time at each depth, increasing downward — which makes age and residence-time diagnostics meaningful from the first timestep instead of only after a multi-century spinup has flushed the column. Under `:zero` the instant of initialization is the epoch and every cell starts at 0, which is what the fidelity paths do regardless (`constant_density` discards the march's density, so its age describes a column that was replaced). Nothing in the physics reads `age`, so this choice is output-only.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.initialize_climate_summary' href='#GEMB.initialize_climate_summary'><span class="jlbinding">GEMB.initialize_climate_summary</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
initialize_climate_summary(cf::ClimateForcing, mp::ModelParameters;
                           max_albedo_passes=20) -> ClimateSummary
```


Reduce a climate forcing to the scalars that determine a steady-state firn column, in a few passes over the series (one to settle each of the coupled quantities; see the albedo note below).

Three quantities are computed simultaneously in the pass:

**Snowfall vs. rainfall.** Precipitation is partitioned on `mp.rain_temperature_threshold`, the same test [`calculate_accumulation`](/internals_physics#GEMB.calculate_accumulation-Tuple{Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20ClimateForcingStep,%20ModelParameters,%20Bool}) applies per timestep. This matters: the forcing's `precipitation_mean` is _total_ precipitation, so using it directly counts rain as accumulating mass.

**Melt, from a surface energy balance.** At each step the skin temperature that closes

```julia
SW(1−α) + LW↓ − εσT⁴ + QH + QE = 0
```


is found by Newton iteration, reusing [`turbulent_heat_flux`](/internals_physics#GEMB.turbulent_heat_flux-Tuple{Float64,%20Float64,%20Float64,%20Float64,%20Float64,%20ClimateForcingStep}) for the turbulent terms. Where that temperature would exceed the melt point, the surface is held at `273.15 K` and the residual energy becomes melt, `Δmelt = Q·Δt/LF`. This replaces a degree-day factor with the radiation, wind and humidity the forcing already carries.

**The annual temperature cycle.** Three scalar accumulators give the least-squares harmonic fit `T(t) = T̄ + A·cos(ωt − φ)`, which the marcher uses for the damped thermal wave.

Melt depends on albedo, and albedo depends on the surface state melt helps determine. That fixed point is _solved_, not relaxed: `g(α) = α_target(α) − α` is monotone in `α`, so a secant iteration bracketed by `mp.albedo_ice` and `mp.albedo_snow` reaches `ALBEDO_FIXED_POINT_TOLERANCE` in about five passes, and `max_albedo_passes` is a backstop rather than the normal exit. (Plain under-relaxed iteration needed 16–25 passes on the same sites — enough to exceed any bounded budget and return an unconverged albedo without saying so.) Failure to converge warns.

The albedo blend ([`_snow_cover_albedo`](/internals_grid#GEMB._snow_cover_albedo-Tuple{Real,%20Real,%20ModelParameters})) is continuous in the melt rate, so a site near `b = 0` does not flip between a snow and an ice albedo. Cost is a few passes over a vector, with no column integration.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.steady_state_profile' href='#GEMB.steady_state_profile'><span class="jlbinding">GEMB.steady_state_profile</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
steady_state_profile(dz, cs::ClimateSummary, mp::ModelParameters; n_age=2000)
    -> NamedTuple
```


Build a steady-state estimate of the full firn column state by marching a single parcel of snow forward in age and recording every state variable as it is buried.

This is the initial-guess engine behind [`initialize_profile`](/api#GEMB.initialize_profile). It replaces a family of regime-specific special cases with one continuous calculation: the burial rate is the net annual mass balance `cs.balance`, so an accumulation site buries snow and grows a firn column, while an ablation site buries nothing and yields the solid ice it exposes.

**There is no threshold to cross** — the sign of `cs.balance` is the only discriminator, and the profile varies continuously as a site approaches `balance = 0` from above (the surface density rises smoothly toward `mp.density_ice` as the burial rate falls). There _is_ a branch at `balance <= 0`, because the age step `dt = dz_step·ρ/b` is singular at `b = 0`; it returns the `b → 0⁺` limit of the march directly rather than stepping to it.

Each state variable is carried along the age coordinate using the same physics the model integrates transiently, so the guess is consistent with the run it seeds:
- **density** relaxes as `ρ(t+dt) = ρᵢ − (ρᵢ − ρ)·exp(−c·dt)` (Arthern et al. 2010 eq. 1) with `c` from [`_densification_rate`](/internals_physics#GEMB._densification_rate-Tuple{GEMB.DensificationCoeffs,%20Float64,%20Float64}) for the configured `mp.densification_method` — _not_ always Herron–Langway.
  
- **depth** follows from mass balance, `z(t) = b·∫₀ᵗ ds/ρ(s)`, where `b = cs.balance` [kg m-2 yr-1] is the mass buried per year.
  
- **temperature** is a damped annual wave about the latent-heat-warmed mean; see [`_steady_state_temperature`](/internals_grid#GEMB._steady_state_temperature-Tuple{Real,%20Real,%20ClimateSummary,%20ModelParameters}).
  
- **grain size** is evolved by [`calculate_grain_size`](/internals_physics#GEMB.calculate_grain_size-Tuple{Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20ClimateForcingStep,%20ModelParameters}) itself, so Marbouty/Brun metamorphism is not duplicated here.
  
- **water** is the irreducible content at the melt point, using the same formula as [`calculate_melt`](/internals_physics#GEMB.calculate_melt-Tuple{Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Float64,%20ModelParameters,%20Bool}).
  
- **age** is the march's own time coordinate — the accumulated `dt` at each depth, in days, which is the residence time a parcel buried at `cs.balance` would have. Below the firn/ice transition the march has stopped, so it is continued at the ice burial rate `ρᵢ/b`; it is zero everywhere for an ablation column, where nothing is buried at all.
  

Returns a `NamedTuple` of `Vector{Float64}`s at each depth in `z_center`, with keys `density`, `temperature`, `water`, `grain_radius`, `grain_dendricity`, `grain_sphericity`, `age`, plus the scalar `ice_depth` used to size the grid: the depth at which the column reaches `mp.density_ice` (to within `DENSITY_FIRN_TOLERANCE`), below which it carries no further firn information. `ice_depth` is `0.0` for an ablation column (ice at the surface) and `max_depth` if the march never reached ice.

**Arguments**
- `dz`: cell thicknesses [m], surface first, as [`initialize_grid`](/internals_grid#GEMB.initialize_grid-Tuple{ModelParameters}) returns. Cell centers are derived with [`dz2z`](/api#GEMB.dz2z); the thicknesses themselves are needed because water content is a per-cell mass, not a density.
  
- `cs`: the climate summary from [`initialize_climate_summary`](/api#GEMB.initialize_climate_summary).
  
- `mp`: model parameters; supplies the densification method, ice density and irreducible-water saturation.
  
- `n_age`: number of march steps spanning the column. The march is stepped in _depth_ (`max_depth/n_age` per step) with the age step derived from the local density, so resolution is uniform in depth and the step count does not depend on the burial rate.
  

**References**
- Arthern, R. J., et al. (2010). J. Geophys. Res., 115, F03011.
  
- Herron, M. and Langway, C. (1980). J. Glaciol., 25, 373-385.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.steady_state_density' href='#GEMB.steady_state_density'><span class="jlbinding">GEMB.steady_state_density</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
steady_state_density(z_center, T_mean, accumulation, ρ0, ρ_ice, mp; n_age=2000)
```


Steady-state firn density profile for `mp.densification_method`, as a standalone vector.

This is the density-only entry point to the same march [`steady_state_profile`](/api#GEMB.steady_state_profile) performs — both drive [`_march_density_step`](/internals_grid#GEMB._march_density_step-Tuple{Any,%20Vararg{Float64,%206}}), so there is one stepping scheme, not two — kept for direct use and for testing the density law in isolation. With `mp.densification_method = :HerronLangway` it reproduces the classic Herron & Langway (1980) analytic profile.

**Arguments**
- `z_center`: cell-center heights [m], negative below surface.
  
- `T_mean`: mean annual temperature [K].
  
- `accumulation`: mean annual accumulation [kg m-2 yr-1].
  
- `ρ0`: surface (fresh-snow) density [kg m-3].
  
- `ρ_ice`: density of ice [kg m-3], the asymptote and clamp.
  
- `mp`: supplies the densification method and its coefficients.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.fresh_snow_density' href='#GEMB.fresh_snow_density'><span class="jlbinding">GEMB.fresh_snow_density</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
fresh_snow_density(mp::ModelParameters, T_air_mean, precip_mean, wind_speed_mean,
                   T_air=T_air_mean, wind_speed=wind_speed_mean)
```


Return the density of fresh snow [kg m-3] for the configured `mp.new_snow_method`.

The first three forcing arguments are climatological means (temperature [K], accumulation [kg m-2 yr-1], wind speed [m s-1]), so this is shared by `calculate_accumulation` (per timestep) and `initialize_profile` (as the surface density ρ₀ of the steady-state firn profile).

`T_air` and `wind_speed` are the _instantaneous_ air temperature [K] and wind speed [m s-1], read only by the methods whose fits are meant to apply at the moment of snowfall rather than to a climatology: `T_air` by `:FaustoFit` and `:Pahaut`, `wind_speed` by `:Pahaut`. They default to their climatological counterparts so that the steady-state march — which has no instantaneous forcing to offer — gets a well-defined climatological ρ₀ from the same call, and so that callers of the older methods need not pass them.

Methods:
- `:Constant150` → 150, `:Constant315` → 315, `:Constant350` → 350: constants.
  
- `:Fausto` → 315: the Fausto et al. (2018) Greenland fit frozen at one temperature (it is `:FaustoFit` evaluated at T ≈ 256.2 K). Inherited from MATLAB and kept as the default for reference fidelity. Numerically equal to `:Constant315`, but `:Fausto` additionally selects the Crocus wind-dependent fresh-grain properties below, where `:Constant315` does not — use `:Constant315` for a bare 315 kg m-3 with the default grain properties.
  
- `:FaustoFit` → `362.1 + 2.78·(T_air - CtoK)`: that same fit carrying its temperature dependence, as implemented in IMAU-FDM (`initialise_model.f90`) for its Greenland domain. Unbounded below as published, so callers clamp — see `calculate_accumulation` and `steady_state_profile`.
  
- `:Pahaut` → `max(50, 109 + 6·(T_air - CtoK) + 26·√wind_speed)`: Pahaut (1975) as implemented in Crocus, via Lafaysse et al. (2026) eq. 35 (SURFEX/Crocus v3.0.2, their `:V12` fresh-snow option). The only method here fitted to _alpine seasonal snow_ rather than to a polar ice sheet, and the only one carrying an explicit wind dependence at the instantaneous timestep. Both dependencies are physical: warmer snowfall gives denser crystals, and wind fragments them before and during deposition. It is much lighter than the polar fits over their common range — at T = -10 °C and 5 m s-1 it gives 107 kg m-3 against `:FaustoFit`'s 334 — because alpine snowfall is warmer, wetter, and much less wind-packed than the katabatic-scoured surfaces the Greenland and Antarctic fits were regressed on. Prefer it for temperate and mid-latitude glaciers; the polar fits will overestimate fresh-snow density there, which suppresses the albedo of new snow and speeds its burial. The published floor of 50 kg m-3 is retained (it binds below about -10 °C in calm air).
  
- `:Kaspers`, `:KuipersMunneke`: temperature/accumulation/wind-dependent fits.
  

**References**
- Pahaut, E. (1975). _La métamorphose des cristaux de neige (Snow crystal metamorphosis)_. Monographies de la Météorologie Nationale 96, Météo-France.
  
- Lafaysse, M., et al. (2026). The SURFEX/Crocus v3.0.2 snowpack model. _Geosci. Model Dev._ 19, 6273-6320.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>


## Running the Model {#Running-the-Model}
<details class='jldocstring custom-block' open>
<summary><a id='GEMB.gemb' href='#GEMB.gemb'><span class="jlbinding">GEMB.gemb</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
gemb(profile::DimStack, climate_forcing::ClimateForcing, mp::ModelParameters; verbose=false)
```


Run the Glacier Energy and Mass Balance (GEMB) model. Matches MATLAB's `gemb.m`.

Returns a DimStack containing time series of surface flux (monolevel) and vertical profiles at the specified output frequency.

**Provenance**

The output metadata records where the forcing came from (`dataset`, `latitude`, `longitude`, `elevation_offset`) and how the initial column was prepared. If `profile` was produced by [`gemb_spinup`](/api#GEMB.gemb_spinup), its `spinup_*` / `climatology_*` provenance is copied onto the output and `spinup_performed => true` is set. If `profile` came straight from [`initialize_profile`](/api#GEMB.initialize_profile) (no spinup), the output records `spinup_performed => false`, so the run is unambiguous about having started from an un-spun-up column.

**Metadata**

Every output layer carries CF-style attributes — `units`, `long_name`, the interval reduction as `cell_methods` (`time: sum` / `time: mean` / `time: point`), and a `standard_name` where a CF standard name matches the quantity exactly — and the stack carries `Conventions = "CF-1.11"` alongside the provenance above. See [`GEMB_CF_ATTRIBUTES`](/api#GEMB.GEMB_CF_ATTRIBUTES) for the table and [`cf_attributes`](/api#GEMB.cf_attributes) to read one layer's attributes. Conformance is at the attribute level: GEMB.jl ships no NetCDF writer, so the attributes are there for whichever writer or plotting code consumes them.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.gemb_spinup' href='#GEMB.gemb_spinup'><span class="jlbinding">GEMB.gemb_spinup</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
gemb_spinup(profile, cf, mp; max_iterations=100,
            convergence_delta_density=nothing, convergence_drift_density=nothing,
            drift_window=10, verbose=false)
```


Run GEMB for multiple spinup cycles to reach quasi-steady state.

Forces `output_frequency=:last` internally to minimize memory usage during spinup. Returns the spun-up profile DimStack.

Convergence is judged on the **column-mean density** — averaged over the whole column, since the column depth is fixed for the run (chosen by [`initialize_profile`](/api#GEMB.initialize_profile), held by [`trim_bottom!`](/internals_grid#GEMB.trim_bottom!-Tuple{NamedTuple,%20Float64,%20ModelParameters})) and so the same domain is compared at every cycle.

**Provenance**

The returned profile carries a metadata `NamedTuple` (accessible via `DimensionalData.metadata(profile)`) recording how the spinup ran and which climatology it used: `spinup_cycles`, `spinup_converged`, `spinup_final_delta_density`, `spinup_final_drift_density`, the convergence parameters (`spinup_max_iterations`, `spinup_convergence_delta_density`, `spinup_convergence_drift_density`, `spinup_drift_window`), and the climatology fields copied forward from `cf` (`climatology_window_start`, `climatology_window_stop`, `climatology_n_years`, `climatology_steps_per_year`). This provenance is propagated onto the [`gemb`](/api#GEMB.gemb) output when the spun-up profile is used to start a transient run.

**Convergence criteria**

Two independent tests, either or both of which may be requested. When **both** are given, **both** must hold — a small per-cycle step and no systematic trend.
- `convergence_delta_density`: exit when the absolute change in column-mean density between consecutive cycles is below this value [kg/m³]. Cheap, but a _step_ test: a column creeping steadily at just under the tolerance passes it while still drifting.
  
- `convergence_drift_density`: exit when the magnitude of the least-squares slope of column-mean density against cycle number, over the trailing `drift_window` cycles, is below this value [kg/m³ per cycle]. This is the trend test the step test cannot make: it distinguishes a column oscillating about a settled mean (slope ≈ 0) from one still densifying monotonically, and by using every point in the window it is not fooled by a single anomalous cycle. Inactive until `drift_window` cycles have run.
  

When neither is given the spinup always runs `max_iterations` cycles.

**Keyword arguments**
- `max_iterations`: maximum number of spinup cycles (default 100). The spinup always exits after this many cycles even if convergence has not been reached.
  
- `convergence_delta_density`, `convergence_drift_density`: see above.
  
- `drift_window`: trailing cycles the drift slope is fitted over (default 10). Must be ≥ 2.
  
- `verbose`: print a convergence message when early exit occurs.
  

Derived from MATLAB's `gemb_spinup.m`; the drift criterion and the whole-column mean have no MATLAB counterpart.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>


## Post-processing {#Post-processing}
<details class='jldocstring custom-block' open>
<summary><a id='GEMB.gemb_profile' href='#GEMB.gemb_profile'><span class="jlbinding">GEMB.gemb_profile</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
gemb_profile(out::DimStack)
gemb_profile(out::DimStack, time_extract::DateTime)
```


Extract a column state from GEMB output as a Profile DimStack.

If `time_extract` is not provided, the last time step is used. If `time_extract` does not exactly match any output time, the nearest time step is used.

Matches MATLAB's `gemb_profile.m`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.gemb_interp' href='#GEMB.gemb_interp'><span class="jlbinding">GEMB.gemb_interp</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
gemb_interp(z_center::AbstractMatrix, A::AbstractMatrix, z_target::AbstractVector; interp_method=:linear)
```


Regularize GEMB Lagrangian output onto a consistent vertical grid.

This function is necessary because the vertical spacing of GEMB output evolves with every timestep.

**Arguments**
- `z_center`: MxN matrix of grid cell center heights (from `dz2z(OutData.dz)`)
  
- `A`: MxN matrix (or `DimArray` with `Z` and `Ti` dimensions) of data to regrid
  
- `z_target`: Vector of target depth coordinates defining the output grid
  
- `interp_method`: Interpolation method (`:linear` default, `:nearest`)
  

**Returns**

A `DimArray` of size `length(z_target) × N` with `Z(z_target)` and `Ti` dimensions. When `A` is a `DimArray` with a `Ti` dimension, the time coordinates are preserved.

Target heights outside the column — above the top cell centre or below the bottom cell centre — are **not extrapolated**; they are left as the missing value `NaN`. The column depth and the surface height both migrate over a run (the surface cell is centimetres thick and the base is pinned only to `z_target`'s own tolerance), so a fixed `z_target` range like `(-20, 0)` inevitably reaches past the column. Extrapolating there continued the near-basal temperature gradient into rock and produced firn temperatures above the melt point — values the model itself never holds. `NaN` says "no column here" instead, and Makie renders it as blank.

Differs from MATLAB's `gemb_interp.m`, which extrapolates.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.surface_timeseries' href='#GEMB.surface_timeseries'><span class="jlbinding">GEMB.surface_timeseries</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
surface_timeseries(A::AbstractMatrix)
```


Return the surface (row 1) value in each column of matrix `A`.

GEMB profile output is top-justified, so the surface cell is always row 1. Retained as a named function because it expresses intent at call sites and matches MATLAB's `surface_timeseries.m`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.firn_air_content' href='#GEMB.firn_air_content'><span class="jlbinding">GEMB.firn_air_content</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
firn_air_content(dz, density, density_ice, z_max=Inf) -> m of air
```


Height of air in the column above depth `z_max`: `Σ dz (1 − ρ/ρ_ice)`, with densities above `density_ice` contributing nothing rather than negative air.

The cell straddling `z_max` contributes in proportion to the fraction of it that lies above the cutoff, so the result is a true integral to `z_max` and does not step as cells cross the boundary during a run. `z_max = Inf` (the default) integrates the whole column.

Depth-limited totals — conventionally to 10 m or 20 m — are what the firn-core literature reports (e.g. Vandecrux et al., 2019), so they are what a whole-column value cannot be compared against.

Allocation-free: called every timestep from the `gemb` time loop.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.close_off_age' href='#GEMB.close_off_age'><span class="jlbinding">GEMB.close_off_age</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
close_off_age(density, age, closeoff_density=DENSITY_PORE_CLOSEOFF) -> d
```


Age [d] of the shallowest cell at or above pore close-off density, or `NaN` when the column never closes off.

This is the model's analogue of the gas age at close-off that firn-air and ice-core work reports: the residence time of the firn at the depth where the pore space seals and the column stops exchanging with the atmosphere. `NaN` (rather than the column's total depth or its deepest age) is the correct answer for an open column, following the same convention as `ice_slab_depth` — there is no close-off depth to report an age at.

Meaningful from the first timestep only under `mp.initialize_age = :steady_state`; under `:zero` it reads 0 until a spinup long enough to bury the epoch has run.

Allocation-free: called every timestep from the `gemb` time loop.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>


## Plotting {#Plotting}

`gemb_plot_output` is exported from GEMB, but its methods live in a package extension that loads only once a Makie backend (`CairoMakie`, `GLMakie`, or `WGLMakie`) is present.
<details class='jldocstring custom-block' open>
<summary><a id='GEMB.gemb_plot_output' href='#GEMB.gemb_plot_output'><span class="jlbinding">GEMB.gemb_plot_output</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
gemb_plot_output(output::DimStack; datelims=nothing, depthlims=(-10, 0),
                 variables=nothing, title="")
```


Create a diagnostic dashboard for a GEMB model run (`output = gemb(profile, cf, mp)`).

The figure has two columns sharing a common (linked) time axis:
- **Left — profile fields.** Each 2-D variable (depth × time, e.g. `temperature`, `density`) is interpolated from the evolving Lagrangian grid onto a fixed vertical grid with [`gemb_interp`](/api#GEMB.gemb_interp) and drawn as a heatmap. Colormaps are chosen per field and color limits are clipped to robust (2–98%) percentiles for maximum contrast. Long runs are strided in time so the raster stays a sensible width. Panels are ordered `temperature`, `age`, `density`, `grain_radius`, `water` — physical order, not the order the fields happen to be stored in.
  
- **Right — scalar fields.** 1-D variables are grouped by physical theme and shared unit (air temperature, energy fluxes, mass fluxes, …) and overplotted as time series on a common axis, with a legend when a panel holds more than one series. Air temperature is placed first (top), aligned with the temperature profile heatmap.
  

Axis, colorbar, and legend labels take their units from each layer's CF metadata (see [`GEMB_CF_ATTRIBUTES`](/api#GEMB.GEMB_CF_ATTRIBUTES)) rather than from a hard-coded table, so a panel cannot disagree with the data it draws. Fields stored in kelvin are converted to °C for display; everything else is shown in its stored units.

Redundant decorations are removed: x ticks and the "year" label appear only on the bottom panel of each column, and a header banner records provenance (GEMB version, forcing source, time span, sampling cadence, spinup climatology window + convergence — or "no spinup" — mass budget, generation time).

**Arguments**
- `output::DimStack`: The `DimStack` returned by [`gemb`](/api#GEMB.gemb).
  

**Keywords**
- `datelims`: Optional `(start, stop)` limits for the time axis. Accepts a tuple of `DateTime`s or of decimal years.
  
- `depthlims`: `(low, high)` depth limits (metres, negative down). Sets both the regridding range and the depth axis of the heatmaps. Defaults to `(-10, 0)` — the near-surface zone where most of the dynamics live. Pass a wider range (or derive it from the data) to see the full column.
  
- `variables`: Iterable of variable name `Symbol`s to include. Defaults to every variable in `output` except the profile fields `dz`, `grain_dendricity`, and `grain_sphericity`, and the densification scalar group (a grid diagnostic, fields already represented by the grain-radius panel, and a lower-priority group); pass any of them explicitly here to include them.
  
- `title`: Optional figure title. Defaults to `"GEMB.jl glacier firn model output"`.
  

**Returns**

A Makie `Figure`.

::: tip Requires Makie

This function lives in a package extension. Load a Makie backend (e.g. `using CairoMakie` or `using GLMakie`) before calling it.

:::

**Example**

```julia
using GEMB, CairoMakie
output = gemb(profile_spunup, cf, mp)
fig = gemb_plot_output(output; depthlims=(-5, 0))
save("gemb_diagnostics.png", fig)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>


## GEMB_ClimateForcing Companion {#GEMB_ClimateForcing-Companion}

A `DimStack` is GEMB's neutral, producer-agnostic forcing interface. The companion package [GEMB_ClimateForcing.jl](https://github.com/alex-s-gardner/GEMB_ClimateForcing.jl) is one producer of a conforming DimStack — it downloads ERA5-Land reanalysis data and generates synthetic forcing.

First, install GEMB_ClimateForcing.jl from GitHub (not yet in the General registry):

```julia
using Pkg
Pkg.add(url="https://github.com/alex-s-gardner/GEMB_ClimateForcing.jl")
```


GEMB_ClimateForcing produces a `DimStack`, which the core method `initialize_forcing(::DimStack)`converts to a`ClimateForcing`:

```julia
using GEMB
using GEMB_ClimateForcing

# Download ERA5-Land data
forcing_data = climate_forcing(:era5land, lat, lon; 
                                time_range=..., 
                                token=ENV["CDS_API_KEY"])

# Convert to ClimateForcing
cf = GEMB.initialize_forcing(forcing_data)

# Or generate synthetic forcing
ds = simulate_climate_forcing("test_1", 3)   # returns DimStack
cf = GEMB.initialize_forcing(ds)

# Compute climatological average (ClimateForcing → ClimateForcing)
cf_clim = forcing_climatology(cf)
```


`initialize_forcing(::DimStack)` validates the required fields and metadata, then forwards to the vector method of `initialize_forcing`. It is a core method (always available); `GEMB_ClimateForcing` is one optional producer of a conforming DimStack.

Humidity conversion utilities (`dewpoint_to_vapor_pressure`, `vapor_pressure_to_relative_humidity`, `relative_humidity_to_vapor_pressure`) and climate fitting functions (`fit_air_temperature`, `fit_precipitation`, `fit_longwave_irradiance_delta`, `fit_seasonal_daily_noise`) are provided by `GEMB_ClimateForcing`. See the [GEMB_ClimateForcing documentation](https://github.com/alex-s-gardner/GEMB_ClimateForcing.jl) for details.

## Output Metadata {#Output-Metadata}

The `gemb` output stack is self-describing: each layer carries CF-style attributes and the stack carries CF global attributes. These are the table and accessors behind that.
<details class='jldocstring custom-block' open>
<summary><a id='GEMB.GEMB_CF_ATTRIBUTES' href='#GEMB.GEMB_CF_ATTRIBUTES'><span class="jlbinding">GEMB.GEMB_CF_ATTRIBUTES</span></a> <Badge type="info" class="jlObjectType jlConstant" text="Constant" /></summary>



```julia
GEMB_CF_ATTRIBUTES::Dict{Symbol,CFAttrs}
```


CF attributes for every layer of the [`gemb`](/api#GEMB.gemb) output stack, keyed by layer name. `cell_methods` records how each layer was reduced over the output interval: `"time: sum"` for mass/thickness accumulations, `"time: mean"` for interval averages, `"time: point"` for instantaneous snapshots taken at the output time.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.CFAttrs' href='#GEMB.CFAttrs'><span class="jlbinding">GEMB.CFAttrs</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
CFAttrs(units, long_name, cell_methods; standard_name="", comment="")
```


CF attributes for one output variable. `cell_methods`, `standard_name`, and `comment` are empty strings when they do not apply, and [`cf_attributes`](/api#GEMB.cf_attributes) omits empty fields from the emitted attribute dictionary rather than writing blanks.

`cell_methods` is held separately from the rest so it can be dropped for products with no time coordinate — a `cell_methods` of `"time: point"` must not appear on a variable whose `time` coordinate does not exist.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.cf_attributes' href='#GEMB.cf_attributes'><span class="jlbinding">GEMB.cf_attributes</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
cf_attributes(key::Symbol; time_axis=true)
```


CF attributes for output layer `key` as a `Dict{String,String}` ready to attach as `DimArray` metadata (or to write as NetCDF variable attributes). Attributes that do not apply to the variable are absent rather than empty.

Pass `time_axis=false` for a product with no time coordinate (e.g. a single column from [`gemb_profile`](/api#GEMB.gemb_profile)) to drop `cell_methods`, which would otherwise reference a `time` coordinate that is not there.

A fresh `Dict` is built on each call, so [`GEMB_CF_ATTRIBUTES`](/api#GEMB.GEMB_CF_ATTRIBUTES) cannot be mutated through the returned value.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.GEMB_CF_GLOBAL_ATTRIBUTES' href='#GEMB.GEMB_CF_GLOBAL_ATTRIBUTES'><span class="jlbinding">GEMB.GEMB_CF_GLOBAL_ATTRIBUTES</span></a> <Badge type="info" class="jlObjectType jlConstant" text="Constant" /></summary>



```julia
GEMB_CF_GLOBAL_ATTRIBUTES::Dict{String,Any}
```


CF global attributes stamped onto the [`gemb`](/api#GEMB.gemb) output stack alongside the run's own provenance. A `history` attribute is deliberately omitted: a generation timestamp would make the output non-deterministic and break the numerical regression fingerprint.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>


The helpers that turn that table into `DimStack` keyword arguments — useful if you are building your own stack of GEMB fields, or writing the output to NetCDF:
<details class='jldocstring custom-block' open>
<summary><a id='GEMB.cf_layermetadata' href='#GEMB.cf_layermetadata'><span class="jlbinding">GEMB.cf_layermetadata</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
cf_layermetadata(layers::NamedTuple; time_axis=true)
```


Build the `layermetadata` NamedTuple for a `DimStack` whose layers are `layers`.

Mapping over `keys(layers)` rather than a fixed order is deliberate: DimensionalData requires the `layermetadata` keys to match the layer keys both in membership and in order, and the extracted-profile stack orders its layers differently from the `gemb` output stack.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.cf_time_attributes' href='#GEMB.cf_time_attributes'><span class="jlbinding">GEMB.cf_time_attributes</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
cf_time_attributes()
```


CF attributes for the `Ti` dimension of GEMB output.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.cf_layer_index_attributes' href='#GEMB.cf_layer_index_attributes'><span class="jlbinding">GEMB.cf_layer_index_attributes</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
cf_layer_index_attributes()
```


CF attributes for the `Z` dimension of GEMB output, which is a bare 1-based cell index rather than a coordinate. It therefore carries no `standard_name` and no `positive`: claiming either would present an index as a physical height. Depths come from the `dz` layer via [`dz2z`](/api#GEMB.dz2z).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.cf_height_attributes' href='#GEMB.cf_height_attributes'><span class="jlbinding">GEMB.cf_height_attributes</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
cf_height_attributes()
```


CF attributes for a genuine vertical coordinate, as produced by [`gemb_interp`](/api#GEMB.gemb_interp) when it regrids onto fixed heights. [`dz2z`](/api#GEMB.dz2z) returns heights measured negative-downward from the surface, hence `positive = "up"`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>


## Utilities {#Utilities}
<details class='jldocstring custom-block' open>
<summary><a id='GEMB.dz2z' href='#GEMB.dz2z'><span class="jlbinding">GEMB.dz2z</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
dz2z(dz::AbstractVector)
```


Convert layer thicknesses `dz` to cell center heights (negative below surface). The surface is at z=0; centers are at negative depths.

Matches MATLAB's `dz2z.m` for vector input.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>



```julia
dz2z(dz::AbstractMatrix)
```


Convert 2D layer thickness matrix to cell center heights. Each column is an independent profile, surface first (row 1).

GEMB profile output is top-justified with a fixed row count, so this is a plain cumulative sum with no padding to work around. Any NaN present in the input (e.g. an output array not yet fully written) propagates through its own column from that row down, as `cumsum` implies.

Matches MATLAB's `dz2z.m` for matrix input.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.fast_divisors' href='#GEMB.fast_divisors'><span class="jlbinding">GEMB.fast_divisors</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
fast_divisors(n::Integer)
```


Find all positive divisors of integer `n`, returned sorted. Matches MATLAB's `fast_divisors.m`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.datetime2decyear' href='#GEMB.datetime2decyear'><span class="jlbinding">GEMB.datetime2decyear</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
datetime2decyear(datetimes::AbstractVector{DateTime})
```


Convert Julia `DateTime` objects to decimal year values.

Automatically accounts for leap years. Inverse of the year-fraction convention used by [`decyear2datenum`](/api#GEMB.decyear2datenum).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.decyear2datenum' href='#GEMB.decyear2datenum'><span class="jlbinding">GEMB.decyear2datenum</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
decyear2datenum(decyear)
```


Convert decimal year to MATLAB datenum format (days since 0000-01-01).

Automatically accounts for leap years.

Matches MATLAB's `decyear2datenum.m`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>


## Index {#Index}
- [`GEMB.AQUIFER_TOLERANCE`](#GEMB.AQUIFER_TOLERANCE)
- [`GEMB.BARNOLA_CLOSEOFF`](#GEMB.BARNOLA_CLOSEOFF)
- [`GEMB.BARNOLA_MIN_DENSITY_ICE`](#GEMB.BARNOLA_MIN_DENSITY_ICE)
- [`GEMB.BARNOLA_N`](#GEMB.BARNOLA_N)
- [`GEMB.CROCUS_HYBRID_DENSITY`](#GEMB.CROCUS_HYBRID_DENSITY)
- [`GEMB.DENSIFICATION_COEFFS_M01`](#GEMB.DENSIFICATION_COEFFS_M01)
- [`GEMB.DT_MIN_WARN`](#GEMB.DT_MIN_WARN)
- [`GEMB.GEMB_CF_ATTRIBUTES`](#GEMB.GEMB_CF_ATTRIBUTES)
- [`GEMB.GEMB_CF_GLOBAL_ATTRIBUTES`](#GEMB.GEMB_CF_GLOBAL_ATTRIBUTES)
- [`GEMB.HEAT_CAPACITY_WATER`](#GEMB.HEAT_CAPACITY_WATER)
- [`GEMB.THERMAL_BC_DERIVATIVE_STEP`](#GEMB.THERMAL_BC_DERIVATIVE_STEP)
- [`GEMB.THERMAL_EXPLICIT_SAFETY_FACTOR`](#GEMB.THERMAL_EXPLICIT_SAFETY_FACTOR)
- [`GEMB.THERMAL_IMPLICIT_DT_TARGET`](#GEMB.THERMAL_IMPLICIT_DT_TARGET)
- [`GEMB.AbstractThermalSolver`](#GEMB.AbstractThermalSolver)
- [`GEMB.CFAttrs`](#GEMB.CFAttrs)
- [`GEMB.ClimateForcing`](#GEMB.ClimateForcing)
- [`GEMB.ClimateForcingStep`](#GEMB.ClimateForcingStep)
- [`GEMB.ClimateSummary`](#GEMB.ClimateSummary)
- [`GEMB.DensificationCoeffs`](#GEMB.DensificationCoeffs)
- [`GEMB.ExplicitThermal`](#GEMB.ExplicitThermal)
- [`GEMB.ImplicitThermal`](#GEMB.ImplicitThermal)
- [`GEMB.ModelParameters`](#GEMB.ModelParameters)
- [`GEMB._ThermalSurface`](#GEMB._ThermalSurface)
- [`GEMB._albedo_gardner`](#GEMB._albedo_gardner-Tuple{Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vararg{Float64,%204}})
- [`GEMB._albedo_residual`](#GEMB._albedo_residual-Tuple{Float64,%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Float64,%20Float64,%20Float64,%20Float64,%20Float64,%20Float64,%20Float64,%20ModelParameters})
- [`GEMB._assert_grid_feasible`](#GEMB._assert_grid_feasible-Tuple{Vector{Float64},%20Float64,%20ModelParameters})
- [`GEMB._barnola_f`](#GEMB._barnola_f-Tuple{Float64,%20Float64})
- [`GEMB._cold_content_refreeze`](#GEMB._cold_content_refreeze-Tuple{Real,%20Real,%20Real,%20Real,%20ModelParameters})
- [`GEMB._compute_output_times`](#GEMB._compute_output_times-Tuple{Vector{Dates.DateTime},%20Symbol})
- [`GEMB._crocus_viscosity`](#GEMB._crocus_viscosity-NTuple{5,%20Float64})
- [`GEMB._densification_rate`](#GEMB._densification_rate-Tuple{GEMB.DensificationCoeffs,%20Float64,%20Float64})
- [`GEMB._densify_cell!`](#GEMB._densify_cell!-Tuple{Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Int64,%20Float64,%20Float64,%20Float64})
- [`GEMB._derive_column_depth`](#GEMB._derive_column_depth-Tuple{ModelParameters,%20ClimateSummary})
- [`GEMB._emissivity_initialize`](#GEMB._emissivity_initialize-Tuple{Float64,%20ModelParameters})
- [`GEMB._extract_profile_at_index`](#GEMB._extract_profile_at_index-Tuple{DimStack,%20Int64})
- [`GEMB._find_dt_divisor`](#GEMB._find_dt_divisor-Tuple{Float64,%20Vector{Float64}})
- [`GEMB._gemb_time_loop!`](#GEMB._gemb_time_loop!-Tuple{Any,%20Any,%20Any,%20Any,%20Bool,%20Vector{Dates.DateTime},%20Any,%20Int64,%20Float64,%20Float64,%20AbstractVector,%20AbstractVector,%20AbstractVector,%20AbstractVector,%20AbstractVector,%20AbstractVector,%20AbstractVector,%20AbstractVector,%20AbstractVector,%20AbstractVector,%20AbstractVector,%20AbstractVector,%20AbstractVector,%20Vararg{Float64,%207}})
- [`GEMB._grain_forcing_step`](#GEMB._grain_forcing_step-Tuple{Float64,%20ClimateSummary})
- [`GEMB._grain_gradient`](#GEMB._grain_gradient-Tuple{Vector{Float64},%20Vector{Float64},%20Int64,%20Int64})
- [`GEMB._grain_lwc`](#GEMB._grain_lwc-Tuple{Float64,%20Float64,%20Float64})
- [`GEMB._interp1`](#GEMB._interp1-Tuple{AbstractVector,%20AbstractVector,%20AbstractVector})
- [`GEMB._interp1_nearest`](#GEMB._interp1_nearest-Tuple{AbstractVector,%20AbstractVector,%20AbstractVector})
- [`GEMB._interp_curve`](#GEMB._interp_curve-Tuple{Vector{Float64},%20Vector{Float64},%20Vector{Float64}})
- [`GEMB._irreducible_water`](#GEMB._irreducible_water-Tuple{Real,%20Real,%20Real,%20ModelParameters})
- [`GEMB._marbouty_Q`](#GEMB._marbouty_Q-Tuple{Float64,%20Float64,%20Float64})
- [`GEMB._march_density_step`](#GEMB._march_density_step-Tuple{Any,%20Vararg{Float64,%206}})
- [`GEMB._march_step_depth`](#GEMB._march_step_depth-Tuple{Real,%20Integer})
- [`GEMB._max_safe_dt`](#GEMB._max_safe_dt-Tuple{Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20ModelParameters})
- [`GEMB._seb_annual_melt`](#GEMB._seb_annual_melt-Tuple{Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Float64,%20Float64,%20Float64,%20Float64,%20Float64,%20ModelParameters})
- [`GEMB._seb_forcing_step`](#GEMB._seb_forcing_step-NTuple{9,%20Float64})
- [`GEMB._seb_residual`](#GEMB._seb_residual-Tuple{Float64,%20Float64,%20Float64,%20Float64,%20Float64,%20Float64,%20Float64,%20Float64,%20ClimateForcingStep})
- [`GEMB._seb_skin_temperature`](#GEMB._seb_skin_temperature-Tuple{Float64,%20Float64,%20Float64,%20Float64,%20Float64,%20Float64,%20Float64,%20ClimateForcingStep})
- [`GEMB._select_merge_pair`](#GEMB._select_merge_pair-Tuple{Vector{Float64},%20Vector{Float64},%20Int64,%20ModelParameters})
- [`GEMB._select_split_cell`](#GEMB._select_split_cell-Tuple{Vector{Float64},%20Vector{Float64},%20Int64,%20ModelParameters})
- [`GEMB._snow_cover_albedo`](#GEMB._snow_cover_albedo-Tuple{Real,%20Real,%20ModelParameters})
- [`GEMB._snow_fraction`](#GEMB._snow_fraction-Tuple{Real,%20Real})
- [`GEMB._snow_fraction_from_albedo`](#GEMB._snow_fraction_from_albedo-Tuple{Real,%20ModelParameters})
- [`GEMB._steady_state_temperature`](#GEMB._steady_state_temperature-Tuple{Real,%20Real,%20ClimateSummary,%20ModelParameters})
- [`GEMB._surface_energy_balance`](#GEMB._surface_energy_balance-Tuple{Float64,%20Float64,%20GEMB._ThermalSurface,%20ClimateForcingStep})
- [`GEMB._surface_energy_balance_slope`](#GEMB._surface_energy_balance_slope-Tuple{Float64,%20Float64,%20GEMB._ThermalSurface,%20ClimateForcingStep,%20Float64,%20Float64,%20Float64})
- [`GEMB._thermal_conductivity_air`](#GEMB._thermal_conductivity_air-Tuple{Real})
- [`GEMB._thermal_conductivity_calonne2019`](#GEMB._thermal_conductivity_calonne2019-Tuple{Real,%20Real})
- [`GEMB._thermal_conductivity_ice`](#GEMB._thermal_conductivity_ice-Tuple{Real})
- [`GEMB._thermal_conductivity_marchenko2019`](#GEMB._thermal_conductivity_marchenko2019-Tuple{Real})
- [`GEMB._thermal_solve!`](#GEMB._thermal_solve!-Tuple{ImplicitThermal,%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Float64,%20Vector{Float64},%20GEMB._ThermalSurface,%20ClimateForcingStep,%20ModelParameters,%20Bool})
- [`GEMB._thermal_solve!`](#GEMB._thermal_solve!-Tuple{ExplicitThermal,%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Float64,%20Vector{Float64},%20GEMB._ThermalSurface,%20ClimateForcingStep,%20ModelParameters,%20Bool})
- [`GEMB._turbulent_heat_flux`](#GEMB._turbulent_heat_flux-Tuple{Float64,%20Float64,%20Float64,%20Float64,%20Float64,%20ClimateForcingStep,%20Vararg{Float64,%207}})
- [`GEMB._uniform_ice_profile`](#GEMB._uniform_ice_profile-Tuple{ModelParameters,%20Real})
- [`GEMB.add_energy_temperature`](#GEMB.add_energy_temperature-Tuple{ModelParameters,%20Real,%20Real,%20Real})
- [`GEMB.air_density`](#GEMB.air_density-Tuple{Real,%20Real})
- [`GEMB.apply_horizontal_strain!`](#GEMB.apply_horizontal_strain!-Tuple{NamedTuple,%20Float64,%20ModelParameters})
- [`GEMB.apply_lateral_drainage!`](#GEMB.apply_lateral_drainage!-Tuple{NamedTuple,%20Float64,%20ModelParameters})
- [`GEMB.aquifer_diagnostics`](#GEMB.aquifer_diagnostics-Tuple{Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20ModelParameters})
- [`GEMB.calculate_accumulation`](#GEMB.calculate_accumulation-Tuple{Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20ClimateForcingStep,%20ModelParameters,%20Bool})
- [`GEMB.calculate_albedo`](#GEMB.calculate_albedo-Tuple{Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Float64,%20ClimateForcingStep,%20ModelParameters})
- [`GEMB.calculate_density`](#GEMB.calculate_density-Tuple{Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20ClimateForcingStep,%20ModelParameters})
- [`GEMB.calculate_grain_size`](#GEMB.calculate_grain_size-Tuple{Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20ClimateForcingStep,%20ModelParameters})
- [`GEMB.calculate_melt`](#GEMB.calculate_melt-Tuple{Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Float64,%20ModelParameters,%20Bool})
- [`GEMB.calculate_shortwave_radiation`](#GEMB.calculate_shortwave_radiation-Tuple{Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Float64,%20Float64,%20ClimateForcingStep,%20ModelParameters})
- [`GEMB.calculate_temperature`](#GEMB.calculate_temperature-Tuple{Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Float64,%20Vector{Float64},%20Vector{Float64},%20ClimateForcingStep,%20ModelParameters,%20Bool})
- [`GEMB.cf_attributes`](#GEMB.cf_attributes)
- [`GEMB.cf_height_attributes`](#GEMB.cf_height_attributes-Tuple{})
- [`GEMB.cf_height_attributes`](#GEMB.cf_height_attributes)
- [`GEMB.cf_layer_index_attributes`](#GEMB.cf_layer_index_attributes)
- [`GEMB.cf_layer_index_attributes`](#GEMB.cf_layer_index_attributes-Tuple{})
- [`GEMB.cf_layermetadata`](#GEMB.cf_layermetadata)
- [`GEMB.cf_layermetadata`](#GEMB.cf_layermetadata-Tuple{NamedTuple})
- [`GEMB.cf_time_attributes`](#GEMB.cf_time_attributes-Tuple{})
- [`GEMB.cf_time_attributes`](#GEMB.cf_time_attributes)
- [`GEMB.chord_heat_capacity`](#GEMB.chord_heat_capacity-Tuple{ModelParameters,%20Real,%20Real})
- [`GEMB.close_off_age`](#GEMB.close_off_age)
- [`GEMB.close_slot!`](#GEMB.close_slot!-Tuple{NamedTuple,%20Int64})
- [`GEMB.cold_content_mass`](#GEMB.cold_content_mass-Tuple{ModelParameters,%20Real,%20Real})
- [`GEMB.cold_content_mass`](#GEMB.cold_content_mass-Tuple{ModelParameters,%20Real,%20Real,%20Real})
- [`GEMB.column_bands!`](#GEMB.column_bands!-Tuple{Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20ModelParameters})
- [`GEMB.column_enthalpy`](#GEMB.column_enthalpy-Tuple{ModelParameters,%20AbstractVector{Float64},%20AbstractVector{Float64}})
- [`GEMB.column_mass`](#GEMB.column_mass-Tuple{NamedTuple})
- [`GEMB.column_mass_energy`](#GEMB.column_mass_energy-Tuple{NamedTuple,%20ModelParameters})
- [`GEMB.column_state`](#GEMB.column_state-NTuple{8,%20Any})
- [`GEMB.datetime2decyear`](#GEMB.datetime2decyear)
- [`GEMB.decyear2datenum`](#GEMB.decyear2datenum)
- [`GEMB.densification_lookup_M01`](#GEMB.densification_lookup_M01-Tuple{Symbol})
- [`GEMB.dilute_age`](#GEMB.dilute_age-Tuple{Float64,%20Float64,%20Float64})
- [`GEMB.dz2z`](#GEMB.dz2z)
- [`GEMB.energy_tolerance`](#GEMB.energy_tolerance-Tuple{Real})
- [`GEMB.enforce_column_length!`](#GEMB.enforce_column_length!-Tuple{NamedTuple,%20Int64,%20ModelParameters})
- [`GEMB.enthalpy_temperature_scale`](#GEMB.enthalpy_temperature_scale-Tuple{ModelParameters})
- [`GEMB.excess_specific_enthalpy`](#GEMB.excess_specific_enthalpy-Tuple{ModelParameters,%20Real})
- [`GEMB.excess_temperature_from_specific_enthalpy`](#GEMB.excess_temperature_from_specific_enthalpy-Tuple{ModelParameters,%20Real})
- [`GEMB.fast_divisors`](#GEMB.fast_divisors)
- [`GEMB.firn_air_content`](#GEMB.firn_air_content)
- [`GEMB.forcing_climatology`](#GEMB.forcing_climatology)
- [`GEMB.fresh_snow_density`](#GEMB.fresh_snow_density)
- [`GEMB.full_melt_excess_temperature`](#GEMB.full_melt_excess_temperature-Tuple{ModelParameters})
- [`GEMB.gemb`](#GEMB.gemb)
- [`GEMB.gemb_core`](#GEMB.gemb_core-Tuple{Any,%20ClimateForcingStep,%20ModelParameters,%20Bool})
- [`GEMB.gemb_interp`](#GEMB.gemb_interp)
- [`GEMB.gemb_plot_output`](#GEMB.gemb_plot_output)
- [`GEMB.gemb_profile`](#GEMB.gemb_profile)
- [`GEMB.gemb_spinup`](#GEMB.gemb_spinup)
- [`GEMB.heat_capacity`](#GEMB.heat_capacity-Tuple{ModelParameters,%20Real})
- [`GEMB.hydraulic_conductivity_saturated`](#GEMB.hydraulic_conductivity_saturated-Tuple{Float64,%20Float64})
- [`GEMB.ice_slab_diagnostics`](#GEMB.ice_slab_diagnostics-Tuple{Vector{Float64},%20Vector{Float64},%20ModelParameters})
- [`GEMB.initialize_climate_summary`](#GEMB.initialize_climate_summary)
- [`GEMB.initialize_forcing`](#GEMB.initialize_forcing)
- [`GEMB.initialize_grid`](#GEMB.initialize_grid-Tuple{ModelParameters})
- [`GEMB.initialize_parameters`](#GEMB.initialize_parameters)
- [`GEMB.initialize_profile`](#GEMB.initialize_profile)
- [`GEMB.irreducible_saturation`](#GEMB.irreducible_saturation-Tuple{ModelParameters,%20Float64})
- [`GEMB.manage_layer_thickness`](#GEMB.manage_layer_thickness-Tuple{Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20ModelParameters,%20Bool})
- [`GEMB.melt_mass_from_excess`](#GEMB.melt_mass_from_excess-Tuple{ModelParameters,%20Real,%20Real,%20Real})
- [`GEMB.merge_pair!`](#GEMB.merge_pair!-Tuple{NamedTuple,%20Int64,%20Int64,%20Float64,%20Float64,%20ModelParameters})
- [`GEMB.mix_age`](#GEMB.mix_age-NTuple{4,%20Float64})
- [`GEMB.mix_temperature`](#GEMB.mix_temperature-Tuple{ModelParameters,%20Vararg{Real,%204}})
- [`GEMB.mix_temperature_liquid`](#GEMB.mix_temperature_liquid-Tuple{ModelParameters,%20Vararg{Real,%204}})
- [`GEMB.open_slot!`](#GEMB.open_slot!-Tuple{NamedTuple,%20Int64})
- [`GEMB.pond_blocked_water!`](#GEMB.pond_blocked_water!-Tuple{Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Int64,%20ModelParameters})
- [`GEMB.pore_capacity`](#GEMB.pore_capacity-Tuple{ModelParameters,%20Float64,%20Float64})
- [`GEMB.refreeze_temperature`](#GEMB.refreeze_temperature-Tuple{ModelParameters,%20Real,%20Real,%20Real})
- [`GEMB.relative_permeability`](#GEMB.relative_permeability-Tuple{Float64,%20Float64,%20Float64})
- [`GEMB.remove_melt!`](#GEMB.remove_melt!-Tuple{ModelParameters,%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Int64,%20Float64})
- [`GEMB.runoff_timescale`](#GEMB.runoff_timescale-Tuple{ModelParameters})
- [`GEMB.specific_enthalpy`](#GEMB.specific_enthalpy-Tuple{ModelParameters,%20Real})
- [`GEMB.specific_enthalpy_water`](#GEMB.specific_enthalpy_water-Tuple{ModelParameters})
- [`GEMB.specific_enthalpy_water`](#GEMB.specific_enthalpy_water-Tuple{ModelParameters,%20Real})
- [`GEMB.split_cell!`](#GEMB.split_cell!-Tuple{NamedTuple,%20Int64})
- [`GEMB.steady_state_density`](#GEMB.steady_state_density)
- [`GEMB.steady_state_profile`](#GEMB.steady_state_profile)
- [`GEMB.surface_timeseries`](#GEMB.surface_timeseries)
- [`GEMB.surplus_energy`](#GEMB.surplus_energy-Tuple{ModelParameters,%20Real,%20Real,%20Real})
- [`GEMB.temperature_from_scaled_enthalpy`](#GEMB.temperature_from_scaled_enthalpy-Tuple{ModelParameters,%20Real})
- [`GEMB.temperature_from_specific_enthalpy`](#GEMB.temperature_from_specific_enthalpy-Tuple{ModelParameters,%20Real})
- [`GEMB.thermal_conductivity`](#GEMB.thermal_conductivity-Tuple{Real,%20Real,%20ModelParameters})
- [`GEMB.thermal_conductivity`](#GEMB.thermal_conductivity-Tuple{AbstractVector,%20AbstractVector,%20ModelParameters})
- [`GEMB.thermal_damping_depth`](#GEMB.thermal_damping_depth-Tuple{Real,%20Real,%20ModelParameters})
- [`GEMB.trim_bottom!`](#GEMB.trim_bottom!-Tuple{NamedTuple,%20Float64,%20ModelParameters})
- [`GEMB.turbulent_heat_flux`](#GEMB.turbulent_heat_flux-Tuple{Float64,%20Float64,%20Float64,%20Float64,%20Float64,%20ClimateForcingStep})

