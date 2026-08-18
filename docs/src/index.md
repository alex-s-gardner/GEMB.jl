```@raw html
---
layout: home

hero:
  name: "GEMB.jl"
  text: "Glacier Energy and Mass Balance"
  tagline: A column model of firn processes for cryosphere research — surface energy balance, densification, meltwater percolation, and grain metamorphism in one vertical column.
  image:
    src: /logo.png
    alt: GEMB.jl
  actions:
    - theme: brand
      text: Quick Start
      link: /#Quick-Start
    - theme: alt
      text: Model Architecture
      link: /architecture
    - theme: alt
      text: API Reference
      link: /api
    - theme: alt
      text: View on GitHub
      link: https://github.com/alex-s-gardner/GEMB.jl

features:
  - title: Surface energy balance
    details: Radiative and turbulent fluxes via Monin-Obukhov similarity theory, coupled to subsurface conduction and solved together.
  - title: Firn densification
    details: Compaction from snow through firn to close-off ice, with a choice of empirical and semi-empirical schemes.
  - title: Meltwater and refreezing
    details: Tipping-bucket percolation with irreducible retention, ice-slab blocking, and optional delayed runoff for firn aquifers.
  - title: Grain metamorphism and age
    details: Effective grain radius, dendricity, and sphericity evolve with the snow; every cell carries a mass-weighted age.
---
```

# GEMB.jl

```@meta
CurrentModule = GEMB
```

GEMB.jl is a Julia implementation of the **Glacier Energy and Mass Balance** model (GEMB —
the "B" is silent), a one-dimensional model of the surface energy balance and vertical firn
evolution of glaciers and ice sheets. It couples atmospheric forcing to subsurface
thermodynamics and densification physics to resolve temperature, density, water content,
grain properties, and layer age through the column.

GEMB is a column model of intermediate complexity with no horizontal communication,
prioritizing computational efficiency to accommodate the multi-millennial spin-ups required
to initialize deep firn columns. It is used for interpreting satellite altimetry, firn and
ice-core studies, surface mass balance inversion, and uncertainty quantification.

A complete description of version 1.0 of the model is given in
[Gardner et al., 2023](https://doi.org/10.5194/gmd-16-2277-2023):

> Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass Balance (GEMB): a
> model of firn processes for cryosphere research, *Geoscientific Model Development*, 16,
> 2277–2302, [https://doi.org/10.5194/gmd-16-2277-2023](https://doi.org/10.5194/gmd-16-2277-2023), 2023.

GEMB.jl is the reference implementation and is developed independently. It began as a
translation of an earlier MATLAB version; since v2.0.0 its physics, defaults, and numerics
are set on their own merits — against the published literature and against the
[Community Firn Model](https://github.com/UWGlaciology/CommunityFirnModel) and
[IMAU-FDM](https://github.com/IMAU-ice-and-climate/IMAU-FDM). The
[MATLAB version](https://github.com/alex-s-gardner/GEMB) is maintained separately and the
two are no longer expected to agree numerically.

## Installation

GEMB.jl is **not yet registered** in the Julia General registry. Install it from GitHub by URL:

```julia
using Pkg
Pkg.add(url="https://github.com/alex-s-gardner/GEMB.jl")
```

Climate forcing generation lives in the companion package
[GEMB_ClimateForcing.jl](https://github.com/alex-s-gardner/GEMB_ClimateForcing.jl), also
unregistered and also installed by URL. It is optional — GEMB consumes any conforming
`DimStack` — but it is what the examples below use:

```julia
Pkg.add(url="https://github.com/alex-s-gardner/GEMB_ClimateForcing.jl")
```

## Quick Start

Running GEMB takes four steps:

1. **Define climate forcing** — [`initialize_forcing`](@ref) converts a `DimStack` of
   meteorological time series into a [`ClimateForcing`](@ref).
2. **Define model parameters** — [`initialize_parameters`](@ref) builds and validates a
   [`ModelParameters`](@ref).
3. **Initialize the firn/ice column** — [`initialize_profile`](@ref) sets up the vertical
   profile of temperature, density, grain properties, and age.
4. **Run GEMB** — [`gemb`](@ref) loops over the forcing and returns a `DimStack`.

```julia
using GEMB
using GEMB_ClimateForcing

# 1. Climate forcing: synthetic 3-hourly series (a DimStack) → ClimateForcing
ds = simulate_climate_forcing("test_1", 3)
cf = initialize_forcing(ds)

# 2. Model parameters
mp = initialize_parameters(output_frequency=:daily)

# 3. Initialize the firn/ice column
profile = initialize_profile(mp, cf)

# 4. Run
output = gemb(profile, cf, mp)
```

## What the output looks like

[`gemb`](@ref) returns a
[DimensionalData.jl](https://github.com/rafaqz/DimensionalData.jl) `DimStack`: 26 monolevel
time series over `Ti`, and 8 profile fields over `Z × Ti`. In the REPL:

```julia-repl
julia> output
┌ 11688×256 DimStack ┐
├────────────────────┴─────────────────────────────────────────────────────────────── dims ┐
  ↓ Ti Sampled{Dates.DateTime} [DateTime("1994-01-01T21:00:00"), …, DateTime("2025-12-31T21:00:00")] ForwardOrdered Irregular Points,
  → Z  Sampled{Int64} 1:256 ForwardOrdered Regular Points
├────────────────────────────────────────────────────────────────────────────────── layers ┤
  :melt                          eltype: Float64 dims: Ti size: 11688
  :runoff                        eltype: Float64 dims: Ti size: 11688
  :refreeze                      eltype: Float64 dims: Ti size: 11688
  :evaporation_condensation      eltype: Float64 dims: Ti size: 11688
  :shortwave_net                 eltype: Float64 dims: Ti size: 11688
  :longwave_net                  eltype: Float64 dims: Ti size: 11688
  :heat_flux_sensible            eltype: Float64 dims: Ti size: 11688
  :heat_flux_latent              eltype: Float64 dims: Ti size: 11688
  :heat_flux_basal               eltype: Float64 dims: Ti size: 11688
  :albedo_broadband              eltype: Float64 dims: Ti size: 11688
  :densification_from_compaction eltype: Float64 dims: Ti size: 11688
  :densification_from_melt       eltype: Float64 dims: Ti size: 11688
  :strain_thinning               eltype: Float64 dims: Ti size: 11688
  :ice_flux                      eltype: Float64 dims: Ti size: 11688
  :firn_air_content              eltype: Float64 dims: Ti size: 11688
  :firn_air_content_10m          eltype: Float64 dims: Ti size: 11688
  :firn_air_content_20m          eltype: Float64 dims: Ti size: 11688
  :close_off_age                 eltype: Float64 dims: Ti size: 11688
  :percolation_depth             eltype: Float64 dims: Ti size: 11688
  :ice_slab_thickness            eltype: Float64 dims: Ti size: 11688
  :ice_slab_depth                eltype: Float64 dims: Ti size: 11688
  :aquifer_thickness             eltype: Float64 dims: Ti size: 11688
  :aquifer_depth                 eltype: Float64 dims: Ti size: 11688
  :temperature_air               eltype: Float64 dims: Ti size: 11688
  :precipitation                 eltype: Float64 dims: Ti size: 11688
  :rain                          eltype: Float64 dims: Ti size: 11688
  :temperature                   eltype: Float64 dims: Z, Ti size: 256×11688
  :dz                            eltype: Float64 dims: Z, Ti size: 256×11688
  :density                       eltype: Float64 dims: Z, Ti size: 256×11688
  :water                         eltype: Float64 dims: Z, Ti size: 256×11688
  :grain_radius                  eltype: Float64 dims: Z, Ti size: 256×11688
  :grain_dendricity              eltype: Float64 dims: Z, Ti size: 256×11688
  :grain_sphericity              eltype: Float64 dims: Z, Ti size: 256×11688
  :age                           eltype: Float64 dims: Z, Ti size: 256×11688
├──────────────────────────────────────────────────────────────────────────────── metadata ┤
  Dict{String, Any} with 10 entries:
  "latitude"         => -73.3307
  "longitude"        => 290.625
  "source"           => "GEMB.jl"
  "dataset"          => "synthetic"
  "spinup_performed" => false
  "title"            => "GEMB (Glacier Energy and Mass Balance) point simulation output"
  "Conventions"      => "CF-1.11"
  ⋮
└──────────────────────────────────────────────────────────────────────────────────────────┘
```

Every layer is self-describing: it carries CF-style `units`, `long_name`, `cell_methods`, and
a `standard_name` where a CF standard name applies exactly.

```julia-repl
julia> output[:melt]
┌ 11688-element DimArray{Float64, 1} melt ┐
├─────────────────────────────────────────┴────────────────────────────────────────── dims ┐
  ↓ Ti Sampled{Dates.DateTime} [DateTime("1994-01-01T21:00:00"), …, DateTime("2025-12-31T21:00:00")] ForwardOrdered Irregular Points
├──────────────────────────────────────────────────────────────────────────────── metadata ┤
  Dict{String, String} with 5 entries:
  "units"         => "kg m-2"
  "cell_methods"  => "time: sum"
  "long_name"     => "meltwater produced"
  "standard_name" => "surface_snow_melt_amount"
  "comment"       => "Column total, including bare-ice melt in the ablation regime."
└──────────────────────────────────────────────────────────────────────────────────────────┘
 1994-01-01T21:00:00  11.4082
 ⋮
 2025-12-31T21:00:00   0.0
```

Index with `DimensionalData` keyword syntax. Profile output is top-justified — the surface
cell is always row 1 — so the surface row of any profile field is `[Z=1]`:

```julia-repl
julia> output[:temperature][Z=1]
┌ 11688-element DimArray{Float64, 1} ┐
├────────────────────────────────────┴─────────────────────────────────────────────── dims ┐
  ↓ Ti Sampled{Dates.DateTime} [DateTime("1994-01-01T21:00:00"), …, DateTime("2025-12-31T21:00:00")] ForwardOrdered Irregular Points
└──────────────────────────────────────────────────────────────────────────────────────────┘
 1994-01-01T21:00:00  272.615
 1994-01-02T21:00:00  272.934
 1994-01-03T21:00:00  270.271
 ⋮
 2025-12-31T21:00:00  261.531

julia> output[:density][Z=1:5, Ti=At(DateTime(2020, 7, 1, 21))]   # top five cells
julia> sum(output[:melt])                                          # column total [kg m-2]
```

`Z` is a **cell index, not a depth** — the grid is Lagrangian, so cells move with the firn.
Convert with [`dz2z`](@ref), or regrid onto fixed depths with [`gemb_interp`](@ref).

## Diagnostic plot

With a Makie backend loaded, [`plot_output`](@ref) draws the whole run — profile fields
as heatmaps on the left, scalar time series on the right, on a shared time axis:

```julia
using CairoMakie
fig = plot_output(output; depthlims=(-10, 0))
save("gemb_diagnostics.png", fig)
```

![GEMB.jl diagnostic output](assets/gemb_output_example.png)

This is the output of `examples/synthetic_example.jl` — 32 years of synthetic forcing run from
a column spun up to convergence first (see [Spinup](#Spinup) below), which is why firn air
content settles into a repeating seasonal cycle instead of drifting through a cold-start
transient. The banner records that provenance: version, forcing source, time span, cadence,
spinup window and convergence, and the closed mass budget. Panels take their units from each
layer's CF metadata, so a label cannot disagree with the data it draws.

## Spinup

For research applications the column should be spun up to a quasi-steady state before running
transient forcing. [`forcing_climatology`](@ref) averages complete years into a one-year
cycle, and [`gemb_spinup`](@ref) repeats it:

```julia
mp = initialize_parameters()                 # gemb_spinup forces output_frequency=:last itself
cf_clim = forcing_climatology(cf)            # ClimateForcing → one-year climatology

profile = initialize_profile(mp, cf_clim)
spun_up = gemb_spinup(profile, cf_clim, mp; max_iterations=200,
                      convergence_delta_density=0.01,
                      convergence_drift_density=0.005)
```

Convergence is judged on column-mean density. Two tests are available, and when both are
given both must hold: `convergence_delta_density` bounds the change between consecutive
cycles, while `convergence_drift_density` bounds the *trend* — the least-squares slope of
column-mean density against cycle over the trailing `drift_window` cycles. The trend test is
the stricter claim: a column creeping steadily at just under the delta tolerance passes the
step test while still densifying.

The spun-up profile carries its provenance, and that provenance propagates onto the output of
the transient run:

```julia-repl
julia> metadata(spun_up)
(spinup_cycles = …, spinup_converged = …, climatology_n_years = …, …)

julia> output = gemb(spun_up, cf, initialize_parameters(output_frequency=:daily));

julia> metadata(output)["spinup_performed"]
true
```

## Using real climate data

For production runs, `GEMB_ClimateForcing` downloads and formats ERA5-Land reanalysis:

```julia
using GEMB, GEMB_ClimateForcing

# Summit Station, Greenland
forcing_data = climate_forcing(:era5land, 72.58, -38.48;
                               time_range=(DateTime(2020,1,1), DateTime(2020,12,31)),
                               token=ENV["CDS_API_KEY"])

cf = initialize_forcing(forcing_data)   # DimStack → ClimateForcing

mp = initialize_parameters(output_frequency=:daily)
profile = initialize_profile(mp, cf)
output = gemb(profile, cf, mp)
```

A `DimStack` is GEMB's neutral, producer-agnostic forcing interface, and
`initialize_forcing(::DimStack)` is a **core** method — always available, no extension
required. It validates the required layers, a `DateTime`-indexed `Ti` dimension, and the
required metadata. Any source that produces a conforming `DimStack` works;
`GEMB_ClimateForcing` is one such producer.

## Output variables

The `Ti` monolevel series and `Z × Ti` profile fields are documented field by field, with
units and `cell_methods`, in the [Variable Reference](@ref "GEMB variables"). The table
itself lives in [`GEMB_CF_ATTRIBUTES`](@ref); read one layer's attributes with
[`cf_attributes`](@ref).

Profile fields are instantaneous snapshots (`cell_methods = "time: point"`), top-justified
with the surface at row 1: `temperature`, `dz`, `density`, `water`, `grain_radius`,
`grain_dendricity`, `grain_sphericity`, and `age` — the mass-weighted mean age in decimal
days of all mass in the cell, matrix plus pore water, measured from column initialization.
Snowfall, rain, and vapour deposition enter at age 0; meltwater carries the age of the firn it
melted from.

## Examples

Runnable scripts are in the `examples/` directory:

- **`examples/synthetic_example.jl`** — complete workflow on synthetic forcing (spinup + run)
- **`examples/era5_example.jl`** — the same against ERA5 reanalysis, with download instructions

## Author information

GEMB was created by Alex Gardner, with contributions from Nicole-Jeanne Schlegel and Chad
Greene. The Julia implementation is at
[github.com/alex-s-gardner/GEMB.jl](https://github.com/alex-s-gardner/GEMB.jl).
