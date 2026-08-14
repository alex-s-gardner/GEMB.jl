using DimensionalData
using DimensionalData: AbstractDimStack, data_eltype
import DimensionalData as DD
using Dates

"""
    ModelParameters

All GEMB model configuration parameters with validation.
Construct with keyword arguments; unspecified fields use defaults.

Matches the 38 fields in MATLAB's `model_initialize_parameters.m`.
"""
Base.@kwdef struct ModelParameters
    # --- General ---
    run_prefix::String = "default"

    # --- Density & Densification ---
    densification_method::Symbol = :Arthern
    densification_coeffs_M01::Symbol = :Gre_RACMO_GS_SW0
    new_snow_method::Symbol = Symbol("350kgm2")
    density_ice::Float64 = 910.0
    rain_temperature_threshold::Float64 = 273.15

    # --- Longwave Emissivity ---
    emissivity_method::Symbol = :uniform
    emissivity::Float64 = 0.97
    emissivity_grain_radius_large::Float64 = 0.97
    emissivity_grain_radius_threshold::Float64 = 10.0
    surface_roughness_effective_ratio::Float64 = 0.10

    # --- Thermal Conductivity ---
    thermal_conductivity_method::Symbol = :Sturm

    # --- Melt & Water ---
    water_irreducible_saturation::Float64 = 0.07

    # --- Albedo & Radiation ---
    albedo_method::Symbol = :GardnerSharp
    albedo_density_threshold::Float64 = Inf
    shortwave_subsurface_absorption::Bool = false
    albedo_snow::Float64 = 0.85
    albedo_ice::Float64 = 0.48
    albedo_fixed::Float64 = 0.85
    shortwave_downward_diffuse::Float64 = 0.0
    solar_zenith_angle::Float64 = 0.0
    cloud_optical_thickness::Float64 = 0.0
    black_carbon_snow::Float64 = 0.0
    black_carbon_ice::Float64 = 0.0
    cloud_fraction::Float64 = 0.1
    albedo_wet_snow_t0::Float64 = 15.0
    albedo_dry_snow_t0::Float64 = 30.0
    albedo_K::Float64 = 7.0

    # --- Output Controls ---
    output_frequency::Symbol = :all

    # --- Grid Geometry ---
    column_ztop::Float64 = 10.0
    column_dztop::Float64 = 0.05
    column_dzmin::Float64 = 0.025
    column_dzmax::Float64 = 0.075
    # `column_depth` sets the depth of the constructed grid (`initialize_grid`) and hence the
    # fixed total column depth held for the whole run.
    column_depth::Float64 = 250.0
    column_zy::Float64 = 1.10

    # --- Thermal Time Stepping ---
    dt_divisors::Vector{Float64} = Float64[]  # pre-computed divisors for thermo sub-stepping; set by gemb driver
end

"""
    ClimateForcing <: DimensionalData.AbstractDimStack

Time-series surface meteorological forcing for GEMB, implemented as an
`AbstractDimStack`. It therefore supports the full DimensionalData stack API —
`keys`, `length`, `dims`, `map`, `layers`, and DimStack-style indexing such as
date-range subsetting `cf[Ti(a .. b)]` and single-time selection `cf[Ti(At(t))]`
— all returning a `ClimateForcing` (a sub-stack), never a scalar step.

# Layers (13 time-series `DimArray`s sharing a common `Ti` dimension)
Required forcing: `temperature_air`, `pressure_air`, `precipitation`,
`wind_speed`, `shortwave_downward`, `longwave_downward`, `vapor_pressure`.
Time-varying model parameters (typically `Fill`-backed): `black_carbon_snow`,
`black_carbon_ice`, `cloud_optical_thickness`, `solar_zenith_angle`,
`shortwave_downward_diffuse`, `cloud_fraction`.

# Metadata (scalars, carried in the stack `metadata` as a `NamedTuple`)
`time_step::Int` [s], plus `temperature_air_mean`, `wind_speed_mean`,
`precipitation_mean`, `temperature_observation_height`, `wind_observation_height`,
source provenance (`dataset`, `latitude`, `longitude`, `elevation_offset`) and
climatology provenance (`climatology_window_start`, `climatology_window_stop`,
`climatology_n_years`, `climatology_steps_per_year`; set by
[`forcing_climatology`](@ref)).

Both layers and scalar metadata are reachable by name via property access
(`cf.temperature_air` returns the layer `DimArray`; `cf.time_step` returns the
concrete scalar), so existing `cf.<field>` code continues to work unchanged.
Construct via the 19-argument positional constructor (13 `DimArray`s then the 6
scalars, in the order listed above) or, more commonly, via [`initialize_forcing`](@ref).
"""
struct ClimateForcing{K,T,N,L,D<:Tuple,R<:Tuple,LD,M,LM} <: AbstractDimStack{K,T,N,L,D}
    data::L
    dims::D
    refdims::R
    layerdims::NamedTuple{K,LD}
    metadata::M
    layermetadata::NamedTuple{K,LM}
    # Inner constructors mirror `DimStack` so the generic AbstractDimStack
    # `rebuild`/`stacktype` machinery (which reconstructs via
    # `basetypeof(s){K,...}(data, dims, refdims, layerdims, metadata, layermetadata)`)
    # finds a matching signature.
    function ClimateForcing(
        data, dims, refdims, layerdims::LD, metadata, layermetadata::NamedTuple{K}
    ) where LD <: NamedTuple{K} where K
        T = data_eltype(data)
        N = length(dims)
        ClimateForcing{K,T,N}(data, dims, refdims, layerdims, metadata, layermetadata)
    end
    function ClimateForcing{K,T,N}(
        data::L, dims::D, refdims::R, layerdims::NamedTuple{K,LD}, metadata::M, layermetadata::NamedTuple{K,LM}
    ) where {K,T,N,L,D,R,LD,M,LM}
        new{K,T,N,L,D,R,LD,M,LM}(data, dims, refdims, layerdims, metadata, layermetadata)
    end
end

# The 13 forcing layers, in positional-constructor order.
const CLIMATE_FORCING_LAYER_KEYS = (
    :temperature_air, :pressure_air, :precipitation, :wind_speed,
    :shortwave_downward, :longwave_downward, :vapor_pressure,
    :black_carbon_snow, :black_carbon_ice, :cloud_optical_thickness,
    :solar_zenith_angle, :shortwave_downward_diffuse, :cloud_fraction,
)
# The scalar metadata carried alongside the layers. The first six are the physics
# scalars (positional in the constructor); the trailing four are provenance
# (passed as keyword args) that trace where the forcing came from and how it was
# adjusted, so they can be surfaced in `gemb` output and diagnostic plots.
const CLIMATE_FORCING_META_KEYS = (
    :time_step, :temperature_air_mean, :wind_speed_mean, :precipitation_mean,
    :temperature_observation_height, :wind_observation_height,
    :dataset, :latitude, :longitude, :elevation, :elevation_native, :elevation_offset,
    :climatology_window_start, :climatology_window_stop,
    :climatology_n_years, :climatology_steps_per_year,
)

"""
    ClimateForcing(temperature_air, pressure_air, precipitation, wind_speed,
                   shortwave_downward, longwave_downward, vapor_pressure,
                   black_carbon_snow, black_carbon_ice, cloud_optical_thickness,
                   solar_zenith_angle, shortwave_downward_diffuse, cloud_fraction,
                   time_step, temperature_air_mean, wind_speed_mean,
                   precipitation_mean, temperature_observation_height,
                   wind_observation_height)

Positional constructor: 13 forcing `DimArray`s (sharing a `Ti` dimension) followed
by the 6 scalar metadata values. Builds the stack layers and stores the scalars in
the stack `metadata` as a `NamedTuple` (keeping their concrete types).

Optional provenance keywords — `dataset` (source name, e.g. `"ERA5-Land"`),
`latitude`, `longitude`, `elevation` (absolute target elevation, m; the surface the
forcing represents after any elevation adjustment), `elevation_native` (the source
dataset's own surface elevation, m, before any adjustment), and `elevation_offset`
(m the forcing was elevation-adjusted by) — are stored alongside the physics scalars
and default to `""`/`NaN`/`NaN`/`NaN`/`NaN`/`0.0`.

`elevation_offset` is forced to `0.0` when `elevation` is not finite: with no known
target elevation there is nothing an offset could be measured against, and carrying a
bare offset invites downstream code to derive a native elevation from `NaN - offset`.

Climatology provenance keywords — `climatology_window_start`,
`climatology_window_stop` (the requested averaging window, `DateTime` or
`nothing`), `climatology_n_years` (number of complete years averaged), and
`climatology_steps_per_year` — are set by [`forcing_climatology`](@ref) and
default to `nothing`/`0` for non-climatological forcing.
"""
function ClimateForcing(
    temperature_air, pressure_air, precipitation, wind_speed,
    shortwave_downward, longwave_downward, vapor_pressure,
    black_carbon_snow, black_carbon_ice, cloud_optical_thickness,
    solar_zenith_angle, shortwave_downward_diffuse, cloud_fraction,
    time_step, temperature_air_mean, wind_speed_mean, precipitation_mean,
    temperature_observation_height, wind_observation_height;
    dataset::AbstractString="", latitude::Real=NaN, longitude::Real=NaN,
    elevation::Real=NaN, elevation_native::Real=NaN, elevation_offset::Real=0.0,
    climatology_window_start=nothing, climatology_window_stop=nothing,
    climatology_n_years::Integer=0, climatology_steps_per_year::Integer=0,
)
    das = NamedTuple{CLIMATE_FORCING_LAYER_KEYS}((
        temperature_air, pressure_air, precipitation, wind_speed,
        shortwave_downward, longwave_downward, vapor_pressure,
        black_carbon_snow, black_carbon_ice, cloud_optical_thickness,
        solar_zenith_angle, shortwave_downward_diffuse, cloud_fraction,
    ))
    meta = NamedTuple{CLIMATE_FORCING_META_KEYS}((
        time_step, temperature_air_mean, wind_speed_mean, precipitation_mean,
        temperature_observation_height, wind_observation_height,
        String(dataset), Float64(latitude), Float64(longitude),
        Float64(elevation), Float64(elevation_native),
        # No target elevation ⇒ no meaningful offset (see docstring).
        isfinite(elevation) ? Float64(elevation_offset) : 0.0,
        climatology_window_start, climatology_window_stop,
        Int(climatology_n_years), Int(climatology_steps_per_year),
    ))
    stackdims = DD.combinedims(collect(das))
    data = map(parent, das)
    layerdims = map(DD.basedims, das)
    layermetadata = map(DD.metadata, das)
    ClimateForcing(data, stackdims, (), layerdims, meta, layermetadata)
end

# Preserve the `cf.<field>` interface on top of the stack layout: forcing fields
# resolve to their layer `DimArray`, scalar metadata to the concrete stored value,
# and everything else (`data`, `dims`, `refdims`, ...) to the underlying struct field.
Base.@constprop :aggressive function Base.getproperty(cf::ClimateForcing, k::Symbol)
    if k in CLIMATE_FORCING_LAYER_KEYS
        return cf[k]
    elseif k in CLIMATE_FORCING_META_KEYS
        return getfield(cf, :metadata)[k]
    else
        return getfield(cf, k)
    end
end

"""
    ClimateForcingStep

Single time-step forcing values extracted from ClimateForcing.
Plain struct of scalars for efficient use in the physics loop.
"""
struct ClimateForcingStep
    dt::Float64
    temperature_air::Float64
    pressure_air::Float64
    precipitation::Float64
    wind_speed::Float64
    shortwave_downward::Float64
    longwave_downward::Float64
    vapor_pressure::Float64
    # Metadata
    temperature_air_mean::Float64
    wind_speed_mean::Float64
    precipitation_mean::Float64
    temperature_observation_height::Float64
    wind_observation_height::Float64
    # Time-varying model params (from forcing or ModelParam defaults)
    black_carbon_snow::Float64
    black_carbon_ice::Float64
    cloud_optical_thickness::Float64
    solar_zenith_angle::Float64
    shortwave_downward_diffuse::Float64
    cloud_fraction::Float64
end
