using DimensionalData
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
    densification_method::String = "Arthern"
    densification_coeffs_M01::String = "Gre_RACMO_GS_SW0"
    new_snow_method::String = "350kgm2"
    density_ice::Float64 = 910.0
    rain_temperature_threshold::Float64 = 273.15

    # --- Longwave Emissivity ---
    emissivity_method::String = "uniform"
    emissivity::Float64 = 0.97
    emissivity_grain_radius_large::Float64 = 0.97
    emissivity_grain_radius_threshold::Float64 = 10.0
    surface_roughness_effective_ratio::Float64 = 0.10

    # --- Thermal Conductivity ---
    thermal_conductivity_method::String = "Sturm"

    # --- Melt & Water ---
    water_irreducible_saturation::Float64 = 0.07

    # --- Albedo & Radiation ---
    albedo_method::String = "GardnerSharp"
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
    output_frequency::String = "all"
    output_padding::Int = 1000

    # --- Grid Geometry ---
    column_ztop::Float64 = 10.0
    column_dztop::Float64 = 0.05
    column_dzmin::Float64 = 0.025
    column_dzmax::Float64 = 0.075
    column_zmax::Float64 = 250.0
    column_zmin::Float64 = 130.0
    column_zy::Float64 = 1.10

    # --- Thermal Time Stepping ---
    dt_divisors::Vector{Float64} = Float64[]  # pre-computed divisors for thermo sub-stepping; set by gemb driver
end

"""
    ClimateForcing

Time-series surface meteorological forcing for GEMB.
All forcing arrays share a common `Ti` (time) dimension.

Required fields: temperature_air, pressure_air, precipitation, wind_speed,
shortwave_downward, longwave_downward, vapor_pressure.

Metadata: temperature_air_mean, wind_speed_mean, precipitation_mean,
temperature_observation_height, wind_observation_height.
"""
struct ClimateForcing
    # Time-series fields (all DimArray with Ti dimension)
    temperature_air::DimArray
    pressure_air::DimArray
    precipitation::DimArray
    wind_speed::DimArray
    shortwave_downward::DimArray
    longwave_downward::DimArray
    vapor_pressure::DimArray

    # Metadata (scalars)
    temperature_air_mean::Float64
    wind_speed_mean::Float64
    precipitation_mean::Float64
    temperature_observation_height::Float64
    wind_observation_height::Float64
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
