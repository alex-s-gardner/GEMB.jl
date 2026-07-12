module GEMB

using DimensionalData
using Dates
using Statistics

# Physical constants
include("constants.jl")

# Type definitions
include("types.jl")

# Utility functions
include("utilities.jl")

# Initialization
include("initialize_parameters.jl")
include("initialize_forcing.jl")
include("initialize_profile.jl")

# Leaf physics
include("thermal_conductivity.jl")
include("turbulent_heat_flux.jl")
include("densification_lookup.jl")

# Core physics modules
include("calculate_grain_size.jl")
include("calculate_albedo.jl")
include("calculate_shortwave_radiation.jl")
include("calculate_temperature.jl")
include("calculate_accumulation.jl")
include("calculate_melt.jl")
include("calculate_density.jl")
include("manage_layers.jl")

# Integration
include("gemb_core.jl")
include("gemb_driver.jl")

# Utilities
include("spinup.jl")
include("profile_extract.jl")
include("interpolation.jl")
include("forcing_climatology.jl")

# Synthetic forcing
include("simulate/simulate_climate_forcing.jl")

# Climate fitting functions
include("fit_climate/fit_air_temperature.jl")
include("fit_climate/fit_precipitation.jl")
include("fit_climate/fit_longwave_irradiance_delta.jl")
include("fit_climate/fit_seasonal_daily_noise.jl")
include("fit_climate/varname2longname.jl")
include("fit_climate/simulate_coeffs_disp.jl")

# Re-export DimensionalData essentials
using DimensionalData: DimArray, DimStack, Ti, Z, dims
export DimArray, DimStack, Ti, Z

# Exports
export ModelParameters, ClimateForcing, ClimateForcingStep
export initialize_parameters, initialize_forcing, initialize_profile
export gemb, gemb_spinup, gemb_profile, gemb_interp
export forcing_climatology, simulate_climate_forcing
export dz2z, surface_timeseries, fast_divisors, decyear2datenum
export dewpoint_to_vapor_pressure, vapor_pressure_to_relative_humidity
export relative_humidity_to_vapor_pressure
export fit_air_temperature, fit_precipitation, fit_longwave_irradiance_delta
export fit_seasonal_daily_noise, varname2longname, simulate_coeffs_disp

end
