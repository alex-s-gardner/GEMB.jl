module GEMB

using DimensionalData
using Dates
using Statistics
using FillArrays: Fill
import DataInterpolations

# Physical constants
include("constants.jl")

# Type definitions
include("types.jl")

# Utility functions
include("utilities.jl")

# Initialization
include("initialize_parameters.jl")
include("initialize_forcing.jl")
include("forcing_climatology.jl")
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

# Plotting (implementation provided by the GEMBMakieExt extension)
include("plotting.jl")

# Re-export DimensionalData essentials
using DimensionalData: DimArray, DimStack, Ti, Z, dims
export DimArray, DimStack, Ti, Z

# Exports
export ModelParameters, ClimateForcing, ClimateForcingStep
export initialize_parameters, initialize_forcing, forcing_climatology, initialize_profile
export gemb, gemb_spinup, gemb_profile, gemb_interp
export gemb_plot_output
export dz2z, surface_timeseries, fast_divisors, decyear2datenum, datetime2decyear
export herron_langway_steady_state, annual_pdd_melt, fresh_snow_density

end
