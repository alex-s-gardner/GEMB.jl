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

# CF attribute table for output variables (consulted when building output stacks)
include("cf_metadata.jl")

# Initialization
include("initialize_parameters.jl")
include("initialize_forcing.jl")
include("forcing_climatology.jl")
include("initialize_climate_summary.jl")
include("steady_state_profile.jl")
include("initialize_profile.jl")

# Leaf physics
include("heat_capacity.jl")
include("thermal_conductivity.jl")
include("turbulent_heat_flux.jl")
include("densification_lookup.jl")

# Shared vertical-grid primitives (used by the accumulation / melt / layer-management
# modules below, so it must be included before them)
include("grid_ops.jl")

# Core physics modules
include("calculate_grain_size.jl")
include("calculate_albedo.jl")
include("calculate_shortwave_radiation.jl")
include("calculate_temperature.jl")
include("calculate_accumulation.jl")
include("calculate_melt.jl")
include("hydrology.jl")
include("calculate_density.jl")
include("manage_layer_thickness.jl")

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
using DimensionalData: DimArray, DimStack, Ti, Z, dims, metadata
export DimArray, DimStack, Ti, Z, metadata

# Exports
export ModelParameters, ClimateForcing, ClimateForcingStep
export AbstractThermalSolver, ExplicitThermal, ImplicitThermal, ThermalWorkspace
export initialize_parameters, initialize_forcing, forcing_climatology, initialize_profile
export gemb, gemb_spinup, gemb_profile, gemb_interp
export plot_output
export dz2z, firn_air_content, close_off_age, fast_divisors, decyear2datenum, datetime2decyear
export CFAttrs, GEMB_CF_ATTRIBUTES, GEMB_CF_GLOBAL_ATTRIBUTES, cf_attributes
export fresh_snow_density
export ClimateSummary, initialize_climate_summary
export steady_state_profile, steady_state_density

end
