# Example of running GEMB with ERA5 reanalysis data using GEMB_ClimateForcing.jl
#
# This example uses the GEMB_ClimateForcing.jl package to automatically download
# and format ERA5-Land climate data.
#
# Setup (required before running):
#   1. Install GEMB_ClimateForcing from GitHub:
#      using Pkg
#      Pkg.add(url="https://github.com/alex-s-gardner/GEMB_ClimateForcing.jl")
#   2. Get a CDS API key from: https://cds.climate.copernicus.eu/api-how-to
#   3. Set environment variable: ENV["CDS_API_KEY"] = "your-key-here"
#
# NOTE: This example requires GEMB_ClimateForcing.jl to be installed.

using GEMB
using DimensionalData
using Dates
using Statistics

# Check if GEMB_ClimateForcing is available
try
    using GEMB_ClimateForcing
catch e
    @error """
    GEMB_ClimateForcing.jl not found!

    To run this example, install GEMB_ClimateForcing.jl from GitHub:
        using Pkg
        Pkg.add(url="https://github.com/alex-s-gardner/GEMB_ClimateForcing.jl")

    Then get a CDS API key from: https://cds.climate.copernicus.eu/api-how-to
    And set: ENV["CDS_API_KEY"] = "your-key-here"
    """
    rethrow(e)
end

## Download ERA5-Land forcing data

# Download data for Summit Station, Greenland (72.58°N, 38.48°W)
# This will download data from the Copernicus Climate Data Store
forcing_data = climate_forcing(:era5land, 72.58, -38.48;
                                time_range=(DateTime(2020,1,1), DateTime(2020,12,31)),
                                token=ENV["CDS_API_KEY"])

# Convert to GEMB ClimateForcing (automatic via package extension)
cf = GEMB.ClimateForcing(forcing_data)

## Run GEMB

# Initialize model parameters
mp = ModelParameters(output_frequency="daily")

# Initialize the firn column
profile = initialize_profile(mp, cf)

# Create climatological forcing for spinup
cf_spinup = forcing_climatology(cf)

# Spin up for 100 years to reach quasi-steady state
mp_spinup = ModelParameters(output_frequency="last")
profile_spunup = gemb_spinup(profile, cf_spinup, mp_spinup, 100)

# Run GEMB with transient forcing and the spun-up profile
output = gemb(profile_spunup, cf, mp)

## Post-processing examples

# Get grid cell centers for plotting
z_center = dz2z(parent(output[:dz]))

# Get surface temperature time series
temp_surface = surface_timeseries(parent(output[:temperature]))

# Regrid to fixed vertical coordinate for plotting
temp_gridded = gemb_interp(z_center, parent(output[:temperature]), profile_spunup)

# Convert vapor pressure back to relative humidity
rh = vapor_pressure_to_relative_humidity(
    parent(cf.vapor_pressure), parent(cf.temperature_air))

println("Simulation complete!")
println("  Mean surface albedo: ", round(mean(parent(output[:albedo_surface])), digits=3))
println("  Mean firn air content: ", round(mean(parent(output[:firn_air_content])), digits=3), " m")
