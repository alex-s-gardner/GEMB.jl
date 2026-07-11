# Example of running GEMB with ERA5 reanalysis data.
# Equivalent to MATLAB's GEMB_example_ERA5.m
#
# This example requires downloading ERA5 data first.
# See: https://github.com/alex-s-gardner/GEMB/blob/main/docs/ERA5_time_series_data.md
#
# Required packages: NetCDF (for reading .nc files)

using GEMB
using Dates
using Statistics

## Load ERA5 data from NetCDF files
# Note: You must download ERA5 data first. Uncomment and modify paths as needed.
#
# using NetCDF
#
# # Define climate data filenames:
# filename_temperature = "reanalysis-era5-land-timeseries-sfc-2m-temperature_summit.nc"
# filename_pressure    = "reanalysis-era5-land-timeseries-sfc-pressure-precipitation_summit.nc"
# filename_radiation   = "reanalysis-era5-land-timeseries-sfc-radiation-heat_summit.nc"
# filename_wind        = "reanalysis-era5-land-timeseries-sfc-wind_summit.nc"
#
# # Read time vector:
# valid_time = ncread(filename_temperature, "valid_time")
# time_vector = DateTime(1970, 1, 1) .+ Second.(Int.(valid_time))
#
# # Read forcing variables:
# temperature_air = Float64.(ncread(filename_temperature, "t2m"))
# pressure_air = Float64.(ncread(filename_pressure, "sp"))
# precipitation = Float64.(ncread(filename_pressure, "tp")) .* 1000  # m to kg/m²
# precipitation[precipitation .< 0] .= 0  # Fix numerical noise
#
# # Wind speed from vector components:
# u10 = Float64.(ncread(filename_wind, "u10"))
# v10 = Float64.(ncread(filename_wind, "v10"))
# wind_speed = hypot.(u10, v10)
#
# # Radiation (convert accumulated J/m² to average W/m²):
# shortwave_downward = Float64.(ncread(filename_radiation, "ssrd")) ./ 3600
# longwave_downward = Float64.(ncread(filename_radiation, "strd")) ./ 3600
#
# # Convert dewpoint temperature to vapor pressure:
# temperature_dewpoint = Float64.(ncread(filename_temperature, "d2m"))
# vapor_pressure = dewpoint_to_vapor_pressure(temperature_dewpoint)
#
# # Build ClimateForcing:
# cf = initialize_forcing(
#     time_vector,
#     temperature_air,
#     pressure_air,
#     precipitation,
#     wind_speed,
#     shortwave_downward,
#     longwave_downward,
#     vapor_pressure;
#     temperature_air_mean=mean(temperature_air),
#     wind_speed_mean=mean(wind_speed),
#     precipitation_mean=mean(precipitation) * 24 * 365.25,  # annual mean
#     temperature_observation_height=2.0,
#     wind_observation_height=10.0
# )
#
# ## Run GEMB
#
# # Initialize model parameters:
# mp = ModelParameters(output_frequency="daily")
#
# # Initialize grid:
# profile = initialize_profile(mp, cf)
#
# # Create climatological forcing for spinup:
# cf_spinup = forcing_climatology(cf)
#
# # Spin up for 100 years:
# mp_spinup = ModelParameters(output_frequency="last")
# profile_spunup = gemb_spinup(profile, cf_spinup, mp_spinup, 100)
#
# # Run GEMB with spun-up profile:
# output = gemb(profile_spunup, cf, mp)
#
# ## Post-processing examples:
#
# # Get grid cell centers for plotting:
# z_center = dz2z(parent(output[:dz]))
#
# # Get surface temperature time series:
# temp_surface = surface_timeseries(parent(output[:temperature]))
#
# # Regrid to fixed vertical coordinate for plotting:
# temp_gridded = gemb_interp(z_center, parent(output[:temperature]), profile_spunup)
#
# # Convert vapor pressure back to relative humidity:
# rh = vapor_pressure_to_relative_humidity(
#     parent(cf.vapor_pressure), parent(cf.temperature_air))
#
# println("Simulation complete!")
# println("  Mean surface albedo: ", round(mean(parent(output[:albedo_surface])), digits=3))
# println("  Mean firn air content: ", round(mean(parent(output[:firn_air_content])), digits=3), " m")
