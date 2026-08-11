# API Reference

```@meta
CurrentModule = GEMB
```

## Types

```@docs
ModelParameters
ClimateForcing
ClimateForcingStep
```

## Initialization

```@docs
initialize_forcing
initialize_profile
```

## Running the Model

```@docs
gemb
gemb_spinup
```

## Post-processing

```@docs
gemb_profile
gemb_interp
surface_timeseries
```

## GEMB_ClimateForcing Extension

GEMB.jl includes a package extension that provides seamless integration with [GEMB_ClimateForcing.jl](https://github.com/alex-s-gardner/GEMB_ClimateForcing.jl) for downloading ERA5, ERA5-Land, and MERRA-2 reanalysis data.

First, install GEMB_ClimateForcing.jl from GitHub (not yet in the General registry):

```julia
using Pkg
Pkg.add(url="https://github.com/alex-s-gardner/GEMB_ClimateForcing.jl")
```

When both packages are loaded, a conversion method `ClimateForcing(::DimStack)` becomes available:

```julia
using GEMB
using GEMB_ClimateForcing

# Download ERA5-Land data
forcing_data = climate_forcing(:era5land, lat, lon; 
                                time_range=..., 
                                token=ENV["CDS_API_KEY"])

# Convert to ClimateForcing (extension method)
cf = GEMB.ClimateForcing(forcing_data)

# Or generate synthetic forcing
ds = simulate_climate_forcing("test_1", 3)   # returns DimStack
cf = GEMB.ClimateForcing(ds)

# Compute climatological average (DimStack → DimStack)
ds_clim = forcing_climatology(ds)
cf_clim = GEMB.ClimateForcing(ds_clim)
```

The extension automatically validates required fields and metadata, then calls `initialize_forcing` internally. See the extension source at `ext/GEMBClimateForcing.jl` for details.

Humidity conversion utilities (`dewpoint_to_vapor_pressure`, `vapor_pressure_to_relative_humidity`,
`relative_humidity_to_vapor_pressure`) and climate fitting functions (`fit_air_temperature`,
`fit_precipitation`, `fit_longwave_irradiance_delta`, `fit_seasonal_daily_noise`) are provided
by `GEMB_ClimateForcing`. See the
[GEMB_ClimateForcing documentation](https://github.com/alex-s-gardner/GEMB_ClimateForcing.jl)
for details.

## Utilities

```@docs
dz2z
fast_divisors
```

## Index

```@index
```
