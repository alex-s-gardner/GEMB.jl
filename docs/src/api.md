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
initialize_parameters
initialize_forcing
forcing_climatology
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

## GEMB_ClimateForcing Companion

A `DimStack` is GEMB's neutral, producer-agnostic forcing interface. The companion package [GEMB_ClimateForcing.jl](https://github.com/alex-s-gardner/GEMB_ClimateForcing.jl) is one producer of a conforming DimStack — it downloads ERA5-Land reanalysis data and generates synthetic forcing.

First, install GEMB_ClimateForcing.jl from GitHub (not yet in the General registry):

```julia
using Pkg
Pkg.add(url="https://github.com/alex-s-gardner/GEMB_ClimateForcing.jl")
```

GEMB_ClimateForcing produces a `DimStack`, which the core method
`initialize_forcing(::DimStack)` converts to a `ClimateForcing`:

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

`initialize_forcing(::DimStack)` validates the required fields and metadata, then
forwards to the vector method of `initialize_forcing`. It is a core method (always
available); `GEMB_ClimateForcing` is one optional producer of a conforming DimStack.

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
