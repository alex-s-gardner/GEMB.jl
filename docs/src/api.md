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
forcing_climatology
surface_timeseries
```

## Synthetic Forcing

```@docs
simulate_climate_forcing
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

# Download data
forcing_data = climate_forcing(:era5land, lat, lon; 
                                time_range=..., 
                                token=ENV["CDS_API_KEY"])

# Convert to ClimateForcing (extension method)
cf = GEMB.ClimateForcing(forcing_data)
```

The extension automatically validates required fields and metadata, then calls `initialize_forcing` internally. See the extension source at `ext/GEMBClimateForcing.jl` for details.

## Humidity Conversions

```@docs
dewpoint_to_vapor_pressure
relative_humidity_to_vapor_pressure
vapor_pressure_to_relative_humidity
```

## Utilities

```@docs
dz2z
fast_divisors
```

## Index

```@index
```
