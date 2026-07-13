"""
Extension for GEMB.jl to convert GEMB_ClimateForcing DimStack to ClimateForcing struct.

This extension is automatically loaded when both GEMB and GEMB_ClimateForcing are loaded.
It provides a conversion constructor that validates the DimStack structure and converts
it to GEMB's ClimateForcing type.
"""
module GEMBClimateForcing

using GEMB
using GEMB_ClimateForcing
using DimensionalData
using Dates

"""
    GEMB.ClimateForcing(stack::DimStack) -> ClimateForcing

Convert a GEMB_ClimateForcing DimStack to a GEMB ClimateForcing struct.

# Arguments
- `stack::DimStack`: DimStack from `climate_forcing()` with required variables

# Required Variables in DimStack
- `temperature_air`: Air temperature (K)
- `pressure_air`: Surface pressure (Pa)
- `vapor_pressure`: Vapor pressure (Pa)
- `wind_speed`: Wind speed (m/s)
- `precipitation`: Precipitation rate (kg/m²/hr or kg/m²/timestep)
- `shortwave_downward`: Downward shortwave radiation (W/m²)
- `longwave_downward`: Downward longwave radiation (W/m²)

# Required Metadata
- `temperature_air_mean`: Mean air temperature (K)
- `wind_speed_mean`: Mean wind speed (m/s)
- `precipitation_mean`: Annual precipitation (kg/m²/year)
- `temperature_observation_height`: Height of temperature observations (m)
- `wind_observation_height`: Height of wind observations (m)

# Optional Variables (defaults used if not present)
- `black_carbon_snow`: Black carbon concentration in snow (default: 0.0)
- `black_carbon_ice`: Black carbon concentration in ice (default: 0.0)
- `cloud_optical_thickness`: Cloud optical thickness (default: 0.0)
- `solar_zenith_angle`: Solar zenith angle (degrees, default: 0.0)
- `shortwave_downward_diffuse`: Diffuse shortwave radiation (W/m², default: 0.0)
- `cloud_fraction`: Cloud fraction (default: 0.1)

# Examples
```julia
using GEMB
using GEMB_ClimateForcing

# Load forcing data
forcing_data = climate_forcing(:era5land, 67.0, -50.0;
                                time_range=(DateTime(2020,1,1), DateTime(2020,12,31)),
                                token=ENV["CDS_API_KEY"])

# Convert to ClimateForcing (extension method)
cf = GEMB.ClimateForcing(forcing_data)

# Use with GEMB
mp = GEMB.ModelParameters()
profile = GEMB.initialize_profile(mp, cf)
output = GEMB.gemb(profile, cf, mp)
```
"""
function GEMB.ClimateForcing(stack::DimStack)
    # Validate required fields
    required_fields = [
        :temperature_air, :pressure_air, :vapor_pressure,
        :wind_speed, :precipitation,
        :shortwave_downward, :longwave_downward
    ]

    missing_fields = filter(f -> !haskey(stack, f), required_fields)
    if !isempty(missing_fields)
        throw(ArgumentError(
            "DimStack missing required fields: $(join(missing_fields, ", ")). " *
            "Required: $(join(required_fields, ", "))"
        ))
    end

    # Validate time dimension
    if !hasdim(stack, Ti)
        throw(ArgumentError(
            "DimStack must have a Ti (time) dimension. " *
            "Found dimensions: $(join(dims(stack), ", "))"
        ))
    end

    # Extract time coordinate
    time_dim = dims(stack, Ti)
    time = lookup(time_dim)

    if !isa(time, AbstractVector{DateTime})
        throw(ArgumentError(
            "Ti dimension must be indexed by DateTime values, got $(typeof(time))"
        ))
    end

    # Extract data vectors from DimStack
    temperature_air = parent(stack[:temperature_air])
    pressure_air = parent(stack[:pressure_air])
    vapor_pressure = parent(stack[:vapor_pressure])
    wind_speed = parent(stack[:wind_speed])
    precipitation = parent(stack[:precipitation])
    shortwave_downward = parent(stack[:shortwave_downward])
    longwave_downward = parent(stack[:longwave_downward])

    # Extract required metadata
    meta = metadata(stack)

    required_meta = [
        "temperature_air_mean", "wind_speed_mean", "precipitation_mean",
        "temperature_observation_height", "wind_observation_height"
    ]

    missing_meta = filter(m -> !haskey(meta, m), required_meta)
    if !isempty(missing_meta)
        throw(ArgumentError(
            "DimStack metadata missing required fields: $(join(missing_meta, ", ")). " *
            "Required: $(join(required_meta, ", "))"
        ))
    end

    temperature_air_mean = Float64(meta["temperature_air_mean"])
    wind_speed_mean = Float64(meta["wind_speed_mean"])
    precipitation_mean = Float64(meta["precipitation_mean"])
    temperature_observation_height = Float64(meta["temperature_observation_height"])
    wind_observation_height = Float64(meta["wind_observation_height"])

    # Optional variables with defaults
    black_carbon_snow = get(meta, "black_carbon_snow", 0.0)
    black_carbon_ice = get(meta, "black_carbon_ice", 0.0)
    cloud_optical_thickness = get(meta, "cloud_optical_thickness", 0.0)
    solar_zenith_angle = get(meta, "solar_zenith_angle", 0.0)
    shortwave_downward_diffuse = get(meta, "shortwave_downward_diffuse", 0.0)
    cloud_fraction = get(meta, "cloud_fraction", 0.1)

    # Call GEMB.initialize_forcing to create ClimateForcing struct
    # This includes GEMB's validation logic
    return GEMB.initialize_forcing(
        time,
        temperature_air,
        pressure_air,
        precipitation,
        wind_speed,
        shortwave_downward,
        longwave_downward,
        vapor_pressure;
        temperature_air_mean = temperature_air_mean,
        wind_speed_mean = wind_speed_mean,
        precipitation_mean = precipitation_mean,
        temperature_observation_height = temperature_observation_height,
        wind_observation_height = wind_observation_height,
        black_carbon_snow = black_carbon_snow,
        black_carbon_ice = black_carbon_ice,
        cloud_optical_thickness = cloud_optical_thickness,
        solar_zenith_angle = solar_zenith_angle,
        shortwave_downward_diffuse = shortwave_downward_diffuse,
        cloud_fraction = cloud_fraction
    )
end

end # module
