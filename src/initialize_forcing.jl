"""
    initialize_forcing(time, temperature_air, pressure_air, precipitation,
                       wind_speed, shortwave_downward, longwave_downward,
                       vapor_pressure; kwargs...)

Create a `ClimateForcing` struct from time-series vectors.
Matches MATLAB's `model_initialize_forcing.m`.

# Arguments
- `time::Vector{DateTime}`: Time stamps
- `temperature_air::Vector{Float64}`: Air temperature [K]
- `pressure_air::Vector{Float64}`: Air pressure [Pa]
- `precipitation::Vector{Float64}`: Precipitation per timestep [kg m-2]
- `wind_speed::Vector{Float64}`: Wind speed [m s-1]
- `shortwave_downward::Vector{Float64}`: Incoming shortwave [W m-2]
- `longwave_downward::Vector{Float64}`: Incoming longwave [W m-2]
- `vapor_pressure::Vector{Float64}`: Vapor pressure [Pa]

# Keyword Arguments
- `temperature_air_mean::Float64`: Climatological mean temperature [K]
- `wind_speed_mean::Float64`: Climatological mean wind speed [m s-1]
- `precipitation_mean::Float64`: Climatological mean precipitation [kg m-2 yr-1]
- `temperature_observation_height::Float64`: Height of temperature observation [m]
- `wind_observation_height::Float64`: Height of wind observation [m]
- `accumulation_mean::Float64`: Climatological mean *snowfall* [kg m-2 yr-1], i.e.
  `precipitation_mean` less the rain fraction. Derived from the record by partitioning
  precipitation on `rain_temperature_threshold` when not supplied. Read by
  `calculate_density` under `densification_accumulation = :accumulation`.
- `temperature_air_effective::Float64`: Arrhenius-weighted mean air temperature [K],
  `Eg/(R·log(mean(exp(Eg/(R·T)))))`. Derived from the record when not supplied. Read by
  `calculate_density` under `mean_temperature_method = :arrhenius`.
- `rain_temperature_threshold::Float64`: Rain/snow partition temperature [K] used only to
  derive `accumulation_mean`; keep it equal to the `ModelParameters` field of the same name.
"""
function initialize_forcing(
    time::AbstractVector{DateTime},
    temperature_air::AbstractVector{<:Real},
    pressure_air::AbstractVector{<:Real},
    precipitation::AbstractVector{<:Real},
    wind_speed::AbstractVector{<:Real},
    shortwave_downward::AbstractVector{<:Real},
    longwave_downward::AbstractVector{<:Real},
    vapor_pressure::AbstractVector{<:Real};
    temperature_air_mean::Real=NaN,
    wind_speed_mean::Real=NaN,
    precipitation_mean::Real=NaN,
    temperature_observation_height::Real=NaN,
    wind_observation_height::Real=NaN,
    accumulation_mean::Real=NaN,
    temperature_air_effective::Real=NaN,
    rain_temperature_threshold::Real=273.15,
    black_carbon_snow::Real=0.0,
    black_carbon_ice::Real=0.0,
    cloud_optical_thickness::Real=0.0,
    solar_zenith_angle::Real=0.0,
    shortwave_downward_diffuse::Real=0.0,
    cloud_fraction::Real=0.1,
    dataset::AbstractString="",
    latitude::Real=NaN,
    longitude::Real=NaN,
    elevation::Real=NaN,
    elevation_native::Real=NaN,
    elevation_offset::Real=0.0
)
    # Validate input sizes
    n = length(time)
    @assert length(temperature_air) == n "All input variables must have the same size."
    @assert length(pressure_air) == n "All input variables must have the same size."
    @assert length(precipitation) == n "All input variables must have the same size."
    @assert length(wind_speed) == n "All input variables must have the same size."
    @assert length(shortwave_downward) == n "All input variables must have the same size."
    @assert length(longwave_downward) == n "All input variables must have the same size."
    @assert length(vapor_pressure) == n "All input variables must have the same size."

    # Validate physical ranges
    @assert all(temperature_air .> 100) "temperature_air values unrealistic. Ensure units are kelvin."
    @assert all(pressure_air .>= 0) && all(pressure_air .< 150000) "pressure_air values unrealistic. Ensure units are pascals."
    @assert all(precipitation .>= 0) && all(precipitation .< 20000) "precipitation values unrealistic. Ensure units are kg/m^2 per timestep."
    @assert all(wind_speed .>= 0) && all(wind_speed .< 1000) "wind_speed values unrealistic. Ensure units are m/s."
    @assert all(shortwave_downward .< 10000) "shortwave_downward values unrealistic. Ensure units are W/m^2."
    @assert all(longwave_downward .< 10000) "longwave_downward values unrealistic. Ensure units are W/m^2."
    @assert all(vapor_pressure .>= 0) && all(vapor_pressure .< 150000) "vapor_pressure values unrealistic. Ensure units are Pa."

    # Set defaults for metadata with warnings
    if isnan(temperature_air_mean)
        @warn "Undeclared temperature_air_mean. Assuming mean(temperature_air) represents the climatological mean temperature."
        temperature_air_mean = Statistics.mean(temperature_air)
    end

    if isnan(wind_speed_mean)
        @warn "Undeclared wind_speed_mean. Assuming mean(wind_speed) represents the climatological mean wind speed."
        wind_speed_mean = Statistics.mean(wind_speed)
    end

    if isnan(precipitation_mean)
        @warn "Undeclared precipitation_mean."
        dt_days = Dates.value(time[2] - time[1]) / (1000 * 86400)  # milliseconds to days
        precipitation_mean = Statistics.mean(precipitation) * 365.25 / dt_days
    end

    # Mean *snowfall*: `precipitation_mean` counts rain, which does not bury the column, so
    # densification driven by it overstates the burial rate at any raining site. Partitioned on
    # the same threshold `calculate_accumulation` applies per timestep, and on the same one
    # `initialize_climate_summary` uses, so the transient run and the initializer agree by
    # construction. Derived as a *fraction* of `precipitation_mean` rather than summed
    # independently, so a caller-supplied climatological `precipitation_mean` (which may be
    # scaled differently from this record's own mean) carries its scale over. No warning when
    # undeclared: this is a derived refinement, not a scalar the caller was expected to supply.
    if isnan(accumulation_mean)
        precip_total = 0.0
        snow_total = 0.0
        @inbounds for i in eachindex(precipitation)
            precip_total += precipitation[i]
            if temperature_air[i] <= rain_temperature_threshold
                snow_total += precipitation[i]
            end
        end
        snow_fraction = precip_total > 0 ? snow_total / precip_total : 1.0
        accumulation_mean = precipitation_mean * snow_fraction
    end

    # Arrhenius-weighted mean temperature. `exp(Eg/RT)` is convex in 1/T, so evaluating an
    # Arrhenius factor at the arithmetic mean is not the mean of the factor; this inverts the
    # mean of the exponential instead, which is the temperature at which a rate factor
    # reproduces the record's true mean rate. Exact for every consumer only because all three
    # grain-growth Arrhenius factors share `GRAIN_GROWTH_EG` — a scheme with a different
    # activation energy would need its own effective temperature.
    if isnan(temperature_air_effective)
        acc = 0.0
        @inbounds for T in temperature_air
            acc += exp(GRAIN_GROWTH_EG / (R_GAS * T))
        end
        temperature_air_effective =
            GRAIN_GROWTH_EG / (R_GAS * log(acc / length(temperature_air)))
    end

    if isnan(temperature_observation_height)
        @warn "Undeclared temperature_observation_height. Assuming 2 m above surface."
        temperature_observation_height = 2.0
    end

    if isnan(wind_observation_height)
        @warn "Undeclared wind_observation_height. Assuming 10 m above surface."
        wind_observation_height = 10.0
    end

    # Compute the forcing time step [s] from the (assumed uniform) time axis
    dt_seconds = Dates.value(time[2] - time[1]) / 1000.0  # milliseconds to seconds
    time_step = round(Int, dt_seconds)

    # Create DimArrays with Ti dimension
    tdim = Ti(time)

    return ClimateForcing(
        DimArray(Float64.(temperature_air), (tdim,)),
        DimArray(Float64.(pressure_air), (tdim,)),
        DimArray(Float64.(precipitation), (tdim,)),
        DimArray(Float64.(wind_speed), (tdim,)),
        DimArray(Float64.(shortwave_downward), (tdim,)),
        DimArray(Float64.(longwave_downward), (tdim,)),
        DimArray(Float64.(vapor_pressure), (tdim,)),
        DimArray(Fill(Float64(black_carbon_snow), n), (tdim,)),
        DimArray(Fill(Float64(black_carbon_ice), n), (tdim,)),
        DimArray(Fill(Float64(cloud_optical_thickness), n), (tdim,)),
        DimArray(Fill(Float64(solar_zenith_angle), n), (tdim,)),
        DimArray(Fill(Float64(shortwave_downward_diffuse), n), (tdim,)),
        DimArray(Fill(Float64(cloud_fraction), n), (tdim,)),
        time_step,
        Float64(temperature_air_mean),
        Float64(wind_speed_mean),
        Float64(precipitation_mean),
        Float64(temperature_observation_height),
        Float64(wind_observation_height);
        accumulation_mean=Float64(accumulation_mean),
        temperature_air_effective=Float64(temperature_air_effective),
        dataset=dataset,
        latitude=latitude,
        longitude=longitude,
        elevation=elevation,
        elevation_native=elevation_native,
        elevation_offset=elevation_offset,
    )
end

"""
    initialize_forcing(stack::DimStack) -> ClimateForcing

Create a `ClimateForcing` struct from a forcing `DimStack`.

A `DimStack` is GEMB's neutral, producer-agnostic forcing interface: any source
(e.g. the companion `GEMB_ClimateForcing.jl` package, or a hand-built stack) can
supply forcing as long as it carries the required layers, a `DateTime`-indexed
`Ti` dimension, and the required metadata. This method extracts the forcing
vectors and metadata from `stack` and forwards them to the vector method of
[`initialize_forcing`](@ref), which applies the same validation and normalization.

# Required Variables in DimStack
- `temperature_air`: Air temperature (K)
- `pressure_air`: Surface pressure (Pa)
- `vapor_pressure`: Vapor pressure (Pa)
- `wind_speed`: Wind speed (m/s)
- `precipitation`: Precipitation rate (kg/m²/timestep)
- `shortwave_downward`: Downward shortwave radiation (W/m²)
- `longwave_downward`: Downward longwave radiation (W/m²)

# Required Metadata
- `temperature_air_mean`: Mean air temperature (K)
- `wind_speed_mean`: Mean wind speed (m/s)
- `precipitation_mean`: Annual precipitation (kg/m²/year)
- `temperature_observation_height`: Height of temperature observations (m)
- `wind_observation_height`: Height of wind observations (m)

# Optional Metadata (defaults used if not present)
- `accumulation_mean` (kg/m²/year; derived from the record by rain/snow partitioning)
- `temperature_air_effective` (K; derived from the record as the Arrhenius-weighted mean)
- `black_carbon_snow` (default: 0.0)
- `black_carbon_ice` (default: 0.0)
- `cloud_optical_thickness` (default: 0.0)
- `solar_zenith_angle` (degrees, default: 0.0)
- `shortwave_downward_diffuse` (W/m², default: 0.0)
- `cloud_fraction` (default: 0.1)

# Examples
```julia
using GEMB
using GEMB_ClimateForcing  # one producer of a conforming DimStack

ds = simulate_climate_forcing("test_1", 3)
cf = initialize_forcing(ds)

mp = initialize_parameters()
profile = initialize_profile(mp, cf)
output = gemb(profile, cf, mp)
```
"""
function initialize_forcing(stack::DimStack)
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

    # Optional refinements of the two scalars above. NaN is the "derive it from the record"
    # sentinel the vector method understands, not a missing-data error.
    accumulation_mean = get(meta, "accumulation_mean", NaN)
    temperature_air_effective = get(meta, "temperature_air_effective", NaN)

    # Optional variables with defaults
    black_carbon_snow = get(meta, "black_carbon_snow", 0.0)
    black_carbon_ice = get(meta, "black_carbon_ice", 0.0)
    cloud_optical_thickness = get(meta, "cloud_optical_thickness", 0.0)
    solar_zenith_angle = get(meta, "solar_zenith_angle", 0.0)
    shortwave_downward_diffuse = get(meta, "shortwave_downward_diffuse", 0.0)
    cloud_fraction = get(meta, "cloud_fraction", 0.1)

    # Provenance metadata (optional): source name, cell coordinates, the absolute
    # target `elevation` (m), the source dataset's own `elevation_native` (m, written as
    # `elevation_reanalysis` by `climate_adjust_for_elevation`), and any elevation
    # adjustment applied upstream. Absent for hand-built stacks, so default gracefully.
    # `elevation_offset` supersedes the older `delta_elevation` key that
    # `climate_adjust_for_elevation` used to write.
    dataset = get(meta, "dataset", "")
    latitude = get(meta, "latitude", NaN)
    longitude = get(meta, "longitude", NaN)
    elevation = get(meta, "elevation", NaN)
    elevation_offset = get(meta, "elevation_offset", get(meta, "delta_elevation", 0.0))
    # An unadjusted stack's `elevation` IS its native elevation (no offset applied yet).
    elevation_native = get(meta, "elevation_reanalysis", elevation)

    # Delegate to the vector method, which applies GEMB's validation logic
    return initialize_forcing(
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
        accumulation_mean = accumulation_mean,
        temperature_air_effective = temperature_air_effective,
        black_carbon_snow = black_carbon_snow,
        black_carbon_ice = black_carbon_ice,
        cloud_optical_thickness = cloud_optical_thickness,
        solar_zenith_angle = solar_zenith_angle,
        shortwave_downward_diffuse = shortwave_downward_diffuse,
        cloud_fraction = cloud_fraction,
        dataset = dataset,
        latitude = latitude,
        longitude = longitude,
        elevation = elevation,
        elevation_native = elevation_native,
        elevation_offset = elevation_offset
    )
end
