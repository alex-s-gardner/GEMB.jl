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
    wind_observation_height::Real=NaN
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

    if isnan(temperature_observation_height)
        @warn "Undeclared temperature_observation_height. Assuming 2 m above surface."
        temperature_observation_height = 2.0
    end

    if isnan(wind_observation_height)
        @warn "Undeclared wind_observation_height. Assuming 10 m above surface."
        wind_observation_height = 10.0
    end

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
        Float64(temperature_air_mean),
        Float64(wind_speed_mean),
        Float64(precipitation_mean),
        Float64(temperature_observation_height),
        Float64(wind_observation_height)
    )
end
