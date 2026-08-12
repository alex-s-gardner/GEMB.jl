"""
    forcing_climatology(cf::ClimateForcing)
    forcing_climatology(cf::ClimateForcing, datetime_range::Tuple{DateTime,DateTime})

Compute climatological average forcing from a [`ClimateForcing`](@ref).

Creates a single-year average forcing by:
1. Optionally subsetting to `datetime_range`
2. Eliminating leap days (day 366 of the year)
3. Excluding partial years
4. Averaging every time-series field across the complete years

Returns a new `ClimateForcing` holding one year of climatological forcing. The
time step and scalar metadata (`temperature_air_mean`, `wind_speed_mean`,
`precipitation_mean`, `temperature_observation_height`, `wind_observation_height`)
are carried forward unchanged.

Typically used to build a repeating forcing cycle for [`gemb_spinup`](@ref).

Matches MATLAB's `forcing_climatology.m`.
"""
function forcing_climatology(cf::ClimateForcing, datetime_range::Tuple{DateTime,DateTime})
    times = collect(lookup(dims(cf.temperature_air, Ti)))
    keep = (times .>= datetime_range[1]) .& (times .<= datetime_range[2])

    tkeep = Ti(times[keep])
    _subset(a) = DimArray(parent(a)[keep], (tkeep,))
    cf_subset = ClimateForcing(
        _subset(cf.temperature_air),
        _subset(cf.pressure_air),
        _subset(cf.precipitation),
        _subset(cf.wind_speed),
        _subset(cf.shortwave_downward),
        _subset(cf.longwave_downward),
        _subset(cf.vapor_pressure),
        _subset(cf.black_carbon_snow),
        _subset(cf.black_carbon_ice),
        _subset(cf.cloud_optical_thickness),
        _subset(cf.solar_zenith_angle),
        _subset(cf.shortwave_downward_diffuse),
        _subset(cf.cloud_fraction),
        cf.time_step,
        cf.temperature_air_mean,
        cf.wind_speed_mean,
        cf.precipitation_mean,
        cf.temperature_observation_height,
        cf.wind_observation_height,
    )
    return forcing_climatology(cf_subset)
end

function forcing_climatology(cf::ClimateForcing)
    times = collect(lookup(dims(cf.temperature_air, Ti)))

    # Eliminate leap days (day 366 of the year)
    non_leap = [Dates.dayofyear(t) != 366 for t in times]
    times_noleap = times[non_leap]

    # Count timesteps per year
    years_all = Dates.year.(times_noleap)
    unique_years = sort(unique(years_all))
    counts_per_year = [count(==(yr), years_all) for yr in unique_years]

    # Find years with the maximum number of entries (complete years)
    max_count = maximum(counts_per_year)
    complete_years = unique_years[counts_per_year .== max_count]

    # Indices (within the leap-day-filtered series) for complete years only
    forcing_index = [yr in complete_years for yr in years_all]
    n_complete_years = length(complete_years)
    steps_per_year = max_count

    times_complete = times_noleap[forcing_index]

    # Reshape into (steps_per_year × n_years) and average across years
    reshape_avg(arr) = vec(Statistics.mean(
        reshape(arr[non_leap][forcing_index], steps_per_year, n_complete_years), dims=2))

    # Use times from the first complete year
    tdim = Ti(times_complete[1:steps_per_year])
    _avg(a) = DimArray(reshape_avg(parent(a)), (tdim,))

    return ClimateForcing(
        _avg(cf.temperature_air),
        _avg(cf.pressure_air),
        _avg(cf.precipitation),
        _avg(cf.wind_speed),
        _avg(cf.shortwave_downward),
        _avg(cf.longwave_downward),
        _avg(cf.vapor_pressure),
        _avg(cf.black_carbon_snow),
        _avg(cf.black_carbon_ice),
        _avg(cf.cloud_optical_thickness),
        _avg(cf.solar_zenith_angle),
        _avg(cf.shortwave_downward_diffuse),
        _avg(cf.cloud_fraction),
        cf.time_step,
        cf.temperature_air_mean,
        cf.wind_speed_mean,
        cf.precipitation_mean,
        cf.temperature_observation_height,
        cf.wind_observation_height,
    )
end
