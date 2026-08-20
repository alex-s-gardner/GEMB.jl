"""
    forcing_climatology(cf::ClimateForcing)
    forcing_climatology(cf::ClimateForcing, datetime_range::Tuple{DateTime,DateTime})
    forcing_climatology(cf::ClimateForcing, selector)

Compute climatological average forcing from a [`ClimateForcing`](@ref).

Creates a single-year average forcing by:
1. Optionally subsetting to a time window (`datetime_range` or `selector`)
2. Eliminating leap days (day 366 of the year)
3. Excluding partial years
4. Averaging every time-series field across the complete years

Returns a new `ClimateForcing` holding one year of climatological forcing. The
time step and scalar metadata (`temperature_air_mean`, `wind_speed_mean`,
`precipitation_mean`, `temperature_observation_height`, `wind_observation_height`)
are carried forward unchanged.

The window can be given as a `(start, stop)` `DateTime` tuple, or — since
`ClimateForcing` is an `AbstractDimStack` — as any DimensionalData `Ti` selector,
e.g. a closed interval `DateTime(1950,1,1) .. DateTime(1980,12,31)`, `At`, `Near`,
or an index range. Both forms subset via `cf[Ti(...)]` before averaging.

Typically used to build a repeating forcing cycle for [`gemb_spinup`](@ref).

Matches MATLAB's `forcing_climatology.m`.
"""
function forcing_climatology(cf::ClimateForcing, datetime_range::Tuple{DateTime,DateTime})
    # Record the requested window as climatology provenance on the result.
    return forcing_climatology(cf[Ti(datetime_range[1] .. datetime_range[2])];
                               window=datetime_range)
end

# Idiomatic DimensionalData form: forward any `Ti` selector (interval, `At`, `Near`,
# index range, ...) through the stack's own indexing, then average. The tuple method
# above is strictly more specific, so it always wins for `(start, stop)` calls.
function forcing_climatology(cf::ClimateForcing, selector)
    return forcing_climatology(cf[Ti(selector)])
end

# `window` is the requested `(start, stop)` averaging window recorded as
# provenance; when `nothing` (e.g. called on an already-subset stack) it falls
# back to the extent of the complete years actually averaged.
function forcing_climatology(cf::ClimateForcing; window=nothing)
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

    # Climatology provenance: the requested window (or, when not given, the extent
    # of the complete years actually averaged) plus the year / step counts.
    window_start = window === nothing ? first(times_complete) : window[1]
    window_stop  = window === nothing ? last(times_complete)  : window[2]

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
        cf.wind_observation_height;
        snow_drift=_avg(cf.snow_drift),
        accumulation_mean=cf.accumulation_mean,
        temperature_air_effective=cf.temperature_air_effective,
        dataset=cf.dataset,
        latitude=cf.latitude,
        longitude=cf.longitude,
        elevation=cf.elevation,
        elevation_native=cf.elevation_native,
        elevation_offset=cf.elevation_offset,
        climatology_window_start=window_start,
        climatology_window_stop=window_stop,
        climatology_n_years=n_complete_years,
        climatology_steps_per_year=steps_per_year,
    )
end
