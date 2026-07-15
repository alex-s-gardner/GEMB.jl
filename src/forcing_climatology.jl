"""
    forcing_climatology(cf::ClimateForcing)
    forcing_climatology(cf::ClimateForcing, datetime_range::Tuple{DateTime,DateTime})

Compute climatological average forcing from a ClimateForcing struct.

Creates a single-year average forcing by:
1. Optionally subsetting to `datetime_range`
2. Eliminating leap days (day 366)
3. Excluding partial years
4. Averaging across complete years

Returns a new ClimateForcing with one year of climatological forcing.

Matches MATLAB's `forcing_climatology.m`.
"""
function forcing_climatology(cf::ClimateForcing, datetime_range::Tuple{DateTime,DateTime})
    # Subset to the specified date range
    times = collect(dims(cf.temperature_air, Ti))
    keep = (times .>= datetime_range[1]) .& (times .<= datetime_range[2])

    cf_subset = ClimateForcing(
        DimArray(parent(cf.temperature_air)[keep], (Ti(times[keep]),)),
        DimArray(parent(cf.pressure_air)[keep], (Ti(times[keep]),)),
        DimArray(parent(cf.precipitation)[keep], (Ti(times[keep]),)),
        DimArray(parent(cf.wind_speed)[keep], (Ti(times[keep]),)),
        DimArray(parent(cf.shortwave_downward)[keep], (Ti(times[keep]),)),
        DimArray(parent(cf.longwave_downward)[keep], (Ti(times[keep]),)),
        DimArray(parent(cf.vapor_pressure)[keep], (Ti(times[keep]),)),
        cf.temperature_air_mean,
        cf.wind_speed_mean,
        cf.precipitation_mean,
        cf.temperature_observation_height,
        cf.wind_observation_height,
    )
    return forcing_climatology(cf_subset)
end

function forcing_climatology(cf::ClimateForcing)
    times = collect(dims(cf.temperature_air, Ti))

    # Eliminate leap days (day 366 of the year)
    non_leap = [Dates.dayofyear(t) != 366 for t in times]
    times_noleap = times[non_leap]

    # Extract data arrays without leap days
    ta = parent(cf.temperature_air)[non_leap]
    pa = parent(cf.pressure_air)[non_leap]
    pr = parent(cf.precipitation)[non_leap]
    ws = parent(cf.wind_speed)[non_leap]
    sw = parent(cf.shortwave_downward)[non_leap]
    lw = parent(cf.longwave_downward)[non_leap]
    vp = parent(cf.vapor_pressure)[non_leap]

    # Count timesteps per year
    years_all = Dates.year.(times_noleap)
    unique_years = sort(unique(years_all))
    counts_per_year = [count(==(yr), years_all) for yr in unique_years]

    # Find years with the maximum number of entries (complete years)
    max_count = maximum(counts_per_year)
    complete_mask = counts_per_year .== max_count
    complete_years = unique_years[complete_mask]

    # Get indices for complete years only
    forcing_index = [yr in complete_years for yr in years_all]
    n_complete_years = length(complete_years)
    steps_per_year = max_count

    # Subset to complete years
    ta_complete = ta[forcing_index]
    pa_complete = pa[forcing_index]
    pr_complete = pr[forcing_index]
    ws_complete = ws[forcing_index]
    sw_complete = sw[forcing_index]
    lw_complete = lw[forcing_index]
    vp_complete = vp[forcing_index]
    times_complete = times_noleap[forcing_index]

    # Reshape into (steps_per_year x n_years) and average
    reshape_avg(arr) = vec(Statistics.mean(reshape(arr, steps_per_year, n_complete_years), dims=2))

    clim_ta = reshape_avg(ta_complete)
    clim_pa = reshape_avg(pa_complete)
    clim_pr = reshape_avg(pr_complete)
    clim_ws = reshape_avg(ws_complete)
    clim_sw = reshape_avg(sw_complete)
    clim_lw = reshape_avg(lw_complete)
    clim_vp = reshape_avg(vp_complete)

    # Use times from first complete year
    clim_times = times_complete[1:steps_per_year]
    tdim = Ti(clim_times)

    return ClimateForcing(
        DimArray(clim_ta, (tdim,)),
        DimArray(clim_pa, (tdim,)),
        DimArray(clim_pr, (tdim,)),
        DimArray(clim_ws, (tdim,)),
        DimArray(clim_sw, (tdim,)),
        DimArray(clim_lw, (tdim,)),
        DimArray(clim_vp, (tdim,)),
        cf.temperature_air_mean,
        cf.wind_speed_mean,
        cf.precipitation_mean,
        cf.temperature_observation_height,
        cf.wind_observation_height,
    )
end
