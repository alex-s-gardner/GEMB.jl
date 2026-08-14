"""
    fast_divisors(n::Integer)

Find all positive divisors of integer `n`, returned sorted.
Matches MATLAB's `fast_divisors.m`.
"""
function fast_divisors(n::Integer)
    k = 1:ceil(Int, sqrt(n))
    d = k[rem.(n, k) .== 0]
    return sort(unique(vcat(d, div.(n, d))))
end

"""
    dz2z(dz::AbstractVector)

Convert layer thicknesses `dz` to cell center heights (negative below surface).
The surface is at z=0; centers are at negative depths.

Matches MATLAB's `dz2z.m` for vector input.
"""
function dz2z(dz::AbstractVector)
    z_center = -cumsum(dz) .+ dz[1] / 2
    # Note: For vector input, surface_timeseries equivalent is just dz[1]
    # since there's only one column
    return z_center
end

"""
    dz2z(dz::AbstractMatrix)

Convert 2D layer thickness matrix to cell center heights.
Each column is an independent profile, surface first (row 1).

GEMB profile output is top-justified with a fixed row count, so this is a plain cumulative
sum with no padding to work around. Any NaN present in the input (e.g. an output array not
yet fully written) propagates through its own column from that row down, as `cumsum`
implies.

Matches MATLAB's `dz2z.m` for matrix input.
"""
function dz2z(dz::AbstractMatrix)
    z_center = similar(dz, Float64)
    @inbounds for j in axes(dz, 2)
        half_top = dz[1, j] / 2
        cs = 0.0
        for i in axes(dz, 1)
            cs += dz[i, j]
            z_center[i, j] = -cs + half_top
        end
    end
    return z_center
end

"""
    surface_timeseries(A::AbstractMatrix)

Return the surface (row 1) value in each column of matrix `A`.

GEMB profile output is top-justified, so the surface cell is always row 1. Retained as a
named function because it expresses intent at call sites and matches MATLAB's
`surface_timeseries.m`.
"""
function surface_timeseries(A::AbstractMatrix)
    return Float64[A[1, j] for j in axes(A, 2)]
end

"""
    datetime2decyear(datetimes::AbstractVector{DateTime})

Convert Julia `DateTime` objects to decimal year values.

Automatically accounts for leap years. Inverse of the year-fraction convention
used by [`decyear2datenum`](@ref).
"""
function datetime2decyear(datetimes::AbstractVector{DateTime})
    dec_year = Vector{Float64}(undef, length(datetimes))
    for i in eachindex(datetimes)
        yr = year(datetimes[i])
        start_of_year = DateTime(yr, 1, 1)
        start_of_next_year = DateTime(yr + 1, 1, 1)
        days_in_year = Dates.value(start_of_next_year - start_of_year) / (1000.0 * 86400.0)
        day_offset = Dates.value(datetimes[i] - start_of_year) / (1000.0 * 86400.0)
        dec_year[i] = yr + day_offset / days_in_year
    end
    return dec_year
end

"""
    decyear2datenum(decyear)

Convert decimal year to MATLAB datenum format (days since 0000-01-01).

Automatically accounts for leap years.

Matches MATLAB's `decyear2datenum.m`.
"""
function decyear2datenum(decyear)
    year_part = floor.(Int, decyear)
    fractional_year = decyear .- year_part

    # Calculate start of year and next year
    start_of_year = DateTime.(year_part, 1, 1)
    start_of_next_year = DateTime.(year_part .+ 1, 1, 1)

    # Days in this year
    days_in_year = (start_of_next_year .- start_of_year) ./ Millisecond(1000 * 60 * 60 * 24)

    # Convert to MATLAB datenum (rata die + 1)
    rata_start = datetime2rata.(start_of_year)
    datenum_out = rata_start .+ 1 .+ (fractional_year .* days_in_year)

    return datenum_out
end

"""
    annual_pdd_melt(cf::ClimateForcing; ddf_snow=3.0)

Estimate the annual potential melt [kg m-2 yr-1] from a positive-degree-day (PDD)
scheme applied to the air-temperature forcing.

`PDD = Σ max(T_air − 273.15, 0) · Δt_days` is summed over the whole forcing series
and annualized by `365.25 / total_days` so forcing that is not exactly one year
still yields an annual rate. The potential melt is `ddf_snow · PDD`, where
`ddf_snow` is the snow degree-day factor [mm w.e. °C-1 d-1] (default 3.0, a
standard value for snow).

Used by [`initialize_profile`](@ref) to decide whether a site accumulates firn
(build a Herron–Langway steady-state profile) or ablates (initialize pure ice).
"""
function annual_pdd_melt(cf::ClimateForcing; ddf_snow::Real=3.0)
    T = Float64.(parent(cf.temperature_air))         # [K]
    dt_days = cf.time_step / 86400.0                 # timestep length [d]

    pdd = 0.0
    @inbounds for t in T
        pd = t - CtoK
        if pd > 0.0
            pdd += pd * dt_days                       # [°C · d]
        end
    end

    total_days = length(T) * dt_days
    total_days <= 0 && return 0.0
    pdd_annual = pdd * (365.25 / total_days)          # [°C · d yr-1]

    return Float64(ddf_snow) * pdd_annual             # [kg m-2 yr-1]
end
