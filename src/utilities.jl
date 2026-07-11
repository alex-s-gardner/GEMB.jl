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
Each column is an independent profile. NaN values are preserved.

Matches MATLAB's `dz2z.m` for matrix input.
"""
function dz2z(dz::AbstractMatrix)
    nrows, ncols = size(dz)
    z_center = similar(dz)
    isn = isnan.(dz)

    # Replace NaN with 0 for cumsum, then restore
    dz_clean = copy(dz)
    dz_clean[isn] .= 0.0

    # cumsum along columns (dim=1)
    cs = cumsum(dz_clean, dims=1)

    # Get top finite value for each column (surface_timeseries equivalent)
    top_vals = surface_timeseries(dz)

    for j in axes(dz, 2)
        for i in axes(dz, 1)
            z_center[i, j] = -cs[i, j] + top_vals[j] / 2
        end
    end
    z_center[isn] .= NaN

    return z_center
end

"""
    surface_timeseries(A::AbstractMatrix)

Return the uppermost finite value in each column of matrix `A`.
Matches MATLAB's `surface_timeseries.m`.
"""
function surface_timeseries(A::AbstractMatrix)
    ncols = size(A, 2)
    result = fill(NaN, ncols)
    for j in axes(A, 2)
        for i in axes(A, 1)
            if isfinite(A[i, j])
                result[j] = A[i, j]
                break
            end
        end
    end
    return result
end

"""
    dewpoint_to_vapor_pressure(temperature_dewpoint)

Convert dewpoint temperature [K] to actual vapor pressure [Pa].

Uses the NOAA formula: vapor_pressure = 611 * 10^(7.5 * Td / (237.3 + Td))
where Td is dewpoint in Celsius.

Matches MATLAB's `dewpoint_to_vapor_pressure.m`.
"""
function dewpoint_to_vapor_pressure(temperature_dewpoint)
    Td = temperature_dewpoint .- 273.15
    return 611.0 .* (10.0 .^ (7.5 .* Td ./ (237.3 .+ Td)))
end

"""
    vapor_pressure_to_relative_humidity(vapor_pressure, temperature_air)

Calculate relative humidity [%] from vapor pressure [Pa] and air temperature [K].

Uses Tetens' formula for saturation vapor pressure and clamps the result to [0, 100].

Matches MATLAB's `vapor_pressure_to_relative_humidity.m`.
"""
function vapor_pressure_to_relative_humidity(vapor_pressure, temperature_air)
    Tc = temperature_air .- 273.15
    A = 610.78
    B = 17.27
    C = 237.3
    es = A .* exp.((B .* Tc) ./ (Tc .+ C))
    relative_humidity = (vapor_pressure ./ es) .* 100.0
    return clamp.(relative_humidity, 0.0, 100.0)
end
