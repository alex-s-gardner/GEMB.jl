"""
    air_density(pressure, temperature) -> ρ_air [kg m-3]

Density of air from the ideal gas law, `ρ = 0.029·P/(R·T)`, with the 0.029 kg mol-1
mean molar mass of dry air the MATLAB implementation uses.

Shared by the transient surface energy balance ([`calculate_temperature`](@ref))
and the initial-guess one (`_seb_annual_melt`), so both use one expression.
"""
@inline air_density(pressure::Real, temperature::Real) =
    0.029 * pressure / (R_GAS * temperature)

"""
    saturation_vapor_pressure_ice(T) -> e_s [Pa]

Saturation vapour pressure over ice, Bolton (1980), `610.78·exp(21.8745584·(T-273.16)/(T-7.66))`
with `T` in kelvin.

Shared by the latent-heat flux in [`turbulent_heat_flux`](@ref) (evaluated at the surface
temperature, for a sub-freezing surface) and by [`blowing_snow_sublimation_rate`](@ref)
(evaluated at the air temperature). Companion of
[`saturation_vapor_pressure_water`](@ref), which the flux uses above freezing.

# References
- Bolton, D. (1980). The computation of equivalent potential temperature.
  *Mon. Weather Rev.* 108, 1046-1053.
"""
@inline saturation_vapor_pressure_ice(T::Real) =
    610.78 * exp(21.8745584 * (T - CtoK - 0.01) / (T - 7.66))

"""
    saturation_vapor_pressure_water(T) -> e_s [Pa]

Saturation vapour pressure over liquid water, Murray (1967),
`610.78·exp(17.2693882·(T-273.16)/(T-35.86))` with `T` in kelvin. The above-freezing
counterpart of [`saturation_vapor_pressure_ice`](@ref).

# References
- Murray, F. W. (1967). On the computation of saturation vapor pressure.
  *J. Appl. Meteorol.* 6, 203-204.
"""
@inline saturation_vapor_pressure_water(T::Real) =
    610.78 * exp(17.2693882 * (T - CtoK - 0.01) / (T - 35.86))

"""
    specific_humidity(vapor_pressure, pressure_air) -> q [kg kg-1]

Specific humidity from vapour and total pressure, `q = εe/(P - (1-ε)e)` with
`ε = MOLAR_MASS_RATIO_VAPOR_AIR`.

Used by [`blowing_snow_sublimation_rate`](@ref), which needs `q` rather than `e` because the
Gordon et al. (2006) sublimation rate is written against the saturation *specific humidity*.
The turbulent latent-heat flux works in vapour pressure directly and does not call this.
"""
@inline specific_humidity(vapor_pressure::Real, pressure_air::Real) =
    MOLAR_MASS_RATIO_VAPOR_AIR * vapor_pressure /
    (pressure_air - (1.0 - MOLAR_MASS_RATIO_VAPOR_AIR) * vapor_pressure)

"""
    surface_roughness(mp, density_surface, water_surface) -> z0 [m]

Aerodynamic roughness length of the surface, after Bougamont et al. (2005): dry snow, wet
snow, and bare ice each get their own value (`Z0_SNOW_DRY`, `Z0_SNOW_WET`, `Z0_ICE`), with
"ice" meaning a surface cell at `mp.density_ice` and "wet" meaning one holding liquid water.

Read by the turbulent fluxes in [`calculate_temperature`](@ref) and by the 5 m gust speed in
[`drift_wind_speed`](@ref), so both see the same surface.
"""
@inline function surface_roughness(mp::ModelParameters, density_surface::Real,
    water_surface::Real)
    if (density_surface < (mp.density_ice - D_TOLERANCE)) && (water_surface < W_TOLERANCE)
        return Z0_SNOW_DRY   # 0.12 mm for dry snow
    elseif density_surface >= (mp.density_ice - D_TOLERANCE)
        return Z0_ICE        # 3.2 mm for ice
    else
        return Z0_SNOW_WET   # 1.3 mm for wet snow
    end
end

"""
    dendritic_grain_radius(dendricity, sphericity) -> r [mm]

Effective grain radius of dendritic snow diagnosed from its dendricity and sphericity, Brun
et al. (1992) as GEMB carries it: grain *size* (diameter) is not independent state while
`dendricity > 0`, so any scheme that changes `(δ, S)` must re-derive the radius.

Three inline copies of this expression predate the helper —
[`calculate_grain_size`](@ref) (the dendritic branch), and two in
[`calculate_accumulation`](@ref) (fresh snow). They are left in place because they floor at
`GDN_TOLERANCE` on the radius where this floors at `GDN_TOLERANCE*2` on the diameter, and the
two differ in the degenerate branch only; unifying them is not bit-identical. New callers
should use this.
"""
@inline dendritic_grain_radius(dendricity::Real, sphericity::Real) =
    max(0.1 * (dendricity / 0.99 + (1.0 - 1.0 * dendricity / 0.99) *
               (sphericity / 0.99 * 3.0 + (1.0 - sphericity / 0.99) * 4.0)),
        GDN_TOLERANCE * 2) / 2

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
    firn_air_content(dz, density, density_ice, z_max=Inf) -> m of air

Height of air in the column above depth `z_max`: `Σ dz (1 − ρ/ρ_ice)`, with densities above
`density_ice` contributing nothing rather than negative air.

The cell straddling `z_max` contributes in proportion to the fraction of it that lies above
the cutoff, so the result is a true integral to `z_max` and does not step as cells cross the
boundary during a run. `z_max = Inf` (the default) integrates the whole column.

Depth-limited totals — conventionally to 10 m or 20 m — are what the firn-core literature
reports (e.g. Vandecrux et al., 2019), so they are what a whole-column value cannot be
compared against.

Allocation-free: called every timestep from the `gemb` time loop.
"""
function firn_air_content(dz::AbstractVector, density::AbstractVector,
    density_ice::Real, z_max::Real=Inf)

    fac = 0.0
    z = 0.0
    @inbounds for i in eachindex(dz)
        z >= z_max && break
        # Thickness of this cell lying above the cutoff — the whole cell for every cell when
        # `z_max` is Inf. Written as a clipped thickness rather than a fraction so a
        # zero-thickness cell cannot produce a 0/0.
        dz_in = min(dz[i], z_max - z)
        fac += dz_in * (density_ice - min(density[i], density_ice))
        z += dz[i]
    end
    return fac / density_ice
end

"""
    close_off_age(density, age, closeoff_density=DENSITY_PORE_CLOSEOFF) -> d

Age [d] of the shallowest cell at or above pore close-off density, or `NaN` when the column
never closes off.

This is the model's analogue of the gas age at close-off that firn-air and ice-core work
reports: the residence time of the firn at the depth where the pore space seals and the
column stops exchanging with the atmosphere. `NaN` (rather than the column's total depth or
its deepest age) is the correct answer for an open column, following the same convention as
`ice_slab_depth` — there is no close-off depth to report an age at.

Meaningful from the first timestep only under `mp.initialize_age = :steady_state`; under
`:zero` it reads 0 until a spinup long enough to bury the epoch has run.

Allocation-free: called every timestep from the `gemb` time loop.
"""
function close_off_age(density::AbstractVector, age::AbstractVector,
    closeoff_density::Real=DENSITY_PORE_CLOSEOFF)

    @inbounds for i in eachindex(density)
        if density[i] >= closeoff_density - D_TOLERANCE
            return age[i]
        end
    end
    return NaN
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

