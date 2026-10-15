"""
    air_density(pressure, temperature) -> ρ_air [kg m-3]

Density of air from the ideal gas law, `ρ = 0.029·P/(R·T)`, with the 0.029 kg mol-1
mean molar mass of dry air.

Shared by the transient surface energy balance ([`calculate_temperature`](@ref))
and the initial-guess one (`_seb_annual_melt`), so both use one expression.
"""
@inline air_density(pressure::Real, temperature::Real) =
    0.029 * pressure / (R_GAS * temperature)

"""
    saturation_vapor_pressure_ice(T) -> e_s [Pa]

Saturation vapor pressure over ice, Murphy & Koop (2005) eq. 7, with `T` in kelvin:

    ln(e_s) = 9.550426 − 5723.265/T + 3.53068·ln(T) − 0.00728332·T

valid above 110 K; 611.657 Pa at the triple point.

Shared by the latent-heat flux in [`turbulent_heat_flux`](@ref) (evaluated at the surface
temperature, for a sub-freezing surface) and by [`blowing_snow_sublimation_rate`](@ref)
(evaluated at the air temperature). Companion of
[`saturation_vapor_pressure_water`](@ref), which the flux uses at the melting point.

# References
- Murphy, D. M. and Koop, T. (2005). Review of the vapor pressures of ice and supercooled
  water for atmospheric applications. *Q. J. R. Meteorol. Soc.* 131, 1539-1565. Eq. 7.
"""
@inline saturation_vapor_pressure_ice(T::Real) =
    exp(9.550426 - 5723.265 / T + 3.53068 * log(T) - 0.00728332 * T)

"""
    saturation_vapor_pressure_water(T) -> e_s [Pa]

Saturation vapor pressure over liquid water, supercooled below the triple point, Murphy &
Koop (2005) eq. 10 (p. 1552), with `T` in kelvin:

    ln(e_s) = 54.842763 − 6763.22/T − 4.210·ln(T) + 0.000367·T
              + tanh(0.0415·(T − 218.8))·(53.878 − 1331.22/T − 9.44523·ln(T) + 0.014025·T)

valid for 123 < T < 332 K; outside that range it throws a `DomainError`. Its Table C1
(p. 1561) gives 611.213 Pa at 273.15 K and 611.657 Pa at the triple point, where it meets
[`saturation_vapor_pressure_ice`](@ref) (eq. 7 of the same paper); below the triple point it
lies above the ice curve.

The latent-heat flux ([`turbulent_heat_flux`](@ref)) uses it only for a surface at the melting
point, where the exchanging phase is liquid water; every colder surface exchanges vapor with
ice.

# References
- Murphy, D. M. and Koop, T. (2005). Review of the vapor pressures of ice and supercooled
  water for atmospheric applications. *Q. J. R. Meteorol. Soc.* 131, 1539-1565. Eq. 10.
"""
@inline function saturation_vapor_pressure_water(T::Real)
    123 < T < 332 || _saturation_vapor_pressure_water_domain(T)
    return exp(54.842763 - 6763.22 / T - 4.210 * log(T) + 0.000367 * T +
               tanh(0.0415 * (T - 218.8)) *
               (53.878 - 1331.22 / T - 9.44523 * log(T) + 0.014025 * T))
end

@noinline _saturation_vapor_pressure_water_domain(T) = throw(DomainError(T,
    "Murphy & Koop (2005) eq. 10 for the vapor pressure over liquid water is valid for " *
    "123 < T < 332 K"))

"""
    specific_humidity(vapor_pressure, pressure_air) -> q [kg kg-1]

Specific humidity from vapour and total pressure, `q = εe/(P - (1-ε)e)` with
`ε = MOLAR_MASS_RATIO_VAPOR_AIR`.

Used by [`blowing_snow_sublimation_rate`](@ref), which needs `q` rather than `e` because the
Gordon et al. (2006) sublimation rate is written against the saturation *specific humidity*.
The turbulent latent-heat flux works in vapour pressure directly and does not call this.

# References
- Gordon, M., Simon, K. and Taylor, P. A. (2006). On snow depth predictions with the Canadian
  land surface scheme including a parametrization of blowing snow sublimation.
  *Atmos.-Ocean* 44, 239-255.
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

# References
- Bougamont, M., Bamber, J. L. and Greuell, W. (2005). A surface mass balance model for the
  Greenland Ice Sheet. *J. Geophys. Res.* 110, F04018.
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

Two inline copies of this expression predate the helper —
[`calculate_grain_size`](@ref) (the dendritic branch), and the wind-dependent fresh-snow radius
in [`calculate_accumulation`](@ref). They are left in place because they floor at
`GDN_TOLERANCE` on the radius where this floors at `GDN_TOLERANCE*2` on the diameter, and the
two differ in the degenerate branch only; unifying them is not bit-identical. New callers
should use this.

# References
- Brun, E., David, P., Sudul, M. and Brunot, G. (1992). A numerical model to simulate
  snow-cover stratigraphy for operational avalanche forecasting. *J. Glaciol.* 38, 13-22.
"""
@inline dendritic_grain_radius(dendricity::Real, sphericity::Real) =
    max(0.1 * (dendricity / 0.99 + (1.0 - 1.0 * dendricity / 0.99) *
               (sphericity / 0.99 * 3.0 + (1.0 - sphericity / 0.99) * 4.0)),
        GDN_TOLERANCE * 2) / 2

"""
    fast_divisors(n::Integer)

Find all positive divisors of integer `n`, returned sorted.
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

# References
- Vandecrux, B., et al. (2019). Firn data compilation reveals widespread decrease of firn air
  content in western Greenland. *The Cryosphere* 13, 845-859.
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

Automatically accounts for leap years. The output format is MATLAB's by definition — this is
an interop convenience for callers exchanging dates with MATLAB tooling, not a statement
about the model agreeing with the MATLAB implementation.
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
    saturation_vapor_pressure_ice_derivative(T) -> de_s/dT [Pa K-1]

Temperature derivative of [`saturation_vapor_pressure_ice`](@ref),
`e_s·(5723.265/T² + 3.53068/T − 0.00728332)`.
"""
@inline saturation_vapor_pressure_ice_derivative(T::Real) =
    saturation_vapor_pressure_ice(T) * (5723.265 / T^2 + 3.53068 / T - 0.00728332)

"""
    saturation_vapor_pressure_water_derivative(T) -> de_s/dT [Pa K-1]

Temperature derivative of [`saturation_vapor_pressure_water`](@ref), the analytic derivative
of Murphy & Koop (2005) eq. 10: with `s = 0.0415·(T − 218.8)` and
`g(T) = 53.878 − 1331.22/T − 9.44523·ln(T) + 0.014025·T`,

    de_s/dT = e_s·(6763.22/T² − 4.210/T + 0.000367 + 0.0415·sech²(s)·g(T)
                   + tanh(s)·(1331.22/T² − 9.44523/T + 0.014025))

Valid over the same 123 < T < 332 K, and throws the same `DomainError` outside it.
"""
@inline function saturation_vapor_pressure_water_derivative(T::Real)
    s = 0.0415 * (T - 218.8)
    g = 53.878 - 1331.22 / T - 9.44523 * log(T) + 0.014025 * T
    dlog_es = 6763.22 / T^2 - 4.210 / T + 0.000367 + 0.0415 * g / cosh(s)^2 +
              tanh(s) * (1331.22 / T^2 - 9.44523 / T + 0.014025)
    return saturation_vapor_pressure_water(T) * dlog_es
end
