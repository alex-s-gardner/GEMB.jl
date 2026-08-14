"""
    calculate_albedo(dz, density, water, grain_radius, melt_surface, cfs::ClimateForcingStep, mp::ModelParameters)

Calculate snow, firn, and ice albedo using one of several methods:
1. "GardnerSharp": function of effective grain radius (Gardner & Sharp, 2010)
2. "BrunLefebre": function of effective grain radius (Lefebre et al., 2003)
3. "GreuellKonzelmann": function of density and cloud amount (Greuell & Konzelmann, 1994)

Returns `(albedo, albedo_diffuse)` as scalars — the broadband surface albedo for the
direct-beam and for the diffuse stream.

Every method here is **diagnostic**: the albedo is a function of the current column
state alone, so nothing is carried across a timestep boundary and nothing is stored per
cell. `albedo_diffuse` differs from `albedo` only under `:GardnerSharp` (recomputed at a
fixed 50° effective zenith angle), which is also the only method whose consumer,
[`calculate_shortwave_radiation`](@ref), splits the direct and diffuse streams; the
other methods return it equal to `albedo`.

# References
- Gardner, A. S. and Sharp, M. J. (2010). J. Geophys. Res., 115, F01009.
- Lefebre, F., et al. (2003). J. Geophys. Res., 108, 4231.
- Greuell, W. and Konzelmann, T. (1994). Global Planet. Change, 9, 91-114.
"""
function calculate_albedo(dz::Vector{Float64},
    density::Vector{Float64}, water::Vector{Float64},
    grain_radius::Vector{Float64}, melt_surface::Float64,
    cfs::ClimateForcingStep, mp::ModelParameters)

    T_tolerance = 1e-10
    d_tolerance = 1e-11

    # constants
    density_fresh_snow = 300.0         # density of fresh snow [kg m-3]
    density_phc = DENSITY_PORE_CLOSEOFF  # Pore closeoff density
    albedo_ice_max = 0.58              # maximum ice albedo, from Lefebre, 2003
    albedo_ice_min = mp.albedo_ice     # minimum ice albedo
    albedo_snow_min = 0.65             # minimum snow albedo, from Alexander 2014

    albedo = 0.0
    albedo_diffuse = 0.0

    if (mp.albedo_method == :None) || ((mp.albedo_density_threshold - density[1]) < d_tolerance)
        albedo = mp.albedo_fixed
        albedo_diffuse = albedo
    else
        if mp.albedo_method == :GardnerSharp
            albedo = _albedo_gardner(grain_radius, dz, density,
                cfs.black_carbon_snow, cfs.black_carbon_ice,
                cfs.solar_zenith_angle, cfs.cloud_optical_thickness)
            albedo_diffuse = _albedo_gardner(grain_radius, dz, density,
                cfs.black_carbon_snow, cfs.black_carbon_ice,
                50.0, cfs.cloud_optical_thickness)

        elseif mp.albedo_method == :BrunLefebre
            # Spectral fractions (Lefebre et al., 2003)
            # [0.3-0.8um 0.8-1.5um 1.5-2.8um]
            sF = [0.606, 0.301, 0.093]

            # convert effective radius to grain size in meters
            gsz = (grain_radius[1] * 2.0) / 1000.0

            # spectral range:
            a1 = min(0.98, 0.95 - 1.58 * gsz^0.5)
            a2 = max(0.0, 0.95 - 15.4 * gsz^0.5)
            a3 = max(0.127, 0.88 + 346.3 * gsz - 32.31 * gsz^0.5)

            # broadband surface albedo
            albedo = sF[1] * a1 + sF[2] * a2 + sF[3] * a3
            albedo_diffuse = albedo

        elseif mp.albedo_method == :GreuellKonzelmann
            albedo = mp.albedo_ice + (density[1] - mp.density_ice) *
                (mp.albedo_snow - mp.albedo_ice) /
                (density_fresh_snow - mp.density_ice) +
                (0.05 * (cfs.cloud_fraction - 0.5))
            albedo_diffuse = albedo
        end

        # If we do not have fresh snow
        if (mp.albedo_method == :GardnerSharp || mp.albedo_method == :BrunLefebre) &&
           ((mp.albedo_density_threshold - density[1]) >= d_tolerance)

            # In a snow layer < 10cm, account for mix of ice and snow
            lice_first = something(findfirst(d -> d >= density_phc - d_tolerance, density), length(density) + 1)
            depthsnow = 0.0
            @inbounds @simd for i in 1:(lice_first-1)
                depthsnow += dz[i]
            end

            if (depthsnow <= (0.1 + d_tolerance)) && (lice_first <= length(density)) &&
               (density[lice_first] >= (density_phc - d_tolerance))

                aice = albedo_ice_max + (albedo_snow_min - albedo_ice_max) *
                    (density[lice_first] - mp.density_ice) /
                    (density_phc - mp.density_ice)

                albedo = aice + max(albedo - aice, 0.0) * (depthsnow / 0.1)
            end

            if (density[1] >= density_phc - d_tolerance)
                if (density[1] < mp.density_ice - d_tolerance)
                    albedo = albedo_ice_max + (albedo_snow_min - albedo_ice_max) *
                        (density[1] - mp.density_ice) / (density_phc - mp.density_ice)
                else
                    M = melt_surface + water[1]
                    albedo = max(albedo_ice_min + (albedo_ice_max - albedo_ice_min) * exp(-1.0 * (M / 200.0)), albedo_ice_min)
                end
            end
        end
    end

    # Check for erroneous values
    if albedo > (1 + T_tolerance)
        @warn "albedo > 1.0"
    elseif albedo < (0 - d_tolerance)
        @warn "albedo is negative"
    elseif isnan(albedo)
        error("albedo == NAN")
    end

    return albedo, albedo_diffuse
end

"""
    _albedo_gardner(grain_radius, dz, density, c1, c2, SZA, t)

Broadband albedo parameterization from Gardner and Sharp (2010).
Accounts for grain size, soot loading, solar zenith angle, and cloud optical thickness.
Two-layer parameterization is applied when an ice layer exists below snow.
"""
function _albedo_gardner(grain_radius::Vector{Float64}, dz::Vector{Float64},
    density::Vector{Float64}, c1::Float64, c2::Float64,
    SZA::Float64, t::Float64)

    d_tolerance = 1e-11

    # convert effective radius to specific surface area [cm2 g-1]
    S1 = 3.0 / (0.091 * grain_radius[1])

    # effective solar zenith angle
    x = min((t / (3 * cos(pi * SZA / 180)))^0.5, 1.0)
    u = 0.64 * x + (1 - x) * cos(pi * SZA / 180)

    # pure snow albedo
    as = 1.48 - S1^(-0.07)

    # change in pure snow albedo due to soot loading
    dac = max(0.04 - as,
        -(c1^0.55) / (0.16 + 0.6 * S1^0.5 + (1.8 * c1^0.6) * (S1^(-0.25))))

    # Two layer albedo parameterization
    lice_first = something(findfirst(d -> d >= DENSITY_PORE_CLOSEOFF - d_tolerance, density),
        length(density) + 1)
    z1 = 0.0
    @inbounds @simd for i in 1:(lice_first-1)
        z1 += dz[i] * density[i]
    end

    m = length(density)
    if (m > 0) && (lice_first <= m) && (z1 > d_tolerance)
        # determine albedo values for bottom layer
        S2 = 3.0 / (0.091 * grain_radius[lice_first])

        # pure snow albedo
        as2 = 1.48 - S2^(-0.07)

        # change in pure snow albedo due to soot loading
        dac2 = max(0.04 - as2,
            -(c2^0.55) / (0.16 + 0.6 * S2^0.5 + (1.8 * c2^0.6) * (S2^(-0.25))))

        # determine the effective change due to finite depth and soot loading
        A = min(1.0, (2.1 * z1^(1.35 * (1 - as) - 0.1 * c1 - 0.13)))

        dac = (as2 + dac2 - as) + A * ((as + dac) - (as2 + dac2))
    end

    # change in albedo due to solar zenith angle
    dasz = 0.53 * as * (1 - (as + dac)) * (1 - u)^1.2

    # change in albedo due to cloud (apart from change in diffuse fraction)
    dat = (0.1 * t * (as + dac)^1.3) / ((1 + 1.5 * t)^as)

    # Broadband albedo
    albedo = as + dac + dasz + dat

    return albedo
end
