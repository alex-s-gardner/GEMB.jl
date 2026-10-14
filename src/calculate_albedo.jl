"""
    calculate_albedo(dz, density, grain_radius, cfs::ClimateForcingStep, mp::ModelParameters)

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

# Firn, ice, and thin snow over them

Under `:GardnerSharp` and `:BrunLefebre` the grain-size albedo describes semi-infinite snow.
Where the column has a surface of firn or ice, or a thin snow layer over one, both returned
albedos are replaced as follows ([`firn_ice_albedo`](@ref)):

- **Exposed ice** (`density[1] ≥ mp.density_ice`): `mp.albedo_ice` exactly.
- **Exposed firn** (`DENSITY_PORE_CLOSEOFF ≤ density[1] < mp.density_ice`): linear in density
  from `mp.albedo_ice` at `mp.density_ice` to `ALBEDO_FIRN_MAX` (0.65) at pore close-off.
- **Thin snow** (snow of total thickness `d ≤ ALBEDO_SNOW_DEPTH_THICK` = 0.1 m above the first
  cell at or above pore close-off): linear in depth between the albedo of that underlying firn
  or ice surface at `d = 0` and the semi-infinite snow albedo at `d = 0.1 m`.

The depth blend is the only thin-snow mechanism. Gardner & Sharp's (2010) two-layer
finite-depth term is not applied.

# References
- Gardner, A. S. and Sharp, M. J. (2010). A review of snow and ice albedo and the development
  of a new physically based broadband albedo parameterization. *J. Geophys. Res.* 115, F01009.
- Lefebre, F., Gallée, H., van Ypersele, J.-P., and Greuell, W. (2003). Modeling of snow and
  ice melt at ETH Camp (West Greenland): a study of surface albedo. *J. Geophys. Res.* 108,
  4231.
- Greuell, W. and Konzelmann, T. (1994). Numerical modelling of the energy balance and the
  englacial temperature of the Greenland ice sheet. *Global Planet. Change* 9, 91-114.
- Alexander, P. M., Tedesco, M., Fettweis, X., van de Wal, R. S. W., Smeets, C. J. P. P., and
  van den Broeke, M. R. (2014). Assessing spatio-temporal variability and trends in modelled
  and measured Greenland Ice Sheet albedo (2000-2013). *The Cryosphere* 8, 2293-2312. (Source
  of the 0.65 minimum snow albedo used as `ALBEDO_FIRN_MAX`.)
"""
function calculate_albedo(dz::Vector{Float64},
    density::Vector{Float64}, grain_radius::Vector{Float64},
    cfs::ClimateForcingStep, mp::ModelParameters)

    density_fresh_snow = 300.0         # density of fresh snow [kg m-3]

    albedo = 0.0
    albedo_diffuse = 0.0

    if (mp.albedo_method == :None) || ((mp.albedo_density_threshold - density[1]) < D_TOLERANCE)
        albedo = mp.albedo_fixed
        albedo_diffuse = albedo
    else
        if mp.albedo_method == :GardnerSharp
            albedo = _albedo_gardner(grain_radius[1], cfs.black_carbon_snow,
                cfs.solar_zenith_angle, cfs.cloud_optical_thickness)
            albedo_diffuse = _albedo_gardner(grain_radius[1], cfs.black_carbon_snow,
                50.0, cfs.cloud_optical_thickness)

        elseif mp.albedo_method == :BrunLefebre
            albedo = _albedo_brun_lefebre(grain_radius[1])
            albedo_diffuse = albedo

        elseif mp.albedo_method == :GreuellKonzelmann
            albedo = mp.albedo_ice + (density[1] - mp.density_ice) *
                (mp.albedo_snow - mp.albedo_ice) /
                (density_fresh_snow - mp.density_ice) +
                (0.05 * (cfs.cloud_fraction - 0.5))
            albedo_diffuse = albedo
        end

        # The grain-size methods describe semi-infinite snow; replace them over firn, ice, or
        # thin snow on either (see the docstring).
        if (mp.albedo_method == :GardnerSharp) || (mp.albedo_method == :BrunLefebre)
            albedo = _underlying_surface_albedo(albedo, dz, density, mp)
            albedo_diffuse = _underlying_surface_albedo(albedo_diffuse, dz, density, mp)
        end
    end

    # Check for erroneous values
    if albedo > (1 + T_TOLERANCE)
        @warn "albedo > 1.0"
    elseif albedo < (0 - D_TOLERANCE)
        @warn "albedo is negative"
    elseif isnan(albedo)
        error("albedo == NAN")
    end

    return albedo, albedo_diffuse
end

"""
    firn_ice_albedo(density, mp) -> α

Albedo of an exposed firn or ice surface of the given `density` [kg m-3]: `mp.albedo_ice` at
or above `mp.density_ice`, and linear in density from there up to `ALBEDO_FIRN_MAX` at
`DENSITY_PORE_CLOSEOFF`. Densities below pore close-off are snow and throw.
"""
@inline function firn_ice_albedo(density::Float64, mp::ModelParameters)
    density >= DENSITY_PORE_CLOSEOFF - D_TOLERANCE ||
        error("firn_ice_albedo: density $(density) kg m-3 is below pore close-off " *
              "($(DENSITY_PORE_CLOSEOFF) kg m-3), so the surface is snow, not firn or ice")
    density >= mp.density_ice - D_TOLERANCE && return mp.albedo_ice
    return mp.albedo_ice + (ALBEDO_FIRN_MAX - mp.albedo_ice) *
        (mp.density_ice - density) / (mp.density_ice - DENSITY_PORE_CLOSEOFF)
end

"""
    _underlying_surface_albedo(albedo_snow, dz, density, mp) -> α

Apply the firn, ice and thin-snow rules of [`calculate_albedo`](@ref) to the semi-infinite
snow albedo `albedo_snow`.
"""
@inline function _underlying_surface_albedo(albedo_snow::Float64, dz::Vector{Float64},
    density::Vector{Float64}, mp::ModelParameters)

    density[1] >= DENSITY_PORE_CLOSEOFF - D_TOLERANCE && return firn_ice_albedo(density[1], mp)

    i_firn = findfirst(d -> d >= DENSITY_PORE_CLOSEOFF - D_TOLERANCE, density)
    i_firn === nothing && return albedo_snow

    depth_snow = 0.0
    for i in firstindex(dz):(i_firn - 1)
        depth_snow += dz[i]
    end
    depth_snow >= ALBEDO_SNOW_DEPTH_THICK && return albedo_snow

    albedo_under = firn_ice_albedo(density[i_firn], mp)
    return albedo_under + (albedo_snow - albedo_under) * (depth_snow / ALBEDO_SNOW_DEPTH_THICK)
end

"""
    _albedo_brun_lefebre(grain_radius_surface) -> α

Broadband albedo of semi-infinite snow from three spectral-band albedos weighted by the
spectral fractions of Lefebre et al. (2003), `[0.606, 0.301, 0.093]` for 0.3–0.8, 0.8–1.5 and
1.5–2.8 µm. The fractions sum to 1.
"""
@inline function _albedo_brun_lefebre(grain_radius_surface::Float64)
    # grain size (diameter) in metres
    gsz = (grain_radius_surface * 2.0) / 1000.0
    a1 = min(0.98, 0.95 - 1.58 * gsz^0.5)
    a2 = max(0.0, 0.95 - 15.4 * gsz^0.5)
    a3 = max(0.127, 0.88 + 346.3 * gsz - 32.31 * gsz^0.5)
    return 0.606 * a1 + 0.301 * a2 + 0.093 * a3
end

"""
    _albedo_gardner(grain_radius_surface, c, SZA, τ) -> α

Broadband albedo of semi-infinite snow from Gardner and Sharp (2010): pure-snow albedo from
the specific surface area, plus the changes due to soot loading `c` [ppmw], solar zenith angle
`SZA` [degrees] and cloud optical thickness `τ`.
"""
function _albedo_gardner(grain_radius_surface::Float64, c::Float64,
    SZA::Float64, t::Float64)

    # convert effective radius to specific surface area [cm2 g-1]
    S1 = 3.0 / (0.091 * grain_radius_surface)

    # effective solar zenith angle
    x = min((t / (3 * cos(pi * SZA / 180)))^0.5, 1.0)
    u = 0.64 * x + (1 - x) * cos(pi * SZA / 180)

    # pure snow albedo
    as = 1.48 - S1^(-0.07)

    # change in pure snow albedo due to soot loading
    dac = max(0.04 - as,
        -(c^0.55) / (0.16 + 0.6 * S1^0.5 + (1.8 * c^0.6) * (S1^(-0.25))))

    # change in albedo due to solar zenith angle
    dasz = 0.53 * as * (1 - (as + dac)) * (1 - u)^1.2

    # change in albedo due to cloud (apart from change in diffuse fraction)
    dat = (0.1 * t * (as + dac)^1.3) / ((1 + 1.5 * t)^as)

    # Broadband albedo
    albedo = as + dac + dasz + dat

    return albedo
end
