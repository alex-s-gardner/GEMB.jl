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
fixed 50° effective zenith angle); the other methods return it equal to `albedo`.
[`calculate_shortwave_radiation`](@ref) applies `albedo` to the direct beam
`SW↓ − SW↓_diffuse` and `albedo_diffuse` to `SW↓_diffuse`.

Under `:GardnerSharp` the direct-beam albedo uses `cfs.solar_zenith_angle` [degrees] and
`cfs.cloud_optical_thickness`. When the sun is at or below the horizon
(`solar_zenith_angle ≥ 90°`) there is no direct beam, and `albedo` is returned equal to
`albedo_diffuse`. A step with `SW↓ = 0` and the sun above the horizon is evaluated normally.
Neither case changes the energy budget, since no shortwave is absorbed when `SW↓ = 0`.

# Firn, ice, and snow over them

Under `:GardnerSharp` and `:BrunLefebre` the grain-size albedo describes semi-infinite snow.
Where the column has a surface of firn or ice, or snow over one, the albedos are
([`firn_ice_albedo`](@ref)):

- **Exposed ice** (`density[1] ≥ mp.density_ice`): `mp.albedo_ice` exactly, for both streams.
- **Exposed firn** (`DENSITY_PORE_CLOSEOFF ≤ density[1] < mp.density_ice`): linear in density
  from `mp.albedo_ice` at `mp.density_ice` to `ALBEDO_FIRN_MAX` (0.65) at pore close-off, for
  both streams.
- **Snow over firn or ice, `:GardnerSharp`**: Gardner & Sharp (2010) eq. 13 (p. 11) for a snow
  layer of depth `z` [m w.e.] over a semi-infinite layer. `z` is the matrix mass of the snow
  cells above the first cell at or above pore close-off, divided by `DENSITY_WATER`; the lower
  layer's albedo is `firn_ice_albedo` of that cell. See [`_albedo_gardner`](@ref). Eq. 13
  keeps the zenith and cloud terms of eqs. 9–11, while exposed firn and ice carry none, so
  with `SZA > 0` or `τ > 0` the albedo drops from the `z → 0` limit of eq. 13 to
  `firn_ice_albedo` when the last snow cell goes.
- **Thin snow over firn or ice, `:BrunLefebre`**: the MAR integrated albedo of Lefebre et al.
  (2003, p. 3 ¶14), linear in the snow thickness `d` between the albedo of the firn or ice
  surface at `d = 0` and the semi-infinite snow albedo at `d ≥ ALBEDO_SNOW_DEPTH_THICK`
  (0.1 m).

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
    elseif mp.albedo_method == :GreuellKonzelmann
        albedo = mp.albedo_ice + (density[1] - mp.density_ice) *
            (mp.albedo_snow - mp.albedo_ice) /
            (density_fresh_snow - mp.density_ice) +
            (0.05 * (cfs.cloud_fraction - 0.5))
        albedo_diffuse = albedo

    elseif (mp.albedo_method == :GardnerSharp) || (mp.albedo_method == :BrunLefebre)
        albedo_under, depth_snow, mass_snow = _snow_over_firn_ice(dz, density, mp)

        if depth_snow == 0.0
            # Exposed firn or ice.
            albedo = albedo_under
            albedo_diffuse = albedo_under

        elseif mp.albedo_method == :GardnerSharp
            # `z_we` is Inf when no firn or ice lies under the snow (semi-infinite snow).
            z_we = mass_snow / DENSITY_WATER
            albedo_diffuse = _albedo_gardner(grain_radius[1], cfs.black_carbon_snow,
                ALBEDO_DIFFUSE_ZENITH_ANGLE, cfs.cloud_optical_thickness, albedo_under, z_we)
            # With the sun at or below the horizon there is no direct beam; see the docstring.
            albedo = cfs.solar_zenith_angle < 90.0 ?
                     _albedo_gardner(grain_radius[1], cfs.black_carbon_snow,
                         cfs.solar_zenith_angle, cfs.cloud_optical_thickness, albedo_under,
                         z_we) :
                     albedo_diffuse

        else
            albedo = _albedo_brun_lefebre(grain_radius[1])
            if depth_snow < ALBEDO_SNOW_DEPTH_THICK
                albedo = albedo_under +
                         (albedo - albedo_under) * (depth_snow / ALBEDO_SNOW_DEPTH_THICK)
            end
            albedo_diffuse = albedo
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
    _snow_over_firn_ice(dz, density, mp) -> (albedo_under, depth_snow, mass_snow)

The snow above the first cell at or above `DENSITY_PORE_CLOSEOFF`: its thickness `depth_snow`
[m] and matrix mass `mass_snow = Σ dz·density` [kg m-2], and `albedo_under`, the
[`firn_ice_albedo`](@ref) of that firn or ice cell. Both amounts are 0 for exposed firn or ice.
For a column with no firn or ice both are `Inf` and `albedo_under` is 0 (unused).
"""
@inline function _snow_over_firn_ice(dz::Vector{Float64}, density::Vector{Float64},
    mp::ModelParameters)

    depth_snow = 0.0
    mass_snow = 0.0
    for i in eachindex(dz, density)
        if density[i] >= DENSITY_PORE_CLOSEOFF - D_TOLERANCE
            return firn_ice_albedo(density[i], mp), depth_snow, mass_snow
        end
        depth_snow += dz[i]
        mass_snow += dz[i] * density[i]
    end
    return 0.0, Inf, Inf
end

# Spectral fractions of Lefebre et al. (2003) for 0.3–0.8, 0.8–1.5 and 1.5–2.8 µm; they sum to 1.
const BRUN_LEFEBRE_SPECTRAL_FRACTIONS = (0.606, 0.301, 0.093)

"""
    _brun_lefebre_band_albedos(grain_radius_surface) -> (a1, a2, a3)

Albedos of semi-infinite snow in the three spectral bands of Lefebre et al. (2003)
(`BRUN_LEFEBRE_SPECTRAL_FRACTIONS`), from the surface effective grain radius [mm].
"""
@inline function _brun_lefebre_band_albedos(grain_radius_surface::Float64)
    # grain size (diameter) in metres
    gsz = (grain_radius_surface * 2.0) / 1000.0
    a1 = min(0.98, 0.95 - 1.58 * gsz^0.5)
    a2 = max(0.0, 0.95 - 15.4 * gsz^0.5)
    a3 = max(0.127, 0.88 + 346.3 * gsz - 32.31 * gsz^0.5)
    return a1, a2, a3
end

"""
    _albedo_brun_lefebre(grain_radius_surface) -> α

Broadband albedo of semi-infinite snow, the band albedos of
[`_brun_lefebre_band_albedos`](@ref) weighted by `BRUN_LEFEBRE_SPECTRAL_FRACTIONS`.
"""
@inline function _albedo_brun_lefebre(grain_radius_surface::Float64)
    sF = BRUN_LEFEBRE_SPECTRAL_FRACTIONS
    a1, a2, a3 = _brun_lefebre_band_albedos(grain_radius_surface)
    return sF[1] * a1 + sF[2] * a2 + sF[3] * a3
end

"""
    _albedo_gardner(grain_radius_surface, c, SZA, τ) -> α
    _albedo_gardner(grain_radius_surface, c, SZA, τ, albedo_under, z_we) -> α

Broadband albedo from Gardner and Sharp (2010, eqs. 7–12, pp. 10–11): the pure-snow albedo
`α_S` from the specific surface area (eq. 7), plus the changes due to soot loading `c` [ppmw]
(eq. 8, giving `α_c = α_S + dα_c`), solar zenith angle `SZA` [degrees, 0 ≤ SZA < 90] (eqs. 9–10)
and cloud optical thickness `τ ≥ 0` (eq. 11), summed by eq. 12. A zenith angle outside that
range throws: the zenith term is defined only for the sun above the horizon.

The four-argument form is semi-infinite snow. The six-argument form is a snow layer `z_we`
[m w.e.] deep over a semi-infinite layer of albedo `albedo_under`, by eq. 13 (p. 11 ¶51–52):
`dα_c` is replaced in eqs. 8–12 by

    dα′_c = (α_c^btm − α_S^top) + A·(α_c^top − α_c^btm),
    A = min(1, 2.1·z^p),   p = 1.35·(1 − α_S^top) − 0.1·c^top − 0.13,

with `α_c^btm = albedo_under` and the snow's `α_S`, `α_c` and `c` as the top layer. The zenith
and cloud terms (eqs. 9 and 11) are then evaluated with `α_c′ = α_S + dα′_c`. As `z_we → 0`
the result therefore tends to eqs. 9–12 with `α_c = albedo_under`, which equals
`albedo_under` at `SZA = 0` and `τ = 0`.

Where `A = 1`, from `z_we = 2.1^(−1/p)` upward, the result is the semi-infinite value exactly;
`z_we = Inf` is semi-infinite snow. Gardner and Sharp fitted eq. 13 for snow over ice or old
snow (p. 11 ¶52) with `z` 0.01–0.5 m w.e., `r_e` 0.05–1 mm, `c` 0–0.3 ppmw, `SZA` 0–75° and
`τ` 0–20 (Table B2, p. 13). Below `z = 0.01` m w.e. the formula is extrapolated toward
`A = 0`. Over that fitted range `p > 0`; for `p ≤ 0` (very fine or sooty snow, outside it) `A`
is 1, which is the formula's value for every `z_we ≤ 1` m w.e.
"""
_albedo_gardner(grain_radius_surface::Float64, c::Float64, SZA::Float64, t::Float64) =
    _albedo_gardner(grain_radius_surface, c, SZA, t, 0.0, Inf)

function _albedo_gardner(grain_radius_surface::Float64, c::Float64,
    SZA::Float64, t::Float64, albedo_under::Float64, z_we::Float64)

    (0.0 <= SZA < 90.0) ||
        throw(DomainError(SZA, "solar zenith angle must be in [0, 90) degrees for the " *
                               "Gardner-Sharp zenith term"))
    t >= 0.0 || throw(DomainError(t, "cloud optical thickness must be non-negative"))
    z_we >= 0.0 || throw(DomainError(z_we, "snow depth in m w.e. must be non-negative"))

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

    # finite snow depth over a lower layer of albedo `albedo_under` (eq. 13); at A = 1 the
    # semi-infinite `dac` is kept as is
    if z_we < Inf
        (0.0 <= albedo_under <= 1.0) ||
            throw(DomainError(albedo_under, "albedo of the layer under the snow must be in [0, 1]"))
        p = 1.35 * (1 - as) - 0.1 * c - 0.13
        A = p > 0.0 ? min(1.0, 2.1 * z_we^p) : 1.0
        if A < 1.0
            dac = (albedo_under - as) + A * ((as + dac) - albedo_under)
        end
    end

    # change in albedo due to solar zenith angle
    dasz = 0.53 * as * (1 - (as + dac)) * (1 - u)^1.2

    # change in albedo due to cloud (apart from change in diffuse fraction)
    dat = (0.1 * t * (as + dac)^1.3) / ((1 + 1.5 * t)^as)

    # Broadband albedo
    albedo = as + dac + dasz + dat

    return albedo
end
