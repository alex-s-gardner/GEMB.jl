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
  layer of depth `z` [m w.e.] over a semi-infinite layer whose albedo is `firn_ice_albedo` of
  the first cell at or above pore close-off, with `z` the matrix mass of the snow cells above
  that cell divided by `DENSITY_WATER`. The zenith and cloud terms of eqs. 9–11 are scaled by
  eq. 13's thin-snow weight `A`, so they vanish with the snow, as they do on exposed firn and
  ice: the albedo is continuous at snow exposure, its `z → 0` limit is exactly the exposed
  albedo, and wherever `A = 1` it is the semi-infinite value of eqs. 7–12. See
  [`_albedo_gardner`](@ref).
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
        albedo, albedo_diffuse = _albedo_snow_over(albedo_under, depth_snow, mass_snow,
            grain_radius[1], cfs, mp)
    end

    _check_albedo(albedo)
    return albedo, albedo_diffuse
end

# `:GardnerSharp` and `:BrunLefebre` albedos of snow `depth_snow` [m] thick, of matrix mass
# `mass_snow` [kg m-2], over a surface of albedo `albedo_under`; see `calculate_albedo`.
@inline function _albedo_snow_over(albedo_under::Float64, depth_snow::Float64,
    mass_snow::Float64, grain_radius_surface::Float64, cfs::ClimateForcingStep,
    mp::ModelParameters)

    if depth_snow == 0.0
        # Exposed firn or ice.
        albedo = albedo_under
        albedo_diffuse = albedo_under

    elseif mp.albedo_method == :GardnerSharp
        # `z_we` is Inf when no firn or ice lies under the snow (semi-infinite snow).
        z_we = mass_snow / DENSITY_WATER
        albedo_diffuse = _albedo_gardner(grain_radius_surface, cfs.black_carbon_snow,
            ALBEDO_DIFFUSE_ZENITH_ANGLE, cfs.cloud_optical_thickness, albedo_under, z_we)
        # With the sun at or below the horizon there is no direct beam; see the docstring.
        albedo = cfs.solar_zenith_angle < 90.0 ?
                 _albedo_gardner(grain_radius_surface, cfs.black_carbon_snow,
                     cfs.solar_zenith_angle, cfs.cloud_optical_thickness, albedo_under,
                     z_we) :
                 albedo_diffuse

    else
        albedo = _albedo_brun_lefebre(grain_radius_surface)
        if depth_snow < ALBEDO_SNOW_DEPTH_THICK
            albedo = albedo_under +
                     (albedo - albedo_under) * (depth_snow / ALBEDO_SNOW_DEPTH_THICK)
        end
        albedo_diffuse = albedo
    end
    return albedo, albedo_diffuse
end

"""
    calculate_albedo_ground(dz, density, grain_radius, cfs, mp) -> (albedo, albedo_diffuse)

[`calculate_albedo`](@ref) for the snow over the ground of `lower_boundary = :ground`. `dz`,
`density` and `grain_radius` are the snow cells only and may be empty.

- Snow-free ground takes `mp.ground_albedo` for both streams.
- Under `:GardnerSharp` and `:BrunLefebre`, snow over the ground blends to `mp.ground_albedo` by
  the same thin-snow rules that blend snow over firn or ice to [`firn_ice_albedo`](@ref): Gardner
  and Sharp (2010) eq. 13, and the linear depth weighting of Lefebre et al. (2003). A firn or ice
  layer within the snow (a refrozen lens at or above `DENSITY_PORE_CLOSEOFF`) is the surface
  underneath, as over a glacier.
- `:None` and `:GreuellKonzelmann` treat snow over the ground as `calculate_albedo` does.
"""
function calculate_albedo_ground(dz::Vector{Float64}, density::Vector{Float64},
    grain_radius::Vector{Float64}, cfs::ClimateForcingStep, mp::ModelParameters)
    isempty(dz) && return mp.ground_albedo, mp.ground_albedo
    if (mp.albedo_method == :GardnerSharp) || (mp.albedo_method == :BrunLefebre)
        if !((mp.albedo_density_threshold - density[1]) < D_TOLERANCE)
            albedo_under, depth_snow, mass_snow = _snow_over_firn_ice(dz, density, mp)
            if isinf(depth_snow)
                albedo_under = mp.ground_albedo
                depth_snow = sum(dz)
                mass_snow = 0.0
                for i in eachindex(dz, density)
                    mass_snow += dz[i] * density[i]
                end
            end
            albedo, albedo_diffuse = _albedo_snow_over(albedo_under, depth_snow, mass_snow,
                grain_radius[1], cfs, mp)
            _check_albedo(albedo)
            return albedo, albedo_diffuse
        end
    end
    return calculate_albedo(dz, density, grain_radius, cfs, mp)
end

@inline function _check_albedo(albedo::Float64)
    if albedo > (1 + T_TOLERANCE)
        @warn "albedo > 1.0"
    elseif albedo < (0 - D_TOLERANCE)
        @warn "albedo is negative"
    elseif isnan(albedo)
        error("albedo == NAN")
    end
    return nothing
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
and cloud optical thickness `τ ≥ 0` (eq. 11), summed by eq. 12:

    α = α_c + dα_u(α_c) + dα_τ(α_c),
    dα_u(α_c) = 0.53·α_S·(1 − α_c)·(1 − u′)^1.2,   dα_τ(α_c) = 0.1·τ·α_c^1.3 / (1 + 1.5τ)^α_S,

with `u′` the cosine of the effective zenith angle (eq. 10). A zenith angle outside `[0, 90)`
throws: the zenith term is defined only for the sun above the horizon.

The four-argument form is semi-infinite snow. The six-argument form is a snow layer `z_we`
[m w.e.] deep over firn or ice whose exposed albedo is `albedo_under`, which carries no zenith
or cloud term ([`firn_ice_albedo`](@ref)). Eq. 13 (p. 11 ¶51–52) replaces `dα_c` in eqs. 8–12 by

    dα′_c = (α_c^btm − α_S^top) + A·(α_c^top − α_c^btm),
    A = min(1, 2.1·z^p),   p = 1.35·(1 − α_S^top) − 0.1·c^top − 0.13,

so `α′_c = α_c^btm + A·(α_c^top − α_c^btm)`, with `α_c^btm = albedo_under` and the snow's
`α_S`, `α_c` and `c` as the top layer. The zenith and cloud increments are evaluated at `α′_c`,
as in the paper, and multiplied by `A`:

    α = α′_c + A·(dα_u(α′_c) + dα_τ(α′_c)).

**Departure from the paper:** the factor `A` on the two increments. Eq. 12 adds them unscaled,
so as `z → 0` the paper's albedo tends to eqs. 9–12 evaluated with `α_c = α_c^btm`, which lies
above an exposed albedo that has no zenith or cloud term (by up to 0.45 for fresh snow at low
sun under clear sky). With the factor, the `z = 0` value is `albedo_under` exactly at every
geometry and the albedo is continuous when the last snow goes. Where `A = 1`, from
`z_we = 2.1^(−1/p)` upward, the factor is 1 and the result is the semi-infinite value of
eqs. 7–12 bit for bit; `z_we = Inf` is semi-infinite snow. Inside the fitted depth range
(`z ≥ 0.01` m w.e.) `A < 1` only where `p > ln 2.1 / ln 100 ≈ 0.161`, which for clean snow
means `r_e` above about 0.18 mm.

Writing `F = dα_u + dα_τ`, `F(α) ≤ k_u + k_τ` on `[0, 1]` with `k_u = 0.53·α_S·(1 − u′)^1.2 < 1`
and `k_τ = 0.1·τ/(1 + 1.5τ)^α_S`, so `|α − albedo_under| ≤ A·(|α_c^top − albedo_under| + k_u +
k_τ)`: the albedo approaches the exposed value as `A`, i.e. as `z^p`. `dα/dA = (α_c^top −
albedo_under)·(1 + A·F′(α′_c)) + F(α′_c)` with `1 + A·F′ ≥ 1 − k_u > 0`, so where the snow's
`α_c` is at least `albedo_under` the albedo increases with `z`. The approach is slow: `p` is
0.08–0.35 over the fitted grain sizes.

Gardner and Sharp fitted eq. 13 for snow over ice or old snow (p. 11 ¶52) with `z`
0.01–0.5 m w.e., `r_e` 0.05–1 mm, `c` 0–0.3 ppmw, `SZA` 0–75° and `τ` 0–20 (Table B2, p. 13).
Below `z = 0.01` m w.e. the formula is extrapolated toward `A = 0`. Over that fitted range
`p > 0`; for `p ≤ 0` (very fine or sooty snow, outside it) `A` is 1, which is the formula's
value for every `z_we ≤ 1` m w.e., so such snow keeps its semi-infinite albedo down to `z = 0`.
"""
_albedo_gardner(grain_radius_surface::Float64, c::Float64, SZA::Float64, t::Float64) =
    _albedo_gardner(grain_radius_surface, c, SZA, t, 0.0, Inf)

function _albedo_gardner(grain_radius_surface::Float64, c::Float64,
    SZA::Float64, t::Float64, albedo_under::Float64, z_we::Float64)

    as, dac, u = _gardner_terms(grain_radius_surface, c, SZA, t, z_we)
    z_we < Inf || return _gardner_sum(as, dac, u, t)

    (0.0 <= albedo_under <= 1.0) ||
        throw(DomainError(albedo_under, "albedo of the layer under the snow must be in [0, 1]"))
    A = _gardner_thin_snow_weight(as, c, z_we)
    A < 1.0 || return _gardner_sum(as, dac, u, t)

    # eq. 13: α′_c, and the zenith and cloud increments scaled by A
    ac = albedo_under + A * ((as + dac) - albedo_under)
    return ac + A * (_gardner_zenith_term(as, ac, u) + _gardner_cloud_term(as, ac, t))
end

# Pure-snow albedo `as` (eq. 7), its change by soot `dac` (eq. 8) and the cosine `u` of the
# effective zenith angle (eq. 10) of Gardner and Sharp (2010), with the argument checks shared by
# every Gardner-Sharp albedo.
@inline function _gardner_terms(grain_radius_surface::Float64, c::Float64, SZA::Float64,
    t::Float64, z_we::Float64)

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

    return as, dac, u
end

# The weight `A` of eq. 13 for snow `z_we` [m w.e.] deep; see `_albedo_gardner` for the `p ≤ 0`
# case.
@inline function _gardner_thin_snow_weight(as::Float64, c::Float64, z_we::Float64)
    p = 1.35 * (1 - as) - 0.1 * c - 0.13
    return p > 0.0 ? min(1.0, 2.1 * z_we^p) : 1.0
end

# Eq. 9: change in albedo due to solar zenith angle, at loaded-snow albedo `ac`.
@inline _gardner_zenith_term(as::Float64, ac::Float64, u::Float64) =
    0.53 * as * (1 - ac) * (1 - u)^1.2

# Eq. 11: change in albedo due to cloud (apart from the change in diffuse fraction), at
# loaded-snow albedo `ac`.
@inline _gardner_cloud_term(as::Float64, ac::Float64, t::Float64) =
    (0.1 * t * ac^1.3) / ((1 + 1.5 * t)^as)

# Eq. 12 for semi-infinite snow: the sum of `as`, `dac`, and the zenith and cloud terms evaluated
# with α_c = as + dac.
@inline function _gardner_sum(as::Float64, dac::Float64, u::Float64, t::Float64)
    return as + dac + _gardner_zenith_term(as, as + dac, u) + _gardner_cloud_term(as, as + dac, t)
end
