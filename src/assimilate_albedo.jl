"""
Updating the near-surface snow from an observed broadband albedo.

A broadband albedo constrains a combination of the snow's specific surface area (SSA) and its
light-absorbing impurity loading, not each separately: a line of `(SSA, c)` pairs gives the same
Gardner & Sharp (2010) albedo. [`albedo_update`](@ref) therefore takes the smallest change from the
model's own state that reproduces the observation, measured in `ln SSA` and `ln(c + c₀)` and weighted
by the prior standard deviations `mp.albedo_update_sigma_log_ssa` and
`mp.albedo_update_sigma_log_impurity`. The split between the two is set by those weights, not by
the observation.

The loading so obtained is an *effective* one, not a measured or inferred impurity concentration: it also absorbs errors of the albedo model, of
wet-snow albedo and of mixing in the observation footprint. The near-surface snow a satellite sees
during melt holds the insoluble residue meltwater did not carry away, so the update sets the
retained class `impurity_effective_retained` of the near-surface layers, which meltwater scavenges at
`mp.impurity_scavenging_efficiency_retained`, and leaves the mobile `impurity_effective` from snowfall
as it is. The retained loading then moves with the snow, buried by new snow and re-exposed by
melt, until the next observation; the grain radius is scaled and evolves by metamorphism.
"""

# Offset c₀ [ppmw] in the update's impurity coordinate ln(c + c₀), so clean snow (c = 0) is finite.
const IMPURITY_UPDATE_OFFSET = 0.01
# Largest effective loading [ppmw] an update may set.
const IMPURITY_EFFECTIVE_MAX = 20.0
# Grain-radius range [mm] an update may set: new snow to the non-spherical grain cap.
const ALBEDO_UPDATE_RADIUS_RANGE = (RE_NEW_SNOW, GRAIN_RADIUS_ICE)
const ALBEDO_UPDATE_TOLERANCE = 1e-6
# Range an observation of exposed ice may set the column's ice albedo to: dark, impurity-rich ice
# to the albedo of firn at pore close-off.
const ALBEDO_ICE_OBSERVED_RANGE = (0.1, ALBEDO_FIRN_MAX)
# Largest fraction of the way from ice to pore close-off density at which an observation of
# exposed firn sets the ice albedo (see `_ice_albedo_update`).
const ALBEDO_ICE_UPDATE_FIRN_MAX = 0.5
# Points of the grain-size scan in `albedo_update`, and the albedo residual the loading solve stops at.
const ALBEDO_UPDATE_SCAN_POINTS = 40
const ALBEDO_LOADING_TOLERANCE = 1e-13

# SSA [cm2 g-1] of an effective grain radius [mm], and back, as in `_gardner_terms`.
@inline _ssa_of_radius(r) = 3.0 / (0.091 * r)
@inline _radius_of_ssa(S) = 3.0 / (0.091 * S)

"""
    albedo_update(albedo_observed, r, c, α, mp) -> (r_new, c_new, converged)

The near-surface grain radius `r_new` [mm] and effective loading `c_new` [ppmw] whose albedo `α(r, c)`
equals `albedo_observed`, by the smallest change from the model's `r` and `c` in
`y = (ln SSA, ln(c + c₀))` weighted by
`σ = (mp.albedo_update_sigma_log_ssa, mp.albedo_update_sigma_log_impurity)`:

    minimize ((y₁ − y₁⁰)/σ₁)² + ((y₂ − y₂⁰)/σ₂)²  subject to  α(y) = albedo_observed,

with `r` in `ALBEDO_UPDATE_RADIUS_RANGE` and `c` in `[0, IMPURITY_EFFECTIVE_MAX]`.

The albedo falls monotonically with the loading at fixed grain size, so for each `y₁` the
constraint has at most one `y₂` ([`_loading_for_albedo`](@ref)). The update is the minimum of the
cost along that contour, a one-dimensional problem in `y₁`: a scan of `ALBEDO_UPDATE_SCAN_POINTS`
points over the bounds of `y₁` finds the best point at which the contour exists, and a
golden-section search refines it within the neighbouring interval. With `σ_ssa = 0` the grain size is
held and only the loading changes. When no state within the bounds reaches the observation, the
update is the reachable state nearest to it (clean snow at the smallest grain for an observation
brighter than any, the largest loading at the largest grain for one darker than any), with
`converged = false`.

`α` must fall with `c` at fixed `r`; [`assimilate_albedo!`](@ref) passes the albedo GEMB applies to the
observation step.
"""
function albedo_update(albedo_observed::Real, r::Real, c::Real, α::F, mp::ModelParameters) where {F}
    σ1, σ2 = mp.albedo_update_sigma_log_ssa, mp.albedo_update_sigma_log_impurity
    (σ1 >= 0 && σ2 > 0) || throw(ArgumentError(
        "albedo_update_sigma_log_ssa must be ≥ 0 and albedo_update_sigma_log_impurity > 0, got " *
        "$(σ1) and $(σ2)"))
    target = Float64(albedo_observed)
    rlo, rhi = ALBEDO_UPDATE_RADIUS_RANGE
    r0 = clamp(Float64(r), rlo, rhi); c0 = clamp(Float64(c), 0.0, IMPURITY_EFFECTIVE_MAX)
    y10 = log(_ssa_of_radius(r0)); y20 = log(c0 + IMPURITY_UPDATE_OFFSET)
    radius(y1) = _radius_of_ssa(exp(y1))
    cost(y1, cc) = (σ1 > 0 ? ((y1 - y10) / σ1)^2 : 0.0) + ((log(cc + IMPURITY_UPDATE_OFFSET) - y20) / σ2)^2

    # Grain size held: a one-dimensional solve in the loading.
    if σ1 == 0
        cc, ok = _loading_for_albedo(target, r0, α)
        return r0, cc, ok
    end
    # Out of reach at every grain size within the bounds.
    brightest = α(rlo, 0.0)
    darkest = α(rhi, IMPURITY_EFFECTIVE_MAX)
    target > brightest && return rlo, 0.0, false
    target < darkest && return rhi, IMPURITY_EFFECTIVE_MAX, false
    # The contour exists where clean snow of that grain is at least as bright as the target and the
    # largest loading at least as dark. Scan y₁ for the best feasible point, then refine it by
    # golden-section search over the neighbouring interval.
    y1lo, y1hi = log(_ssa_of_radius(rhi)), log(_ssa_of_radius(rlo))
    # Inf where the contour does not exist at this grain size: clean snow darker than the target, or
    # the largest loading brighter.
    function f(y1)
        r1 = radius(y1)
        a_clean = α(r1, 0.0)
        a_dark = α(r1, IMPURITY_EFFECTIVE_MAX)
        (a_clean >= target && a_dark <= target) || return Inf
        return cost(y1, _loading_in_bracket(target, r1, α, a_clean, a_dark))
    end
    n = ALBEDO_UPDATE_SCAN_POINTS
    grid = range(y1lo, y1hi; length = n)
    vals = f.(grid)
    k = argmin(vals)
    lo = grid[max(k - 1, 1)]; hi = grid[min(k + 1, n)]
    φ = (sqrt(5.0) - 1) / 2
    x1 = hi - φ * (hi - lo); x2 = lo + φ * (hi - lo); f1 = f(x1); f2 = f(x2)
    for _ in 1:80
        hi - lo < 1e-10 && break
        if f1 < f2
            hi = x2; x2 = x1; f2 = f1; x1 = hi - φ * (hi - lo); f1 = f(x1)
        else
            lo = x1; x1 = x2; f1 = f2; x2 = lo + φ * (hi - lo); f2 = f(x2)
        end
    end
    y1 = (lo + hi) / 2
    rr = radius(y1)
    cc, ok = _loading_for_albedo(target, rr, α)
    return rr, cc, ok
end

"""
    _loading_for_albedo(target, r, α) -> (c, converged)

The effective loading `c` [ppmw] in `[0, IMPURITY_EFFECTIVE_MAX]` at which the albedo `α(r, c)` of grain
radius `r` equals `target` ([`_loading_in_bracket`](@ref)). Clamped to a bound when
the target lies outside that range, with `converged` true only if the bound's albedo is within
`ALBEDO_UPDATE_TOLERANCE` of it.
"""
function _loading_for_albedo(target::Float64, r::Float64, α::F) where {F}
    a_clean = α(r, 0.0)
    a_clean <= target && return 0.0, abs(a_clean - target) < ALBEDO_UPDATE_TOLERANCE
    a_dark = α(r, IMPURITY_EFFECTIVE_MAX)
    a_dark >= target && return IMPURITY_EFFECTIVE_MAX, abs(a_dark - target) < ALBEDO_UPDATE_TOLERANCE
    c = _loading_in_bracket(target, r, α, a_clean, a_dark)
    return c, abs(α(r, c) - target) < ALBEDO_UPDATE_TOLERANCE
end

"""
    _loading_in_bracket(target, r, α, a_clean, a_dark) -> c

The loading `c` [ppmw] at which the albedo `α(r, c)` of grain radius `r` equals `target`, given
that clean snow (`a_clean`) is at least as bright as `target` and the largest loading (`a_dark`) at
most as bright. The albedo falls monotonically with `c`; the root is found in `ln(c + c₀)`, where
the albedo is close to linear, by regula falsi with the Illinois modification (Dowell and Jarratt,
1971), which keeps the bracket and converges superlinearly. It stops when the albedo is within
`ALBEDO_LOADING_TOLERANCE` of `target` or the bracket has closed.

# References
- Dowell, M. and Jarratt, P. (1971). A modified regula falsi method for computing the root of an
  equation. *BIT* 11, 168-174.
"""
function _loading_in_bracket(target::Float64, r::Float64, α::F, a_clean::Float64,
    a_dark::Float64) where {F}
    c0 = IMPURITY_UPDATE_OFFSET
    loading(y) = clamp(exp(y) - c0, 0.0, IMPURITY_EFFECTIVE_MAX)
    ylo, yhi = log(c0), log(IMPURITY_EFFECTIVE_MAX + c0)
    glo, ghi = a_clean - target, a_dark - target
    glo == 0.0 && return 0.0
    ghi == 0.0 && return IMPURITY_EFFECTIVE_MAX
    side = 0
    y = ylo
    for _ in 1:100
        y = (ylo * ghi - yhi * glo) / (ghi - glo)
        g = α(r, loading(y)) - target
        (abs(g) < ALBEDO_LOADING_TOLERANCE || yhi - ylo < 1e-14) && break
        if g > 0
            ylo, glo = y, g
            side == -1 && (ghi /= 2)
            side = -1
        else
            yhi, ghi = y, g
            side == 1 && (glo /= 2)
            side = 1
        end
    end
    return loading(y)
end

"""
    assimilate_albedo!(dz, density, grain_radius, age, impurity_effective, impurity_effective_retained, cfs, mp;
                       albedo_ice = mp.albedo_ice) -> (dlog_ssa, dlog_impurity, applied, albedo_ice)

Update the near-surface snow from `cfs.albedo_observed` when it is finite, the albedo method is
`:GardnerSharp` and the surface is snow (its top layer below pore close-off density). The update
is [`albedo_update`](@ref) from the surface grain radius and the total [`surface_impurity_effective`](@ref).
It sets the surface layer's retained film so that the loading the albedo reads is the new one: the
film supplies `max(0, c_new − c_mobile)`, `c_mobile` being the mobile loading over the top
`mp.albedo_impurity_depth`, held as retained mass in the surface layer
([`surface_impurity_effective`](@ref)); the retained loading of the layers below is left as it is. It
scales `grain_radius` of every snow layer within the top `mp.albedo_impurity_depth` by
`r_new / r_surface`, leaving the profile's shape. Where the mobile loading alone exceeds `c_new` the
film is 0 and the surface stays darker than the update.

No update is made while the surface layer is younger than `mp.albedo_update_snow_age_min` [d], 0 by
default, so every observation of a snow surface is used.

The state matched is the one whose albedo GEMB applies to the step's shortwave, `1 − SW_net / SW↓`
([`effective_albedo`](@ref)), equals the observation, which is `Σ SW↑ / Σ SW↓` under the sky it was measured
in.

An observation of exposed firn or ice sets the column's ice albedo instead
([`_ice_albedo_update`](@ref)). `albedo_ice` is the column's current ice albedo, read for the snow
over firn or ice; the caller carries it as column state (`albedo_ice` of the [`gemb_core`](@ref)
state), so an observed ice albedo holds until the next observation of exposed firn or ice and is
also the ice albedo under snow. The `:ground` and `:prescribed` columns, whose firn or ice is a
refrozen lens in the snow, do not carry it and keep `mp.albedo_ice`.

Returns the change of the surface state, `Δ ln SSA` and `Δ ln(c + c₀)`, whether an update was
applied, and the column's ice albedo after the step's update; `(0, 0, false, albedo_ice)` when
none was.
"""
function assimilate_albedo!(dz::Vector{Float64}, density::Vector{Float64},
    grain_radius::Vector{Float64}, age::Vector{Float64}, impurity_effective::Vector{Float64},
    impurity_effective_retained::Vector{Float64}, cfs::ClimateForcingStep, mp::ModelParameters;
    albedo_ice::Float64=mp.albedo_ice)
    (isfinite(cfs.albedo_observed) && mp.albedo_method === :GardnerSharp && !isempty(dz)) ||
        return (0.0, 0.0, false, albedo_ice)
    if density[1] >= DENSITY_PORE_CLOSEOFF - D_TOLERANCE
        applied, albedo_ice = _ice_albedo_update(cfs.albedo_observed, density[1], albedo_ice, mp)
        return (0.0, 0.0, applied, albedo_ice)
    end
    age[1] < mp.albedo_update_snow_age_min && return (0.0, 0.0, false, albedo_ice)
    albedo_under, depth_snow, mass_snow = _snow_over_firn_ice(dz, density, mp, albedo_ice)
    r0 = grain_radius[1]
    c0 = surface_impurity_effective(dz, density, impurity_effective, impurity_effective_retained, mp)
    # The albedo GEMB applies to this step's shortwave from a surface state, 1 − SW_net / SW↓.
    α(rr, cc) = effective_albedo(_albedo_snow_over(albedo_under, depth_snow, mass_snow, rr, cc, cfs, mp)..., cfs)
    r, c, _ = albedo_update(cfs.albedo_observed, r0, c0, α, mp)
    scale = r / r0
    film = impurity_effective_retained[1] / mp.albedo_impurity_depth
    impurity_effective_retained[1] = max(0.0, c - (c0 - film)) * mp.albedo_impurity_depth
    remaining = mp.albedo_impurity_depth
    for i in eachindex(dz, density, grain_radius)
        density[i] >= DENSITY_PORE_CLOSEOFF - D_TOLERANCE && break
        grain_radius[i] = clamp(grain_radius[i] * scale, ALBEDO_UPDATE_RADIUS_RANGE...)
        remaining -= dz[i]
        remaining <= 0.0 && break
    end
    return (log(_ssa_of_radius(r)) - log(_ssa_of_radius(r0)),
            log(c + IMPURITY_UPDATE_OFFSET) - log(c0 + IMPURITY_UPDATE_OFFSET), true, albedo_ice)
end

"""
    _ice_albedo_update(albedo_observed, density_surface, albedo_ice, mp) -> (applied, albedo_ice_new)

The ice albedo whose [`firn_ice_albedo`](@ref) at the exposed surface density `density_surface`
[kg m-3] equals `albedo_observed`, clamped to `ALBEDO_ICE_OBSERVED_RANGE`. Exposed firn and ice
reflect that albedo in both streams, so it is also the albedo GEMB applies. Firn whose
`_firn_fraction` exceeds `ALBEDO_ICE_UPDATE_FIRN_MAX`, whose albedo hardly depends on the ice
albedo, is not updated.
"""
function _ice_albedo_update(albedo_observed::Float64, density_surface::Float64,
    albedo_ice::Float64, mp::ModelParameters)
    f = _firn_fraction(density_surface, mp)
    f > ALBEDO_ICE_UPDATE_FIRN_MAX && return (false, albedo_ice)
    # Inverse of `firn_ice_albedo`: α = (1 − f)·albedo_ice + f·ALBEDO_FIRN_MAX.
    albedo_ice_new = (albedo_observed - f * ALBEDO_FIRN_MAX) / (1.0 - f)
    return (true, clamp(albedo_ice_new, ALBEDO_ICE_OBSERVED_RANGE...))
end
