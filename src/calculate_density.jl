"""
    calculate_density(temperature, dz, density, grain_radius, water, cfs::ClimateForcingStep, mp::ModelParameters)

Compute the densification of snow/firn using one of several models:
- `:HerronLangway`: Herron and Langway (1980)
- `:Arthern`: semi-empirical model of Arthern et al. (2010) [default]
- `:ArthernB`: physical model from Appendix B of Arthern et al. (2010)
- `:GSFC2020`: recalibrated Arthern of Medley et al. (2022), GSFC-FDM v1.2.1
- `:Simonsen2013`: Arthern retuned for Greenland by Simonsen et al. (2013)
- `:LiZwally`: empirical model of Li and Zwally (2004)
- `:Helsen`: modified empirical model by Helsen et al. (2008)
- `:Ligtenberg`: semi-empirical model of Ligtenberg et al. (2011)

Returns `(dz, density)`. `density` is updated in place; `dz` is returned as a
new array (recomputed from the conserved cell mass).

`water` [kg m-2] is only read by `:ArthernB`, the one stress-driven scheme here;
the accumulation-driven schemes proxy overburden with mean accumulation and
ignore it.

The accumulation-driven branches are scalar-loop implementations that are
numerically identical, element by element, to the reference vectorized MATLAB
translation, but avoid the mask / gather / broadcast temporaries (`mass_init`,
`idx`, `H`, `c0`, `c1`, ...) the vectorized form allocated per call. `:ArthernB`
deviates from the reference — see the comment on that branch.

# References
- Arthern, R. J., et al. (2010). J. Geophys. Res., 115, F03011.
- Herron, M. and Langway, C. (1980). J. Glaciol., 25, 373-385.
- Li, J. and Zwally, H. (2004). Ann. Glaciol., 38, 309-313.
- Helsen, M. M., et al. (2008). Science, 320, 1626-1629.
- Ligtenberg, S. R. M., et al. (2011). The Cryosphere, 5, 809-819.
- Medley, B., et al. (2022). The Cryosphere, 16, 3971-4011.
- Simonsen, S. B., et al. (2013). J. Glaciol., 59, 545-558.
- Lundin, J. M. D., et al. (2017). J. Glaciol., 63, 401-422. (FirnMICE; eqs. A36-A37)
"""
function calculate_density(temperature::Vector{Float64}, dz::Vector{Float64},
    density::Vector{Float64}, grain_radius::Vector{Float64}, water::Vector{Float64},
    cfs::ClimateForcingStep, mp::ModelParameters)

    d_tolerance = D_TOLERANCE

    # specify constants
    dt = cfs.dt / 86400.0   # convert from [s] to [d]

    m = length(density)
    density_ice = mp.density_ice
    pm = cfs.precipitation_mean
    tam = cfs.temperature_air_mean

    # New grid-cell lengths (fresh array; density is updated in place).
    dz_out = similar(dz)

    method = mp.densification_method

    if method == :HerronLangway
        # Rate coefficients `_hl_c0`/`_hl_c1` are shared with `steady_state_density`.
        @inbounds for i in 1:m
            T = temperature[i]
            c = density[i] <= DENSITY_STAGE_TRANSITION + d_tolerance ?
                _hl_c0(T, pm) : _hl_c1(T, pm)
            _densify_cell!(density, dz, dz_out, i, c, dt, density_ice, d_tolerance)
        end

    elseif method == :Arthern || method == :Simonsen2013 || method == :Ligtenberg
        # The scaled-Arthern family: one kernel, differing only in the hoisted per-stage
        # scalars `k0`/`k1` that `DensificationCoeffs` computes (1.0 for `:Arthern` itself,
        # Ligtenberg's M0/M1, Simonsen's F0 and F1·γ).
        p = DensificationCoeffs(mp, pm, tam)
        @inbounds for i in 1:m
            c = _arthern_scaled_c(density[i], temperature[i], tam, p.precip_force, p.k0, p.k1)
            _densify_cell!(density, dz, dz_out, i, c, dt, density_ice, d_tolerance)
        end

    elseif method == :ArthernB
        # Arthern et al. (2010) eq. B1, the one stress-driven scheme here.
        #
        # Deviates from the MATLAB reference in two ways (upstream GEMB issue #200):
        #
        #  1. Overburden is a true integral of the overlying load, `Σ (ρⱼdzⱼ + waterⱼ) g`
        #     [Pa]. The reference computes `cumsum(dz[1:end-1]) .* density[1:end-1]`,
        #     which applies the single immediately-overlying cell's density to the whole
        #     overlying depth, and omits both gravity and the pore water.
        #  2. `kc1`/`kc2` are the paper's SI coefficients: per second, stress in Pa,
        #     grain radius in m. `_densify_cell!` consumes `c` as yr-1 (it applies
        #     `c/DAYS_PER_YEAR_DENSIFICATION*dt` with `dt` in days), so they are converted
        #     against that same year length so the two cancel exactly.
        #
        # Cross-checked against the Community Firn Model's `Arthern2010T`, which uses the
        # same `kc·exp(-Ec/RT)·σ/r²` form with `σ = cumsum((mass + LWC·ρ_w)·g)`.
        kc_per_year = DAYS_PER_YEAR_DENSIFICATION * 86400.0
        load = 0.0       # Σ (ρⱼdzⱼ + waterⱼ) over cells above the current one [kg m-2]
        @inbounds for i in 1:m
            d0 = density[i]
            T = temperature[i]
            gr = grain_radius[i] / 1000          # [mm] -> [m]
            obp = load * GRAVITY                 # overburden stress [Pa]
            H = exp((-60000.0 / (T * R_GAS))) * obp / gr^2
            c = (d0 <= DENSITY_STAGE_TRANSITION + d_tolerance ? 9.2e-9 : 3.7e-9) *
                kc_per_year * H
            _densify_cell!(density, dz, dz_out, i, c, dt, density_ice, d_tolerance)
            # accumulate the load using this cell's original density
            load += d0 * dz[i] + water[i]
        end

    elseif method == :GSFC2020
        # Medley et al. (2022). Both accumulation powers are loop-invariant, so they are
        # hoisted with the rest via DensificationCoeffs.
        p = DensificationCoeffs(mp, pm, tam)
        @inbounds for i in 1:m
            c = _gsfc2020_c(density[i], temperature[i], tam, p.k0, p.k1)
            _densify_cell!(density, dz, dz_out, i, c, dt, density_ice, d_tolerance)
        end

    elseif method == :LiZwally || method == :Helsen
        # Li and Zwally (2004) and Helsen et al. (2008) share a form and differ only in the
        # two fitted coefficients of the mean-temperature term. Both remain validator-gated.
        base = (pm / density_ice) * 8.36 * (method == :LiZwally ?
            max(139.21 - 0.542 * tam, 1.0) : max(76.138 - 0.28965 * tam, 1.0))
        @inbounds for i in 1:m
            c = base * max(CtoK - temperature[i], 1.0)^(-2.061)
            _densify_cell!(density, dz, dz_out, i, c, dt, density_ice, d_tolerance)
        end

    else
        error("unrecognized densification method")
    end

    return dz_out, density
end

"""
    _densify_cell!(density, dz_in, dz_out, i, c, dt, density_ice, d_tolerance)

Apply the densification increment for cell `i` given rate coefficient `c`:
update `density[i]` in place, clamp it to the density of ice, and write the
mass-conserving new grid-cell length to `dz_out[i]`.
"""
@inline function _densify_cell!(density::Vector{Float64}, dz_in::Vector{Float64},
    dz_out::Vector{Float64}, i::Int, c::Float64, dt::Float64,
    density_ice::Float64, d_tolerance::Float64)
    @inbounds begin
        d0 = density[i]
        mass = d0 * dz_in[i]
        d = d0 + (c * (density_ice - d0) / DAYS_PER_YEAR_DENSIFICATION * dt)
        if d > density_ice - d_tolerance
            d = density_ice
        end
        density[i] = d
        dz_out[i] = mass / d
    end
end

# ---------------------------------------------------------------------------
# Densification rate coefficients `c` [yr-1]
#
# These are the single source of truth for the three methods `validate_parameters`
# permits: they are called both by the transient branches of `calculate_density`
# above and by the steady-state marcher (`steady_state_density`), so the two can
# never drift. Each is a plain `@inline` function of the layer state, and the
# loop-invariant terms are hoisted by the caller into `DensificationCoeffs`, so
# routing the hot loop through them costs nothing relative to the inline form.
#
# `T` is layer temperature [K], `ρ` layer density [kg m-3], `pm` mean annual
# accumulation [kg m-2 yr-1], `tam` mean annual air temperature [K].
# `R_GAS = 8.314 J mol-1 K-1`.
# ---------------------------------------------------------------------------

# Herron & Langway (1980) two-stage coefficients, equivalent to Arthern et al.
# (2010) eq. (2).
@inline _hl_c0(T::Float64, pm::Float64) = (11.0 * exp(-10160.0 / (T * R_GAS))) * (pm / 1000.0)          # ρ <= 550
@inline _hl_c1(T::Float64, pm::Float64) = (575.0 * exp(-21400.0 / (T * R_GAS))) * sqrt(pm / 1000.0)     # ρ > 550

# Arthern et al. (2010) semi-empirical, and the family of schemes that scale it.
#
# `:Arthern`, `:Ligtenberg` and `:Simonsen2013` are the same law with a per-stage scalar
# applied: 1 for Arthern itself, Ligtenberg's accumulation-dependent M0/M1, Simonsen's F0 and
# F1·γ. They therefore share one kernel rather than each restating `H` — the shared form is
# what makes the `Simonsen ≡ F·Arthern` and `Ligtenberg ≡ M·Arthern` identities structural
# rather than merely tested. `k0`/`k1` are the hoisted stage scalars (see
# `DensificationCoeffs`) and `precip_force = pm * GRAVITY`.
#
# `:GSFC2020` is *not* in this family: it refits both the accumulation exponent and the
# activation energies, so it keeps its own kernel below.
@inline function _arthern_scaled_c(ρ::Float64, T::Float64, tam::Float64,
    precip_force::Float64, k0::Float64, k1::Float64)
    H = exp((-60000.0 / (T * R_GAS)) + (42400.0 / (tam * R_GAS))) * precip_force
    return ρ <= DENSITY_STAGE_TRANSITION + D_TOLERANCE ? k0 * (0.07 * H) : k1 * (0.03 * H)
end

@inline _arthern_c(ρ::Float64, T::Float64, tam::Float64, precip_force::Float64) =
    _arthern_scaled_c(ρ, T, tam, precip_force, 1.0, 1.0)

# Medley et al. (2022) GSFC-FDM v1.2.1 eqs. 9-10 with the calibrated parameters of eq. 18.
# Arthern's form with the accumulation raised to a fitted exponent and a separate activation
# energy per stage. `A_pow = pm^α` is stage-dependent, so unlike `_arthern_c` the caller
# hoists two of them; `GRAVITY` is folded in there as well.
#
# The refit replaced Arthern's Ec = 60000 with 59500 (ρ <= 550) and 56870 (ρ > 550), and
# absorbed the site calibration factors R0/R1 of eqs. 16-17 — those do not appear here.
#
# `pm` must be kg m-2 yr-1. For `:Arthern` the accumulation enters linearly and its units
# only rescale a fitted prefactor, but α != 1 here, so the unit convention is load-bearing:
# feeding m i.e. yr-1 would not merely rescale the rate, it would change its accumulation
# dependence. CFM converts to kg m-2 yr-1 for both schemes, which is the convention followed.
@inline function _gsfc2020_c(ρ::Float64, T::Float64, tam::Float64,
    A_pow_g_0::Float64, A_pow_g_1::Float64)
    grain_growth = 42400.0 / (tam * R_GAS)
    return if ρ <= DENSITY_STAGE_TRANSITION + D_TOLERANCE
        0.07 * A_pow_g_0 * exp(-59500.0 / (T * R_GAS) + grain_growth)
    else
        0.03 * A_pow_g_1 * exp(-56870.0 / (T * R_GAS) + grain_growth)
    end
end

# Simonsen et al. (2013): Arthern's form and activation energies (Ec = 60000, Eg = 42400)
# retuned for Greenland, with a constant factor F0 = 0.8 below 550 kg m-3 and F1 = 1.25 above.
# The second stage carries an additional accumulation- and climate-dependent factor
# `γ = 61.7 / sqrt(pm) · exp(-3800/(R·tam))`, which applies to stage 2 only — stage 1 gets F0
# alone. `γ` and both F factors are climatological, so the caller hoists `F0` and `F1·γ`.
#
# The structure is Lundin et al. (2017) FirnMICE eqs. A36-A37, which give
# `c0(SIM) = f0·c0(ART)` and `c1(SIM) = f1·(61.7/ḃ^0.5)·exp(-3800/(R·Ta))·c1(ART)` — i.e. the
# γ factor multiplies the second stage only, and `Ta` is the mean annual temperature. FirnMICE
# states the form but publishes no numeric f0/f1.
#
# The values F0 = 0.8, F1 = 1.25 are the ones the Community Firn Model applies; its in-source
# comment attributes them to correspondence with the author. CFM also carries commented-out
# alternatives F0 = 0.68, F1 = 1.03 labelled "firnmice value?" — those are not in the FirnMICE
# paper, and are not used here.
#
# `pm` must be kg m-2 yr-1: `γ ∝ 1/sqrt(pm)` makes the unit convention load-bearing, as for
# `_gsfc2020_c`. FirnMICE eq. A14 states this unit explicitly for the Arthern base coefficients.
#
# Both `:Simonsen2013` and `:Ligtenberg` (Arthern scaled by the accumulation-dependent M0/M1)
# are evaluated by `_arthern_scaled_c` above; only their hoisted stage scalars differ.

"""
    DensificationCoeffs(mp::ModelParameters, pm, tam)

Loop-invariant densification terms for `mp.densification_method`, hoisted out of
the per-cell / per-age-step loop. Shared by [`calculate_density`](@ref) and
[`steady_state_density`](@ref) so both derive `c` from identical inputs.

`k0` and `k1` are the hoisted per-stage scalars (below and above
`DENSITY_STAGE_TRANSITION`). Every scheme that reaches them needs exactly two, so
their *meaning* is per-method while their arity is not:

| method | `k0`, `k1` |
|---|---|
| `:Arthern`, `:ArthernB`, `:HerronLangway` | unused (1.0) |
| `:Ligtenberg` | the accumulation-dependent `M0`, `M1` |
| `:Simonsen2013` | `F0 = 0.8` and `F1·γ` |
| `:GSFC2020` | `pm^α · g` with α0 = 0.91, α1 = 0.644 |

`pm` is mean annual accumulation [kg m-2 yr-1] and `tam` mean annual air
temperature [K].
"""
struct DensificationCoeffs
    method::Symbol
    pm::Float64
    tam::Float64
    precip_force::Float64
    k0::Float64
    k1::Float64
end

function DensificationCoeffs(mp::ModelParameters, pm::Real, tam::Real)
    pm = Float64(pm)
    tam = Float64(tam)
    method = mp.densification_method

    # Non-positive mean accumulation reaches here from the steady-state marcher on an ablation
    # column. It is not a valid climatology for any of these schemes, so the policy is uniform:
    # the *fractional powers of `pm`* are guarded to keep them real and finite, while
    # `precip_force` is left signed. Guarding `precip_force` too would break the identity that
    # the scaled-Arthern schemes are exactly `k · :Arthern`, since `:Arthern` itself is signed.
    pm_nonneg = max(pm, 0.0)

    k0 = k1 = 1.0
    if method == :GSFC2020
        # `pm` is kg m-2 yr-1 (see `_gsfc2020_c`); `^0.91` of a negative would be NaN.
        k0 = pm_nonneg^0.91 * GRAVITY
        k1 = pm_nonneg^0.644 * GRAVITY
    elseif method == :Simonsen2013
        k0 = 0.8
        # `γ ∝ 1/sqrt(pm)` is Inf at pm = 0, and the rate also carries `precip_force ∝ pm`, so
        # the stage-2 product would be 0·Inf = NaN rather than the 0 the accumulation implies.
        k1 = pm_nonneg > 0.0 ?
            1.25 * (61.7 / sqrt(pm_nonneg) * exp(-3800.0 / (tam * R_GAS))) : 0.0
    elseif method == :Ligtenberg
        M01 = densification_lookup_M01(mp.densification_coeffs_M01)
        if length(M01) == 4
            k0 = max(M01[1] - (M01[2] * log(pm)), 0.25)
            k1 = max(M01[3] - (M01[4] * log(pm)), 0.25)
        elseif abs(mp.density_ice - 820.0) < D_TOLERANCE
            k0 = max(M01[1, 1] - (M01[1, 2] * log(pm)), 0.25)
            k1 = max(M01[1, 3] - (M01[1, 4] * log(pm)), 0.25)
        else
            k0 = max(M01[2, 1] - (M01[2, 2] * log(pm)), 0.25)
            k1 = max(M01[2, 3] - (M01[2, 4] * log(pm)), 0.25)
        end
    end
    return DensificationCoeffs(method, pm, tam, pm * GRAVITY, k0, k1)
end

"""
    _densification_rate(p::DensificationCoeffs, ρ, T) -> c [yr-1]

Densification rate coefficient for the method in `p`, dispatching on
`p.method`. Used by [`steady_state_density`](@ref), where the branch is taken
once per age step rather than once per cell. `calculate_density`'s hot loop calls
the per-method functions directly, so it takes the branch once per call.

`:ArthernB` falls back to `:Arthern` here. The steady-state march carries neither
overburden stress nor grain radius, which eq. B1 needs, and it only produces an
initial guess that the spinup then relaxes — so the accumulation-driven form of
the same paper's eq. 1 is the closest available stand-in. The transient run uses
the real `:ArthernB` law.
"""
@inline function _densification_rate(p::DensificationCoeffs, ρ::Float64, T::Float64)
    if p.method == :HerronLangway
        return ρ <= DENSITY_STAGE_TRANSITION + D_TOLERANCE ? _hl_c0(T, p.pm) : _hl_c1(T, p.pm)
    elseif p.method == :GSFC2020
        return _gsfc2020_c(ρ, T, p.tam, p.k0, p.k1)
    elseif p.method == :Arthern || p.method == :ArthernB ||
           p.method == :Simonsen2013 || p.method == :Ligtenberg
        # One kernel for the whole scaled-Arthern family; `k0`/`k1` are 1.0 for `:Arthern`.
        return _arthern_scaled_c(ρ, T, p.tam, p.precip_force, p.k0, p.k1)
    end
    error("unrecognized densification method: $(p.method)")
end
