"""
    calculate_density(temperature, dz, density, grain_radius, water, cfs::ClimateForcingStep, mp::ModelParameters;
                      grain_dendricity, grain_sphericity, viscosity)

Compute the densification of snow/firn using one of several models:
- `:HerronLangway`: Herron and Langway (1980)
- `:Arthern`: semi-empirical model of Arthern et al. (2010) [default]
- `:ArthernB`: physical model from Appendix B of Arthern et al. (2010)
- `:Crocus`: viscous settling of Vionnet et al. (2012) below `CROCUS_HYBRID_DENSITY`,
  `:GSFC2020` above; the only scheme here in which liquid water weakens the matrix
- `:CrocusPure`: Vionnet et al. (2012) at every density, without the firn handover
- `:GSFC2020`: recalibrated Arthern of Medley et al. (2022), GSFC-FDM v1.2.1
- `:Simonsen2013`: Arthern retuned for Greenland by Simonsen et al. (2013)
- `:Barnola1991`: Herron-Langway stage 1 below `DENSITY_STAGE_TRANSITION`, pressure
  sintering above; the only scheme here that models the firn-ice transition mechanistically
- `:LiZwally`: empirical model of Li and Zwally (2004)
- `:Helsen`: modified empirical model by Helsen et al. (2008)
- `:Ligtenberg`: semi-empirical model of Ligtenberg et al. (2011)

Returns `(dz, density)`. `density` is updated in place; `dz` is returned as a
new array (recomputed from the conserved cell mass).

Pass a length-`m` `viscosity` vector to have the per-cell effective viscosity [Pa s]
written into it. It is filled with `NaN` first and written only by the
`:Crocus`/`:CrocusPure` settling branch, the one place a viscosity is genuinely formed;
under every other scheme, and for the cells `:Crocus` hands to `:GSFC2020`, it stays
`NaN`. Omitting the keyword leaves the hot path untouched: the default is the empty
`NO_VISCOSITY`, and an empty vector means "not requested" (a concrete type on both
settings — see the note in `gemb_core.jl`).

With `mp.destructive_metamorphism` set, the `:Crocus`/`:CrocusPure` viscous branch adds the
destructive-metamorphism settling of fresh snow ([`destructive_metamorphism_rate`](@ref)) to
its overburden strain rate. Both are strain rates `(1/D)·dD/dt`, and van Kampenhout et al.
(2017, Sect. 3.4) and Schneider et al. (2022, p. 3) sum them, so one step is
`D·exp(-(σ/η + r_dm)·dt)` with both rates taken at the start-of-step state. The cells
`:Crocus` hands to `:GSFC2020` (ρ ≥ `CROCUS_HYBRID_DENSITY`) do not get the term; there it is
below `exp(-0.046·275) = 3.2e-6` of its fresh-snow value. Requesting it with any other
densification scheme throws.

`water` [kg m-2] is read only by the stress-driven schemes: `:ArthernB` and `:Barnola1991`
(as part of the overburden) and `:Crocus`/`:CrocusPure` (overburden, eq. 8's viscosity
reduction, and the grain-type gate of eq. 9); the accumulation-driven schemes proxy overburden
with mean accumulation and ignore it.

`grain_dendricity` and `grain_sphericity` (keywords, per cell) are read only by
`:Crocus`/`:CrocusPure`, which apply eq. 9's grain-type stiffening only to cells that
`crocus_f2_applies` admits. Those two schemes throw an `ArgumentError` when either is
missing or of the wrong length; every other scheme ignores them.

The accumulation-driven branches are scalar-loop implementations, numerically identical
element by element to the vectorized form they replaced but without its mask / gather /
broadcast temporaries (`mass_init`, `idx`, `H`, `c0`, `c1`, ...) allocated per call.

# References
- Arthern, R. J., Vaughan, D. G., Rankin, A. M., Mulvaney, R., and Thomas, E. R. (2010).
  In situ measurements of Antarctic snow compaction compared with predictions of models.
  *J. Geophys. Res.* 115, F03011.
- Herron, M. M. and Langway, C. C. (1980). Firn densification: an empirical model.
  *J. Glaciol.* 25, 373-385.
- Li, J. and Zwally, H. J. (2004). Modeling the density variation in the shallow firn layer.
  *Ann. Glaciol.* 38, 309-313.
- Helsen, M. M., van den Broeke, M. R., van de Wal, R. S. W., et al. (2008). Elevation changes
  in Antarctica mainly determined by accumulation variability. *Science* 320, 1626-1629.
- Ligtenberg, S. R. M., Helsen, M. M., and van den Broeke, M. R. (2011). An improved
  semi-empirical model for the densification of Antarctic firn. *The Cryosphere* 5, 809-819.
- Medley, B., Neumann, T. A., Zwally, H. J., Smith, B. E., and Stevens, C. M. (2022).
  Simulations of firn processes over the Greenland and Antarctic ice sheets: 1980-2021.
  *The Cryosphere* 16, 3971-4011.
- Simonsen, S. B., Stenseng, L., Aðalgeirsdóttir, G., et al. (2013). Assessing a multilayered
  dynamic firn-compaction model for Greenland with ASIRAS radar measurements.
  *J. Glaciol.* 59, 545-558.
- Vionnet, V., Brun, E., Morin, S., et al. (2012). The detailed snowpack scheme Crocus and its
  implementation in SURFEX v7.2. *Geosci. Model Dev.* 5, 773-791.
- Lundin, J. M. D., Stevens, C. M., Arthern, R., et al. (2017). Firn Model Intercomparison
  Experiment (FirnMICE). *J. Glaciol.* 63, 401-422. (eqs. A36-A37.)
- Pimienta, P. and Duval, P. (1987). Rate controlling processes in the creep of polar glacier
  ice. *J. Phys. Colloques* 48, C1-243-C1-248.
- Barnola, J.-M., Pimienta, P., Raynaud, D., and Korotkevich, Y. S. (1991). CO2-climate
  relationship as deduced from the Vostok ice core: a re-examination based on new measurements
  and on a re-evaluation of the air dating. *Tellus B* 43, 83-90.
- Goujon, C., Barnola, J.-M., and Ritz, C. (2003). Modeling the densification of polar firn
  including heat diffusion: application to close-off characteristics and gas isotopic
  fractionation for Antarctica and Greenland sites. *J. Geophys. Res.* 108, 4792. (eqs.
  A11-A12, the closed-form pressure term `:Barnola1991` omits — see the note on that branch.)
- van Kampenhout, L., Lenaerts, J. T. M., Lipscomb, W. H., Sacks, W. J., Lawrence, D. M.,
  Slater, A. G., and van den Broeke, M. R. (2017). Improving the representation of polar snow
  and firn in the Community Earth System Model. *J. Adv. Model. Earth Syst.* 9, 2583-2600.
  (eq. 8, the grain-size treatment the Crocus branch adopts; eq. 5, destructive
  metamorphism.)
- Schneider, A. M., Zender, C. S., and Price, S. F. (2022). More realistic intermediate depth
  dry firn densification in the Energy Exascale Earth System Model (E3SM). *J. Adv. Model.
  Earth Syst.* 14, e2021MS002542. (eq. 3, the destructive-metamorphism density taper.)
- Anderson, E. A. (1976). A point energy and mass balance model of a snow cover. NOAA
  Technical Report NWS 19. (Destructive metamorphism.)
- Wilkinson, D. S. and Ashby, M. F. (1975). Pressure sintering by power law creep.
  *Acta Metall.* 23, 1277-1285. (The closed-pore `fs(ρ)` of Barnola's eq. 4.)
- Doake, C. S. M. and Wolff, E. W. (1985). Flow law for ice in polar ice sheets. *Nature* 314,
  255-257. (Cited by Barnola et al. for the *n* = 1 creep exponent — see `_barnola_n`.)
"""
function calculate_density(temperature::Vector{Float64}, dz::Vector{Float64},
    density::Vector{Float64}, grain_radius::Vector{Float64}, water::Vector{Float64},
    cfs::ClimateForcingStep, mp::ModelParameters;
    grain_dendricity::Vector{Float64}=NO_GRAIN_STATE,
    grain_sphericity::Vector{Float64}=NO_GRAIN_STATE,
    viscosity::Vector{Float64}=NO_VISCOSITY)

    # specify constants
    dt = cfs.dt / 86400.0   # convert from [s] to [d]

    m = length(density)
    density_ice = mp.density_ice
    # The two climatological scalars every accumulation- or Arrhenius-driven scheme below
    # compacts against. `:accumulation`/`:arrhenius` select the refined forms; the
    # `:precipitation`/`:arithmetic` settings reproduce the pre-refinement behaviour exactly.
    pm = mp.densification_accumulation === :accumulation ?
         cfs.accumulation_mean : cfs.precipitation_mean
    tam = mp.mean_temperature_method === :arrhenius ?
          cfs.temperature_air_effective : cfs.temperature_air_mean

    # New grid-cell lengths (fresh array; density is updated in place).
    dz_out = similar(dz)

    # Effective viscosity is reported only where a scheme genuinely forms one — the
    # `:Crocus`/`:CrocusPure` `dz = dz₀exp(-σ/η·dt)` branch. Every other scheme is a
    # `ρ̇ = c(ρᵢ - ρ)` relaxation with no stress in it, so `η = σ/(ρ̇/ρ)` would be a
    # reconstruction from a σ the scheme never used, reported in units that invite
    # comparison with a quantity it is not. NaN says "not computed here", which is
    # honest; a number would not be.
    isempty(viscosity) || fill!(viscosity, NaN)

    method = mp.densification_method
    if mp.destructive_metamorphism && !(method == :Crocus || method == :CrocusPure)
        error("destructive_metamorphism requires densification_method :Crocus or " *
              ":CrocusPure, got :$(method)")
    end

    if method == :HerronLangway
        # Rate coefficients `_hl_c0`/`_hl_c1` are shared with `steady_state_density`.
        @inbounds for i in 1:m
            T = temperature[i]
            c = density[i] <= DENSITY_STAGE_TRANSITION + D_TOLERANCE ?
                _hl_c0(T, pm) : _hl_c1(T, pm)
            _densify_cell!(density, dz, dz_out, i, c, dt, density_ice)
        end

    elseif method == :Arthern || method == :Simonsen2013 || method == :Ligtenberg
        # The scaled-Arthern family: one kernel, differing only in the hoisted per-stage
        # scalars `k0`/`k1` that `DensificationCoeffs` computes (1.0 for `:Arthern` itself,
        # Ligtenberg's M0/M1, Simonsen's F0 and F1·γ).
        p = DensificationCoeffs(mp, pm, tam)
        @inbounds for i in 1:m
            c = _arthern_scaled_c(density[i], temperature[i], tam, p.precip_force, p.k0, p.k1)
            _densify_cell!(density, dz, dz_out, i, c, dt, density_ice)
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
            c = (d0 <= DENSITY_STAGE_TRANSITION + D_TOLERANCE ? 9.2e-9 : 3.7e-9) *
                kc_per_year * H
            _densify_cell!(density, dz, dz_out, i, c, dt, density_ice)
            # accumulate the load using this cell's original density
            load += d0 * dz[i] + water[i]
        end

    elseif method == :Barnola1991
        # Barnola et al. (1991) pressure sintering, with Herron & Langway (1980) stage 1
        # below `DENSITY_STAGE_TRANSITION`. The only scheme here that treats the firn-ice
        # transition mechanistically rather than by extrapolating a snow/firn fit.
        #
        # Not in the MATLAB reference. Ported from the published law; see `_barnola_f`.
        #
        # Calibration range: the paper fits eq. 2 to Antarctic and Greenland profiles at
        # -14 to -57 °C and 2.2 to 65 g cm-2 yr-1. Warmer or wetter sites are extrapolation,
        # and the law has no liquid-water term (`:Crocus` is the scheme that does).
        #
        # Integrated as a rate on ρ rather than through `_densify_cell!`: this is a genuine
        # `dρ/dt`, not the `c·(ρᵢ − ρ)` relaxation the accumulation-driven schemes share, and
        # it does not vanish as ρ → ρᵢ. Forward Euler on `dρ/dt` is used directly, with the
        # same ice-density clamp the other branches apply. The paper's own note that the
        # model overshoots ρᵢ by ~1e-3 is what that clamp absorbs.
        #
        # Below 550 the law is exactly `_hl_c0` — see the identity noted on `_barnola_f`.
        #
        # The 550 handover is a rate discontinuity, and a large one. Herron-Langway stage 1 is
        # accumulation-driven and depth-independent, while pressure sintering scales as σ³, so
        # the ratio across 550 depends entirely on the overburden the cell carries. On a 250 K
        # column at 300 kg m-2 yr-1 the sintering rate is ~1.3e4x slower than Herron-Langway
        # just above 550 at 1 m depth, ~17x slower at 5 m, and overtakes it at ~13 m. This is
        # the paper's own construction — "From the surface to ρ = 0.55 g cm-3 the Herron and
        # Langway (1980) model was used" — not a defect here: their columns reach 550 at a depth
        # where the two rates are comparable. On a column whose 550 horizon sits much shallower
        # than that (a high-accumulation or wind-slab site) densification effectively stalls at
        # 550 until enough load accumulates. `:Barnola1991` is not the scheme to use there.
        load = 0.0      # Σ (ρⱼdzⱼ + waterⱼ) over cells above the current one [kg m-2]
        @inbounds for i in 1:m
            d0 = density[i]
            dz0 = dz[i]
            self_load = d0 * dz0 + water[i]
            if d0 <= DENSITY_STAGE_TRANSITION + D_TOLERANCE
                # Herron & Langway stage 1, shared with `:HerronLangway` so the two can
                # never drift.
                _densify_cell!(density, dz, dz_out, i, _hl_c0(temperature[i], pm), dt, density_ice)
            else
                # Cell-midpoint overburden, the same convention `:Crocus` uses: the load of
                # everything above plus half this cell's own weight.
                #
                # The paper's σ is `ΔP`, "the pressure due to the overburden load, minus the
                # bubble pressure". Only the overburden is applied here. In open-pore firn the
                # pore space is atmospherically connected, so the bubble pressure is the
                # surface pressure and ΔP is the overburden as measured against it; the term
                # only becomes a real back-pressure once pores close. This branch therefore
                # overstates the effective stress above `DENSITY_PORE_CLOSEOFF`.
                #
                # Not implemented because this paper gives no expression for it. Goujon et al.
                # (2003) eqs. A11-A12 do, in closed form and without new state:
                # `P_eff = P + P_atm - P_b` with `P_b = P_c·D(1-D_c)/(D_c(1-D))`, needing only
                # the close-off density and the overburden at close-off. Adding it here would
                # mean grafting another paper's law onto this one's calibration, so it belongs
                # in a `:Goujon2003` scheme rather than as a correction to this branch.
                σ = (load + 0.5 * self_load) * GRAVITY          # [Pa]
                # `σ*σ*σ` rather than `σ^BARNOLA_N`: the latter routes through `pow` for a
                # non-`Integer` exponent. `BARNOLA_N` documents the exponent and is pinned by
                # the test suite; the two must be kept in step.
                dρ_dt = d0 * BARNOLA_A0 * exp(-BARNOLA_Q / (temperature[i] * R_GAS)) *
                        _barnola_f(d0, density_ice) * σ * σ * σ        # [kg m-3 s-1]
                d = d0 + dρ_dt * cfs.dt
                if d > density_ice - D_TOLERANCE
                    d = density_ice
                end
                density[i] = d
                dz_out[i] = d0 * dz0 / d
            end
            load += self_load
        end

    elseif method == :Crocus || method == :CrocusPure
        # Crocus viscous settling, Vionnet et al. (2012) eqs. 5-9. The second stress-driven
        # scheme here, and the only one in which liquid water weakens the matrix.
        #
        # Not in the MATLAB reference. Ported from the paper; see `_crocus_viscosity`.
        #
        # Integrated in thickness rather than through `_densify_cell!`: eq. 5 is a strain rate
        # on `D`, not the `c·(ρᵢ − ρ)` relaxation the other schemes share. Over a step of
        # constant σ and η its exact integral is `D exp(-σ/η·dt)`, which is used directly —
        # unconditionally positive, where a forward-Euler step on the same law is not.
        #
        # `:Crocus` hands cells at or above `CROCUS_HYBRID_DENSITY` to `:GSFC2020`, because
        # eq. 7 is fitted to a 1-2 m alpine snowpack and `exp(bη·ρ)` saturates in firn —
        # see `CROCUS_HYBRID_DENSITY`. `:CrocusPure` applies eq. 5 at every density.
        (length(grain_dendricity) == m && length(grain_sphericity) == m) ||
            throw(ArgumentError(":$(method) densification gates eq. 9 on grain shape and needs " *
                                "grain_dendricity and grain_sphericity for all $m cells; got " *
                                "$(length(grain_dendricity)) and $(length(grain_sphericity))"))
        hybrid = method == :Crocus
        p = hybrid ? DensificationCoeffs(mp, pm, tam) : nothing
        dt_seconds = cfs.dt
        settling = mp.destructive_metamorphism
        density_dm = mp.destructive_metamorphism_density
        load = 0.0       # Σ (ρⱼdzⱼ + waterⱼ) over cells above the current one [kg m-2]
        @inbounds for i in 1:m
            d0 = density[i]
            dz0 = dz[i]
            self_load = d0 * dz0 + water[i]
            if hybrid && d0 >= CROCUS_HYBRID_DENSITY - D_TOLERANCE
                c = _gsfc2020_c(d0, temperature[i], tam, p.k0, p.k1)
                _densify_cell!(density, dz, dz_out, i, c, dt, density_ice)
            else
                # Eq. 6 with cos(Θ) = 1 (GEMB carries no local slope), evaluated at the cell
                # midpoint: every cell carries the overlying load plus half its own weight.
                # Crocus adds the half-own-weight term to the uppermost layer only (the paper's
                # note under eq. 6; `snowcro.F90` SNOWCROCOMPACTN), so a lower cell there is
                # loaded at its top. The midpoint stress converges to the continuous profile
                # from closer as cells coarsen; with fine layering the two agree. For the first
                # cell the self term is the difference between compacting and not.
                σ = (load + 0.5 * self_load) * GRAVITY          # [Pa]
                η = _crocus_viscosity(d0, temperature[i], water[i], dz0, grain_radius[i],
                    crocus_f2_applies(grain_dendricity[i], grain_sphericity[i], water[i],
                        grain_radius[i]))
                isempty(viscosity) || (viscosity[i] = η)
                strain_rate = σ / η                              # [s-1]
                if settling
                    strain_rate += destructive_metamorphism_rate(d0, temperature[i],
                        water[i] > WATER_TOLERANCE, density_dm)
                end
                dz_new = dz0 * exp(-strain_rate * dt_seconds)
                d = d0 * dz0 / dz_new
                if d > density_ice - D_TOLERANCE
                    d = density_ice
                    dz_new = d0 * dz0 / d
                end
                density[i] = d
                dz_out[i] = dz_new
            end
            load += self_load
        end

    elseif method == :GSFC2020
        # Medley et al. (2022). Both accumulation powers are loop-invariant, so they are
        # hoisted with the rest via DensificationCoeffs.
        p = DensificationCoeffs(mp, pm, tam)
        @inbounds for i in 1:m
            c = _gsfc2020_c(density[i], temperature[i], tam, p.k0, p.k1)
            _densify_cell!(density, dz, dz_out, i, c, dt, density_ice)
        end

    elseif method == :LiZwally || method == :Helsen
        # Li and Zwally (2004) and Helsen et al. (2008) share a form and differ only in the
        # two fitted coefficients of the mean-temperature term. Both remain validator-gated.
        base = (pm / density_ice) * 8.36 * (method == :LiZwally ?
            max(139.21 - 0.542 * tam, 1.0) : max(76.138 - 0.28965 * tam, 1.0))
        @inbounds for i in 1:m
            c = base * max(CtoK - temperature[i], 1.0)^(-2.061)
            _densify_cell!(density, dz, dz_out, i, c, dt, density_ice)
        end

    else
        error("unrecognized densification method")
    end

    return dz_out, density
end

"""
    _densify_cell!(density, dz_in, dz_out, i, c, dt, density_ice)

Apply the densification increment for cell `i` given rate coefficient `c`:
update `density[i]` in place, clamp it to the density of ice, and write the
mass-conserving new grid-cell length to `dz_out[i]`.
"""
@inline function _densify_cell!(density::Vector{Float64}, dz_in::Vector{Float64},
    dz_out::Vector{Float64}, i::Int, c::Float64, dt::Float64, density_ice::Float64)
    @inbounds begin
        d0 = density[i]
        mass = d0 * dz_in[i]
        d = d0 + (c * (density_ice - d0) / DAYS_PER_YEAR_DENSIFICATION * dt)
        if d > density_ice - D_TOLERANCE
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
# `_hl_k0`/`_hl_k1` are the temperature-only rate constants k₀ and k₁ of Herron & Langway's eq. 7.
@inline _hl_k0(T::Float64) = 11.0 * exp(-10160.0 / (T * R_GAS))
@inline _hl_k1(T::Float64) = 575.0 * exp(-21400.0 / (T * R_GAS))
@inline _hl_c0(T::Float64, pm::Float64) = _hl_k0(T) * (pm / 1000.0)          # ρ <= 550
@inline _hl_c1(T::Float64, pm::Float64) = _hl_k1(T) * sqrt(pm / 1000.0)      # ρ > 550

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
    H = exp((-60000.0 / (T * R_GAS)) + (GRAIN_GROWTH_EG / (tam * R_GAS))) * precip_force
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
    grain_growth = GRAIN_GROWTH_EG / (tam * R_GAS)
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

# ---------------------------------------------------------------------------
# Barnola et al. (1991) pressure sintering
# ---------------------------------------------------------------------------

# Creep parameters of Barnola et al. (1991) eq. 2. The paper publishes
# `A0 = 2.54e4 MPa-3 s-1`; the `/1e18` below is the MPa-3 -> Pa-3 conversion (1 MPa-3 =
# 1e-18 Pa-3), so with σ in Pa and `n = 3`, `ρ·A0·exp(-Q/RT)·f·σ³` is kg m-3 s-1 directly —
# no year conversion, unlike `:ArthernB`.
const BARNOLA_A0 = 2.54e4 / 1.0e18   # [Pa-3 s-1], from 2.54e4 MPa-3 s-1
const BARNOLA_Q = 60.0e3             # ice lattice diffusion activation energy [J mol-1]

"""
    BARNOLA_N

Stress exponent in `dρ/dt = ρ·A0·exp(-Q/RT)·f(ρ)·σⁿ`.

Fixed at 3, matching Barnola et al. (1991) for the whole firn range this scheme covers:
"*n* taken to be equal to 3, as the effective stress in the firn is rapidly higher than
0.1 MPa". Their `A0` was fitted against `n = 3` over 0.55-0.8 g cm-3, so the exponent and
the prefactor are one calibration.

The paper's adjacent remark that "the exponent *n* is 1 when the effective stress is lower
than 0.1 MPa" is *not* an unimplemented firn branch, though it reads like one. It scopes
itself to bulk ice: "**below the close-off**, a good fit to the experimental data is obtained
by taking successively *n* = 3 and *n* = 1 (Pimienta, 1987)", citing Doake and Wolff (1985)
and Pimienta and Duval (1987). Those are shear-creep studies of dense polar ice — Pimienta
and Duval torsion-tested 0.85 g cm-3 samples and analysed Dye 3 inclinometry at 1000-1784 m,
and report no densification rate at all — so they supply no volumetric prefactor for porous
firn, and none is published here to pair with `n = 1`.

Three reasons not to add the branch, beyond the missing prefactor:

  - Pimienta and Duval's own exponent is not 1. Their abstract and conclusion both give
    "smaller than 2"; regressing their Table 1 (South Pole ice, -15 °C) gives n ≈ 1.55.
    Barnola's "*n* is 1" rounds a sub-2 shear result.
  - Matching an `n = 1` branch continuously at 0.1 MPa forces `A1 = A0·(1e5)²`, so the rate
    would rise by `(1e5/σ)²` below it — ×2 at 70 kPa, ×41 at 15 kPa, ×100 at 10 kPa.
  - Barnola reproduced observed 0.55-0.83 g cm-3 profiles at Vostok and other Antarctic and
    Greenland sites with `n = 3` throughout that range. A 40× shallow acceleration would
    break the agreement the coefficients exist to produce.

The Community Firn Model reaches the same conclusion, carrying the switch as dead code
(`# nBa[sigmaEff<1.0e5]=1.0`).
"""
const BARNOLA_N = 3.0

# Coefficients of the log10 f(ρ) polynomial in `_barnola_f`, in the paper's g cm-3.
const BARNOLA_ALPHA = -37.455
const BARNOLA_BETA = 99.743
const BARNOLA_DELTA = -95.027
const BARNOLA_GAMMA = 30.673

"""
    BARNOLA_CLOSEOFF

Density [kg m-3] at which `_barnola_f` switches from the fitted polynomial to the
analytic closed-pore expression.

800, not GEMB's `DENSITY_PORE_CLOSEOFF = 830`: this is the density at which Barnola's
*polynomial fit* is superseded by the closed-pore geometry, which is a property of that fit
rather than of pore close-off as GEMB's other code means it. The polynomial is fitted to
Pimienta and Duval's data over roughly 0.55-0.8 g cm-3 and diverges from the analytic form
beyond it — at ρ = 900 it gives 0.0432 against the analytic 0.0087, a factor of 5 — which is
why the handover exists rather than running the polynomial to ρᵢ.

The handover is meant to be smooth. The paper states that the polynomial "was calculated in
order to make the two functions, `fe(ρ)` and `fs(ρ)`, and their first derivatives equal for
ρ = 0.8 g cm-3" — so both value and slope should match there, and any step is a mismatch
rather than an intended feature of the law.

The size of that step depends on `mp.density_ice`, because only the closed-pore branch
carries it. The two branches cross at ρᵢ = 919.96, and their first derivatives cross at
920.06 — the paper's C¹ statement is what fixes both, and it identifies the fit's own ice
density as ~920 kg m-3. The branches agree to 0.06% at 920 against 4.3% at 917 and 14.1% at
GEMB's default 910, so the gap at any other `density_ice` measures the mismatch between the
configured value and the fit's. At the default 910 the rate steps *down* by 14% crossing 800,
which is a discontinuity in `dρ/dt` but not in ρ, and is small against the factor-5 error the
alternative would introduce deeper in the column.

The sign of that step is not fixed: it is downward below the fit's ~920 and upward above it
(+7.6% at 925, +49.5% at 950). See [`BARNOLA_MIN_DENSITY_ICE`](@ref) for why the low side is gated
and the high side is not.
"""
const BARNOLA_CLOSEOFF = 800.0

"""
    BARNOLA_MIN_DENSITY_ICE

Lowest `mp.density_ice` [kg m-3] that `:Barnola1991` accepts, enforced in
`validate_parameters`. The rest of the model permits `[800, 950]`.

The gate exists because `BARNOLA_CLOSEOFF` is an absolute density while the polynomial
branch below it carries no `density_ice`. As the configured ice density falls toward 800 the
handover moves into firn whose true porosity is nearly zero, and the polynomial — which
cannot be rescaled, its coefficients being fitted against the fit's own ρᵢ ≈ 920 — overstates
`f` there: by 5.9x at ρ = 700 and 15.9x at ρ = 799.9 for `density_ice = 820`. At exactly 800
the polynomial covers the entire column, `f` never reaches 0, and the scheme loses the
self-limiting behaviour that motivates it, leaving only the ice-density clamp.

900 admits GEMB's default 910 (a -14% step at the handover) and the common 917 (-4%) while
refusing the range where the law degrades: the step is -27% at 900 and -94% at 820.

Deliberately one-sided. Above ~920 the step inverts rather than growing without bound
(+7.6% at 925, +49.5% at 950) and the polynomial stays inside its fitted porosity range, so
those configurations are merely inconsistent with the fit rather than unphysical. They are
documented on [`BARNOLA_CLOSEOFF`](@ref) instead of rejected.
"""
const BARNOLA_MIN_DENSITY_ICE = 900.0

"""
    _barnola_f(ρ, density_ice) -> f [-]

Densification factor of Barnola et al. (1991), the geometric term in their eq. 2,

    dρ/dt = ρ · A0 · exp(-Q/RT) · f(ρ) · ΔPⁿ

Two regimes, meeting at [`BARNOLA_CLOSEOFF`](@ref):

  - Open pores, `ρ <= 800`: the paper's eq. 3, `fe(ρ)`, an empirical fit
    `log10 f = α(ρ/1000)³ + β(ρ/1000)² + δ(ρ/1000) + γ` "empirically deduced for the
    0.55-0.8 g cm-3 density range", parameterizing the Pimienta and Duval (1987)
    pressure-sintering data. Note the polynomial is in g cm-3, so the coefficients are only
    meaningful against `ρ/1000`.
  - Closed pores, `ρ > 800`: the paper's eq. 4, `fs(ρ)`, the Wilkinson and Ashby (1975)
    spherical pore model, `f = (3/16)·(1 − ρ/ρᵢ) / (1 − (1 − ρ/ρᵢ)^(1/3))³`, which the paper
    describes as "already valid below ρ = 0.8 g cm-3".

`ρᵢ` is `mp.density_ice`, so the closed-pore branch follows the configured ice density.
The polynomial branch cannot: its coefficients were fitted with a fixed ice density (~920
kg m-3; see [`BARNOLA_CLOSEOFF`](@ref)), and rescaling them is not something the paper licenses.

Below `DENSITY_STAGE_TRANSITION` the scheme does not use this factor at all — it uses
Herron & Langway stage 1, matching the paper ("from the surface to ρ = 0.55 g cm-3 the
Herron and Langway (1980) model was used"), and that is exactly GEMB's `_hl_c0`. The
identity was verified
against the Community Firn Model's `Barnola1991` zone 1 to a ratio of 1.0 (12 significant
figures) across three (T, ρ, accumulation) triples, so `:Barnola1991` and `:HerronLangway`
share the one kernel rather than restating it.

Cross-checked against the Community Firn Model's `Barnola1991`, which agrees on both
branches and on all four polynomial coefficients. CFM applies the full overburden `σ`
where this uses the cell-midpoint value, matching the convention `:Crocus` already uses
here; and CFM's zone 1 carries an `A^aHL` accumulation term with `aHL = 1`, which is the
linear accumulation dependence already inside `_hl_c0`.
"""
@inline function _barnola_f(ρ::Float64, density_ice::Float64)
    if ρ <= BARNOLA_CLOSEOFF
        ρ_cgs = ρ / 1000.0
        return exp10(BARNOLA_ALPHA * ρ_cgs^3 + BARNOLA_BETA * ρ_cgs^2 +
                     BARNOLA_DELTA * ρ_cgs + BARNOLA_GAMMA)
    else
        porosity = 1.0 - ρ / density_ice
        # `porosity -> 0` as ρ -> ρᵢ; the cube of the cube-root difference vanishes with it,
        # and the ratio tends to 9/16 · porosity^(1/3) · ... -> 0. Guard the exact endpoint,
        # where both numerator and denominator are 0.
        #
        # Reaching this branch at all requires `density_ice > BARNOLA_CLOSEOFF`, which
        # `BARNOLA_MIN_DENSITY_ICE` guarantees for `:Barnola1991`. Without that gate a
        # configured ice density at or below 800 would send every cell down the polynomial
        # branch instead, and `f` would never vanish.
        porosity <= 0.0 && return 0.0
        return (3.0 / 16.0) * porosity / (1.0 - cbrt(porosity))^3
    end
end

# ---------------------------------------------------------------------------
# Crocus viscous settling (Vionnet et al., 2012, GMD 5, 773-791, eqs. 7-9)
# ---------------------------------------------------------------------------

# Eq. 7 coefficients. `η0` in kg s-1 m-1, so `η` is too, and eq. 5's `σ/η` is s-1.
const CROCUS_ETA0 = 7.62237e6   # [kg s-1 m-1]
const CROCUS_A_ETA = 0.1        # [K-1]
const CROCUS_B_ETA = 0.023      # [m3 kg-1]
const CROCUS_C_ETA = 250.0      # [kg m-3]

# Eq. 9 grain-size thresholds, in metres (the paper quotes mm).
const CROCUS_G1 = 4.0e-4        # cap on the grain-size excess entering the exponential
const CROCUS_G2 = 2.0e-4        # grain size at which the exponential is 1
const CROCUS_G3 = 1.0e-4        # e-folding grain size
# Eq. 9's explicit ceiling, reached at gs − g2 = ln(4)·g3, i.e. gs = 0.339 mm. It therefore
# binds before `CROCUS_G1` ever does (that would need gs >= 0.6 mm), so `g1` is redundant in
# the published form; it is kept to keep the expression recognizable against the paper.
const CROCUS_F2_MAX = 4.0
# Floor, not in the paper. See `_crocus_viscosity`.
const CROCUS_F2_MIN = 1.0
# Sphericity below which a non-dendritic cell counts as angular for eq. 9: Crocus's `XVGRAN6`
# (`modd_snow_metamo.F90:128` of the Crocus 3.0.1 release), the test at `snowcro.F90:2114`.
const CROCUS_F2_SPHERICITY_MAX = 51 / 99

"""
    crocus_f2_applies(dendricity, sphericity, water, grain_radius) -> Bool

Whether eq. 9's grain-type factor `f2` of the Crocus viscosity applies to a cell. It applies when
the cell is

- non-dendritic: `dendricity` at most `GDN_TOLERANCE` (the classification
  [`calculate_grain_size`](@ref) uses) *and* an optical diameter `d` of at least
  `1e-4·(4 − S)` m, Crocus's own non-dendritic test, with `d` from [`optical_diameter`](@ref)
  (Lafaysse et al., 2026, eq. D5, GEMB's grain size `2·grain_radius`);
- angular: `sphericity < CROCUS_F2_SPHERICITY_MAX = 51/99`;
- dry: `water` at most `WATER_TOLERANCE` [kg m-2].

Every other cell gets `f2 = 1`. The optical-diameter test admits only grains of at least
0.33 mm (0.40 mm at sphericity 0), where eq. 9 is 3.7 to 4: angular snow finer than that, which
GEMB's `:FlannerZender2006` radius readily produces, is dendritic snow to Crocus and gets no
stiffening there either.

Crocus 3.0.1 (`snowcro.F90:2112-2121`, `SNOWCROCOMPACTN`) applies `f2` only to a layer that is
non-dendritic by that optical-diameter test, has sphericity below 51/99, *and* carries the
historical tracker value 1: faceted crystals at some time and liquid water never (Lafaysse et al.,
2026, Table 2). The tracker is set to 1 when a dry non-dendritic layer first reaches sphericity
below 0.1/99 (`snowcro.F90:2596-2598`), and wetting that rounds the grains moves it to 3, after
which `f2` never applies again (`snowcro.F90:2600-2603`). Neither Vionnet et al. (2012, eq. 9) nor
Lafaysse et al. (2026, eq. 55) prints the history condition; the code is the reference
implementation.

GEMB carries no historical tracker, so this gate reads the cell's present state. It differs
from Crocus in three ways:

- coarse angular snow that never reached sphericity 0 is stiffened here and not in Crocus;
- a layer that was wet and has refrozen and faceted again is stiffened here and not in
  Crocus, which excludes it for good;
- a faceted layer that holds water but has not yet rounded to sphericity 1 is excluded here
  and stiffened in Crocus, where `f1` then dominates its viscosity anyway.

The first two err toward stiffening, the third toward softening.
"""
@inline crocus_f2_applies(dendricity::Float64, sphericity::Float64, water::Float64,
    grain_radius::Float64) =
    dendricity <= GDN_TOLERANCE && sphericity < CROCUS_F2_SPHERICITY_MAX &&
    water <= WATER_TOLERANCE &&
    optical_diameter(dendricity, sphericity, grain_radius) >=
        OPTICAL_DIAMETER_DENDRITIC_UNIT * (4 - sphericity)

"""
    CROCUS_HYBRID_DENSITY

Density [kg m-3] at which `:Crocus` hands a cell over to `:GSFC2020`.

Eq. 7 is fitted to a 1-2 m alpine snowpack, and `exp(bη·ρ)` reaches ~1e9 at ρ = 900: on a
Greenland firn profile the pure law gives 0.7-0.8 of `:Arthern`'s compaction rate in the top
few metres but only 0.02-0.09 below 20 m. That is the published law behaving as fitted, not
a unit error, so `:Crocus` is a composite — Crocus viscosity where liquid water and grain
type govern settling, `:GSFC2020` where creep does.

450 kg m-3 is the threshold the Community Firn Model's `Crocus` scheme uses for the same
purpose (its `RHO_CC`), and it sits between the two laws' calibration ranges rather than
inside either. The blend is a discontinuous handover, not a weighted average: nothing in
either paper prescribes a blending function. `:CrocusPure` applies eq. 5 at every density.
"""
const CROCUS_HYBRID_DENSITY = 450.0

"""
    _crocus_viscosity(ρ, T, W_liq, D, grain_radius, f2_applies=true) -> η [kg s-1 m-1]

Snow viscosity of Vionnet et al. (2012) eq. 7, with the two microstructure correction
factors of eqs. 8 and 9:

- `f1 = 1/(1 + 60·W_liq/(ρ_w·D))` (eq. 8) reduces viscosity in the presence of liquid
  water, bottoming out around 1/61 for a saturated cell. This is the term GEMB has no
  other expression of: every accumulation-driven scheme here is blind to `water`, and
  `:ArthernB` reads it only as overburden.
- `f2 = min(4, exp(min(g1, gs − g2)/g3))` (eq. 9) *raises* viscosity for coarse angular
  grains, suppressing compaction in depth hoar. It applies only when `f2_applies` is true;
  otherwise `f2 = 1`. `calculate_density` sets `f2_applies` from
  `crocus_f2_applies`: non-dendritic by Crocus's optical-diameter test, angular (sphericity
  below 51/99) and dry, the cells Crocus scores with eq. 9.

`W_liq` is cell water [kg m-2], `D` cell thickness [m], `grain_radius` [mm]; `gs` in eq. 9
is a grain size, taken here as the grain diameter `2·grain_radius` in metres.

Where it applies, `f2` is also floored at 1, which the published eq. 9 does not do. Eq. 9 is
bounded above (`min(4.0, ...)`) but not below: it crosses 1 at `gs = g2 = 0.2 mm` and decays
as `exp((gs − g2)/g3)` beneath that, so a fine-grained cell would get `f2 < 1`, a softening from
a factor introduced for "the increase of viscosity with angular grains". The gate admits only
grains of at least 0.33 mm, where eq. 9 is 3.7 to 4, so the floor binds only for a direct
`f2_applies = true` call on a finer grain.

Deviates from the Community Firn Model's `Crocus`, which was used to cross-check eq. 7 and
the `dρ/dt = ρσ/η` form: CFM hardcodes `f2 = 4.0` (eq. 9's ceiling, so it stiffens every
cell including fresh fine-grained snow) and adopts van Kampenhout et al. (2017) eq. 8's
retuned `cη = 358`. The paper's `f2`, Crocus's gate and `cη = 250` are used here.
"""
@inline function _crocus_viscosity(ρ::Float64, T::Float64, W_liq::Float64,
    D::Float64, grain_radius::Float64, f2_applies::Bool=true)
    f1 = 1.0 / (1.0 + 60.0 * W_liq / (DENSITY_WATER * D))
    gs = 2.0 * grain_radius / 1000.0        # radius [mm] -> diameter [m]
    f2 = f2_applies ? clamp(exp(min(CROCUS_G1, gs - CROCUS_G2) / CROCUS_G3),
                            CROCUS_F2_MIN, CROCUS_F2_MAX) : 1.0
    return f1 * f2 * CROCUS_ETA0 * (ρ / CROCUS_C_ETA) *
           exp(CROCUS_A_ETA * (CtoK - T) + CROCUS_B_ETA * ρ)
end

"""
    destructive_metamorphism_rate(density, temperature, wet, density_dm) -> r [s-1]

Compaction rate `-(1/D)·dD/dt` of a snow layer by destructive metamorphism, the settling of
fresh snow as its crystals lose their branches. Anderson (1976) in the form of van Kampenhout
et al. (2017) eq. 5, p. 2587, with the density taper of Schneider et al. (2022) eq. 3, p. 3:

    r = c3·c2·c1·exp(-c4·(T_frz - min(T, T_frz)))
    c1 = 1                           ρ < ρ_dm
         exp(-c5·(ρ - ρ_dm))          ρ ≥ ρ_dm
    c2 = 2 with liquid water, else 1

with `c3 = 2.777e-6 s-1`, `c4 = 0.04 K-1` (van Kampenhout et al., p. 2587) and
`c5 = 0.046 m3 kg-1` (Schneider et al., p. 3). `density` is the ice-matrix density of the
layer, excluding liquid water, as eq. 5 specifies; `ρ_dm` is
`mp.destructive_metamorphism_density`, 175 kg m-3 in CLM5. `wet` selects `c2 = 2`, which the
papers state as "depending on the presence of liquid water" without a threshold; any pore
water above `WATER_TOLERANCE` counts.

The paper prints eq. 5 as `∂ρ/∂t` but its units are s-1: it is the fractional rate that
Schneider et al. eq. 3 writes as a strain rate, and the one the Community Land Model sums with
the overburden strain rate. Temperature is capped at `T_frz`, so wet snow at the melting point
takes the largest rate, `2·c3`. Below `ρ_dm` the rate is `c3·exp(-c4·ΔT)` and does not depend
on density; above it the rate falls by a factor `e` every 21.7 kg m-3, so the term moves fresh
snow to `ρ_dm` within days and then stops mattering.

Throws an `ArgumentError` for a non-positive density or a non-finite or non-positive
temperature.
"""
@inline function destructive_metamorphism_rate(density::Float64, temperature::Float64,
    wet::Bool, density_dm::Float64)
    density > 0 ||
        throw(ArgumentError("density must be positive [kg m-3], got $density"))
    (isfinite(temperature) && temperature > 0) ||
        throw(ArgumentError("temperature must be finite and positive [K], got $temperature"))
    c1 = density < density_dm ? 1.0 :
         exp(-DESTRUCTIVE_METAMORPHISM_C5 * (density - density_dm))
    c2 = wet ? DESTRUCTIVE_METAMORPHISM_WET_FACTOR : 1.0
    return DESTRUCTIVE_METAMORPHISM_C3 * c2 * c1 *
           exp(-DESTRUCTIVE_METAMORPHISM_C4 * (CtoK - min(temperature, CtoK)))
end

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
| `:Arthern`, `:ArthernB`, `:Barnola1991`, `:CrocusPure`, `:HerronLangway` | unused (1.0) |
| `:Ligtenberg` | the accumulation-dependent `M0`, `M1` |
| `:Simonsen2013` | `F0 = 0.8` and `F1·γ` |
| `:GSFC2020`, `:Crocus` | `pm^α · g` with α0 = 0.91, α1 = 0.644 |

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
    if method == :GSFC2020 || method == :Crocus
        # `:Crocus` is a hybrid whose above-threshold branch *is* `:GSFC2020`, so it needs
        # that scheme's hoisted accumulation powers. `:CrocusPure` has no such branch.
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

The stress-driven schemes cannot be evaluated here: the march carries neither overburden
stress, grain radius, nor water, and it only produces an initial guess that the spinup then
relaxes. Each therefore falls back to the nearest accumulation-driven law — `:Crocus` to
`:GSFC2020`, which is already its own above-threshold branch and so governs most of the
column the march builds; `:Barnola1991` to `:HerronLangway`, whose stage 1 *is* its own
below-threshold branch, so the two agree exactly below `DENSITY_STAGE_TRANSITION` and the
fallback is only an approximation above it; `:ArthernB` and `:CrocusPure` to `:Arthern`. The
transient run uses the real law in every case.
"""
@inline function _densification_rate(p::DensificationCoeffs, ρ::Float64, T::Float64)
    if p.method == :HerronLangway || p.method == :Barnola1991
        return ρ <= DENSITY_STAGE_TRANSITION + D_TOLERANCE ? _hl_c0(T, p.pm) : _hl_c1(T, p.pm)
    elseif p.method == :GSFC2020 || p.method == :Crocus
        return _gsfc2020_c(ρ, T, p.tam, p.k0, p.k1)
    elseif p.method == :Arthern || p.method == :ArthernB || p.method == :CrocusPure ||
           p.method == :Simonsen2013 || p.method == :Ligtenberg
        # One kernel for the whole scaled-Arthern family; `k0`/`k1` are 1.0 for `:Arthern`.
        return _arthern_scaled_c(ρ, T, p.tam, p.precip_force, p.k0, p.k1)
    end
    error("unrecognized densification method: $(p.method)")
end
