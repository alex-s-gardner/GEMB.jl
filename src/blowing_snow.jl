"""
Wind rework of snow that has already landed.

Wind enters GEMB elsewhere only through the turbulent fluxes, the `:Pahaut` fresh-snow
density, and the Crocus wind-dependent fresh-grain properties in
[`calculate_accumulation`](@ref) — all of them acting at or before the moment of deposition.
Everything wind does to snow *after* it is on the ground is here: wind-slab compaction, grain
fragmentation, and (optionally) sublimation of the suspended fraction.

This matters on the katabatic-scoured parts of both ice sheets, where blowing snow is a
first-order surface-mass-balance term. Noël et al. (2018) found that halving RACMO2's
saltation coefficient changed drifting-snow sublimation by order 50 mm w.e. yr-1 over
Greenland; Purdie et al. (2011) found wind deflation overriding the orographic precipitation
gradient entirely on maritime glaciers.

# Two independent paths

1. **Computed** ([`apply_blowing_snow!`](@ref), `mp.blowing_snow_method = :Crocus`) — the
   SURFEX/Crocus `SNOWDRIFT` scheme (Vionnet et al., 2012; Lafaysse et al., 2026 eqs. 59-66):
   a mobility index per layer from its microstructure and density, a driftability index from
   mobility and the 5 m wind, an exponential decay of driftability with the overburden it must
   penetrate, then compaction toward `DRIFT_DENSITY_MAX` and grain fragmentation at the
   resulting rate. **Mass-conserving**: `dz` shrinks in exact proportion to the density rise,
   so column mass is untouched and neither grid invariant is at stake. Use this when no drift
   product is available.
2. **Prescribed** ([`apply_prescribed_drift!`](@ref), `mp.drift_rate` or the optional
   `snow_drift` forcing layer) — the IMAU-FDM treatment, which reads `SnowDrif` from RACMO and
   applies it as a surface mass sink (`source/firn_physics.f90`, `Update_Surface`), eroding at
   the surface cell's own density and depositing at the fresh-snow density. All drift physics
   is upstream in the atmospheric model (Lenaerts et al., 2012). Use this when a drift product
   *is* available; it is the only path comparable against RACMO-forced firn models.

Optionally, path 1 also sublimates suspended snow ([`blowing_snow_sublimation_rate`](@ref),
`mp.blowing_snow_sublimation`), after Gordon et al. (2006) as Crocus implements it. That is
the only part of path 1 that removes mass.

The Community Firn Model has no blowing-snow physics at all, so there is no third
implementation to check against; the SURFEX source (`src/SURFEX/snowcro.F90:4685-4950`,
constants in `src/SURFEX/modd_snow_par.F90:290-326`) is what the port was verified line by
line against.

None of this runs at the defaults (`blowing_snow_method = :none`, `drift_rate = 0.0`,
`blowing_snow_sublimation = false`), which leaves output bit-identical to a run without the
term.

# References
- Guyomarc'h, G. & Mérindol, L. (1998). Validation of an application for forecasting blowing
  snow. *Ann. Glaciol.* 26, 138-143.
- Gordon, M., Simon, K. & Taylor, P. A. (2006). On snow depth predictions with the Canadian
  land surface scheme including a parametrization of blowing snow sublimation.
  *Atmos.-Ocean* 44, 239-255. https://doi.org/10.3137/ao.440303
- Vionnet, V., et al. (2012). The detailed snowpack scheme Crocus and its implementation in
  SURFEX v7.2. *Geosci. Model Dev.* 5, 773-791.
- Lenaerts, J. T. M., et al. (2012). Modeling drifting snow in Antarctica with a regional
  climate model: 1. Methods and model evaluation. *J. Geophys. Res.* 117, D05108.
- Lafaysse, M., et al. (2026). Version 3.0.2 of the Crocus snowpack model.
  *Geosci. Model Dev.* 19, 6273-6334.
"""

# Crocus SNOWDRIFT calibration constants, from `src/SURFEX/modd_snow_par.F90:290-326`. The
# Fortran name of each is given so the port can be re-checked against the source; the paper
# symbols (Lafaysse et al., 2026 Appendix K3) are given where the paper names them.
const DRIFT_TAU = 172800.0          # XVTIME, τ_DRIFT: drift compaction timescale [s] (2 d)
const DRIFT_DENSITY_MAX = 350.0     # XVROMAX, ρ_MAX: density a wind slab tends to [kg m-3]
const DRIFT_DENSITY_MIN = 50.0      # XVROMIN, ρ_min: density below which F(ρ) saturates
const DRIFT_MOB1 = 0.295            # XVMOB1: sets the density slope of F(ρ)
const DRIFT_MOB2 = 0.833            # XVMOB2: sphericity weight in the non-dendritic mobility
const DRIFT_MOB3 = 0.583            # XVMOB3: grain-size weight in the non-dendritic mobility
const DRIFT_MOB4 = -0.0583          # XVMOB4: mobility ceiling for snow that has been wet
const DRIFT_D1 = 2.868              # XVDRIFT1, c_SUBL: driftability wind-threshold prefactor
const DRIFT_D2 = 0.085              # XVDRIFT2, d_SUBL: driftability wind-threshold decay
const DRIFT_D3 = 3.25               # XVDRIFT3: overburden-decay offset
const DRIFT_GRAIN_SIZE_MIN = 3.0e-4 # XVSIZEMIN: floor on the fragmented grain size [m]
const DRIFT_GUST_COEF = 1.25        # XCOEF_FF: gust factor on the extrapolated wind
const DRIFT_WIND_EFFECT = 1.0       # XCOEF_EFFECT: weight of the wind term in the drift effect
const DRIFT_QS_REF = 2.0e-5         # XQS_REF: sublimation rate scale in the drift effect
const DRIFT_WIND_HEIGHT = 5.0       # PPHREF_WIND: height the driftability wind is taken at [m]

# Gordon et al. (2006) blowing-snow sublimation coefficients as Crocus writes them
# (`snowcro.F90:4857-4858`), with the paper's Appendix K3 symbol for each.
const DRIFT_SUBLIM_A = 1.8e-3       # a_SUBL: rate prefactor
# NOTE: the paper and the source disagree on which ratio carries which exponent. Lafaysse et
# al. (2026) eq. 67 with Table K3 puts γ_SUBL = 3.6 on the temperature ratio (T₀/T_a) and
# b_SUBL = 4 on the wind ratio (U₅/U_t); `snowcro.F90:4857-4858` writes `**4` on the
# temperature ratio and `**3.6` on the wind ratio — exactly swapped. We follow the source,
# because it is the code that was evaluated, and name the two constants after the ratio they
# are applied to rather than after the paper symbol, so the discrepancy cannot be reintroduced
# by matching names. The choice is not academic: at T_a = 253 K and U₅ = 3·U_t the two
# orderings differ by a factor of ~1.6 in the sublimation rate.
const DRIFT_SUBLIM_EXP_TEMPERATURE = 4.0    # paper calls this b_SUBL, on the wind ratio
const DRIFT_SUBLIM_EXP_WIND = 3.6           # paper calls this γ_SUBL, on the temperature ratio

"""
    drift_wind_speed(cfs, z0) -> U5 [m s-1]

Gust speed at `DRIFT_WIND_HEIGHT` (5 m) above the snow surface, which is the wind the
driftability index is defined against.

`cfs.wind_speed` is given at `cfs.wind_observation_height`, so it is extrapolated to 5 m
through a neutral logarithmic profile over roughness `z0` and multiplied by the
`DRIFT_GUST_COEF` gust factor:

    U₅ = 1.25·U·log(5/z₀)/log(z_u/z₀)

This is `snowcro.F90:4801-4802` exactly, including the three-way floor on the roughness
length, `z₀ = min(z0, z_u/2, 2.5)`. The first two bounds keep the two logs finite and
same-signed when the observation height is low; the 2.5 m bound (half the 5 m reference
height, `PPHREF_MIN` in the source) keeps a very rough surface from inverting the ratio.
Vionnet added the systematic use of a 5 m wind in 2014 precisely so the index does not depend
on the forcing's observation height.

Note the gust factor is *removed* again wherever a rate is formed from `U₅`
(`.../DRIFT_GUST_COEF` in [`apply_blowing_snow!`](@ref) and
[`blowing_snow_sublimation_rate`](@ref)), matching `XCOEF_FF` appearing in both places in the
source: the gust speed sets *whether and how vigorously* snow moves, but the mass flux is
scaled back to the mean wind.

`z0` comes from [`surface_roughness`](@ref), so the drift scheme and the turbulent fluxes see
the same surface.
"""
@inline function drift_wind_speed(cfs::ClimateForcingStep, z0::Float64)
    z_u = cfs.wind_observation_height
    z0_eff = min(z0, 0.5 * z_u, 0.5 * DRIFT_WIND_HEIGHT)
    return DRIFT_GUST_COEF * cfs.wind_speed * log(DRIFT_WIND_HEIGHT / z0_eff) /
           log(z_u / z0_eff)
end

"""
    drift_density_factor(density) -> F(ρ) [-]

Density term of the mobility index, Lafaysse et al. (2026) eq. 60:

    F(ρ) = 1.25 - 0.004·(max(ρ_min, ρ) - ρ_min),   ρ_min = 50 kg m-3

Denser snow is harder to mobilize; the term is flat below `DRIFT_DENSITY_MIN` and goes
negative above ~360 kg m-3, which is what stops a wind slab from being re-eroded indefinitely.

The source writes the slope as `1.25/(1000·XVMOB1)` (`snowcro.F90:4808`), which is
0.0042373, not the 0.004 the paper rounds it to; the source form is used here.
"""
@inline drift_density_factor(density::Float64) =
    DRIFT_GUST_COEF - DRIFT_GUST_COEF *
                      (max(density, DRIFT_DENSITY_MIN) - DRIFT_DENSITY_MIN) / 1000.0 /
                      DRIFT_MOB1

"""
    mobility_index(grain_radius, dendricity, sphericity, density, wet) -> MOB [-]

Mobility index of one layer, Guyomarc'h & Mérindol (1998) as Crocus's `B92` branch implements
it (`snowcro.F90:4812-4820`; Lafaysse et al., 2026 eq. 59, `SNOWMOB = GM98` line combined
with the 0.34/0.66 weighting):

    dendritic:      MOB = 0.34·(0.5 - (0.75δ' + 0.5S')/99) + 0.66·F(ρ)
    non-dendritic:  MOB = 0.34·(0.833 - 0.833·S'/99 - 0.583·gs) + 0.66·F(ρ)

Fresh dendritic crystals and small rounded grains are the most mobile; large faceted grains
and dense snow the least. `F` is [`drift_density_factor`](@ref).

Crocus stores dendricity and sphericity as the Brun et al. (1992) coded integers, negative
dendricity in `PSNOWGRAN1` and sphericity on 0-99 in `PSNOWGRAN2`, where GEMB carries both as
fractions on [0, 1]; the `/99` and sign conventions of the source are undone here so the
formula reads in GEMB's units. `gs` is the grain *size* (diameter) in metres, which is
`2·grain_radius·1e-3`.

`wet` caps the index at `DRIFT_MOB4` (eq. 61): snow that has held liquid water is
bonded and resists drift regardless of its grain shape.

!!! note "Deviation from Crocus: no wet-snow history"
    Crocus applies that cap on `PSNOWHIST ≥ 2`, i.e. whether the layer has held liquid water
    at **any time** in its history. GEMB's column state has no such field, and adding one
    would extend `RESTART_LAYERS` and break the restart contract, so the caller passes the
    *instantaneous* `water[i] > WATER_TOLERANCE` instead. The consequence is one-directional:
    a refrozen layer that is currently dry is not capped here, so it is more driftable than
    Crocus would make it, and a melt-affected surface can be wind-packed a second time. If
    this proves to matter, a `wet_history` layer field is a separate, restart-breaking change.
"""
@inline function mobility_index(grain_radius::Float64, dendricity::Float64,
    sphericity::Float64, density::Float64, wet::Bool)

    F = drift_density_factor(density)
    if dendricity > 0.0 + GDN_TOLERANCE
        # Dendritic. Crocus's `PSNOWGRAN1 = -99·δ` and `PSNOWGRAN2 = 99·S`, so
        # `(0.75·PSNOWGRAN1 + 0.5·PSNOWGRAN2)/99` is `0.5·S - 0.75·δ`.
        mob = 0.34 * (0.5 + 0.75 * dendricity - 0.5 * sphericity) + 0.66 * F
    else
        # Non-dendritic. `PSNOWGRAN1 = 99·S` here (the source reuses the slot for sphericity
        # once dendricity is spent) and `PSNOWGRAN2` is the grain size in metres.
        grain_size = 2.0 * grain_radius * 1e-3      # mm radius -> m diameter
        mob = 0.34 * (DRIFT_MOB2 - DRIFT_MOB2 * sphericity -
                      DRIFT_MOB3 * grain_size * 1000.0) + 0.66 * F
    end
    return wet ? min(mob, DRIFT_MOB4) : mob
end

"""
    driftability_index(mobility, wind_5m) -> D [-]

Driftability of one layer from its mobility and the 5 m gust speed, Lafaysse et al. (2026)
eq. 62:

    D = max(MOB - (2.868·exp(-0.085·U₅) - 1), 0)

The bracket is the mobility a layer must have for the given wind to move it: it falls from
1.868 at zero wind to 0 at `U₅ = log(2.868)/0.085 = 12.4 m s-1`, above which every layer with
non-negative mobility drifts. Inverting it gives the threshold wind speed
[`drift_threshold_wind_speed`](@ref) that the sublimation rate is written against.

The `max(·, 0)` is in the paper but not in the source, which instead exits its layer loop at
the first non-positive value — see [`apply_blowing_snow!`](@ref), which follows the source.
Both are kept: the clamp here makes the function well-defined for a caller asking about a
single layer.
"""
@inline driftability_index(mobility::Float64, wind_5m::Float64) =
    max(mobility - (DRIFT_D1 * exp(-DRIFT_D2 * wind_5m) - 1.0), 0.0)

"""
    drift_threshold_wind_speed(mobility) -> Ut [m s-1]

The 5 m wind speed at which a layer of the given mobility just begins to drift, Lafaysse et
al. (2026) eq. 68 — [`driftability_index`](@ref) solved for `D = 0`:

    U_t = -log((MOB + 1)/2.868)/0.085

Read only by [`blowing_snow_sublimation_rate`](@ref), where it sets both the scale of the
sublimation rate and the wind ratio it is driven by. This is the "modification of a threshold
wind speed U_t to account for the microstructure-related mobility index" the paper describes:
Gordon et al. (2006) used a fixed threshold, Crocus makes it a function of the snow's own
state.

Undefined for `MOB ≥ 1.868` (the log's argument exceeds 1, giving a negative threshold) and
divergent as `MOB → -1`; neither is reachable from the mobility formula at a physical density,
and the caller gates on a positive driftability, which implies a threshold below the current
wind.
"""
@inline drift_threshold_wind_speed(mobility::Float64) =
    -log((mobility + 1.0) / DRIFT_D1) / DRIFT_D2

"""
    blowing_snow_sublimation_rate(mobility, wind_5m, cfs) -> Q_S [kg m-2 s-1]

Sublimation rate of suspended blowing snow, Gordon et al. (2006) as Crocus implements it
(`snowcro.F90:4851-4859`; Lafaysse et al., 2026 eq. 67):

    Q_S = a·(T₀/T_a)^4·U_t·ρ_a·q_sat(T_a,P_s)·(1 - q_a/q_sat)·(U₅/U_t)^3.6

Suspended particles have a far larger surface-to-mass ratio than the snow surface, so they
sublimate much faster than the surface itself; the rate is driven by the saturation deficit
`(1 - RH_i)` and rises steeply with wind through the ratio to the threshold speed `U_t`
([`drift_threshold_wind_speed`](@ref)).

Humidity is with respect to **ice**, via [`saturation_vapor_pressure_ice`](@ref) and
[`specific_humidity`](@ref) — the deficit that matters is the one against the sublimating
solid, and using the liquid saturation instead would understate it below freezing.

Can be negative (the caller clamps at zero) when the air is supersaturated with respect to
ice, which is common in polar winter; the scheme does not represent the deposition that
physically corresponds to that.

See `DRIFT_SUBLIM_EXP_TEMPERATURE` for the paper-versus-source disagreement over which of the
two exponents goes on which ratio, and why the source's ordering is used.
"""
@inline function blowing_snow_sublimation_rate(mobility::Float64, wind_5m::Float64,
    cfs::ClimateForcingStep)

    U_t = drift_threshold_wind_speed(mobility)
    U_t <= 0.0 && return 0.0

    density_air = air_density(cfs.pressure_air, cfs.temperature_air)
    e_sat = saturation_vapor_pressure_ice(cfs.temperature_air)
    q_sat = specific_humidity(e_sat, cfs.pressure_air)
    q_air = specific_humidity(cfs.vapor_pressure, cfs.pressure_air)
    relative_humidity_ice = q_sat > 0.0 ? q_air / q_sat : 1.0

    return DRIFT_SUBLIM_A * (CtoK / cfs.temperature_air)^DRIFT_SUBLIM_EXP_TEMPERATURE *
           U_t * density_air * q_sat * (1.0 - relative_humidity_ice) *
           (wind_5m / U_t)^DRIFT_SUBLIM_EXP_WIND
end

"""
    apply_blowing_snow!(cols, cfs, mp) -> (mass_sublimated, E_sublimated)

Compact and fragment near-surface snow under wind, and optionally sublimate the suspended
fraction. The Crocus `SNOWDRIFT` scheme (`snowcro.F90:4685-4950`; Lafaysse et al., 2026
eqs. 59-66), ported for GEMB's state variables.

Per layer, from the surface down:

1. Mobility ([`mobility_index`](@ref)) and driftability ([`driftability_index`](@ref)) from
   the layer's own microstructure, density, and the 5 m gust speed.
2. **Stop at the first layer that will not drift.** The overburden a layer sits under can only
   grow with depth, so once driftability reaches zero it stays zero. This is the source's
   `EXIT` at `snowcro.F90:4844`; the paper's eq. 63 does not show it, and a literal port of
   the equation would walk the whole column computing zeros.
3. The effective driftability `DEFF` (eq. 63): the layer's own driftability decayed by
   `exp(-100·Σ 0.5·dz_j·0.1·(3.25 - D_j))` over the overburden, counted a half-layer at a time
   so that each layer is decayed by the half of itself above its own centre. Deeper snow is
   reworked less, and snow under *mobile* snow (large `D_j`) is decayed less than snow under a
   crust — the overburden's own driftability enters the exponent.
4. Compaction toward `DRIFT_DENSITY_MAX` (eq. 64) at the rate `DEFF·dt/τ_DRIFT`, exactly
   conservative in mass: `dz` is rescaled by the density ratio. Snow already at or above
   350 kg m-3 is left alone, which is the source's `IF (PSNOWRHO < XVROMAX)` gate.
5. Grain fragmentation (eqs. 65-66): dendricity falls, sphericity rises toward 1, and grain
   size falls toward `DRIFT_GRAIN_SIZE_MIN` — saltating grains break up and round off. This is
   the mechanism by which drifted snow is optically and mechanically distinct from the snow it
   came from, not merely denser.

When `mp.blowing_snow_sublimation` is set, the surface layer also loses mass at
[`blowing_snow_sublimation_rate`](@ref), capped at half the layer's thickness per step
(the source's `MAX(0.5*PSNOWDZ, ...)`, which keeps a long step from emptying the layer), and
that rate feeds back into the compaction rate through the `ZQS_EFFECT` term: a vigorously
sublimating surface is also being vigorously reworked.

Returns the sublimated `(mass, energy)`, both zero unless `mp.blowing_snow_sublimation` — the
compaction and fragmentation path moves no mass and carries no energy across the surface, so
it contributes nothing to either budget. Both are **negative**, signed as a flux into the
column the way [`apply_horizontal_strain!`](@ref)'s are, so the caller adds rather than
subtracts them.

!!! note "The latent heat is not charged to the column"
    The energy returned is the departing snow's own enthalpy only — no `LS`. Blowing-snow
    sublimation happens to snow *already airborne*, in the saltation and suspension layers, so
    the latent heat is drawn from the air the particles are suspended in, not from the pack:
    that cooling and moistening of the near-surface air is precisely the effect Gordon et al.
    (2006) added the parameterization to capture. Crocus is consistent with this — `SNOWDRIFT`
    writes `PSNOWDZ` and reports the flux in `PSNDRIFT`, but applies no temperature or enthalpy
    tendency to the layer it took the mass from.

    The distinction is against *surface* sublimation, which GEMB does charge `LS` to the column
    (through `heat_flux_latent` in [`calculate_temperature`](@ref)) because that phase change
    happens at the ice-air interface. Charging `LS` here as well would double-count it: with a
    2 h step and a windy forcing it fails the verbose energy check by order 1e5 J m-2 per step.

    The consequence is that GEMB, like Crocus, does not close the *atmospheric* energy budget
    for this term. A coupled model would need the flux back; `flux.mass_blowing_snow` and the
    `blowing_snow` output layer are what it would read.

Returns `(0.0, 0.0)` and touches nothing when `mp.blowing_snow_method === :none`, which is
what keeps the default bit-identical to a run without the term.
"""
function apply_blowing_snow!(cols::NamedTuple, cfs::ClimateForcingStep, mp::ModelParameters)
    mp.blowing_snow_method === :none && return 0.0, 0.0

    dz = cols.dz
    density = cols.density
    water = cols.water
    grain_radius = cols.grain_radius
    dendricity = cols.grain_dendricity
    sphericity = cols.grain_sphericity

    # The gust speed sees the same surface the turbulent fluxes do.
    z0 = surface_roughness(mp, density[1], water[1])
    wind_5m = drift_wind_speed(cfs, z0)

    mass_sublimated = 0.0
    E_sublimated = 0.0

    # Cumulative overburden decay coefficient, `ZPROFEQU`. Advanced by half the current
    # layer before its own `DEFF` is formed and by the other half afterwards, so each layer
    # is decayed by the overburden above its centre.
    decay = 0.0

    @inbounds for i in eachindex(dz)
        wet = water[i] > WATER_TOLERANCE
        mobility = mobility_index(grain_radius[i], dendricity[i], sphericity[i],
            density[i], wet)
        drift = mobility - (DRIFT_D1 * exp(-DRIFT_D2 * wind_5m) - 1.0)

        # Nothing below a non-drifting layer can drift either: the overburden term only grows
        # with depth. `snowcro.F90:4844`.
        drift <= 0.0 && break

        decay += 0.5 * dz[i] * 0.1 * (DRIFT_D3 - drift)
        drift_effective = max(0.0, drift * exp(-decay * 100.0))

        # Blowing-snow sublimation, surface layer only: it is the suspended snow above the
        # surface that sublimates, and the scheme charges the loss to the layer it came from.
        sublimation_rate = 0.0
        if mp.blowing_snow_sublimation && i == 1
            sublimation_rate = max(0.0,
                blowing_snow_sublimation_rate(mobility, wind_5m, cfs))
            if sublimation_rate > 0.0
                # `/DRIFT_GUST_COEF` removes the gust factor from the mass flux; see
                # `drift_wind_speed`. Capped at half the layer.
                dz_loss = min(0.5 * dz[i],
                    sublimation_rate * cfs.dt / DRIFT_GUST_COEF / density[i])
                mass_loss = dz_loss * density[i]
                dz[i] -= dz_loss
                # Signed as a flux into the column, hence negative for a loss. The enthalpy is
                # the departing snow's own, with **no** latent heat: see the docstring.
                mass_sublimated -= mass_loss
                E_sublimated -= mass_loss * specific_enthalpy(mp, cols.temperature[i])
            end
        end

        # Drift effect: the wind term always, plus a sublimation term that saturates at 3
        # (`snowcro.F90:4871-4873`). The gust factor is removed here too.
        qs_effect = min(3.0, sublimation_rate / DRIFT_QS_REF) * drift_effective
        wind_effect = DRIFT_WIND_EFFECT * drift_effective
        drift_rate = (qs_effect + wind_effect) * cfs.dt / DRIFT_GUST_COEF / DRIFT_TAU

        # Wind-slab compaction, at constant mass. Skipped once the layer has reached the
        # wind-slab density: eq. 64's `max(ρ_MAX - ρ, 0)` would be zero anyway, but gating
        # also leaves `dz` exactly untouched rather than multiplied by 1.
        if density[i] < DRIFT_DENSITY_MAX
            density_old = density[i]
            density[i] = min(DRIFT_DENSITY_MAX,
                density_old + drift_rate * (DRIFT_DENSITY_MAX - density_old))
            dz[i] *= density_old / density[i]
        end

        # Fragmentation. Dendritic snow loses dendricity and gains sphericity; non-dendritic
        # snow gains sphericity and loses size. Both branches keep their variables in range
        # by construction — the increments are proportional to the distance remaining.
        if dendricity[i] > 0.0 + GDN_TOLERANCE
            d_dendricity = min(drift_rate * dendricity[i] * 0.5, 0.99 * dendricity[i])
            dendricity[i] -= d_dendricity
            sphericity[i] = min(1.0, sphericity[i] + drift_rate * (1.0 - sphericity[i]))
            # Grain size is diagnostic while dendricity is nonzero, so re-derive it rather
            # than evolving it independently.
            grain_radius[i] = dendritic_grain_radius(dendricity[i], sphericity[i])
        else
            sphericity[i] = min(1.0, sphericity[i] + drift_rate * (1.0 - sphericity[i]))
            # `ZDGR2 = ZDRIFT_EFFECT * 5/10000` in metres of diameter; GEMB carries radius in
            # mm, so 5e-4 m of diameter is 0.25 mm of radius.
            grain_radius[i] = max(0.5 * DRIFT_GRAIN_SIZE_MIN * 1e3,
                grain_radius[i] - drift_rate * 0.25)
        end

        decay += 0.5 * dz[i] * 0.1 * (DRIFT_D3 - drift)
    end

    return mass_sublimated, E_sublimated
end

"""
    apply_prescribed_drift!(cols, cfs, mp) -> (mass_drift, E_drift)

Apply a prescribed net drift divergence to the surface cell, and return the mass and enthalpy
it moved.

This is the IMAU-FDM treatment (`source/firn_physics.f90:81-88`, `Update_Surface`), which
takes drift as forcing computed upstream by RACMO (Lenaerts et al., 2012) rather than
computing it: the flux is a mass sink or source at the surface, with a **sign-dependent
density**, because erosion removes the snow that is actually there while deposition brings
snow of whatever density fresh drifted snow has:

- **Erosion** (positive flux) removes mass at the surface cell's own `density[1]`, thinning it
  and taking its pore water in proportion — the cell's density, temperature, and grain
  properties are unchanged, exactly as sublimation of a fraction of a cell would leave them.
  Age-neutral for the same reason (see [`dilute_age`](@ref)).
- **Deposition** (negative flux) adds mass at [`fresh_snow_density`](@ref) — the same density
  fresh snowfall arrives at under `mp.new_snow_method` — mixing it into the surface cell.
  Age-zero mass, so it dilutes `age[1]`, and it arrives at the air temperature, so it mixes
  the cell's temperature the way snowfall does.

Deposition merges into the surface cell rather than creating one, unlike snowfall: a drift
flux is a slow net divergence, not a discrete event, so at any realistic rate the mass added
per step is far below `mp.column_dzmin` and a fresh cell would be merged straight back by
[`manage_layer_thickness`](@ref). Erosion never removes more than half the surface cell in one
step, the same cap [`apply_blowing_snow!`](@ref) applies to sublimation and for the same
reason.

The rate is the optional `snow_drift` forcing layer when the forcing carries one, and
`mp.drift_rate` (a constant [kg m-2 yr-1], positive for erosion) otherwise. Returns the mass
signed as a loss from the column — negative for erosion, positive for deposition — so the
caller adds it to the mass budget the way it adds [`trim_bottom!`](@ref)'s basal flux, and the
matching enthalpy.

Unlike [`apply_blowing_snow!`](@ref)'s sublimation this is a *transport* term: the mass leaves
or arrives as snow, not vapour, so it carries only its own enthalpy and no latent heat.

Returns `(0.0, 0.0)` and touches nothing when the rate is zero, which is what keeps the
`drift_rate = 0.0` default bit-identical to a run without the term.
"""
function apply_prescribed_drift!(cols::NamedTuple, cfs::ClimateForcingStep,
    mp::ModelParameters)

    # kg m-2 yr-1 -> kg m-2 per step. `SECONDS_PER_YEAR` (Julian) rather than
    # `DAYS_PER_YEAR_DENSIFICATION`: this is an observed mass flux, not a densification
    # coefficient fitted against a 365-day year — the same argument as
    # `apply_horizontal_strain!`.
    flux = cfs.snow_drift * cfs.dt / SECONDS_PER_YEAR
    flux == 0.0 && return 0.0, 0.0

    dz = cols.dz
    density = cols.density
    water = cols.water
    temperature = cols.temperature

    @inbounds if flux > 0.0
        # Erosion. Remove a fraction of the surface cell at its own density, water included.
        dz_loss = min(0.5 * dz[1], flux / density[1])
        water_loss = dz[1] > 0.0 ? water[1] * dz_loss / dz[1] : 0.0
        mass_loss = dz_loss * density[1] + water_loss

        dz[1] -= dz_loss
        water[1] -= water_loss

        # Age-neutral: a fraction of the cell leaves, which does not change its mean age.
        return -mass_loss, -(dz_loss * density[1] * specific_enthalpy(mp, temperature[1]) +
                             water_loss * specific_enthalpy_water(mp))
    else
        # Deposition, at the fresh-snow density and the air temperature. Grain properties are
        # left at the surface cell's own: drifted snow is already fragmented and rounded, and
        # `apply_blowing_snow!` is the path that models that — this one has no microstructure
        # information to bring.
        mass_added = -flux
        density_new = clamp(fresh_snow_density(mp, cfs.temperature_air_mean,
                cfs.precipitation_mean, cfs.wind_speed_mean, cfs.temperature_air,
                cfs.wind_speed),
            1.0, mp.density_ice)

        M_surface = dz[1] * density[1]
        S_surface = M_surface + water[1]

        dz[1] += mass_added / density_new
        temperature[1] = mix_temperature(mp, cfs.temperature_air, mass_added,
            temperature[1], M_surface)
        density[1] = (M_surface + mass_added) / dz[1]
        cols.age[1] = dilute_age(cols.age[1], S_surface, mass_added)

        return mass_added, mass_added * specific_enthalpy(mp, cfs.temperature_air)
    end
end
