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

Path 1 has a second scheme, `mp.blowing_snow_method = :Keenan2021`
([`apply_redeposition!`](@ref)): the SNOWPACK erosion-redeposition scheme of Keenan et al.
(2021), in which a saltation mass flux is eroded from the top of the column and redeposited at
a wind-dependent density that reaches 394 kg m-3 at 10 m s-1. It replaces the `:Crocus`
compaction and is also mass-conserving.

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
  *Geosci. Model Dev.* 19, 6273-6334. https://doi.org/10.5194/gmd-19-6273-2026
- Noël, B., van de Berg, W. J., van Wessem, J. M., et al. (2018). Modelling the climate and
  surface mass balance of polar ice sheets using RACMO2 - Part 1: Greenland (1958-2016).
  *The Cryosphere* 12, 811-831.
- Purdie, H., Mackintosh, A., Lawson, W., Anderson, B., Morgenstern, U., Chinn, T., and
  Mayewski, P. (2011). Interannual variability in net accumulation on Tasman Glacier and its
  relationship with climate. *Global Planet. Change* 77, 142-152.
- Keenan, E., Wever, N., Dattler, M., Lenaerts, J. T. M., Medley, B., Kuipers Munneke, P.,
  and Reijmer, C. (2021). Physics-based SNOWPACK model improves representation of
  near-surface Antarctic snow and firn density. *The Cryosphere* 15, 1065-1085.
  https://doi.org/10.5194/tc-15-1065-2021
- Wever, N., Keenan, E., Amory, C., Lehning, M., Sigmund, A., Huwald, H., and Lenaerts,
  J. T. M. (2023). Observations and simulations of new snow density in the drifting
  snow-dominated environment of Antarctica. *J. Glaciol.* 69, 823-840.
  https://doi.org/10.1017/jog.2022.102
- Brun, E., David, P., Sudul, M., and Brunot, G. (1992). A numerical model to simulate
  snow-cover stratigraphy for operational avalanche forecasting. *J. Glaciol.* 38, 13-22.
  (The coded dendricity/sphericity integers Crocus stores.)
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
    wind_speed_10m(wind_speed, wind_observation_height, z0) -> U₁₀ [m s-1]

Wind speed at 10 m above the surface, from a `wind_speed` [m s-1] observed at
`wind_observation_height` [m], through a neutral logarithmic profile over roughness length
`z0` [m]:

    U₁₀ = U·log(10/z₀)/log(z_u/z₀)

The identity when `wind_observation_height = 10`. It is the profile [`drift_wind_speed`](@ref)
uses, without the gust factor. Read by `new_snow_method = :vanKampenhout2017`, whose wind
term is defined on the 10 m wind (van Kampenhout et al., 2017, eq. 4, p. 2587), and by
`blowing_snow_method = :Keenan2021`, whose friction velocity and redeposited density are
defined on it (Keenan et al., 2021, pp. 1067-1068). Callers pass `z0` from [`surface_roughness`](@ref), so
the 10 m wind is formed over the same surface as the turbulent fluxes.

Throws an `ArgumentError` for a negative or non-finite wind speed, a non-positive `z0`, or an
observation height at or below `z0`.
"""
function wind_speed_10m(wind_speed::Real, wind_observation_height::Real, z0::Real)
    (isfinite(wind_speed) && wind_speed >= 0) ||
        throw(ArgumentError("wind speed must be finite and non-negative [m s-1], got $wind_speed"))
    z0 > 0 || throw(ArgumentError("roughness length must be positive [m], got $z0"))
    wind_observation_height > z0 ||
        throw(ArgumentError("wind observation height ($wind_observation_height m) must exceed the roughness length ($z0 m)"))
    return wind_speed * (log(WIND_HEIGHT_10M / z0) /
                         log(wind_observation_height / z0))
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
what keeps the default bit-identical to a run without the term. Under
`mp.blowing_snow_method === :Keenan2021` it runs [`apply_redeposition!`](@ref) instead of the
compaction above; the two are not combined, since both represent the densification of
wind-transported snow.
"""
function apply_blowing_snow!(cols::NamedTuple, cfs::ClimateForcingStep, mp::ModelParameters)
    mp.blowing_snow_method === :none && return 0.0, 0.0
    mp.blowing_snow_method === :Keenan2021 && return apply_redeposition!(cols, cfs, mp)
    mp.blowing_snow_method === :Crocus ||
        error("unrecognized blowing_snow_method: $(mp.blowing_snow_method)")

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
            # than evolving it independently, unless it is state of its own.
            independent_grain_radius(mp) ||
                (grain_radius[i] = dendritic_grain_radius(dendricity[i], sphericity[i]))
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
                cfs.wind_speed; wind_observation_height=cfs.wind_observation_height,
                z0=surface_roughness(mp, density[1], water[1])),
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

# Keenan et al. (2021) erosion-redeposition constants.
const REDEPOSIT_Z0 = 2.0e-3          # roughness length for u* [m], Keenan et al. p. 1067
const REDEPOSIT_FLUX_A = 0.0014      # eq. 3 prefactor [-]
const REDEPOSIT_FLUX_B = 7.6         # eq. 3 threshold weight [-]
const REDEPOSIT_FLUX_C = 205.0       # eq. 3 offset, as printed (see `saltation_mass_flux`)
const REDEPOSIT_DENSITY_A = 361.0    # eq. 4 slope [kg m-3 per decade of U₁₀]
const REDEPOSIT_DENSITY_B = 33.0     # eq. 4 intercept, and its value at U₁₀ <= 1 [kg m-3]
const REDEPOSIT_WIND_MIN = 1.0       # eq. 4 branch point [m s-1]
# Microstructure of new snow above 5 m s-1 in the SNOWPACK "Redeposit" setup, Wever et al.
# (2023) Table 1, p. 827.
const REDEPOSIT_DENDRICITY = 0.325
const REDEPOSIT_SPHERICITY = 0.875

"""
    redeposition_density(wind_speed_10m) -> ρ_rd [kg m-3]

Density of snow redeposited after drifting, Keenan et al. (2021) eq. 4, p. 1068:

    ρ_rd = 361·log10(U₁₀) + 33     U₁₀ > 1 m s-1
           33                      U₁₀ ≤ 1 m s-1

with `U₁₀` the hourly 10 m wind speed. The fit is Groot Zwaaftink et al. (2013), from field
measurements of snow deposited during drifting-snow events; Keenan et al. drive it with hourly
rather than 100 h mean wind. It is continuous at 1 m s-1 and gives 285.3 kg m-3 at 5 m s-1,
394.0 at 10 and 502.7 at 20.

Throws an `ArgumentError` unless `wind_speed_10m` is finite and non-negative.
"""
function redeposition_density(wind_speed_10m::Real)
    (isfinite(wind_speed_10m) && wind_speed_10m >= 0) ||
        throw(ArgumentError("wind speed must be finite and non-negative [m s-1], got $wind_speed_10m"))
    wind_speed_10m > REDEPOSIT_WIND_MIN || return REDEPOSIT_DENSITY_B
    return REDEPOSIT_DENSITY_A * log10(wind_speed_10m) + REDEPOSIT_DENSITY_B
end

"""
    saltation_friction_velocity(wind_speed_10m) -> u* [m s-1]

Friction velocity of the Keenan et al. (2021) saltation flux from the 10 m wind, through a
neutral logarithmic profile over `REDEPOSIT_Z0 = 2 mm` (p. 1067):

    u* = κ·U₁₀/log(10/z₀)

Keenan et al. add stability corrections (Michlmayr et al., 2008); they are omitted here, so
this is the neutral value. `z₀ = 2 mm` is the paper's choice, about ten times the roughness of
the snow surface GEMB's turbulent fluxes see, and is kept because the fetch length `L` of
eq. 3 was chosen with it.
"""
@inline saltation_friction_velocity(wind_speed_10m::Real) =
    VON_KARMAN * wind_speed_10m / log(WIND_HEIGHT_10M / REDEPOSIT_Z0)

"""
    saltation_mass_flux(density_air, u_star, u_star_threshold, fetch_length) -> Φ [kg m-2 s-1]

Saltation mass flux eroded from the snow surface, Keenan et al. (2021) eq. 3, p. 1067:

    Φ = 0.0014·ρ_a·u*·(u* − u*_t)·(u* + 7.6·u*_t + 205)/L

The saltation transport rate `Q` [kg m-1 s-1] of Lehning and Fierz (2008) divided by the
fetch length `L` [m]. Zero when `u* ≤ u*_t`. The additive 205 is as printed in Keenan et al.
and as SNOWPACK evaluates it; Wever et al. (2023, pp. 827-828) note that the coefficients of
this expression carry inconsistent units, and that the simulations agree with observations
regardless. It dominates the last factor, so in practice `Φ ≈ 0.287·ρ_a·u*·(u* − u*_t)/L`.

Throws an `ArgumentError` for a non-positive fetch length or a negative friction velocity.
"""
function saltation_mass_flux(density_air::Real, u_star::Real, u_star_threshold::Real,
    fetch_length::Real)
    fetch_length > 0 ||
        throw(ArgumentError("fetch length must be positive [m], got $fetch_length"))
    (u_star >= 0 && u_star_threshold >= 0) ||
        throw(ArgumentError("friction velocities must be non-negative, got $u_star and $u_star_threshold"))
    u_star > u_star_threshold || return 0.0
    return REDEPOSIT_FLUX_A * density_air * u_star * (u_star - u_star_threshold) *
           (u_star + REDEPOSIT_FLUX_B * u_star_threshold + REDEPOSIT_FLUX_C) / fetch_length
end

"""
    redeposition_threshold_wind_speed(mobility) -> U₅,t [m s-1]

The 5 m gust speed at which a layer of the given Crocus mobility index starts to drift,
[`drift_threshold_wind_speed`](@ref), extended to the whole range of the index: `Inf` for
`mobility ≤ -1`, where no wind moves the layer (dense or wet snow).

Throws for `mobility ≥ 1.868`, where the threshold would be negative; [`mobility_index`](@ref)
does not exceed 1.25.
"""
@inline function redeposition_threshold_wind_speed(mobility::Float64)
    mobility < DRIFT_D1 - 1.0 ||
        error("mobility index $mobility is outside the range of the drift threshold")
    mobility <= -1.0 && return Inf
    return drift_threshold_wind_speed(mobility)
end

"""
    apply_redeposition!(cols, cfs, mp) -> (0.0, 0.0)

Wind packing of near-surface snow by erosion and redeposition, the SNOWPACK scheme of Keenan
et al. (2021) eqs. 3-4, pp. 1067-1068, with the erosion threshold taken from GEMB's Crocus
mobility index. Selected by `mp.blowing_snow_method = :Keenan2021`.

Per step:

1. The 10 m wind `U₁₀` ([`wind_speed_10m`](@ref)) and the 5 m gust speed `U₅`
   ([`drift_wind_speed`](@ref)), both over the roughness of the surface cell.
2. Erosion starts only if `U₅` exceeds the threshold `U₅,t` of the surface cell
   ([`redeposition_threshold_wind_speed`](@ref) of its [`mobility_index`](@ref)) and the
   surface cell holds no liquid water.
3. The friction velocity `u* = κ·U₁₀/log(10/0.002)` ([`saltation_friction_velocity`](@ref))
   and its threshold `u*_t = u*·U₅,t/U₅`, the friction velocity at which `U₅` would equal the
   threshold.
4. The eroded mass `M = Φ·dt`, with `Φ` the saltation flux of eq. 3
   ([`saltation_mass_flux`](@ref)) at the air density of the forcing and fetch length
   `mp.drift_fetch_length`.
5. `M` is taken from the top of the column down. Erosion stops at the first cell that is wet
   or whose own threshold is at or above `U₅`, or when `M` is used up. Of a cell of mass `m`,
   the reworked mass `m_r = min(M, m)` is redeposited in the same cell at
   `ρ_rd = 361·log10(U₁₀) + 33` ([`redeposition_density`](@ref), eq. 4), so the cell becomes
   `dz = (m − m_r)/ρ + m_r/ρ_rd`, `ρ = m/dz`. A cell already at or above `ρ_rd` keeps its
   density: redeposition never lowers it.
6. In a dendritic cell the reworked fraction `f = m_r/m` takes the microstructure of new snow
   above 5 m s-1 in the SNOWPACK "Redeposit" setup, dendricity 0.325 and sphericity 0.875
   (Wever et al., 2023, Table 1, p. 827), mixed by mass and capped so that neither moves away
   from that state: `d ← (1−f)·d + f·min(d, 0.325)`, `s ← (1−f)·s + f·max(s, 0.875)`, with the
   radius re-derived by [`dendritic_grain_radius`](@ref) (kept where
   [`independent_grain_radius`](@ref) holds). Non-dendritic cells keep their grains.

**Mass and energy.** The eroded mass is redeposited in the column (Keenan et al. p. 1068: the
mass flux "is redeposited at the snow surface"; Wever et al. p. 827: deposition "is equal to
the calculated eroded mass"), so nothing crosses the column boundary. Each cell keeps its own
mass, temperature, water and age, so column mass and enthalpy are unchanged to rounding, and
the function returns `(0.0, 0.0)` for the blowing-snow budget terms. The grid is not touched
beyond the `dz` of the reworked cells; the grid controllers that run later in the step restore
the thickness bands.

Departures from SNOWPACK as Keenan et al. describe it:

- **Threshold.** SNOWPACK computes `u*_t` from Lehning and Fierz (2008), Keenan et al. eq. 2,
  which needs a bond radius and a coordination number. GEMB carries neither, so the threshold
  comes from the Crocus mobility index, which is GEMB's own measure of how readily a layer
  drifts given its dendricity, sphericity, grain size and density (Vionnet et al., 2012). For
  fresh dendritic snow (d = 1, s = 0.5, ρ = 180 kg m-3) it gives `U₅,t = 5.47 m s-1`, about
  4.7 m s-1 at 10 m.
- **Placement.** SNOWPACK removes the eroded elements and deposits a new layer at the surface.
  Here the reworked mass stays in the cell it came from. The mass and its redeposited density
  are the same; the vertical order within the eroded depth and the cell count are kept, so the
  scheme creates and removes no cells.
- **Stability and air density.** `u*` is the neutral value, and `ρ_a` is the air density of
  the forcing step where Keenan et al. quote 1.1 kg m-3 (p. 1067, for eq. 2).
- **Wet snow** does not erode.

Throws an `ArgumentError` for a negative or non-finite wind speed.
"""
function apply_redeposition!(cols::NamedTuple, cfs::ClimateForcingStep, mp::ModelParameters)
    dz = cols.dz
    density = cols.density
    water = cols.water
    grain_radius = cols.grain_radius
    dendricity = cols.grain_dendricity
    sphericity = cols.grain_sphericity

    z0 = surface_roughness(mp, density[first(eachindex(density))],
        water[first(eachindex(water))])
    wind_10m = wind_speed_10m(cfs.wind_speed, cfs.wind_observation_height, z0)
    wind_5m = drift_wind_speed(cfs, z0)

    _cell_threshold(i) = redeposition_threshold_wind_speed(mobility_index(grain_radius[i],
        dendricity[i], sphericity[i], density[i], false))

    i1 = first(eachindex(dz))
    (water[i1] > WATER_TOLERANCE || wind_5m <= _cell_threshold(i1)) && return 0.0, 0.0

    u_star = saltation_friction_velocity(wind_10m)
    u_star_threshold = u_star * _cell_threshold(i1) / wind_5m
    density_air = air_density(cfs.pressure_air, cfs.temperature_air)
    eroded = saltation_mass_flux(density_air, u_star, u_star_threshold,
        mp.drift_fetch_length) * cfs.dt
    density_redeposit = redeposition_density(wind_10m)

    for i in eachindex(dz)
        eroded > 0.0 || break
        (water[i] > WATER_TOLERANCE || wind_5m <= _cell_threshold(i)) && break
        mass = density[i] * dz[i]
        reworked = min(eroded, mass)
        if density_redeposit > density[i]
            dz[i] = (mass - reworked) / density[i] + reworked / density_redeposit
            density[i] = mass / dz[i]
        end
        if dendricity[i] > GDN_TOLERANCE
            f = reworked / mass
            dendricity[i] = (1 - f) * dendricity[i] + f * min(dendricity[i], REDEPOSIT_DENDRICITY)
            sphericity[i] = (1 - f) * sphericity[i] + f * max(sphericity[i], REDEPOSIT_SPHERICITY)
            independent_grain_radius(mp) ||
                (grain_radius[i] = dendritic_grain_radius(dendricity[i], sphericity[i]))
        end
        eroded -= reworked
    end
    return 0.0, 0.0
end
