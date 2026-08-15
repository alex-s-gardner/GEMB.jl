"""
Lateral meltwater drainage.

GEMB's column is one-dimensional and per unit area, so water leaving it sideways cannot be
routed to a neighbour — it is simply removed and reported as runoff. What this file adds is
*when*: a rate law, so that water in excess of irreducible saturation leaves over a finite
timescale set by the surface slope rather than instantaneously.

The two laws here are the ones RetMIP (Vandecrux et al., 2020) found at the good end of the
bucket-scheme spread. Their Table 5 firn-temperature mean errors at KAN_U, an ice-slab site
where the runoff treatment dominates: DMIHH −1.6 °C (Zuo and Oerlemans timescale) and GEUS
+0.6 °C (Darcy flux), against a multi-model spread reaching +4.7 °C; DTU, which runs water off
immediately, produced runoff unrealistic enough to be excluded from the paper's multi-model
mean. Both laws being available makes the pair comparable within one model, which RetMIP could
only do across models that differ in everything else.

Neither runs at the `runoff_method = :instantaneous` default, which leaves output identical to
MATLAB. Both require [`calculate_melt`](@ref) to have let water exceed irreducible saturation
in the first place — a rate law has nothing to act on otherwise, which is why the two changes
are one feature.
"""

# Zuo and Oerlemans (1996) eq. 22 timescale coefficients, in seconds and dimensionless.
#
# `c1`/`c2` are quoted here as Lefebre et al. (2003) and Langen et al. (2017) give them
# (0.33 d and 25 d), which is the DMIHH lineage RetMIP evaluated. Note that the Community Firn
# Model's `runoffZuoOerlemans` uses `c1 = 1.5 d` instead, citing the same source; the two
# published lineages disagree on this coefficient and we follow the one whose model was
# evaluated. The difference matters only at steep slopes, where `c2·exp(-c3·S)` is small and
# `c1` sets the floor: at S = 0.01 the timescale is 21 d either way, at S = 0.1 it is 0.33 d
# against 1.5 d.
const ZUO_OERLEMANS_C1 = 0.33 * 86400.0   # [s]
const ZUO_OERLEMANS_C2 = 25.0 * 86400.0   # [s]
const ZUO_OERLEMANS_C3 = 140.0            # [-]

# Dynamic viscosity of water at the melting point [kg m-1 s-1], for the Darcy flux.
const VISCOSITY_WATER = 1.792e-3

"""
    runoff_timescale(mp) -> tau [s]

Characteristic lateral runoff timescale of Zuo and Oerlemans (1996) eq. 22,
`τ = c1 + c2·exp(-c3·S)`, with `S = mp.surface_slope` [m m-1].

Steeper terrain drains faster, asymptotically to `c1`; the timescale saturates at `c1 + c2`
(≈ 25 d) on flat terrain rather than diverging, so a zero slope still drains. Independent of
grain size, density, and saturation — the whole slope dependence of the scheme is here.
"""
@inline runoff_timescale(mp::ModelParameters) =
    ZUO_OERLEMANS_C1 + ZUO_OERLEMANS_C2 * exp(-ZUO_OERLEMANS_C3 * mp.surface_slope)

"""
    hydraulic_conductivity_saturated(grain_radius, density) -> K_sat [m s-1]

Saturated hydraulic conductivity of snow/firn, Calonne et al. (2012) eq. 6:

    K_sat = 3r²·ρ_w·g/μ·exp(-0.013ρ)

with `r` the grain radius **in metres** and `μ` the dynamic viscosity of water. GEMB carries
`grain_radius` in millimetres (see `RE_NEW_SNOW`), so this converts; the conductivity goes as
`r²`, making a missed conversion a factor of 10⁶.

Cross-checked against the Community Firn Model's `hydrconducsat_Calonne`
(`CFM_main/darcy_funcs.py`).
"""
@inline function hydraulic_conductivity_saturated(grain_radius::Float64, density::Float64)
    r = grain_radius * 1e-3    # mm -> m
    return 3.0 * r^2 * DENSITY_WATER * GRAVITY / VISCOSITY_WATER * exp(-0.013 * density)
end

"""
    relative_permeability(grain_radius, density, saturation_effective) -> K_rel [-]

Relative hydraulic conductivity from the van Genuchten (1980) model with the Yamaguchi et al.
(2012) parameterization, as Hirashima et al. (2010) eq. 10 writes it:

    n     = 1 + 2.7e-3·(ρ/2r)^0.61        (Yamaguchi 2012 eq. 7, r in metres)
    m     = 1 - 1/n
    K_rel = θ_e^½·(1 - (1 - θ_e^(1/m))^m)²

`saturation_effective` is `θ_e`, the water content above irreducible as a fraction of the
drainable pore space, clamped away from both endpoints by the caller. `K_rel` rises from 0 at
`θ_e = 0` to 1 at `θ_e = 1`, so a barely-wet cell drains far more slowly than a saturated one
at the same conductivity — the nonlinearity that distinguishes this from a fixed timescale.

`grain_radius` is in millimetres, as everywhere in GEMB; the Yamaguchi fit needs `ρ/2r` in
kg m-4, so the conversion happens here.

Cross-checked against the Community Firn Model's `vG_Yama_params`/`krel_vG`
(`CFM_main/darcy_funcs.py`).
"""
@inline function relative_permeability(grain_radius::Float64, density::Float64,
    saturation_effective::Float64)
    r = grain_radius * 1e-3    # mm -> m
    n = 1.0 + 2.7e-3 * (density / (2.0 * r))^0.61
    m = 1.0 - 1.0 / n
    return sqrt(saturation_effective) *
           (1.0 - (1.0 - saturation_effective^(1.0 / m))^m)^2
end

"""
    apply_lateral_drainage!(cols, dt_seconds, mp) -> runoff [kg m-2]

Drain water held above irreducible saturation out of the column laterally, and return the
mass removed.

Runs after [`calculate_melt`](@ref) in [`gemb_core`](@ref), on the water that
[`pond_blocked_water!`](@ref) left standing. Only `cols.water` is touched — not `dz`,
`density`, or the cell count — so both grid invariants are untouched and the caller needs only
to add the returned mass to its runoff total. No energy term is needed either: the runoff
enthalpy in [`gemb_core`](@ref) is derived from the runoff mass at
[`specific_enthalpy_water`](@ref), the same value the water carried as pore water.

Per cell, with `excess = water - irreducible` and by `mp.runoff_method`:

- `:instantaneous` — returns 0 and touches nothing. Blocked water has already left inside
  [`calculate_melt`](@ref), and no cell holds any excess to drain.
- `:ZuoOerlemans` — `excess·min(1, dt/τ)` with `τ` from [`runoff_timescale`](@ref). Slope
  enters only through `τ`. The `min(1, ·)` matters at long steps: at the 25 d flat-terrain
  timescale a daily step drains 4% of the excess, but a monthly forcing step would otherwise
  ask for 120% of it.
- `:Darcy` — `ρ_w·dt·K_sat·K_rel·S`, capped at the excess. Darcy's law for a lateral flux
  under a hydraulic gradient equal to the surface slope, with conductivity from
  [`hydraulic_conductivity_saturated`](@ref) scaled by [`relative_permeability`](@ref). Slope,
  grain size, density, and saturation all enter, so unlike the timescale law this drains a
  coarse saturated layer orders of magnitude faster than a fine barely-wet one.

Both are strictly capped at the excess, so a cell is never drained below irreducible: the
scheme removes standing water, never capillary-held water.

`age` is untouched. GEMB carries one mean age per cell, over matrix and pore water together,
so the water leaving carries `age[i]` by construction and the mean of what remains is
unchanged — the same argument that makes melt removal age-neutral in
[`calculate_melt`](@ref).
"""
function apply_lateral_drainage!(cols::NamedTuple, dt_seconds::Float64, mp::ModelParameters)
    mp.runoff_method === :instantaneous && return 0.0

    zuo = mp.runoff_method === :ZuoOerlemans
    slope = mp.surface_slope
    # Hoisted: neither depends on the cell under `:ZuoOerlemans`, and the exponential in
    # `runoff_timescale` would otherwise be recomputed per cell.
    drain_fraction = zuo ? min(1.0, dt_seconds / runoff_timescale(mp)) : 0.0

    runoff = 0.0
    @inbounds for i in eachindex(cols.water)
        cols.water[i] <= WATER_TOLERANCE && continue

        density = cols.density[i]
        dz = cols.dz[i]
        capacity = pore_capacity(mp, density, dz)
        irreducible = (mp.density_ice - density) * irreducible_saturation(mp, density) * dz
        excess = cols.water[i] - irreducible
        excess <= WATER_TOLERANCE && continue

        if zuo
            drain = excess * drain_fraction
        else
            # Effective saturation: excess as a fraction of the drainable pore space. Clamped
            # off both endpoints because `K_rel` has an infinite derivative at 0 and the
            # `(1 - θ_e^(1/m))^m` term is evaluated at 0 for θ_e = 1; the same stabilisation
            # the Community Firn Model applies in `thetae_update`.
            drainable = capacity - irreducible
            drainable <= WATER_TOLERANCE && continue
            θ_e = clamp(excess / drainable, 1e-9, 1.0 - 1e-9)
            flux = DENSITY_WATER * dt_seconds * slope *
                   hydraulic_conductivity_saturated(cols.grain_radius[i], density) *
                   relative_permeability(cols.grain_radius[i], density, θ_e)
            drain = min(excess, flux)
        end

        cols.water[i] -= drain
        runoff += drain
    end
    return runoff
end
