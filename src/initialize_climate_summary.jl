"""
    ClimateSummary

Climatological summary of a [`ClimateForcing`](@ref), reduced to the handful of
scalars needed to guess a steady-state firn column. Built by
[`initialize_climate_summary`](@ref).

Every field is derived from a small, bounded number of passes over the forcing
series (one for the sums and the harmonic, plus one per iteration of the
albedo/melt fixed point — typically 1-6 in total), so building one is cheap
relative to any model integration.

# Mass balance [kg m-2 yr-1]
- `accumulation`: snowfall only — precipitation falling below
  `mp.rain_temperature_threshold`. Rain is excluded because it does not add mass
  to the column as snow.
- `rainfall`: the complement, precipitation falling as rain.
- `melt`: annual melt from a surface-energy-balance estimate (see
  [`initialize_climate_summary`](@ref)).
- `refreeze`: the part of `melt + rainfall` the column can refreeze, capped by
  one annual layer's cold content.
- `balance`: net annual mass balance, `accumulation + refreeze - melt`. **Its sign
  is the only regime discriminator**: positive buries snow and grows firn,
  non-positive exposes ice. Rainfall is not subtracted — see the note in
  [`initialize_climate_summary`](@ref).
- `accumulation_effective`: `accumulation + refreeze`, the mass flux that drives
  densification.

# Temperature [K]
- `temperature_air_mean`: mean annual air temperature.
- `temperature_amplitude`, `temperature_phase`: amplitude [K] and phase [rad] of
  the fitted annual cycle, `T(t) = T̄ + A·cos(ωt − φ)`.
- `latent_warming`: `LF·refreeze/(c_p·accumulation_effective)`, the mean
  warming from refreezing latent heat.

# Other means
- `wind_speed_mean` [m s-1], `pressure_air_mean` [Pa], and the observation
  heights, carried through for the physics the marcher reuses.
"""
struct ClimateSummary
    accumulation::Float64
    rainfall::Float64
    melt::Float64
    refreeze::Float64
    balance::Float64
    accumulation_effective::Float64
    temperature_air_mean::Float64
    temperature_amplitude::Float64
    temperature_phase::Float64
    latent_warming::Float64
    wind_speed_mean::Float64
    pressure_air_mean::Float64
    temperature_observation_height::Float64
    wind_observation_height::Float64
end

"""
    initialize_climate_summary(cf::ClimateForcing, mp::ModelParameters;
                               max_albedo_passes=20) -> ClimateSummary

Reduce a climate forcing to the scalars that determine a steady-state firn
column, in a few passes over the series (one to settle each of the coupled
quantities; see the albedo note below).

Three quantities are computed simultaneously in the pass:

**Snowfall vs. rainfall.** Precipitation is partitioned on
`mp.rain_temperature_threshold`, the same test [`calculate_accumulation`](@ref)
applies per timestep. This matters: the forcing's `precipitation_mean` is *total*
precipitation, so using it directly counts rain as accumulating mass.

**Melt, from a surface energy balance.** At each step the skin temperature that
closes

    SW(1−α) + LW↓ − εσT⁴ + QH + QE = 0

is found by Newton iteration, reusing [`turbulent_heat_flux`](@ref) for the
turbulent terms. Where that temperature would exceed the melt point, the surface
is held at `273.15 K` and the residual energy becomes melt, `Δmelt = Q·Δt/LF`.
This replaces a degree-day factor with the radiation, wind and humidity the
forcing already carries.

**The annual temperature cycle.** Three scalar accumulators give the
least-squares harmonic fit `T(t) = T̄ + A·cos(ωt − φ)`, which the marcher uses for
the damped thermal wave.

Melt depends on albedo, and albedo depends on the surface state melt helps
determine. That fixed point is *solved*, not relaxed: `g(α) = α_target(α) − α` is
monotone in `α`, so a secant iteration bracketed by `mp.albedo_ice` and
`mp.albedo_snow` reaches `ALBEDO_FIXED_POINT_TOLERANCE` in about five passes, and
`max_albedo_passes` is a backstop rather than the normal exit. (Plain
under-relaxed iteration needed 16–25 passes on the same sites — enough to exceed
any bounded budget and return an unconverged albedo without saying so.) Failure to
converge warns.

The albedo blend ([`_snow_cover_albedo`](@ref)) is continuous in the melt rate, so
a site near `b = 0` does not flip between a snow and an ice albedo. Cost is a few
passes over a vector, with no column integration.
"""
function initialize_climate_summary(cf::ClimateForcing, mp::ModelParameters;
    max_albedo_passes::Int=20)

    T_air = Float64.(parent(cf.temperature_air))
    precip = Float64.(parent(cf.precipitation))
    wind = Float64.(parent(cf.wind_speed))
    pressure = Float64.(parent(cf.pressure_air))
    sw = Float64.(parent(cf.shortwave_downward))
    lw = Float64.(parent(cf.longwave_downward))
    vapor = Float64.(parent(cf.vapor_pressure))

    n = length(T_air)
    n == 0 && error("initialize_climate_summary: forcing has no timesteps")

    dt = Float64(cf.time_step)                       # [s]
    dt_days = dt / 86400.0
    total_days = n * dt_days
    per_year = 365.25 / total_days                   # series total → annual rate

    T_mean = Float64(cf.temperature_air_mean)
    wind_mean = Float64(cf.wind_speed_mean)
    pressure_mean = sum(pressure) / n
    zT_obs = Float64(cf.temperature_observation_height)
    zW_obs = Float64(cf.wind_observation_height)

    # --- snowfall / rainfall split, and the annual temperature harmonic ---
    snow_total = 0.0
    rain_total = 0.0
    Σcos = 0.0
    Σsin = 0.0
    ω = 2π / SECONDS_PER_YEAR
    @inbounds for i in 1:n
        if T_air[i] <= mp.rain_temperature_threshold
            snow_total += precip[i]
        else
            rain_total += precip[i]
        end
        θ = ω * (i - 1) * dt
        Tp = T_air[i] - T_mean
        Σcos += Tp * cos(θ)
        Σsin += Tp * sin(θ)
    end

    # Least-squares harmonic: T' ≈ A·cos(ωt − φ) ⇒ A = 2/n·hypot(Σcos, Σsin).
    T_amplitude = 2.0 * sqrt(Σcos^2 + Σsin^2) / n
    T_phase = atan(Σsin, Σcos)

    accumulation = snow_total * per_year
    rainfall = rain_total * per_year

    # --- melt: solve the albedo/melt fixed point ---
    #
    # Each pass costs a sweep of the forcing, so this is solved rather than
    # relaxed. `g(α) = target(α) − α` is monotone and smooth, so a secant
    # iteration bracketed by the ice and snow albedos converges in ~5 passes.
    # Under-relaxed fixed-point iteration needed 16–25 to reach the same root,
    # which silently exceeded any bounded pass budget.
    # Returns the residual alongside the melt and refreeze it was computed from, so
    # the accepted iterate's values are carried out explicitly rather than left
    # behind as closure side-effects (which would also box them).
    albedo_residual(α) = _albedo_residual(α, T_air, wind, pressure, sw, lw, vapor,
        dt, per_year, zT_obs, zW_obs, accumulation, rainfall, T_mean, mp)

    # Snow side first, and only bracket with the ice side if it does not already
    # converge: the convergence test reads the *second* point, so evaluating the ice
    # side up front spends a full forcing sweep that a cold, melt-free site never
    # needs. (Measured: 2 sweeps -> 1 at dry-cold and percolation sites; unchanged at
    # melting ones, which need the bracket anyway.)
    α1 = mp.albedo_snow
    g1, melt, refreeze = albedo_residual(α1)
    albedo = α1
    converged = abs(g1) < ALBEDO_FIXED_POINT_TOLERANCE
    α0, g0 = mp.albedo_ice, 0.0
    if !converged
        g0, _, _ = albedo_residual(α0)
    end
    for _ in 1:max_albedo_passes
        converged && break
        denom = g1 - g0
        abs(denom) < 1e-14 && break
        α2 = clamp(α1 - g1 * (α1 - α0) / denom, mp.albedo_ice, mp.albedo_snow)
        α0, g0 = α1, g1
        α1 = α2
        g1, melt, refreeze = albedo_residual(α1)
        albedo = α1
        converged = abs(g1) < ALBEDO_FIXED_POINT_TOLERANCE
    end
    if !converged
        @warn "initialize_climate_summary: albedo/melt fixed point did not converge; \
               the initial guess may be inconsistent with the forcing." albedo residual = g1
    end

    # Net annual balance. Rain does **not** appear as a loss term: it is not part of
    # `accumulation` (which is snowfall only), so the column never gained it, and
    # rain that runs off carries away only its own mass — it does not remove snow.
    # The part of it that stays is already in `refreeze`. Subtracting `rainfall`
    # here would report a site with zero melt and real snowfall as ablating purely
    # because it also rains, which inverts the scheme's only regime discriminator.
    balance = accumulation + refreeze - melt
    A_eff = accumulation + refreeze
    # Mean warming from refreezing one year of liquid water into the effective accumulation.
    # Under `:constant` this is the reference's `LF·R/(c·A_eff)`, evaluated in that grouping.
    # Under `:CuffeyPaterson` the same energy buys more warming in cold firn, so it is applied
    # as an enthalpy increment about the mean temperature rather than divided by a fixed `c`.
    latent_warming = if A_eff <= 0.0
        0.0
    elseif mp.heat_capacity_method === :constant
        LF * refreeze / (mp.heat_capacity_ice * A_eff)
    else
        add_energy_temperature(mp, T_mean, A_eff, LF * refreeze) - T_mean
    end

    return ClimateSummary(accumulation, rainfall, melt, refreeze, balance, A_eff,
        T_mean, T_amplitude, T_phase, latent_warming,
        wind_mean, pressure_mean, zT_obs, zW_obs)
end

"""
    _albedo_residual(α, ...) -> (g, melt, refreeze)

One evaluation of the albedo fixed point: the melt an albedo of `α` produces, the
refreeze that implies, and the residual `g = α_target − α` whose root the secant
iteration in [`initialize_climate_summary`](@ref) seeks.

Returns the melt and refreeze alongside the residual so the caller keeps the pair
consistent with the accepted albedo without re-running a forcing sweep.
"""
function _albedo_residual(α::Float64, T_air::Vector{Float64}, wind::Vector{Float64},
    pressure::Vector{Float64}, sw::Vector{Float64}, lw::Vector{Float64},
    vapor::Vector{Float64}, dt::Float64, per_year::Float64, zT_obs::Float64,
    zW_obs::Float64, accumulation::Float64, rainfall::Float64, T_mean::Float64,
    mp::ModelParameters)

    melt = _seb_annual_melt(T_air, wind, pressure, sw, lw, vapor,
        α, dt, per_year, zT_obs, zW_obs, mp)
    refreeze = _cold_content_refreeze(melt, rainfall, accumulation, T_mean, mp)
    g = _snow_cover_albedo(accumulation + refreeze, melt, mp) - α
    return g, melt, refreeze
end

"""
    _snow_fraction(accumulation_effective, melt) -> f

Fraction of the year the surface is expected to carry snow,
`f = A_eff/(A_eff + melt)`: with no melt the surface is snow all year; with melt
far exceeding accumulation it is bare ice.

This is the single definition of the snow/ice mixing weight, used by
[`_snow_cover_albedo`](@ref) to blend the two albedos.
"""
@inline function _snow_fraction(accumulation_effective::Real, melt::Real)
    A = max(Float64(accumulation_effective), 0.0)
    M = max(Float64(melt), 0.0)
    total = A + M
    total <= 0.0 && return 1.0
    return A / total
end

"""
    _snow_fraction_from_albedo(α, mp) -> f

Recover the snow fraction from a blended albedo — the exact inverse of
[`_snow_cover_albedo`](@ref), kept adjacent to it so the pair is maintained
together.

[`_seb_annual_melt`](@ref) needs the fraction to interpolate surface roughness, but
it is called *with* a candidate albedo and *before* the melt that would determine
the fraction directly, so inverting is the only order that closes. **If the albedo
blend ever stops being linear in `f`, this inverse must change with it** — hence
one named function rather than the arithmetic inlined at the use site.
"""
@inline function _snow_fraction_from_albedo(α::Real, mp::ModelParameters)
    span = mp.albedo_snow - mp.albedo_ice
    return clamp((Float64(α) - mp.albedo_ice) / max(span, eps()), 0.0, 1.0)
end

"""
    _snow_cover_albedo(accumulation_effective, melt, mp) -> α

Annual-mean surface albedo implied by an accumulation and a melt rate, blending
`mp.albedo_snow` and `mp.albedo_ice` by the fraction of the year the surface is
expected to carry snow.

The blend is `f = A_eff / (A_eff + melt)`: with no melt the surface is snow all
year; with melt far exceeding accumulation it is bare ice. This is deliberately
**continuous**. A binary pick on the sign of the mass balance would reintroduce
exactly the threshold this scheme exists to remove — one level up, in the albedo
that determines the melt that determines the balance — and would make the whole
initialization jump discontinuously as a site crosses `b = 0`.
"""
@inline function _snow_cover_albedo(accumulation_effective::Real, melt::Real,
    mp::ModelParameters)
    f = _snow_fraction(accumulation_effective, melt)
    return f * mp.albedo_snow + (1.0 - f) * mp.albedo_ice
end

"""
    _cold_content_refreeze(melt, rainfall, accumulation, T_mean, mp) -> R [kg m-2 yr-1]

Annual refreeze, limited by both the available liquid water and the cold content
of one annual accumulation layer:

    R = min(melt + rainfall, accumulation·(h(273.15) − h(T̄)) / LF)

The second term is the standard Pfeffer-style retention capacity: the energy
needed to warm one year's accumulation to the melt point, expressed as a mass of
refrozen water. Under `:constant` the enthalpy difference is `c_p·(273.15 − T̄)`,
the reference's form. A temperate column (`T̄ ≥ 273.15`) has no cold content and
refreezes nothing.
"""
@inline function _cold_content_refreeze(melt::Real, rainfall::Real,
    accumulation::Real, T_mean::Real, mp::ModelParameters)
    available = Float64(melt) + Float64(rainfall)
    capacity = if mp.heat_capacity_method === :constant
        # Grouped as the reference groups it, so the default path is bit-identical.
        mp.heat_capacity_ice * max(CtoK - Float64(T_mean), 0.0) * Float64(accumulation) / LF
    else
        cold_content_mass(mp, min(Float64(T_mean), CtoK), Float64(accumulation))
    end
    return min(available, capacity)
end

"""
    _seb_annual_melt(T_air, wind, pressure, sw, lw, vapor, albedo, dt, per_year,
                     zT_obs, zW_obs, mp) -> melt [kg m-2 yr-1]

Annual melt from a zero-layer surface energy balance over the forcing series.

For each step, [`_seb_skin_temperature`](@ref) finds the skin temperature closing
the flux balance. If that exceeds the melt point the surface is pinned at
`273.15 K` and the leftover energy flux `Q` is converted to melt via
`Q·Δt/LF`. Sublimation/condensation is deliberately ignored — it is a small term
for an initial guess, and including it would need the mass feedback the marcher
resolves anyway.
"""
function _seb_annual_melt(T_air::Vector{Float64}, wind::Vector{Float64},
    pressure::Vector{Float64}, sw::Vector{Float64}, lw::Vector{Float64},
    vapor::Vector{Float64}, albedo::Float64, dt::Float64, per_year::Float64,
    zT_obs::Float64, zW_obs::Float64, mp::ModelParameters)

    # Roughness interpolated between the snow and ice values `calculate_temperature`
    # uses, on the same snow-cover fraction the albedo blend implies — so a partly
    # snow-covered surface does not have to be called one or the other.
    f_snow = _snow_fraction_from_albedo(albedo, mp)
    z0 = f_snow * Z0_SNOW_DRY + (1.0 - f_snow) * Z0_ICE
    zT = z0 * mp.surface_roughness_effective_ratio
    zQ = z0 * mp.surface_roughness_effective_ratio

    # Emissivity through the same selector the run uses, so a non-`:uniform`
    # `mp.emissivity_method` is not silently initialized under a different one. The
    # surface is fresh snow at initialization, hence `RE_NEW_SNOW`.
    ε, _ = _emissivity_initialize(RE_NEW_SNOW, mp)

    melt_total = 0.0
    @inbounds for i in eachindex(T_air)
        cfs = _seb_forcing_step(dt, T_air[i], pressure[i], wind[i], sw[i], lw[i],
            vapor[i], zT_obs, zW_obs)
        density_air = air_density(pressure[i], T_air[i])
        sw_net = sw[i] * (1.0 - albedo)

        T_skin = _seb_skin_temperature(sw_net, lw[i], ε, density_air, z0, zT, zQ, cfs)
        T_skin <= CtoK && continue

        # Surface cannot exceed the melt point: the residual flux drives melt.
        Q = _seb_residual(CtoK, sw_net, lw[i], ε, density_air, z0, zT, zQ, cfs)
        Q > 0.0 && (melt_total += Q * dt / LF)
    end

    return melt_total * per_year
end

"""
    _seb_residual(T_surface, sw_net, lw_in, ε, density_air, z0, zT, zQ, cfs) -> Q [W m-2]

Net energy flux into the surface at skin temperature `T_surface`:

    Q = SW_net + LW↓ − εσT⁴ + QH + QE

Turbulent terms come from [`turbulent_heat_flux`](@ref). Positive `Q` warms (or
melts) the surface.
"""
@inline function _seb_residual(T_surface::Float64, sw_net::Float64, lw_in::Float64,
    ε::Float64, density_air::Float64, z0::Float64, zT::Float64, zQ::Float64,
    cfs::ClimateForcingStep)
    shf, lhf, _ = turbulent_heat_flux(T_surface, density_air, z0, zT, zQ, cfs)
    return sw_net + ε * lw_in - ε * SB * T_surface^4 + shf + lhf
end

"""
    _seb_skin_temperature(sw_net, lw_in, ε, density_air, z0, zT, zQ, cfs) -> T [K]

Skin temperature closing the surface energy balance, by secant iteration on
[`_seb_residual`](@ref) bracketed around the air temperature.

The balance is monotone decreasing in `T_surface` (both the `−εσT⁴` radiative
term and the turbulent terms oppose warming), so the iteration is well behaved. A
few iterations suffice for an initial guess; the result is only used to decide
whether the surface reaches melting and by how much energy it overshoots.
"""
@inline function _seb_skin_temperature(sw_net::Float64, lw_in::Float64, ε::Float64,
    density_air::Float64, z0::Float64, zT::Float64, zQ::Float64,
    cfs::ClimateForcingStep; max_iter::Int=8, tol::Float64=1e-3)

    T_a = cfs.temperature_air
    T0 = T_a - 5.0
    T1 = T_a + 5.0
    f0 = _seb_residual(T0, sw_net, lw_in, ε, density_air, z0, zT, zQ, cfs)
    f1 = _seb_residual(T1, sw_net, lw_in, ε, density_air, z0, zT, zQ, cfs)

    for _ in 1:max_iter
        denom = f1 - f0
        abs(denom) < 1e-12 && break
        T2 = T1 - f1 * (T1 - T0) / denom
        # Keep the iterate physical; a runaway secant step is not informative.
        T2 = clamp(T2, T_a - 60.0, T_a + 60.0)
        T0, f0 = T1, f1
        T1 = T2
        f1 = _seb_residual(T1, sw_net, lw_in, ε, density_air, z0, zT, zQ, cfs)
        abs(T1 - T0) < tol && break
    end

    return T1
end

"""
    _seb_forcing_step(dt, T_air, pressure, wind, sw, lw, vapor, zT_obs, zW_obs)

Minimal [`ClimateForcingStep`](@ref) for the surface-energy-balance estimate.
Only the fields [`turbulent_heat_flux`](@ref) reads are meaningful; the
climatological-mean and time-varying-parameter slots are unused here and set to
zero.
"""
@inline function _seb_forcing_step(dt::Float64, T_air::Float64, pressure::Float64,
    wind::Float64, sw::Float64, lw::Float64, vapor::Float64,
    zT_obs::Float64, zW_obs::Float64)
    return ClimateForcingStep(; dt,
        temperature_air=T_air, pressure_air=pressure, wind_speed=wind,
        shortwave_downward=sw, longwave_downward=lw, vapor_pressure=vapor,
        temperature_air_mean=T_air, wind_speed_mean=wind,
        temperature_observation_height=zT_obs, wind_observation_height=zW_obs)
end
