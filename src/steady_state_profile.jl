"""
    steady_state_profile(dz, cs::ClimateSummary, mp::ModelParameters; n_age=2000)
        -> NamedTuple

Build a steady-state estimate of the full firn column state by marching a single
parcel of snow forward in age and recording every state variable as it is buried.

This is the initial-guess engine behind [`initialize_profile`](@ref). It replaces a
family of regime-specific special cases with one continuous calculation: the
burial rate is the net annual mass balance `cs.balance`, so an accumulation site
buries snow and grows a firn column, while an ablation site buries nothing and
yields the solid ice it exposes.

**There is no threshold to cross** — the sign of `cs.balance` is the only
discriminator, and the profile varies continuously as a site approaches
`balance = 0` from above (the surface density rises smoothly toward
`mp.density_ice` as the burial rate falls). There *is* a branch at `balance <= 0`,
because the age step `dt = dz_step·ρ/b` is singular at `b = 0`; it returns the
`b → 0⁺` limit of the march directly rather than stepping to it.

Each state variable is carried along the age coordinate using the same physics
the model integrates transiently, so the guess is consistent with the run it
seeds:

- **density** relaxes as `ρ(t+dt) = ρᵢ − (ρᵢ − ρ)·exp(−c·dt)` (Arthern et al.
  2010 eq. 1) with `c` from [`_densification_rate`](@ref) for the configured
  `mp.densification_method` — *not* always Herron–Langway.
- **depth** follows from mass balance, `z(t) = b·∫₀ᵗ ds/ρ(s)`, where
  `b = cs.balance` [kg m-2 yr-1] is the mass buried per year.
- **temperature** is a damped annual wave about the mean surface temperature, warmed
  near the surface by refreezing latent heat; see
  [`_steady_state_temperature`](@ref).
- **grain size** is evolved by [`calculate_grain_size`](@ref) itself, so
  Marbouty/Brun metamorphism is not duplicated here.
- **water** is the irreducible content at the melt point, using the same formula
  as [`calculate_melt`](@ref).
- **age** is the march's own time coordinate — the accumulated `dt` at each depth, in days,
  which is the residence time a parcel buried at `cs.balance` would have. Below the
  firn/ice transition the march has stopped, so it is continued at the ice burial rate
  `ρᵢ/b`; it is zero everywhere for an ablation column, where nothing is buried at all.

Returns a `NamedTuple` of `Vector{Float64}`s at each depth in `z_center`, with
keys `density`, `temperature`, `water`, `grain_radius`, `grain_dendricity`,
`grain_sphericity`, `age`, plus the scalar `ice_depth` used to size the grid: the depth
at which the column reaches `mp.density_ice` (to within
`DENSITY_FIRN_TOLERANCE`), below which it carries no further firn information.
`ice_depth` is `0.0` for an ablation column (ice at the surface) and `max_depth`
if the march never reached ice.

# Arguments
- `dz`: cell thicknesses [m], surface first, as [`initialize_grid`](@ref) returns.
  Cell centers are derived with [`dz2z`](@ref); the thicknesses themselves are
  needed because water content is a per-cell mass, not a density.
- `cs`: the climate summary from [`initialize_climate_summary`](@ref).
- `mp`: model parameters; supplies the densification method, ice density and
  irreducible-water saturation.
- `n_age`: number of march steps spanning the column. The march is stepped in
  *depth* (`max_depth/n_age` per step) with the age step derived from the local
  density, so resolution is uniform in depth and the step count does not depend on
  the burial rate.

# References
- Arthern, R. J., Vaughan, D. G., Rankin, A. M., Mulvaney, R., and Thomas, E. R. (2010).
  In situ measurements of Antarctic snow compaction compared with predictions of models.
  *J. Geophys. Res.* 115, F03011. eq. 1.
- Herron, M. M. and Langway, C. C. (1980). Firn densification: an empirical model.
  *J. Glaciol.* 25, 373-385.
"""
function steady_state_profile(dz::AbstractVector, cs::ClimateSummary,
    mp::ModelParameters; n_age::Int=2000)

    dz = Float64.(dz)
    depth = -dz2z(dz)                               # positive downward [m]
    m = length(depth)
    max_depth = maximum(depth)

    ρi = mp.density_ice
    b = cs.balance                                  # burial rate [kg m-2 yr-1]

    # Surface values: true fresh snow, matching calculate_accumulation.
    ρ0 = clamp(fresh_snow_density(mp, cs.temperature_air_mean, cs.accumulation_effective,
                                  cs.wind_speed_mean), 1.0, ρi)

    # ---------------------------------------------------------------------
    # Net ablation (or no burial at all): nothing is ever buried, so the column
    # is the ice it exposes. This is the old is_ice branch, reached as a limit
    # rather than by a threshold test.
    # ---------------------------------------------------------------------
    if b <= 0.0
        # The surface mean, for the same reason `_steady_state_temperature` uses it: this fills
        # the whole column, including the Dirichlet bottom cell. `latent_warming` is retained
        # here without a depth decay because an ablation column has no annual layer for it to
        # decay over — nothing is buried — and it is zero at any site with no refreeze anyway.
        T_surf = min(cs.temperature_surface_mean + cs.latent_warming, CtoK)
        return (
            density          = fill(ρi, m),
            temperature      = fill(T_surf, m),
            water            = zeros(m),
            grain_radius     = fill(GRAIN_RADIUS_ICE, m),
            grain_dendricity = zeros(m),
            grain_sphericity = zeros(m),
            # Nothing is ever buried, so the march never advances and there is no residence
            # time to report. Zero, not the true (large, unknowable) age of exposed ice.
            age              = zeros(m),
            ice_depth        = 0.0,
        )
    end

    # `accumulation_effective` is already the snowfall-plus-refreeze flux, so the
    # `densification_accumulation` gate has nothing to switch here — the initializer never
    # used total precipitation. Only the temperature gate applies.
    p = DensificationCoeffs(mp, cs.accumulation_effective,
        mp.mean_temperature_method === :arrhenius ? cs.temperature_air_effective :
                                                    cs.temperature_air_mean)

    # Age-marched curves, interpolated onto z_center at the end.
    z_curve   = Float64[0.0]
    age_curve = Float64[0.0]                        # cumulative age along the march [yr]
    ρ_curve   = Float64[ρ0]
    re_curve  = Float64[RE_NEW_SNOW]
    gdn_curve = Float64[GDN_NEW_SNOW]
    gsp_curve = Float64[GSP_NEW_SNOW]

    ρ   = ρ0
    z   = 0.0
    age = 0.0                                       # residence time of the marched parcel [yr]

    # Grain state marched as a 3-cell stencil: `calculate_grain_size` derives the
    # temperature gradient from neighbours (`_grain_gradient` returns 0 for a
    # single cell, which would silently disable Marbouty metamorphism), so carry
    # the parcel in the middle with its neighbours half an age step either side.
    gs_T     = fill(cs.temperature_air_mean, 3)
    gs_dz    = zeros(3)                              # set from the step below
    gs_ρ     = fill(ρ0, 3)
    gs_w     = zeros(3)
    gs_re    = fill(RE_NEW_SNOW, 3)
    gs_gdn   = fill(GDN_NEW_SNOW, 3)
    gs_gsp   = fill(GSP_NEW_SNOW, 3)

    ice_depth = 0.0
    iter = 0
    dz_step = _march_step_depth(max_depth, n_age)
    max_iter = _march_max_iterations(n_age)

    @inbounds while z < max_depth && ρ < ρi - DENSITY_FIRN_TOLERANCE && iter < max_iter
        # Age step sized to advance a fixed depth, then one step of the shared
        # density kernel so this march and `steady_state_density` cannot diverge
        # in the density law.
        dt = _march_age_step(dz_step, ρ, b)
        z, ρ = _march_density_step(p, ρ, z, b, ρi, dt, gs_T[2])
        dt_seconds = dt * SECONDS_PER_YEAR
        # The march's own time coordinate. `dt` is exactly the interval over which the parcel
        # was buried from the previous depth to this one, so accumulating it *is* the residence
        # time at depth `z` — the same quantity `age` counts in the transient run.
        age += dt

        if ice_depth == 0.0 && ρ >= ρi - DENSITY_FIRN_TOLERANCE
            ice_depth = z
        end

        # --- grain growth, delegated to the transient implementation ---
        # Cell thickness first: the water content below is a per-cell mass, so it
        # needs this step's own thickness rather than a carried-over one.
        dz_here = max(b * dt / ρ, D_TOLERANCE)

        # --- temperature and water at this depth ---
        T_here = _steady_state_temperature(z, ρ, cs, mp)
        w_here = _irreducible_water(T_here, ρ, dz_here, mp)

        # Centre cell is the parcel; neighbours bracket it in depth so the
        # gradient stencil sees the real local temperature profile.
        gs_dz[1] = gs_dz[2] = gs_dz[3] = dz_here
        gs_T[1] = _steady_state_temperature(max(z - dz_here, 0.0), ρ, cs, mp)
        gs_T[2] = T_here
        gs_T[3] = _steady_state_temperature(z + dz_here, ρ, cs, mp)
        gs_ρ[1] = gs_ρ[2] = gs_ρ[3] = ρ
        gs_w[1] = gs_w[2] = gs_w[3] = w_here

        cfs_grain = _grain_forcing_step(dt_seconds, cs)
        gs_re, gs_gdn, gs_gsp = calculate_grain_size(gs_T, gs_dz, gs_ρ, gs_w,
            gs_re, gs_gdn, gs_gsp, cfs_grain, mp)

        push!(z_curve, z)
        push!(age_curve, age)
        push!(ρ_curve, ρ)
        push!(re_curve, gs_re[2])
        push!(gdn_curve, gs_gdn[2])
        push!(gsp_curve, gs_gsp[2])
        iter += 1
    end

    # Extend past the deepest requested cell so interpolation is well posed.
    if z_curve[end] < max_depth
        push!(z_curve, max_depth + 1.0)
        # Age is the one curve that must not be extended by a constant: the march stops when
        # the column reaches ice, but burial does not, so a flat tail would report the whole
        # solid-ice base as the same age as its top. Continue at the ice burial rate, `ρi/b`
        # years per metre, which is exactly the `dt` the loop would have taken.
        push!(age_curve, age + (max_depth + 1.0 - z_curve[end-1]) * ρi / b)
        push!(ρ_curve, ρi)
        push!(re_curve, re_curve[end])
        push!(gdn_curve, gdn_curve[end])
        push!(gsp_curve, gsp_curve[end])
    end

    # `ice_depth` was recorded in the loop at the crossing. If the march never got
    # there the column is still firn at its base, so there is no ice depth to
    # report and the full requested depth stands.
    if ice_depth == 0.0
        ice_depth = max_depth
    end

    # Interpolate the marched curves onto the model depths. Constant
    # extrapolation handles the z=0 surface and any overshoot at depth.
    density = _interp_curve(ρ_curve, z_curve, depth; lo=ρ0, hi=ρi)
    grain_radius = _interp_curve(re_curve, z_curve, depth; lo=RE_NEW_SNOW, hi=GRAIN_RADIUS_ICE)
    grain_dendricity = _interp_curve(gdn_curve, z_curve, depth; lo=0.0, hi=1.0)
    grain_sphericity = _interp_curve(gsp_curve, z_curve, depth; lo=0.0, hi=1.0)
    # In days, matching the transient `age` counter (`cfs.dt/86400`). Converted with
    # `SECONDS_PER_YEAR`, the same year length the march already uses for `dt_seconds`, so the
    # marched age and the grain growth driven by it cannot disagree about how long a step was.
    # `hi` is unreachable given the tail pushed above (which extends past `max_depth`); the
    # surface is age 0 by construction.
    days_per_year = SECONDS_PER_YEAR / 86400.0
    age_profile = _interp_curve(age_curve .* days_per_year, z_curve, depth;
        lo=0.0, hi=age_curve[end] * days_per_year)

    temperature = Vector{Float64}(undef, m)
    water = Vector{Float64}(undef, m)
    @inbounds for i in 1:m
        temperature[i] = _steady_state_temperature(depth[i], density[i], cs, mp)
        water[i] = _irreducible_water(temperature[i], density[i], dz[i], mp)
    end

    return (
        density          = density,
        temperature      = temperature,
        water            = water,
        grain_radius     = grain_radius,
        grain_dendricity = grain_dendricity,
        grain_sphericity = grain_sphericity,
        age              = age_profile,
        ice_depth        = ice_depth,
    )
end

# `p` is a `DensificationCoeffs`, left unannotated because `calculate_density.jl`
# (which defines it) is included after this file; `@inline` specializes regardless.
"""
    _march_density_step(p::DensificationCoeffs, ρ, z, b, ρi, dt, T) -> (z, ρ)

Advance one step of a steady-state density march by the age step `dt` [yr]:

    ρ  = ρᵢ − (ρᵢ − ρ)·exp(−c·dt)   (Arthern et al. 2010 eq. 1)
    z += b·dt/ρ_mid                 (mass balance, trapezoidal in ρ)

This is the single definition of the relaxation-and-burial arithmetic, shared by
[`steady_state_profile`](@ref) and [`steady_state_density`](@ref) so the density
law cannot drift between them.

The two differ only in how they choose `dt`, which is deliberate and stays with
each caller: `steady_state_profile` derives it from a target *depth* advance (see
[`_march_step_depth`](@ref)), while `steady_state_density` uses a fixed age step
to reproduce the classic Herron & Langway discretization bit-for-bit.
"""
@inline function _march_density_step(p, ρ::Float64, z::Float64,
    b::Float64, ρi::Float64, dt::Float64, T::Float64)

    c = _densification_rate(p, ρ, T)
    ρ_next = ρi - (ρi - ρ) * exp(-c * dt)
    if ρ_next > ρi
        ρ_next = ρi
    end

    ρ_mid = 0.5 * (ρ + ρ_next)
    return z + b * dt / ρ_mid, ρ_next
end

"""
    _march_step_depth(max_depth, n_age) -> dz_step [m]
    _march_age_step(dz_step, ρ, b) -> dt [yr]
    _march_max_iterations(n_age) -> Int

Step-size policy for [`steady_state_profile`](@ref)'s march, which is stepped in
**depth** rather than age: `dz_step` is the target advance per step and the age
step follows from the local density, `dt = dz_step·ρ/b` (the mass balance
`dz = b·dt/ρ`, inverted).

This keeps resolution uniform in the coordinate the result is sampled on, and
makes the step count independent of the burial rate — with a fixed age step, a
low-accumulation site needed orders of magnitude more steps to reach the same
depth, hit the iteration cap partway down, and had the rest of its column filled
by extrapolating straight to ice.

Using the current density for `dt` is first-order in `dz_step`, the same order as
the density update itself. `_march_max_iterations` is a backstop against a
pathological rate law, not the normal exit: the depth bound is what terminates a
healthy march.
"""
@inline _march_step_depth(max_depth::Real, n_age::Integer) = Float64(max_depth) / n_age
@inline _march_age_step(dz_step::Float64, ρ::Float64, b::Float64) = dz_step * ρ / b
@inline _march_max_iterations(n_age::Integer) = 4 * n_age

"""
    _interp_curve(values, z_curve, depth; lo, hi) -> Vector{Float64}

Linearly interpolate an age-marched curve onto `depth`, with constant
extrapolation beyond either end and the result clamped to `[lo, hi]`.
"""
function _interp_curve(values::Vector{Float64}, z_curve::Vector{Float64},
    depth::Vector{Float64}; lo::Float64, hi::Float64)
    itp = DataInterpolations.LinearInterpolation(values, z_curve;
        extrapolation=DataInterpolations.ExtrapolationType.Constant)
    out = Vector{Float64}(undef, length(depth))
    @inbounds for i in eachindex(depth)
        out[i] = clamp(itp(depth[i]), lo, hi)
    end
    return out
end

"""
    _steady_state_temperature(z, ρ, cs::ClimateSummary, mp::ModelParameters) -> T [K]

Temperature at depth `z` [m, positive down] for a column in thermal steady state:
a damped annual wave superposed on the mean, warmed near the surface by the latent heat that
refreezing releases.

    T(z) = T_surface + ΔT_latent·exp(−z/z_annual) + A_T·exp(−z/d)·cos(φ − z/d)

The damping depth `d = sqrt(2K/(ρ·c_p·ω))` uses the local density's own thermal
conductivity from [`thermal_conductivity`](@ref), so the wave decays through a
physically consistent column rather than a nominal one. `ω = 2π/yr`.

## The mean is the *surface* mean, not the air mean

The column exchanges heat with the atmosphere only through the surface energy balance, so its
deep mean tends to `cs.temperature_surface_mean` — the mean skin temperature that balance
produces — and not to the screen-level `cs.temperature_air_mean`. The two differ by the whole
radiative and turbulent budget, and the difference is neither small nor of predictable sign: at
the synthetic `test_1` site the surface is 5.1 K *colder* than the air.

This matters more than a mean-state error normally would, because the deepest cell is a
**Dirichlet reservoir** — `calculate_temperature` returns `temperature[m]` bit-unchanged by
construction — so whatever this function returns at the base of the column is a boundary
condition for the entire run, which no length of spinup relaxes. Using the air mean pinned the
base ~9 K too warm at `test_1` and biased the converged mass-weighted mean density by
~11 kg m-3. See `bench/calibrate_initial_guess.jl` for the measurement and the self-consistency
criterion (zero temperature jump across that frozen cell) it is checked against.

## Latent warming decays over the annual layer

`ΔT_latent` is derived in [`initialize_climate_summary`](@ref) as the warming of **one annual
accumulation layer** by one year of refreezing, `LF·R/(c_p·A_eff)`. Applying it at every depth
therefore double-counts it: the heat released into this year's layer is not also released into
firn buried a century ago. It is applied with an exponential decay over `z_annual`, the depth
one year of accumulation occupies at the local density (`A_eff/ρ`), so the near-surface warming
the reference intends is preserved while the deep asymptote is the surface mean alone.

It remains self-limiting where it acts: at the cold-content cap it equals
`(273.15 − T̄)·accumulation/A_eff`, strictly less than the gap to the melt point. The
`min(·, CtoK)` clamp is therefore belt-and-braces, and also guards the ablation case where the
mean itself is above freezing.
"""
@inline function _steady_state_temperature(z::Real, ρ::Real, cs::ClimateSummary,
    mp::ModelParameters)
    T_bar = cs.temperature_surface_mean
    d = thermal_damping_depth(T_bar, ρ, mp)

    # Depth one year of accumulation occupies at this density. Floored so a vanishing
    # accumulation cannot divide by zero; at that limit there is no annual layer to warm and
    # the exponential collapses to zero below the surface, which is the correct limit.
    z_annual = max(cs.accumulation_effective / max(Float64(ρ), 1.0), D_TOLERANCE)

    latent = cs.latent_warming * exp(-Float64(z) / z_annual)
    decay = exp(-z / d)
    T = T_bar + latent + cs.temperature_amplitude * decay * cos(cs.temperature_phase - z / d)
    return min(T, CtoK)
end

"""
    thermal_damping_depth(T, ρ, mp) -> d [m]

E-folding depth of the annual thermal wave in a medium of density `ρ` at
temperature `T`:

    d = sqrt(2K/(ρ·c_p·ω)),   ω = 2π/yr

with `K` from [`thermal_conductivity`](@ref). About 3.3 m at ice density.

Used both to decay the wave through the initialized column
([`_steady_state_temperature`](@ref)) and to set the floor on the column depth
that resolves it ([`_derive_column_depth`](@ref)), so the two agree by
construction.
"""
@inline function thermal_damping_depth(T::Real, ρ::Real, mp::ModelParameters)
    K = thermal_conductivity(T, ρ, mp)
    ω = 2π / SECONDS_PER_YEAR                       # [rad s-1]
    return sqrt(2 * K / (Float64(ρ) * heat_capacity(mp, T) * ω))
end

"""
    _irreducible_water(T, ρ, dz, mp) -> water [kg m-2]

Irreducible (capillary-held) water content of a cell at the melt point, zero for
a cold cell. Uses the same expression as [`calculate_melt`](@ref),
`(ρᵢ − ρ)·S·(M/ρ)` with `M = ρ·dz`, so initialization and runtime agree on what
"irreducible" means and the first timestep has no water to redistribute.

Two gates are added here, because this seeds a column rather than routing water through
an existing one: a cold cell holds nothing, and a cell at or past pore close-off
(`DENSITY_PORE_CLOSEOFF`) has no *connected* pore space to hold water in, so it
starts dry rather than with the small amount `(ρᵢ − ρ)` alone would imply. Under
`water_irreducible_method = :ColeouLesaffre` [`calculate_melt`](@ref) applies the
close-off gate too; under `:constant` it does not.
"""
@inline function _irreducible_water(T::Real, ρ::Real, dz::Real, mp::ModelParameters)
    T < CtoK - T_TOLERANCE && return 0.0
    ρ >= DENSITY_PORE_CLOSEOFF - D_TOLERANCE && return 0.0   # no connected pore space
    M = Float64(ρ) * Float64(dz)
    return (mp.density_ice - Float64(ρ)) * irreducible_saturation(mp, Float64(ρ)) * (M / Float64(ρ))
end

"""
    _grain_forcing_step(dt_seconds, cs::ClimateSummary) -> ClimateForcingStep

Minimal `ClimateForcingStep` for the age march's [`calculate_grain_size`](@ref)
call. Grain growth reads only `dt` from the forcing step (metamorphism depends on
the column's own temperature, density and water), so the remaining fields carry
the climatological means and are not otherwise used.
"""
@inline function _grain_forcing_step(dt_seconds::Float64, cs::ClimateSummary)
    return ClimateForcingStep(; dt=dt_seconds,
        temperature_air=cs.temperature_air_mean,
        pressure_air=cs.pressure_air_mean,
        wind_speed=cs.wind_speed_mean,
        temperature_air_mean=cs.temperature_air_mean,
        wind_speed_mean=cs.wind_speed_mean,
        precipitation_mean=cs.accumulation_effective,
        temperature_observation_height=cs.temperature_observation_height,
        wind_observation_height=cs.wind_observation_height)
end

"""
    steady_state_density(z_center, T_mean, accumulation, ρ0, ρ_ice, mp; n_age=2000)

Steady-state firn density profile for `mp.densification_method`, as a standalone
vector.

This is the density-only entry point to the same march
[`steady_state_profile`](@ref) performs — both drive
[`_march_density_step`](@ref), so there is one stepping scheme, not two — kept for
direct use and for testing the density law in isolation. With
`mp.densification_method = :HerronLangway` it reproduces the classic
Herron & Langway (1980) analytic profile.

# Arguments
- `z_center`: cell-center heights [m], negative below surface.
- `T_mean`: mean annual temperature [K].
- `accumulation`: mean annual accumulation [kg m-2 yr-1].
- `ρ0`: surface (fresh-snow) density [kg m-3].
- `ρ_ice`: density of ice [kg m-3], the asymptote and clamp.
- `mp`: supplies the densification method and its coefficients.
"""
function steady_state_density(z_center::AbstractVector, T_mean::Real,
    accumulation::Real, ρ0::Real, ρ_ice::Real, mp::ModelParameters; n_age::Int=2000)

    depth = -Float64.(z_center)
    m = length(depth)

    # No accumulation → no firn column forms; the profile is ice.
    accumulation <= 0.0 && return fill(Float64(ρ_ice), m)

    A = Float64(accumulation)
    ρi = Float64(ρ_ice)
    ρ = clamp(Float64(ρ0), 1.0, ρi)
    T = Float64(T_mean)

    p = DensificationCoeffs(mp, A, T)
    max_depth = maximum(depth)

    # Fixed age step, iteration cap and float-epsilon ice test, all as in the
    # Herron & Langway (1980) discretization this reproduces bit-for-bit. Only the
    # density law itself is shared with `steady_state_profile`, which sizes its own
    # step in depth instead (see `_march_step_depth`).
    z_curve = Float64[0.0]
    ρ_curve = Float64[ρ]
    z = 0.0
    iter = 0
    dt = 0.02                                       # age step [yr]
    max_iter = n_age * 50
    @inbounds while z < max_depth && ρ < ρi - D_TOLERANCE && iter < max_iter
        z, ρ = _march_density_step(p, ρ, z, A, ρi, dt, T)
        push!(z_curve, z)
        push!(ρ_curve, ρ)
        iter += 1
    end
    if z_curve[end] < max_depth
        push!(z_curve, max_depth + 1.0)
        push!(ρ_curve, ρi)
    end

    return _interp_curve(ρ_curve, z_curve, depth; lo=Float64(ρ0), hi=ρi)
end
