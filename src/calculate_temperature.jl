"""
    calculate_temperature(temperature, dz, density, water_surface, grain_radius, shortwave_flux, cfs::ClimateForcingStep, mp::ModelParameters, verbose::Bool)

Compute new temperature profile accounting for energy absorption and thermal diffusion.

Solves the 1D heat transfer equation using a finite-volume explicit scheme (Patankar, 1980).
Accounts for:
- Surface energy balance (turbulent fluxes, radiative fluxes)
- Subsurface thermal diffusion
- Shortwave penetration as a source term
- Thermal conductivity updates (Sturm, 1997)

Sub-time steps are determined by Von Neumann stability analysis.

Returns `(temperature, longwave_upward, heat_flux_sensible, heat_flux_latent, heat_flux_basal, evaporation_condensation)`.

# References
- Bougamont, M., et al. (2005). (Surface roughness).
- Foken, T. (2008). Micrometeorology. (Roughness lengths).
- Patankar, S. V. (1980). Numerical Heat Transfer and Fluid Flow.
- Sturm, M., et al. (1997). (Thermal conductivity).
"""
function calculate_temperature(temperature::Vector{Float64}, dz::Vector{Float64},
    density::Vector{Float64}, water_surface::Float64,
    grain_radius::Vector{Float64}, shortwave_flux::Vector{Float64},
    cfs::ClimateForcingStep, mp::ModelParameters, verbose::Bool)

    # Note: temperature is modified in-place by the thermal solver.

    ds = density[1]      # density of top grid cell

    # calculated air density [kg/m3]
    density_air = air_density(cfs.pressure_air, cfs.temperature_air)

    # determine grid point 'center' vector size
    # At least two cells are required: the solve below indexes `Q_sw[m-1]` and `A_face[m-1]`,
    # so a single-cell column is not just degenerate but out of bounds.
    m = length(density)
    if m < 2
        error("column must have at least 2 gridcells for the thermal solve: " *
              "length(density) = $m")
    end

    ## SURFACE ROUGHNESS (Bougamont, 2005)
    if (ds < (mp.density_ice - D_TOLERANCE)) && (water_surface < W_TOLERANCE)
        z0 = Z0_SNOW_DRY   # 0.12 mm for dry snow
    elseif ds >= (mp.density_ice - D_TOLERANCE)
        z0 = Z0_ICE        # 3.2 mm for ice
    else
        z0 = Z0_SNOW_WET   # 1.3 mm for wet snow
    end

    # determine emissivity
    emissivity, emissivity_melt_switch = _emissivity_initialize(grain_radius[1], mp)

    # zT and zQ are percentage of z0
    zT = z0 * mp.surface_roughness_effective_ratio
    zQ = z0 * mp.surface_roughness_effective_ratio

    # minimum wind speed to avoid division by zero in turbulent flux calculations
    min_wind_speed = 0.01

    # Sub-timestep-invariant turbulent-flux quantities (only T_surface changes
    # across sub-steps). Hoisted here and passed to _turbulent_heat_flux so the
    # bulk coefficient, Exner factor, and roughness logs are computed once per
    # timestep rather than once per sub-step.
    thf_wind_speed = max(cfs.wind_speed, min_wind_speed)
    thf_C = VON_KARMAN^2 * thf_wind_speed
    thf_pressure_factor = (100000 / cfs.pressure_air)^0.286
    thf_logM = log(cfs.wind_observation_height / z0)
    thf_logHT = log(cfs.temperature_observation_height / zT)
    thf_logHQ = log(cfs.temperature_observation_height / zQ)
    thf_wind_ratio_sq = (thf_wind_speed / cfs.wind_observation_height)^2

    ## THERMAL CONDUCTIVITY (Sturm, 1997)
    K = thermal_conductivity(temperature, density, mp)

    sfc = _ThermalSurface(density_air, z0, zT, zQ, emissivity, emissivity_melt_switch,
        thf_wind_speed, thf_C, thf_pressure_factor, thf_logM, thf_logHT, thf_logHQ,
        thf_wind_ratio_sq)

    # Scheme selected by dispatch on `mp.thermal_solver`'s type, so the call is resolved at
    # compile time. See `AbstractThermalSolver`.
    longwave_upward_out, heat_flux_sensible_out, heat_flux_latent_out,
        heat_flux_basal_out, evaporation_condensation_out =
        _thermal_solve!(mp.thermal_solver, temperature, dz, density, K, shortwave_flux,
            water_surface, grain_radius, sfc, cfs, mp, verbose)

    return temperature, longwave_upward_out, heat_flux_sensible_out, heat_flux_latent_out, heat_flux_basal_out, evaporation_condensation_out
end

"""
    _ThermalSurface

Sub-timestep-invariant surface quantities, computed once per `calculate_temperature` call and
handed to the thermal solver. `isbits`, so it is stack-allocated and the field reads fold away.

Carries only the surface energy balance's own inputs — roughness lengths, air density, the
hoisted turbulent-flux factors, and the emissivity state. Diagnostics travel separately.
"""
struct _ThermalSurface
    density_air::Float64
    z0::Float64
    zT::Float64
    zQ::Float64
    emissivity::Float64
    emissivity_melt_switch::Bool
    wind_speed::Float64
    C::Float64
    pressure_factor::Float64
    logM::Float64
    logHT::Float64
    logHQ::Float64
    wind_ratio_sq::Float64
end

"""
    _thermal_solve!(::ExplicitThermal, temperature, dz, density, K, shortwave_flux, water_surface, grain_radius, sfc, cfs, mp, verbose)

Advance the column temperature over one forcing timestep with the explicit finite-volume
scheme, sub-stepped to the Von Neumann stability limit (see [`_max_safe_dt`](@ref)).

`temperature` is updated in place. Returns the forcing-step averages
`(longwave_upward, heat_flux_sensible, heat_flux_latent, heat_flux_basal, evaporation_condensation)`;
the first four are W m-2, the last is kg m-2 accumulated over the step.

`water_surface` and `grain_radius` are carried for the verbose diagnostic only.

This is one implementation of the [`AbstractThermalSolver`](@ref) interface; see there for the
invariants every scheme must honour.

# References
- Patankar, S. V. (1980). *Numerical Heat Transfer and Fluid Flow*, Ch. 3-4.
"""
function _thermal_solve!(::ExplicitThermal,
    temperature::Vector{Float64}, dz::Vector{Float64},
    density::Vector{Float64}, K::Vector{Float64}, shortwave_flux::Vector{Float64},
    water_surface::Float64, grain_radius::Vector{Float64},
    sfc::_ThermalSurface, cfs::ClimateForcingStep, mp::ModelParameters, verbose::Bool)

    m = length(density)

    # Unpack the surface state. `emissivity` is mutated by the melt switch below, so it
    # becomes a local; the rest are read-only aliases of the carrier's fields.
    density_air = sfc.density_air
    z0 = sfc.z0
    zT = sfc.zT
    zQ = sfc.zQ
    emissivity = sfc.emissivity
    emissivity_melt_switch = sfc.emissivity_melt_switch
    thf_wind_speed = sfc.wind_speed
    thf_C = sfc.C
    thf_pressure_factor = sfc.pressure_factor
    thf_logM = sfc.logM
    thf_logHT = sfc.logHT
    thf_logHQ = sfc.logHQ
    thf_wind_ratio_sq = sfc.wind_ratio_sq

    # initialize cumulative quantities
    longwave_upward_cumulative = 0.0
    EC_cumulative = 0.0
    lhf_cumulative = 0.0
    shf_cumulative = 0.0
    heat_flux_basal_cumulative = 0.0

    if verbose
        T_bottom = temperature[end]
    end

    ## FIND STABLE dt (allocation-free)
    # The Von Neumann limit uses the pointwise heat capacity. Under `:CuffeyPaterson` this is
    # conservative: c_p < 2102 everywhere below the melting point, so the limit only shrinks.
    max_safe_dt = _max_safe_dt(temperature, dz, density, K, mp)
    dt = _find_dt_divisor(max_safe_dt * 0.8, mp.dt_divisors)

    ## THERMAL DIFFUSION IN ENTHALPY FORM (Patankar 1980, Ch. 3&4)
    #
    # The prognostic is cell enthalpy H [J m-2], not temperature. The previous form folded
    # the heat capacity into normalized coefficients (`Np`, `Nu`, `Nd`, each a conductance
    # divided by `ρ·dz·c/dt`), which is exact only for constant c_p: with c_p(T) = a + bT,
    # a temperature increment ΔT carries h(T+ΔT) − h(T) = ΔT·(a + b(T + ΔT/2)), so freezing
    # c_p at the cell's start-of-substep value injects (b/2)ΔT² J kg-1 of fictitious energy
    # every sub-step. The surface cell can swing tens of K per sub-step, which is joules per
    # kilogram, not round-off.
    #
    # Diffusion is applied as one flux per interior face, added to the cell below the face
    # and subtracted from the cell above it, so the pairwise cancellation is exact to the
    # last bit and the column total is conserved independently of c_p.
    #
    # The mass is carried as its reciprocal: enthalpy → temperature is on the inner loop
    # twice per cell per sub-step, and a division there costs more than the whole flux
    # evaluation. `H` is *specific* enthalpy scaled by mass, so `H[i] * M_inv[i]` is the
    # J kg-1 that the inverse map wants. The `:constant` inverse is itself a division by
    # `c`, so that constant is folded into the same reciprocal and
    # `temperature_from_scaled_enthalpy` becomes the identity — no division on the inner loop.
    M_inv = Vector{Float64}(undef, m)    # 1 / (cell mass * enthalpy scale)
    H = Vector{Float64}(undef, m)        # cell enthalpy [J m-2] — carried across sub-steps
    Q_sw = Vector{Float64}(undef, m)     # SW energy absorbed per sub-step [J m-2]
    A_face = Vector{Float64}(undef, m - 1)  # face conductance * dt [J m-2 K-1]

    h_scale = enthalpy_temperature_scale(mp)
    @inbounds for i in 1:m
        M_cell_i = density[i] * dz[i]
        M_inv[i] = 1.0 / (M_cell_i * h_scale)
        H[i] = M_cell_i * specific_enthalpy(mp, temperature[i])
        Q_sw[i] = shortwave_flux[i] * dt
    end
    @inbounds for i in 1:m-1
        # `dt` is fixed across sub-steps, so it is folded in here rather than multiplied per
        # cell per sub-step. A_face is therefore joules per kelvin of face temperature drop.
        A_face[i] = dt / (dz[i+1] / (2 * K[i+1]) + dz[i] / (2 * K[i]))
    end

    # Ensure no SW reaches the bottom cell: its share is absorbed by the cell above. This
    # keeps the Dirichlet bottom cell's enthalpy untouched for the whole solve.
    @inbounds begin
        Q_sw[m-1] += shortwave_flux[m] * dt
        Q_sw[m] = 0.0
    end

    # energy supplied by downward longwave radiation to the top grid cell [J]
    longwave_downward = cfs.longwave_downward * dt

    # Total SW absorbed per sub-step, for the verbose budget. Loop-invariant.
    sw_total = verbose ? sum(shortwave_flux) * dt : 0.0

    ## CALCULATE ENERGY SOURCES AND DIFFUSION FOR EVERY TIME STEP [dt]
    n_steps = round(Int, cfs.dt / dt)

    # Local variables for loop outputs
    local longwave_upward::Float64
    local heat_flux_sensible::Float64
    local heat_flux_latent::Float64
    local evaporation_condensation::Float64

    for _ in 1:n_steps
        # Store initial column enthalpy for the energy conservation check
        if verbose
            E_initial = 0.0
            @inbounds for i in 1:m
                E_initial += H[i]
            end
        end

        # calculate temperature of snow surface
        T_surface = min(273.15, temperature[1])

        # TURBULENT HEAT FLUX (invariant quantities hoisted above the loop)
        heat_flux_sensible, heat_flux_latent, latent_heat = _turbulent_heat_flux(
            T_surface, density_air, z0, zT, zQ, cfs,
            thf_wind_speed, thf_C, thf_pressure_factor, thf_logM, thf_logHT, thf_logHQ,
            thf_wind_ratio_sq)

        lhf_cumulative += heat_flux_latent * dt
        shf_cumulative += heat_flux_sensible * dt

        # mass loss (-)/accretion(+) due to evaporation/condensation [kg]
        evaporation_condensation = heat_flux_latent / latent_heat * dt

        # energy supplied by turbulent fluxes [J]
        thf = (heat_flux_sensible + heat_flux_latent) * dt

        # upward longwave radiation
        T2 = T_surface * T_surface
        longwave_upward = -(SB * T2 * T2 * emissivity) * dt
        longwave_upward_cumulative += -longwave_upward

        # net energy delivered to the surface cell this sub-step [J]
        Q_surface = (longwave_downward + longwave_upward) + thf

        # Single fused pass: absorb this cell's source term, map it to its pre-diffusion
        # temperature, then close out the *previous* cell, whose two faces are now both known.
        #
        # Diffusion is conservative and explicit. F_i is the energy crossing face i (between
        # cells i and i+1) into cell i, so cell i's net is F_i − F_{i-1}. Applying one flux per
        # face with opposite signs makes the pairwise cancellation exact to the last bit, so
        # the column total is conserved independently of c_p. Cell m is the Dirichlet
        # reservoir: face m-1 supplies F_{m-1} = heat_flux_basal to cell m-1 and nothing is
        # taken from cell m, whose Q_sw is zero — `temperature[m]` is therefore left exactly
        # as it came in rather than round-tripped through h⁻¹, which matters because the
        # check below compares it with an exact `!=`.
        #
        # The pre-diffusion temperatures are carried in registers (`T_pre_prev`, `T_pre_i`) —
        # each is read only at the two faces that bound its cell, so trailing the enthalpy
        # update one cell behind the source pass keeps it alive exactly long enough. Two
        # `M_inv` reads and one `A_face` read per cell, one pass, no second sweep.
        @inbounds begin
            H[1] += Q_sw[1] + Q_surface
            T_pre_prev = temperature_from_scaled_enthalpy(mp, H[1] * M_inv[1])
            F_prev = 0.0
            for i in 1:m-2
                H[i+1] += Q_sw[i+1]
                T_pre_i = temperature_from_scaled_enthalpy(mp, H[i+1] * M_inv[i+1])
                F = A_face[i] * (T_pre_i - T_pre_prev)
                H[i] += F - F_prev
                F_prev = F
                T_pre_prev = T_pre_i
            end

            # Face m-1 draws on the fixed bottom cell, so its temperature is used directly.
            heat_flux_basal = A_face[m-1] * (temperature[m] - T_pre_prev)
            heat_flux_basal_cumulative += heat_flux_basal
            H[m-1] += heat_flux_basal - F_prev

            # `temperature[1]` is the only cell read again before the next sub-step (surface
            # fluxes, the emissivity switch, and the verbose diagnostic), so it is the only one
            # reconstructed here. The interior is recovered once after the sub-step loop
            # instead of being rewritten `n_steps` times and read by nobody.
            temperature[1] = temperature_from_scaled_enthalpy(mp, H[1] * M_inv[1])
        end

        # calculate cumulative evaporation (+)/condensation(-)
        EC_cumulative += evaporation_condensation

        # Emissivity melt switch check. Offset by `T_MELT_SWITCH_TOLERANCE` (1e-4 K), not
        # the 1e-10 `T_TOLERANCE` used for branch boundaries elsewhere.
        if emissivity_melt_switch
            if temperature[1] < (CtoK - T_MELT_SWITCH_TOLERANCE)
                emissivity = mp.emissivity
            else
                emissivity = mp.emissivity_grain_radius_large
            end
        end

        # CHECK FOR ENERGY CONSERVATION
        if verbose
            E_final = 0.0
            @inbounds for i in 1:m
                E_final += H[i]
            end
            E_used = E_final - E_initial
            E_supplied = sw_total + longwave_downward + longwave_upward + thf + heat_flux_basal
            E_delta = E_used - E_supplied

            E_tolerance = energy_tolerance(E_initial)
            if (abs(E_delta) > E_tolerance) || isnan(E_delta)
                @error "inputs" temperature[1] water_surface grain_radius[1] sum(shortwave_flux) cfs.longwave_downward cfs.temperature_air cfs.wind_speed cfs.vapor_pressure cfs.pressure_air
                @error "internals" sw_total longwave_downward longwave_upward thf heat_flux_basal
                error("energy not conserved in thermodynamics equations: supplied = $(E_supplied) J, used = $(E_used) J")
            end

            if T_bottom != temperature[end]
                error("temperature of bottom grid cell changed inside of thermal function: original = $(T_bottom) K, updated = $(temperature[end]) K")
            end
        end
    end

    # Recover the interior temperatures from the prognostic enthalpy. Cell 1 is already current
    # (the sub-step loop needs it) and cell m is the untouched Dirichlet reservoir.
    @inbounds for i in 2:m-1
        temperature[i] = temperature_from_scaled_enthalpy(mp, H[i] * M_inv[i])
    end

    heat_flux_latent_out = lhf_cumulative / cfs.dt    # J -> W/m2
    heat_flux_sensible_out = shf_cumulative / cfs.dt  # J -> W/m2
    longwave_upward_out = longwave_upward_cumulative / cfs.dt  # J -> W/m2
    heat_flux_basal_out = heat_flux_basal_cumulative / cfs.dt  # J -> W/m2
    evaporation_condensation_out = EC_cumulative

    return longwave_upward_out, heat_flux_sensible_out, heat_flux_latent_out, heat_flux_basal_out, evaporation_condensation_out
end

"""
    _surface_energy_balance(T1, emissivity, sfc, cfs)

Surface energy balance at surface-cell temperature `T1` [K]. Returns
`(longwave_upward, heat_flux_sensible, heat_flux_latent, latent_heat, T_surface)` with the three
fluxes in W m-2 (signed as energy *into* the column, so `longwave_upward` is negative) and
`T_surface = min(CtoK, T1)`.

The clamp is what makes the balance non-smooth at the melting point: above it the surface stops
responding to `T1` entirely, so every flux — and every derivative — goes flat. Shared by the
implicit solver's Newton iteration and its final applied-flux evaluation so the two cannot drift.
"""
@inline function _surface_energy_balance(T1::Float64, emissivity::Float64,
    sfc::_ThermalSurface, cfs::ClimateForcingStep)

    T_surface = min(CtoK, T1)
    heat_flux_sensible, heat_flux_latent, latent_heat = _turbulent_heat_flux(
        T_surface, sfc.density_air, sfc.z0, sfc.zT, sfc.zQ, cfs,
        sfc.wind_speed, sfc.C, sfc.pressure_factor, sfc.logM, sfc.logHT, sfc.logHQ,
        sfc.wind_ratio_sq)
    T2 = T_surface * T_surface
    longwave_upward = -(SB * T2 * T2 * emissivity)
    return longwave_upward, heat_flux_sensible, heat_flux_latent, latent_heat, T_surface
end

"""
    _surface_energy_balance_slope(T1, emissivity, sfc, cfs, longwave_upward, heat_flux_sensible, heat_flux_latent)

Slope `dQ_surface/dT1` [W m-2 K-1] of [`_surface_energy_balance`](@ref) at `T1`, given the fluxes
already evaluated there. Returns a value `≤ 0`.

The longwave term is analytic, `-4σεT³`. The two turbulent terms come from a single extra
`_turbulent_heat_flux` call, differenced with [`THERMAL_BC_DERIVATIVE_STEP`](@ref): the
Beljaars-Holtslag stability branches make hand-differentiation a maintenance liability for no
benefit, because this slope sets only the Newton convergence *rate* and never the converged
answer (see [`_thermal_solve!`](@ref)).

Three deliberate crudenesses, each licensed by that:

- **Above the melting point the slope is exactly zero.** The `min(CtoK, T1)` clamp makes `Q`
  independent of `T1` there, so this is not an approximation.
- **The difference is taken away from the clamp.** A forward step that would cross `CtoK` is
  replaced by a backward one; `Q` is smooth below the melting point, so either is equally good.
- **Each component is clamped at `≤ 0`.** A positive slope would weaken the Newton diagonal
  instead of strengthening it, and it can arise spuriously — the `LV`/`LS` latent-heat switch at
  `turbulent_heat_flux.jl` is a genuine jump in `heat_flux_latent` right at the melting point,
  where wet columns sit. Clamping costs convergence rate, never correctness.

!!! note "An analytic turbulent slope was tried and reverted"
    The extra `_turbulent_heat_flux` call is the dominant cost of the implicit path, so an
    analytic replacement was measured. `dQ_shf/dT` and `dQ_lhf/dT` are both algebraic in
    quantities the flux evaluation already computes — the transfer-coefficient products and the
    saturation vapour pressure — *if* the coefficients themselves are held fixed. They cannot
    be: `coefM`, `coefHT`, `coefHQ` depend on `T_surface` through the bulk Richardson number,
    and the term dropped by freezing them substantially cancels the rest. Measured over 675
    sampled surface states, the frozen-coefficient slope overshot the true derivative by up to
    10x, with 9% of states outside a factor of two. Overshoot under-relaxes Newton: the
    iteration count more than doubled, hit `THERMAL_IMPLICIT_MAX_ITERATIONS`, and whole-model
    runtime *rose* 43% (15.2 s to 21.8 s) while the unconverged output drifted 15 K. Correct in
    principle — `Λ` really does affect only the rate — but the rate is exactly what it cost.
    Differentiating the stability branches too would be sound, and is the open option here.
"""
@inline function _surface_energy_balance_slope(T1::Float64, emissivity::Float64,
    sfc::_ThermalSurface, cfs::ClimateForcingStep, longwave_upward::Float64,
    heat_flux_sensible::Float64, heat_flux_latent::Float64)

    # Clamped: the surface balance no longer sees T1 at all.
    T1 >= CtoK && return 0.0

    T_surface = min(CtoK, T1)
    slope_longwave = -4 * SB * emissivity * T_surface * T_surface * T_surface

    # Step away from the melting-point clamp so the difference stays on the smooth branch.
    step = (T1 + THERMAL_BC_DERIVATIVE_STEP <= CtoK) ? THERMAL_BC_DERIVATIVE_STEP :
           -THERMAL_BC_DERIVATIVE_STEP
    shf_p, lhf_p, _ = _turbulent_heat_flux(
        T1 + step, sfc.density_air, sfc.z0, sfc.zT, sfc.zQ, cfs,
        sfc.wind_speed, sfc.C, sfc.pressure_factor, sfc.logM, sfc.logHT, sfc.logHQ,
        sfc.wind_ratio_sq)
    slope_turbulent = ((shf_p - heat_flux_sensible) + (lhf_p - heat_flux_latent)) / step

    return min(0.0, slope_longwave) + min(0.0, slope_turbulent)
end

"""
    _thermal_solve!(::ImplicitThermal, temperature, dz, density, K, shortwave_flux, water_surface, grain_radius, sfc, cfs, mp, verbose)

Advance the column temperature over one forcing timestep with a backward-Euler scheme on a
tridiagonal system, solved by the Thomas algorithm. Unconditionally stable, so it takes no
stability sub-steps: [`_max_safe_dt`](@ref), [`_find_dt_divisor`](@ref), `mp.dt_divisors` and
[`DT_MIN_WARN`](@ref) are all unused on this path, and a thin refrozen lens costs nothing.

`temperature` is updated in place; the return is the [`AbstractThermalSolver`](@ref) 5-tuple.
`water_surface` and `grain_radius` are carried for the verbose diagnostic only.

## The system

Unknowns are cells `1 … m-1`. Cell `m` is the Dirichlet reservoir and is never an unknown, so its
temperature is returned bit-unchanged by construction; it enters row `m-1` as the known term
`A_face[m-1]·temperature[m]` on the right-hand side. Its shortwave share is absorbed by cell
`m-1`, exactly as on the explicit path.

Row `i` is the backward-Euler enthalpy balance
`M_i·(h(T_i) − h(T_i^old)) = Q_sw,i + F_i − F_{i-1}` with the diffusive fluxes evaluated at the
*new* temperatures, giving diagonal `M_i·c_i + A_face[i-1] + A_face[i]` and off-diagonals
`−A_face[i]`. The face conductances `A_face` are the same harmonic means the explicit path uses,
so the two schemes discretize identical conductances and differ only in time integration.

`c_i` is the [`chord_heat_capacity`](@ref) between the old and new temperatures — an exact
identity for `h = aT + (b/2)T²`, not a linearization. This is what lets enthalpy stay the
conserved quantity while the system solved is linear in temperature. It also makes the system
mildly nonlinear under `:CuffeyPaterson`, resolved by the same iteration as the surface row; under
the default `:constant` it is a true constant and the iteration only sees the surface.

`M_i·c_i > 0` always, so the matrix is an irreducibly diagonally dominant M-matrix: non-singular,
**no pivoting required**, and `A⁻¹ ≥ 0`, which is the discrete maximum principle — the solve
cannot manufacture a new extremum.

## The surface row: Newton, not lagged Picard

The surface flux `Q(T_1)` is strongly nonlinear and must be taken implicitly. Lagging it — a
Picard iteration on the previous iterate's flux — *diverges* here: the fixed-point gain
`|dQ/dT_1|·Δt/(M_1 c_1)` exceeds 1 for a centimetre-scale surface cell at any Δt of interest.
(That is precisely why the explicit path's lagged sub-steps are stable: their tiny `dt` drives
the same gain far below 1.)

So the flux is linearized into the diagonal, `Q(T_1) ≈ Q_k + Λ_k(T_1 − T_{1,k})` with
`Λ_k ≤ 0` from [`_surface_energy_balance_slope`](@ref), which *strengthens* the already dominant
diagonal. Because `Q_k` is the true nonlinear flux at the iterate, the two `Λ` terms cancel at
convergence and the residual satisfied is the true nonlinear surface balance. **`Λ` therefore
affects only the convergence rate, never the answer** — the licence for its finite-difference
turbulent terms and its clamp.

The two discontinuities stay *outside* the iteration: the emissivity melt switch is evaluated
once per sub-step, as on the explicit path, and the `LV`/`LS` latent-heat switch is absorbed by
the `Λ` clamp. Inside, either would chatter and prevent convergence.

## Conservation

A Thomas sweep leaves a round-off residual, so the solved field is used as a *predictor of the
implicit face temperatures* rather than written to `temperature` directly. The enthalpy update is
then applied as one flux per face, `+F` to cell `i` and `−F` to cell `i+1` — structurally the same
pass the explicit path uses. The pairwise cancellation is exact to the last bit, so the column
total conserves independently of the solve residual, the Newton tolerance, and `c_p(T)`, and the
verbose budget check is the explicit path's with no tolerance moved.

The surface flux applied in that pass is the *true* nonlinear flux at the converged iterate, not
its linearization, and the returned averages are those same applied values — so `gemb_core`'s
energy budget and the `evaporation_condensation` mass budget close by construction.

## Sub-stepping

Sub-steps here buy **accuracy, not stability**: backward Euler is first-order and strongly
damping, and the emissivity melt switch only resolves between sub-steps. The count follows
[`THERMAL_IMPLICIT_DT_TARGET`](@ref) and is independent of the cell count, of `dz`, and of the
stiffest cell in the column — the property the explicit path lacks.

## Static condensation: Newton on a scalar, not on the column

The nonlinearity is confined to row 1 — the surface energy balance. Rows `2 … n` are linear in
the iterate under the default `:constant` heat capacity, so re-eliminating them once per Newton
iteration is wasted work. Instead the interior is condensed **bottom-up, once per sub-step**,
each row reduced to an affine function of the cell above it (`T_i = rhs[i] + sup[i]·T_{i-1}`).
Newton then iterates on a single scalar equation in `T_1`, at O(1) per iteration, and one forward
substitution propagates the answer back down the column.

Under `:CuffeyPaterson` the chord capacities depend on the iterate, so the condensation is
repeated in an outer loop until they settle; under `:constant` the outer loop provably runs once.

This is what the exponential-integrator idea in the design notes reduces to in cheap form: the
operator is constant across all Newton iterations of a sub-step, so factorize it once and reuse
it. Measured: **15.20 s to 5.70 s (2.67x)** on the benchmark below, with whole-model output
agreeing to 1.2e-6 K and annual melt to 3.2e-8 kg m-2 — the same converged answer to round-off,
as the `Λ`-independence property requires.

## Performance: still slower than explicit, and why

Measured on a year of 3-hourly synthetic forcing over a 264-cell column, whole-model runtime:

| Solver | `DT_TARGET` | Runtime |
|---|---|---|
| `ExplicitThermal` | — | 2.44 s |
| `ImplicitThermal` | 1800 s (6 sub-steps) | 3.59 s |
| `ImplicitThermal` | 900 s (12 sub-steps, default) | 5.75 s |
| `ImplicitThermal` | 450 s (24 sub-steps) | 9.62 s |

So 1.5x at the coarsest usable accuracy and 2.4x at the calibrated default. Unconditional
stability is the reason this scheme exists; throughput is not, and on a column with no stiff cell
it does not beat the explicit path.

What remains is the Newton iteration itself: **4.19 iterations per sub-step solve**, measured over
1.12M solves. Capping the iteration at 1 takes the run to 4.98 s, so the iteration is most of the
remaining cost, and it is irreducible — surface-only and whole-profile convergence criteria give
identical counts (4.190 vs 4.190), i.e. cell 1 is always the last cell to converge.

Two things are measured *not* to be the lever, recorded so they are not retried:

- **The turbulent-flux calls.** Direct timing puts `_surface_energy_balance` at 65 ns and
  `_surface_energy_balance_slope` at 69 ns, so all flux evaluations together are 0.70 s of the
  pre-condensation 15.20 s. The earlier attribution of the cost to doubled flux calls was wrong;
  the cost was the O(n) sweep per iteration, which is what condensation removed.
- **The convergence tolerance.** Relaxing `THERMAL_IMPLICIT_T_TOLERANCE` from 1e-10 to 1e-3 —
  seven orders of magnitude, far looser than the model's other branch tolerances — cut the
  pre-condensation runtime only from 15.2 s to 10.8 s.

Where this scheme wins is the case the explicit path cannot bound: its cost is independent of the
stiffest cell. A single 1e-4 m refrozen lens drops the explicit stability limit from 903 s to
3.98 s on this column (measured, `test_calculate_temperature.jl`) — a 227x sub-step increase for
one thin cell — while the implicit count does not move. That is the reason to keep it available.

# References
- Patankar, S. V. (1980). *Numerical Heat Transfer and Fluid Flow*, Ch. 4.
- Versteeg, H. K. & Malalasekera, W. (2007). *An Introduction to Computational Fluid Dynamics*,
  Ch. 8.
- Thomas, L. H. (1949). *Elliptic Problems in Linear Difference Equations over a Network*.
- Beljaars, A. C. M. & Holtslag, A. A. M. (1991). Flux parameterization over land surfaces.
  *J. Appl. Meteorol.* 30, 327-341.
"""
function _thermal_solve!(::ImplicitThermal,
    temperature::Vector{Float64}, dz::Vector{Float64},
    density::Vector{Float64}, K::Vector{Float64}, shortwave_flux::Vector{Float64},
    water_surface::Float64, grain_radius::Vector{Float64},
    sfc::_ThermalSurface, cfs::ClimateForcingStep, mp::ModelParameters, verbose::Bool)

    m = length(density)
    n = m - 1                     # unknowns: cells 1 … m-1; cell m is the Dirichlet reservoir
    T_bottom = temperature[m]

    emissivity = sfc.emissivity

    # Accuracy sub-stepping only — see the docstring. No stability limit is consulted.
    n_sub = clamp(ceil(Int, cfs.dt / THERMAL_IMPLICIT_DT_TARGET), 1,
        THERMAL_IMPLICIT_SUBSTEPS_MAX)
    dt = cfs.dt / n_sub

    M = Vector{Float64}(undef, n)         # cell mass [kg m-2]
    H = Vector{Float64}(undef, n)         # cell enthalpy [J m-2] — the prognostic
    Q_sw = Vector{Float64}(undef, n)      # SW absorbed per sub-step [J m-2]
    A_face = Vector{Float64}(undef, n)    # face conductance * dt [J m-2 K-1]
    rhs = Vector{Float64}(undef, n)       # right-hand side, then the Thomas solution
    sup = Vector{Float64}(undef, n)       # eliminated superdiagonal
    T_it = Vector{Float64}(undef, n)      # current Newton iterate

    @inbounds for i in 1:n
        M[i] = density[i] * dz[i]
        H[i] = M[i] * specific_enthalpy(mp, temperature[i])
        Q_sw[i] = shortwave_flux[i] * dt
        # `dt` is fixed across sub-steps, so it is folded into the conductance once here.
        A_face[i] = dt / (dz[i+1] / (2 * K[i+1]) + dz[i] / (2 * K[i]))
    end
    # The bottom cell absorbs no shortwave; its share warms the cell above it, which keeps the
    # Dirichlet reservoir's enthalpy out of the solve entirely.
    @inbounds Q_sw[n] += shortwave_flux[m] * dt

    longwave_downward = cfs.longwave_downward * dt
    sw_total = verbose ? sum(shortwave_flux) * dt : 0.0

    # Fixed reference for the relative part of the energy tolerance, on the same scale as the
    # explicit path's (which sums the untouched bottom cell into its running total).
    E_reference = 0.0
    if verbose
        @inbounds for i in 1:n
            E_reference += H[i]
        end
        E_reference += density[m] * dz[m] * specific_enthalpy(mp, T_bottom)
    end

    longwave_upward_cumulative = 0.0
    shf_cumulative = 0.0
    lhf_cumulative = 0.0
    heat_flux_basal_cumulative = 0.0
    EC_cumulative = 0.0

    for _ in 1:n_sub
        if verbose
            E_initial = 0.0
            @inbounds for i in 1:n
                E_initial += H[i]
            end
        end

        # Start the iteration from the old profile.
        @inbounds for i in 1:n
            T_it[i] = temperature[i]
        end

        # Under a constant heat capacity the interior rows are exactly linear, so a single
        # elimination is exact and this outer pass never repeats. `:CuffeyPaterson` makes the
        # chord capacities functions of the iterate, so the elimination is re-run until they
        # settle.
        n_outer = (mp.heat_capacity_method === :constant) ? 1 :
                  THERMAL_IMPLICIT_MAX_ITERATIONS

        for _ in 1:n_outer
            # Eliminate the interior rows bottom-up, expressing each as an affine function of
            # the cell above it: `T_i = rhs[i] + sup[i]*T_{i-1}`. Row 1 is the only nonlinear
            # row, and it is the only one Newton then has to touch.
            @inbounds if n > 1
                MC = M[n] * chord_heat_capacity(mp, temperature[n], T_it[n])
                d = MC + A_face[n-1] + A_face[n]
                sup[n] = A_face[n-1] / d
                rhs[n] = (MC * temperature[n] + Q_sw[n] + A_face[n] * T_bottom) / d

                for i in n-1:-1:2
                    MC = M[i] * chord_heat_capacity(mp, temperature[i], T_it[i])
                    d = MC + A_face[i-1] + A_face[i] - A_face[i] * sup[i+1]
                    sup[i] = A_face[i-1] / d
                    rhs[i] = (MC * temperature[i] + Q_sw[i] + A_face[i] * rhs[i+1]) / d
                end
            end

            # Newton on the scalar surface row. Each iteration is O(1): the interior is already
            # condensed into `rhs[2]`/`sup[2]`, so no sweep of the column is repeated here.
            T1 = @inbounds T_it[1]
            for _ in 1:THERMAL_IMPLICIT_MAX_ITERATIONS
                longwave_upward, heat_flux_sensible, heat_flux_latent, _ =
                    _surface_energy_balance(T1, emissivity, sfc, cfs)
                slope = _surface_energy_balance_slope(T1, emissivity, sfc, cfs,
                    longwave_upward, heat_flux_sensible, heat_flux_latent)

                # Joules this sub-step, and the slope in J m-2 K-1.
                Q_surface = longwave_downward +
                            (longwave_upward + heat_flux_sensible + heat_flux_latent) * dt
                Lambda = slope * dt

                @inbounds begin
                    MC = M[1] * chord_heat_capacity(mp, temperature[1], T1)
                    # No face above cell 1; the surface slope enters as `-Λ ≥ 0`.
                    d = MC + A_face[1] - Lambda
                    r = MC * temperature[1] + Q_sw[1] + (Q_surface - Lambda * T1)
                    if n > 1
                        # Substitute the condensed row 2.
                        d -= A_face[1] * sup[2]
                        r += A_face[1] * rhs[2]
                    else
                        # Cell 1 is also the last unknown, so the reservoir term lands here.
                        r += A_face[1] * T_bottom
                    end
                    T1_next = r / d
                end

                converged = abs(T1_next - T1) <= THERMAL_IMPLICIT_T_TOLERANCE
                T1 = T1_next
                converged && break
            end

            # Forward-substitute back down the column.
            @inbounds begin
                dT_max = abs(T1 - T_it[1])
                T_it[1] = T1
                for i in 2:n
                    Ti = rhs[i] + sup[i] * T_it[i-1]
                    dT_max = max(dT_max, abs(Ti - T_it[i]))
                    T_it[i] = Ti
                end
            end
            dT_max <= THERMAL_IMPLICIT_T_TOLERANCE && break
        end

        # Applied fluxes: the *true* nonlinear surface balance at the converged iterate, so the
        # returned averages are exactly what the column received.
        longwave_upward, heat_flux_sensible, heat_flux_latent, latent_heat, _ =
            _surface_energy_balance(T_it[1], emissivity, sfc, cfs)
        longwave_upward_step = longwave_upward * dt
        thf = (heat_flux_sensible + heat_flux_latent) * dt
        Q_surface = (longwave_downward + longwave_upward_step) + thf

        longwave_upward_cumulative += -longwave_upward_step
        shf_cumulative += heat_flux_sensible * dt
        lhf_cumulative += heat_flux_latent * dt
        evaporation_condensation = heat_flux_latent / latent_heat * dt
        EC_cumulative += evaporation_condensation

        # Reconciliation: one flux per face, `+F` below and `−F` above, so the column total
        # conserves to the last bit regardless of the solve residual.
        local heat_flux_basal::Float64
        @inbounds begin
            H[1] += Q_surface
            F_prev = 0.0
            for i in 1:n-1
                F = A_face[i] * (T_it[i+1] - T_it[i])
                H[i] += Q_sw[i] + F - F_prev
                F_prev = F
            end
            heat_flux_basal = A_face[n] * (T_bottom - T_it[n])
            H[n] += Q_sw[n] + heat_flux_basal - F_prev
        end
        heat_flux_basal_cumulative += heat_flux_basal

        @inbounds for i in 1:n
            temperature[i] = temperature_from_specific_enthalpy(mp, H[i] / M[i])
        end

        # Emissivity melt switch, held outside the Newton iteration so it cannot chatter.
        if sfc.emissivity_melt_switch
            if temperature[1] < (CtoK - T_MELT_SWITCH_TOLERANCE)
                emissivity = mp.emissivity
            else
                emissivity = mp.emissivity_grain_radius_large
            end
        end

        if verbose
            E_final = 0.0
            @inbounds for i in 1:n
                E_final += H[i]
            end
            E_used = E_final - E_initial
            E_supplied = sw_total + longwave_downward + longwave_upward_step + thf +
                         heat_flux_basal
            E_delta = E_used - E_supplied

            if (abs(E_delta) > energy_tolerance(E_reference)) || isnan(E_delta)
                @error "inputs" temperature[1] water_surface grain_radius[1] sum(shortwave_flux) cfs.longwave_downward cfs.temperature_air cfs.wind_speed cfs.vapor_pressure cfs.pressure_air
                @error "internals" sw_total longwave_downward longwave_upward_step thf heat_flux_basal
                error("energy not conserved in thermodynamics equations: supplied = $(E_supplied) J, used = $(E_used) J")
            end

            if T_bottom != temperature[m]
                error("temperature of bottom grid cell changed inside of thermal function: original = $(T_bottom) K, updated = $(temperature[m]) K")
            end
        end
    end

    heat_flux_latent_out = lhf_cumulative / cfs.dt              # J -> W/m2
    heat_flux_sensible_out = shf_cumulative / cfs.dt            # J -> W/m2
    longwave_upward_out = longwave_upward_cumulative / cfs.dt   # J -> W/m2
    heat_flux_basal_out = heat_flux_basal_cumulative / cfs.dt   # J -> W/m2
    evaporation_condensation_out = EC_cumulative

    return longwave_upward_out, heat_flux_sensible_out, heat_flux_latent_out, heat_flux_basal_out, evaporation_condensation_out
end

"""
    _emissivity_initialize(grain_radius_surface, mp::ModelParameters)

Initialize emissivity based on surface grain radius and model parameters.
Returns `(emissivity, emissivity_melt_switch)`.
"""
function _emissivity_initialize(grain_radius_surface::Float64, mp::ModelParameters)

    if mp.emissivity_method == :uniform
        emissivity = mp.emissivity
        emissivity_melt_switch = false
    elseif mp.emissivity_method == :grain_radius_threshold
        if grain_radius_surface <= (mp.emissivity_grain_radius_threshold + GDN_TOLERANCE)
            emissivity = mp.emissivity
        else
            emissivity = mp.emissivity_grain_radius_large
        end
        emissivity_melt_switch = false
    elseif mp.emissivity_method == :grain_radius_w_threshold
        if grain_radius_surface <= (mp.emissivity_grain_radius_threshold + GDN_TOLERANCE)
            emissivity = mp.emissivity
            emissivity_melt_switch = true
        else
            emissivity = mp.emissivity_grain_radius_large
            emissivity_melt_switch = false
        end
    else
        error("Unrecognized emissivity_method: $(mp.emissivity_method)")
    end

    return emissivity, emissivity_melt_switch
end

"""
    _max_safe_dt(temperature, dz, density, K, mp)

Largest explicit sub-step [s] that keeps the thermal solve stable. Allocation-free.

The scheme in [`calculate_temperature`](@ref) updates cell `i` as
`H[i] += F_i - F_{i-1}` with `F_i = G_i·dt·(T_{i+1} - T_i)`, so the coefficient of `T_i` in
the equivalent temperature update is `1 - dt·(G_i + G_{i-1})/(ρᵢcᵢdzᵢ)`. Requiring it to stay
non-negative — the positivity condition that is Von Neumann stability for this stencil — gives

    dt ≤ ρᵢ cᵢ dzᵢ / (Gᵢ + Gᵢ₋₁),    Gᵢ = 1 / (dz[i+1]/(2K[i+1]) + dz[i]/(2K[i]))

evaluated here by reusing the very harmonic-mean face conductances the solve uses, so the
limit and the scheme cannot drift apart.

This replaced the textbook uniform-grid form `0.5·ρᵢcᵢdzᵢ²/Kᵢ`, which substitutes `2Kᵢ/dzᵢ`
for `Gᵢ + Gᵢ₋₁`. The two agree exactly when `dz` and `K` are uniform, but GEMB's grid is
graded, and the substitution is not conservative there: since `Gᵢ ≤ 2Kᵢ/dzᵢ` with equality
only in the limit of a vanishing neighbour resistance, `Gᵢ + Gᵢ₋₁` can reach `4Kᵢ/dzᵢ`, so
the old form could overestimate the true limit by up to a factor of two. Measured per cell on
GEMB's own column the ratio spanned 0.66 to 1.82 — an error in both directions, the low end
of which the 0.8 safety factor does not cover. It is also *looser* at the cell that binds, so
the correct limit is both sound and cheaper.

Cell `m` is excluded: it is the Dirichlet reservoir, its enthalpy is never updated
(`Q_sw[m] = 0`, no flux is drawn from it), so it imposes no stability constraint. The old
form included it.

# References
- Patankar, S. V. (1980). *Numerical Heat Transfer and Fluid Flow*, Ch. 4 (positivity of the
  explicit finite-volume coefficients).
"""
function _max_safe_dt(temperature::Vector{Float64}, dz::Vector{Float64},
    density::Vector{Float64}, K::Vector{Float64}, mp::ModelParameters)

    m = length(dz)
    max_safe_dt = Inf

    # G_prev is the conductance of the face above cell i; zero above the surface cell.
    G_prev = 0.0
    @inbounds for i in 1:m-1
        G = 1.0 / (dz[i+1] / (2 * K[i+1]) + dz[i] / (2 * K[i]))
        sl = density[i] * heat_capacity(mp, temperature[i]) * dz[i] / (G + G_prev)
        max_safe_dt = min(max_safe_dt, sl)
        G_prev = G
    end

    return max_safe_dt
end

"""
    _find_dt_divisor(dt_target, dt_divisors)

Find the largest dt_divisor that is <= dt_target. Allocation-free.
"""
function _find_dt_divisor(dt_target::Float64, dt_divisors::Vector{Float64})
    if dt_target < DT_MIN_WARN
        @warn "Timestep is extremely small ($dt_target). Check for near-zero dz layers."
    end

    dt = dt_divisors[1]  # fallback to smallest
    @inbounds for i in eachindex(dt_divisors)
        if dt_divisors[i] <= dt_target
            dt = dt_divisors[i]
        else
            break
        end
    end
    return dt
end
