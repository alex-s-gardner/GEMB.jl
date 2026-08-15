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

Returns `(temperature, longwave_upward, heat_flux_sensible, heat_flux_latent, ghf, evaporation_condensation)`.

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

    # initialize cumulative quantities
    longwave_upward_cumulative = 0.0
    EC_cumulative = 0.0
    lhf_cumulative = 0.0
    shf_cumulative = 0.0
    ghf_cumulative = 0.0

    if verbose
        T_bottom = temperature[end]
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

    ## THERMAL CONDUCTIVITY (Sturm, 1997)
    K = thermal_conductivity(temperature, density, mp)

    ## FIND STABLE dt (allocation-free)
    # The Von Neumann limit uses the pointwise heat capacity. Under `:CuffeyPaterson` this is
    # conservative: c_p < 2102 everywhere below the melting point, so the limit only shrinks.
    max_safe_dt = Inf
    @inbounds for i in eachindex(dz)
        sl = 0.5 * density[i] * heat_capacity(mp, temperature[i]) * dz[i]^2 / K[i]
        max_safe_dt = min(max_safe_dt, sl)
    end
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
            thf_wind_speed, thf_C, thf_pressure_factor, thf_logM, thf_logHT, thf_logHQ)

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
        # reservoir: face m-1 supplies F_{m-1} = ghf to cell m-1 and nothing is taken from
        # cell m, whose Q_sw is zero — `temperature[m]` is therefore left exactly as it came
        # in rather than round-tripped through h⁻¹, which matters because the check below
        # compares it with an exact `!=`.
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
            ghf = A_face[m-1] * (temperature[m] - T_pre_prev)
            ghf_cumulative += ghf
            H[m-1] += ghf - F_prev

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
            E_supplied = sw_total + longwave_downward + longwave_upward + thf + ghf
            E_delta = E_used - E_supplied

            E_tolerance = energy_tolerance(E_initial)
            if (abs(E_delta) > E_tolerance) || isnan(E_delta)
                @error "inputs" temperature[1] water_surface grain_radius[1] sum(shortwave_flux) cfs.longwave_downward cfs.temperature_air cfs.wind_speed cfs.vapor_pressure cfs.pressure_air
                @error "internals" sw_total longwave_downward longwave_upward thf ghf
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
    ghf_out = ghf_cumulative / cfs.dt  # J -> W/m2
    evaporation_condensation_out = EC_cumulative

    return temperature, longwave_upward_out, heat_flux_sensible_out, heat_flux_latent_out, ghf_out, evaporation_condensation_out
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
