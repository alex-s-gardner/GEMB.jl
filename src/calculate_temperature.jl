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

    d_tolerance = 1e-11
    T_tolerance = 1e-4
    W_tolerance = 1e-13
    ds = density[1]      # density of top grid cell

    # calculated air density [kg/m3]
    density_air = 0.029 * cfs.pressure_air / (R_GAS * cfs.temperature_air)

    # thermal capacity of top grid cell [J/k]
    TCs = density[1] * dz[1] * C_ICE

    # determine grid point 'center' vector size
    m = length(density)
    if m == 0
        error("column has no gridcells: length(density) = 0")
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
    if (ds < (mp.density_ice - d_tolerance)) && (water_surface < W_tolerance)
        z0 = 0.00012       # 0.12 mm for dry snow
    elseif ds >= (mp.density_ice - d_tolerance)
        z0 = 0.0032        # 3.2 mm for ice
    else
        z0 = 0.0013        # 1.3 mm for wet snow
    end

    # determine emissivity
    emissivity, emissivity_melt_switch = _emissivity_initialize(grain_radius[1], mp)

    # zT and zQ are percentage of z0
    zT = z0 * mp.surface_roughness_effective_ratio
    zQ = z0 * mp.surface_roughness_effective_ratio

    # minimum wind speed to avoid division by zero in turbulent flux calculations
    min_wind_speed = 0.01

    ## THERMAL CONDUCTIVITY (Sturm, 1997)
    K = thermal_conductivity(temperature, density, mp)

    ## FIND STABLE dt (allocation-free)
    max_safe_dt = Inf
    @inbounds for i in eachindex(dz)
        sl = 0.5 * density[i] * C_ICE * dz[i]^2 / K[i]
        max_safe_dt = min(max_safe_dt, sl)
    end
    dt = _find_dt_divisor(max_safe_dt * 0.8, mp.dt_divisors)

    ## THERMAL DIFFUSION COEFFICIENTS (fused loop - Patankar 1980, Ch. 3&4)
    Nu = Vector{Float64}(undef, m)
    Nd = Vector{Float64}(undef, m)
    Np = Vector{Float64}(undef, m)
    T_delta_sw = Vector{Float64}(undef, m)
    Ad_penultimate = 0.0  # Ad[m-1] needed for ground heat flux

    @inbounds for i in 1:m
        ap = density[i] * dz[i] * C_ICE / dt
        if i > 1
            Nu[i] = (1.0 / (dz[i-1] / (2 * K[i-1]) + dz[i] / (2 * K[i]))) / ap
        else
            Nu[i] = 0.0
        end
        if i < m
            ad_i = 1.0 / (dz[i+1] / (2 * K[i+1]) + dz[i] / (2 * K[i]))
            Nd[i] = ad_i / ap
            if i == m - 1
                Ad_penultimate = ad_i
            end
        else
            Nd[i] = 0.0
        end
        Np[i] = 1.0 - Nu[i] - Nd[i]
        T_delta_sw[i] = shortwave_flux[i] * dt / (C_ICE * density[i] * dz[i])
    end

    # Boundary conditions
    Nu[1] = 0.0
    Np[1] = 1.0 - Nd[1]
    Nu[m] = 0.0; Nd[m] = 0.0; Np[m] = 1.0

    # Ensure no SW reaches bottom cell: add bottom cell's SW energy to cell above,
    # dividing by cell m-1's thermal mass (matching original sw[m-1]+=sw[m] before T_delta conversion)
    T_delta_sw[m-1] += shortwave_flux[m] * dt / (C_ICE * density[m-1] * dz[m-1])
    T_delta_sw[m] = 0.0

    # energy supplied by downward longwave radiation to the top grid cell [J]
    longwave_downward = cfs.longwave_downward * dt

    # temperature change due to longwave_downward
    T_delta_longwave_downward = longwave_downward / TCs

    ## CALCULATE ENERGY SOURCES AND DIFFUSION FOR EVERY TIME STEP [dt]
    n_steps = round(Int, cfs.dt / dt)

    # Local variables for loop outputs
    local longwave_upward::Float64
    local heat_flux_sensible::Float64
    local heat_flux_latent::Float64
    local evaporation_condensation::Float64

    for _ in 1:n_steps
        # Store initial temperature for energy conservation check
        if verbose
            E_initial = sum(temperature .* (C_ICE .* density .* dz))
        end

        # calculate temperature of snow surface
        T_surface = min(273.15, temperature[1])

        # TURBULENT HEAT FLUX
        heat_flux_sensible, heat_flux_latent, latent_heat = turbulent_heat_flux(T_surface, density_air, z0, zT, zQ, cfs; min_wind_speed)

        lhf_cumulative += heat_flux_latent * dt
        shf_cumulative += heat_flux_sensible * dt

        # mass loss (-)/accretion(+) due to evaporation/condensation [kg]
        evaporation_condensation = heat_flux_latent / latent_heat * dt

        # temperature change due to turbulent fluxes
        thf = (heat_flux_sensible + heat_flux_latent) * dt
        T_delta_thf = thf / TCs

        # upward longwave radiation
        T2 = T_surface * T_surface
        longwave_upward = -(SB * T2 * T2 * emissivity) * dt
        longwave_upward_cumulative += -longwave_upward
        T_delta_longwave_upward = longwave_upward / TCs

        # new grid point temperature
        # SW penetrates surface
        temperature .+= T_delta_sw
        temperature[1] += T_delta_longwave_downward + T_delta_longwave_upward + T_delta_thf

        # energy flux across lower boundary
        ghf = Ad_penultimate * (temperature[end] - temperature[end-1]) * dt
        ghf_cumulative += ghf

        # temperature diffusion - single in-place pass (Patankar 1980).
        # new T[i] = Np[i]*T[i] + Nu[i]*T[i-1] + Nd[i]*T[i+1], all reading the
        # OLD temperature field. A one-element carry holds the old upstream
        # value so no shifted Tu/Td copies are needed. Edge stencils reduce to
        # T[1]/T[m] (matching the original boundary handling), but Nu[1]=Nd[m]
        # =Nu[m]=0 make those terms vanish anyway.
        @inbounds begin
            Told_prev = temperature[1]   # old T[i-1]; unused at i=1 (Nu[1]=0)
            for i in 1:m
                Tu_i = i == 1 ? temperature[1] : Told_prev
                Td_i = i == m ? temperature[m] : temperature[i+1]
                Told_i = temperature[i]
                temperature[i] = (Np[i] * Told_i) + (Nu[i] * Tu_i) + (Nd[i] * Td_i)
                Told_prev = Told_i
            end
        end

        # calculate cumulative evaporation (+)/condensation(-)
        EC_cumulative += evaporation_condensation

        # emissivity melt switch check
        if emissivity_melt_switch
            if temperature[1] < (CtoK - T_tolerance)
                emissivity = mp.emissivity
            else
                emissivity = mp.emissivity_grain_radius_large
            end
        end

        # CHECK FOR ENERGY CONSERVATION
        if verbose
            E_used = sum(temperature .* (C_ICE .* density .* dz)) - E_initial
            sw_total = sum(shortwave_flux) * dt
            E_supplied = sw_total + longwave_downward + longwave_upward + thf + ghf
            E_delta = E_used - E_supplied

            E_tolerance = 1e-3
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
    gdn_tolerance = 1e-10

    if mp.emissivity_method == :uniform
        emissivity = mp.emissivity
        emissivity_melt_switch = false
    elseif mp.emissivity_method == :grain_radius_threshold
        if grain_radius_surface <= (mp.emissivity_grain_radius_threshold + gdn_tolerance)
            emissivity = mp.emissivity
        else
            emissivity = mp.emissivity_grain_radius_large
        end
        emissivity_melt_switch = false
    elseif mp.emissivity_method == :grain_radius_w_threshold
        if grain_radius_surface <= (mp.emissivity_grain_radius_threshold + gdn_tolerance)
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
    if dt_target < 1e-4
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
