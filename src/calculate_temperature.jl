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

    # Copy inputs to avoid mutation
    temperature = copy(temperature)

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

    # if wind_speed = 0, goes to infinity; set minimum
    wind_speed_local = max(cfs.wind_speed, 0.01)
    # Create a modified ClimateForcingStep with clamped wind speed
    cfs_local = ClimateForcingStep(
        cfs.dt, cfs.temperature_air, cfs.pressure_air, cfs.precipitation,
        wind_speed_local, cfs.shortwave_downward, cfs.longwave_downward,
        cfs.vapor_pressure, cfs.temperature_air_mean, cfs.wind_speed_mean,
        cfs.precipitation_mean, cfs.temperature_observation_height,
        cfs.wind_observation_height, cfs.black_carbon_snow, cfs.black_carbon_ice,
        cfs.cloud_optical_thickness, cfs.solar_zenith_angle,
        cfs.shortwave_downward_diffuse, cfs.cloud_fraction
    )

    ## THERMAL CONDUCTIVITY (Sturm, 1997)
    K = thermal_conductivity(temperature, density, mp)

    ## THERMAL DIFFUSION COEFFICIENTS
    # Patankar 1980, Ch. 3&4

    # u, d, and p conductivities
    KU = vcat([NaN], K[1:m-1])
    KD = vcat(K[2:m], [NaN])
    KP = K

    # determine u, d & p cell widths
    dzU = vcat([NaN], dz[1:m-1])
    dzD = vcat(dz[2:m], [NaN])

    # find stable dt for thermodynamics loop
    dt = _thermo_optimal_dt(dz, density, C_ICE, K, mp.dt_divisors)

    # determine mean (harmonic mean) of K/dz for u, d, & p
    Au = (dzU ./ (2 .* KU) .+ dz ./ (2 .* KP)) .^ (-1)
    Ad = (dzD ./ (2 .* KD) .+ dz ./ (2 .* KP)) .^ (-1)
    Ap = (density .* dz .* C_ICE) ./ dt

    # Create neighbor coefficient arrays
    Nu = Au ./ Ap
    Nd = Ad ./ Ap
    Np = 1 .- Nu .- Nd

    # specify boundary conditions
    # Constant Temperature (Dirichlet) boundary condition at bottom
    Nu[m] = 0.0
    Np[m] = 1.0
    Nd[m] = 0.0

    # zero flux at surface
    Nu[1] = 0.0         # Disconnect from the node above (Air/Ghost node)
    Np[1] = 1 - Nd[1]   # Balance the center node to conserve energy

    ## RADIATIVE FLUXES

    # energy supplied by shortwave radiation [J]
    sw = shortwave_flux .* dt

    # ensure no sw reaches bottom cell, add any flux to bottom cell to the cell above
    sw[end-1] = sw[end-1] + sw[end]
    sw[end] = 0.0

    # temperature change due to SW
    T_delta_sw = sw ./ (C_ICE .* density .* dz)

    # energy supplied by downward longwave radiation to the top grid cell [J]
    longwave_downward = cfs.longwave_downward * dt

    # temperature change due to longwave_downward
    T_delta_longwave_downward = longwave_downward / TCs

    ## PREALLOCATE ARRAYS
    Tu = zeros(m)
    Td = zeros(m)

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
        heat_flux_sensible, heat_flux_latent, latent_heat = turbulent_heat_flux(T_surface, density_air, z0, zT, zQ, cfs_local)

        lhf_cumulative += heat_flux_latent * dt
        shf_cumulative += heat_flux_sensible * dt

        # mass loss (-)/accretion(+) due to evaporation/condensation [kg]
        evaporation_condensation = heat_flux_latent / latent_heat * dt

        # temperature change due to turbulent fluxes
        thf = (heat_flux_sensible + heat_flux_latent) * dt
        T_delta_thf = thf / TCs

        # upward longwave radiation
        longwave_upward = -(SB * T_surface^4.0 * emissivity) * dt
        longwave_upward_cumulative += -longwave_upward
        T_delta_longwave_upward = longwave_upward / TCs

        # new grid point temperature
        # SW penetrates surface
        temperature .+= T_delta_sw
        temperature[1] += T_delta_longwave_downward + T_delta_longwave_upward + T_delta_thf

        # energy flux across lower boundary
        ghf = Ad[end-1] * (temperature[end] - temperature[end-1]) * dt
        ghf_cumulative += ghf

        # temperature diffusion - optimized with in-place operations
        @inbounds begin
            Tu[1] = temperature[1]
            @simd for i in 2:m
                Tu[i] = temperature[i-1]
            end

            @simd for i in 1:m-1
                Td[i] = temperature[i+1]
            end
            Td[m] = temperature[m]

            # In-place fused broadcast (eliminates allocation)
            @. temperature = (Np * temperature) + (Nu * Tu) + (Nd * Td)
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
            E_supplied = sum(sw) + longwave_downward + longwave_upward + thf + ghf
            E_delta = E_used - E_supplied

            E_tolerance = 1e-3
            if (abs(E_delta) > E_tolerance) || isnan(E_delta)
                @error "inputs" temperature[1] water_surface grain_radius[1] sum(shortwave_flux) cfs.longwave_downward cfs.temperature_air cfs.wind_speed cfs.vapor_pressure cfs.pressure_air
                @error "internals" sum(sw) longwave_downward longwave_upward thf ghf
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
    _thermo_optimal_dt(dz, density, C_ice, K, global_dt_or_dt_divisors)

Find optimal time step for numerical stability of the explicit diffusion scheme.
Uses Von Neumann stability analysis with a 0.8 safety factor.
"""
function _thermo_optimal_dt(dz::Vector{Float64}, density::Vector{Float64},
    C_ice::Float64, K::Vector{Float64},
    global_dt_or_dt_divisors)

    # Calculate the theoretical stability limit for every single grid cell
    stability_limit_per_cell = 0.5 .* (density .* C_ice .* dz .^ 2) ./ K

    # Find the bottleneck
    max_safe_dt = minimum(stability_limit_per_cell)

    # Apply a Safety Factor (0.8)
    dt_target = max_safe_dt * 0.8

    # Sanity check
    if dt_target < 1e-4
        @warn "Timestep is extremely small ($dt_target). Check for near-zero dz layers."
    end

    # Fit this target into input data frequency
    if isa(global_dt_or_dt_divisors, Number)
        dt_divisors = fast_divisors(round(Int, global_dt_or_dt_divisors * 10000)) ./ 10000
    else
        dt_divisors = global_dt_or_dt_divisors
    end

    idx = findlast(dt_divisors .<= dt_target)

    if isnothing(idx)
        @warn "thermo dt_target < all dt_divisors, setting thermo == to the smallest dt_divisors... this may make thermo diffusion unstable"
        dt = dt_divisors[1]  # Fallback to smallest possible step
    else
        dt = dt_divisors[idx]
    end

    return dt
end
