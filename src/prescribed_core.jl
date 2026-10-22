# One time step of a `lower_boundary = :prescribed` column: snow, possibly none, on a
# fixed-temperature reservoir whose temperature is a forcing series. Every `:ice` path is in
# `gemb_core`; this file is reached only under `:prescribed`.

"""
    _thermal_solve_prescribed!(temperature, dz, density, water, T_base, shortwave_flux, sfc,
                               emissivity, emissivity_melt_switch, cfs, mp)
        -> (longwave_upward, heat_flux_sensible, heat_flux_latent, heat_flux_basal,
            evaporation_condensation, evaporation_condensation_liquid, refrozen)

Advance `ns ≥ 1` snow cells over one forcing step with the explicit finite-volume enthalpy scheme
of [`_thermal_solve_ground!`](@ref), the snow resting on a reservoir held at `T_base` [K]. The
reservoir is the face below the deepest cell: the conductance there is that cell's half-thickness,
`2K_n/dz_n`, and the snow side of the face is taken at `min(T_n, 273.15)` as at the snow–soil
face of `:ground` (Crocus's lower-boundary flux, `snowcro.F90`, `PGBAS`). `shortwave_flux` [W
m-2] covers the `ns` cells and, last, what the snow transmits, which the reservoir absorbs.

`heat_flux_basal` [W m-2, positive upward] is the step mean of the conductive flux from the
reservoir into the deepest snow cell. The surface balance, the wet-cell plateau and
[`_wet_cell_handoff!`](@ref) are those of the snow cells of `:ground`. The sub-step is the
explicit stability limit of the snow cells, the reservoir face included, and the surface
relaxation bound. `refrozen` [kg m-2] is the pore water frozen by `_wet_cell_handoff!`.
"""
function _thermal_solve_prescribed!(temperature::Vector{Float64}, dz::Vector{Float64},
    density::Vector{Float64}, water::Vector{Float64}, T_base::Float64,
    shortwave_flux::Vector{Float64}, sfc::_ThermalSurface, emissivity::Float64,
    emissivity_melt_switch::Bool, cfs::ClimateForcingStep, mp::ModelParameters)

    ns = length(dz)
    ns >= 1 || throw(ArgumentError("_thermal_solve_prescribed! needs at least one snow cell"))
    length(shortwave_flux) == ns + 1 || throw(DimensionMismatch(
        "shortwave_flux has $(length(shortwave_flux)) entries, the snow $(ns) cells plus the base"))

    h_water = specific_enthalpy_water(mp)
    M = [density[i] * dz[i] for i in 1:ns]
    W = copy(water)
    H = [M[i] * specific_enthalpy(mp, temperature[i]) + water[i] * h_water for i in 1:ns]
    K = thermal_conductivity(temperature, density, mp)
    Cmin = [density[i] * heat_capacity(mp, temperature[i]) * dz[i] for i in 1:ns]
    T = copy(temperature)
    W[1] > 0.0 && (T[1] = _wet_cell_temperature(mp, H[1], M[1], W[1]))

    # Face conductances [W m-2 K-1]: G[i] between cells i and i+1, G[ns] to the reservoir.
    G = Vector{Float64}(undef, ns)
    for i in 1:ns-1
        G[i] = 1.0 / (dz[i+1] / (2 * K[i+1]) + dz[i] / (2 * K[i]))
    end
    G[ns] = 2 * K[ns] / dz[ns]

    max_safe_dt = Inf
    G_prev = 0.0
    for i in 1:ns
        max_safe_dt = min(max_safe_dt, Cmin[i] / (G[i] + G_prev))
        G_prev = G[i]
    end
    Λ0 = _surface_energy_balance_and_slope(T[1], emissivity, sfc, cfs)[6]
    C1_0 = (M[1] + W[1]) * heat_capacity(mp, T[1])
    max_safe_dt = min(max_safe_dt, THERMAL_SURFACE_STEP_FRACTION * C1_0 / (G[1] + abs(Λ0)))
    dt = _find_dt_divisor(max_safe_dt * mp.thermal_explicit_safety_factor, mp.dt_divisors)
    n_steps = round(Int, cfs.dt / dt)

    Q_sw = [shortwave_flux[i] * dt for i in 1:ns]
    A = G .* dt
    longwave_downward = cfs.longwave_downward * dt
    cell_T(i, Hi) = W[i] > 0.0 ? _wet_cell_temperature(mp, Hi, M[i], W[i]) :
                    temperature_from_specific_enthalpy(mp, Hi / M[i])

    lw_cum = 0.0; shf_cum = 0.0; lhf_cum = 0.0; basal_cum = 0.0
    EC_cum = 0.0; EC_liquid_cum = 0.0
    T_pre = Vector{Float64}(undef, ns)

    for _ in 1:n_steps
        T1_0 = T[1]
        liquid = _surface_liquid(mp, H[1], M[1], W[1]) + EC_liquid_cum > 0.0
        lw_up, shf, lhf, latent_heat, _, dlw, dshf, dlhf =
            _surface_fluxes_and_slopes(T1_0, emissivity, sfc, cfs, liquid)
        Λ = dlw + dshf + dlhf
        if Λ < 0.0
            C1 = (M[1] + W[1]) * heat_capacity(mp, T1_0)
            T_below = ns > 1 ? cell_T(2, H[2]) : T_base
            F1 = (ns > 1 ? A[1] : A[ns]) * (T_below - (ns > 1 ? T1_0 : min(T1_0, CtoK)))
            Q0 = longwave_downward + (lw_up + shf + lhf) * dt
            ΔT1 = (Q_sw[1] + Q0 + F1) / (C1 - Λ * dt)
            lw_up += dlw * ΔT1
            shf += dshf * ΔT1
            lhf += dlhf * ΔT1
        end
        lhf_cum += lhf * dt
        shf_cum += shf * dt
        lw_cum -= lw_up * dt
        EC = lhf / latent_heat * dt
        EC_cum += EC
        latent_heat == LV && (EC_liquid_cum += EC)
        Q_surface = longwave_downward + lw_up * dt + (shf + lhf) * dt

        for i in 1:ns
            T_pre[i] = cell_T(i, H[i])
        end
        H[1] += Q_surface
        F_prev = 0.0
        for i in 1:ns-1
            F = A[i] * (T_pre[i+1] - T_pre[i])
            H[i] += Q_sw[i] + F - F_prev
            F_prev = F
        end
        F_basal = A[ns] * (T_base - min(T_pre[ns], CtoK))
        H[ns] += Q_sw[ns] + F_basal - F_prev
        basal_cum += F_basal

        T[1] = cell_T(1, H[1])
        if emissivity_melt_switch
            emissivity = T[1] < (CtoK - T_MELT_SWITCH_TOLERANCE) ? mp.emissivity :
                         mp.emissivity_grain_radius_large
        end
    end

    refrozen = zeros(ns)
    for i in 1:ns
        if water[i] > 0.0
            _wet_cell_handoff!(mp, temperature, dz, density, water, refrozen, i, H[i], M[i])
        else
            temperature[i] = temperature_from_specific_enthalpy(mp, H[i] / M[i])
        end
    end
    return lw_cum / cfs.dt, shf_cum / cfs.dt, lhf_cum / cfs.dt, basal_cum / cfs.dt, EC_cum,
        EC_liquid_cum, sum(refrozen)
end

"""
    _prescribed_bare_surface(T_surface, emissivity, β, cfs, mp, table)
        -> (longwave_upward, heat_flux_sensible, heat_flux_latent, evaporation_condensation,
            evaporation_condensation_liquid)

Surface fluxes [W m-2, positive into the surface] of a surface held at `T_surface` [K], which is
not solved for: snow-free prescribed ground, or trace snow on it. Emitted and reflected longwave
at `emissivity`, and the turbulent fluxes of [`_ground_turbulent_heat_flux`](@ref) (free
convection included, the latent flux scaled by `β`), with the transfer coefficients of `table`.
`evaporation_condensation` [kg m-2] is the step's vapor mass, of which
`evaporation_condensation_liquid` was exchanged at `LV`.
"""
function _prescribed_bare_surface(T_surface::Float64, emissivity::Float64, β::Float64,
    cfs::ClimateForcingStep, mp::ModelParameters, table::StabilityTable)
    U = max(cfs.wind_speed, 0.01)
    sfc = _ThermalSurface(air_density(cfs.pressure_air, cfs.temperature_air), emissivity, false,
        VON_KARMAN^2 * U, GRAVITY * cfs.wind_observation_height / U^2, table)
    lw_up, shf, lhf, latent_heat, _, _, _ = _ground_surface_fluxes_and_slopes(T_surface,
        emissivity, sfc, cfs, β)
    EC = lhf / latent_heat * cfs.dt
    return -lw_up, shf, lhf, EC, latent_heat == LV ? EC : 0.0
end

"""
    _trace_snow_surface(cols, T_surface, cfs, mp, table)
        -> (longwave_upward, heat_flux_sensible, heat_flux_latent, evaporation_condensation,
            evaporation_condensation_liquid)

Surface fluxes of trace snow held at `T_surface` [K]: the snow surface balance of
[`_surface_fluxes_and_slopes`](@ref) at the snow emissivity, evaporating at `LV` only while the
cell (if any is left) has pore water to give. Snow that has melted out leaves the surface at the
melting point, its fluxes those of the snow it was at the start of the step.
"""
function _trace_snow_surface(cols::NamedTuple, T_surface::Float64, cfs::ClimateForcingStep,
    mp::ModelParameters, table::StabilityTable)
    U = max(cfs.wind_speed, 0.01)
    sfc = _ThermalSurface(air_density(cfs.pressure_air, cfs.temperature_air), mp.emissivity,
        false, VON_KARMAN^2 * U, GRAVITY * cfs.wind_observation_height / U^2, table)
    liquid = !isempty(cols.water) && cols.water[1] > 0.0
    lw_up, shf, lhf, latent_heat, _, _, _, _ = _surface_fluxes_and_slopes(T_surface,
        mp.emissivity, sfc, cfs, liquid)
    EC = lhf / latent_heat * cfs.dt
    return -lw_up, shf, lhf, EC, latent_heat == LV ? EC : 0.0
end

"""
    _trace_snow_on_base!(cols, T_base, mp) -> (E_change, melt, runoff, refreeze)

Trace snow (lighter than [`snow_layer_min_mass`](@ref)`(dt)`, merged into one cell) on the
prescribed base at `T_base` [K]. Its thermal resistance and heat capacity are negligible against
the reservoir's, so it takes the base temperature:

- `T_base < 273.15 K`: the cell is set to `T_base`, its pore water frozen into the matrix at
  constant density;
- `T_base = 273.15 K`: the cell is at the melting point, its pore water kept;
- `T_base > 273.15 K`: the snow melts and runs off.

Crocus treats a vanishing snow cover the same way, handing its mass and energy to the ground
below (`snowcro.F90`, `SNOWCROGONE`). Returns the change of the cell's enthalpy [J m-2] (the
melt water counted at the melting-point enthalpy it leaves with), and the melt, runoff and
refreeze [kg m-2].
"""
function _trace_snow_on_base!(cols::NamedTuple, T_base::Float64, mp::ModelParameters)
    _merge_all_snow!(cols, mp)
    (; temperature, dz, density, water) = cols
    h_water = specific_enthalpy_water(mp)
    M = dz[1] * density[1]
    W = water[1]
    H_old = M * specific_enthalpy(mp, temperature[1]) + W * h_water
    if T_base > CtoK
        close_slot!(cols, 1)
        return (M + W) * h_water - H_old, M, M + W, 0.0
    elseif T_base == CtoK
        temperature[1] = CtoK
        return M * specific_enthalpy(mp, CtoK) + W * h_water - H_old, 0.0, 0.0, 0.0
    end
    temperature[1] = T_base
    dz[1] = (M + W) / density[1]
    water[1] = 0.0
    return (M + W) * specific_enthalpy(mp, T_base) - H_old, 0.0, 0.0, W
end

"""
    gemb_core_prescribed(snow, T_base, cfs, mp, verbose;
                         thermal_workspace=ThermalWorkspace(), column_workspace=ColumnWorkspace())
        -> (snow, flux)

One time step of the snow `snow` ([`column_state`](@ref), possibly empty) on a fixed-temperature
reservoir at `T_base` [K], this step's value of the `temperature_ground` forcing layer
(`lower_boundary = :prescribed`). No ground column is solved: the reservoir stands for the ground
surface, whose temperature is prescribed, typically from a measured near-surface soil
temperature. It has no heat capacity in the model and is an argument of the step, not a cell, so
its temperature changes only between steps, to the next prescribed value.

With snow, the step is that of the snow cells of [`gemb_core_ground`](@ref), the thermal solve
being [`_thermal_solve_prescribed!`](@ref):

1. grain metamorphism;
2. albedo, the snow blending to `mp.ground_albedo` as over `:ground`
   ([`calculate_albedo_ground`](@ref));
3. shortwave, the part the snow transmits absorbed by the reservoir;
4. the thermal solve, the deepest cell exchanging heat with the reservoir;
5. to 8. vapor exchange, snowfall and rain, blowing snow, melt and percolation (meltwater leaving
   the base of the snow runs off), thickness bands and densification
   ([`_snow_after_thermal!`](@ref)).

Snow lighter than [`snow_layer_min_mass`](@ref)`(dt)` is not solved thermally: it takes the base
temperature ([`_trace_snow_on_base!`](@ref)), and its surface balance is that of snow at that
temperature capped at 273.15 K ([`_trace_snow_surface`](@ref)).

With no snow the surface is the prescribed ground: albedo `mp.ground_albedo`, its energy balance
evaluated at `T_base` with `mp.ground_emissivity` and `mp.ground_evaporation_factor`
([`_prescribed_bare_surface`](@ref)), and the reservoir absorbs the net. Snowfall on it starts a
new snow cell with [`fresh_snow_properties`](@ref); rain runs off at the melting point, its
sensible heat above it going to the reservoir.

`heat_flux_basal` [W m-2, positive upward] is the energy the reservoir supplies to the snow
column over the step: the conductive flux into the deepest snow cell, less the melt energy the
snow passes down ([`calculate_melt`](@ref)'s basal surplus) and the sensible heat of rain on
snow-free ground; with no or trace snow, minus the energy the surface delivers to the ground, plus
the trace snow's enthalpy change. Shortwave the snow transmits is absorbed by the reservoir and
is not part of it. `heat_flux_ground` is the same quantity. With `verbose` the step's mass and
energy budgets are checked with the energy through the base counted as `heat_flux_basal`.
"""
function gemb_core_prescribed(snow::NamedTuple, T_base::Float64, cfs::ClimateForcingStep,
    mp::ModelParameters, verbose::Bool;
    thermal_workspace::ThermalWorkspace=ThermalWorkspace(),
    column_workspace::ColumnWorkspace=ColumnWorkspace())

    isfinite(T_base) || error("gemb_core_prescribed: the prescribed base temperature " *
        "(forcing layer temperature_ground) is $(T_base) K")
    (; temperature, dz, density, water, grain_radius, grain_dendricity,
       grain_sphericity, age) = snow
    age .+= cfs.dt / 86400.0
    h_water = specific_enthalpy_water(mp)

    if verbose
        M_initial = column_mass(snow)
        E_initial = column_enthalpy(mp, dz .* density, temperature, water)
    end

    melt = 0.0; runoff = 0.0; refreeze = 0.0

    # 1. Grain metamorphism
    if !isempty(dz)
        grain_radius, grain_dendricity, grain_sphericity =
            calculate_grain_size(temperature, dz, density, water, grain_radius,
                grain_dendricity, grain_sphericity, cfs, mp)
    end
    cols = column_state(temperature, dz, density, water, grain_radius, grain_dendricity,
        grain_sphericity, age)

    # 2. Albedo and 3. shortwave
    albedo_broadband, albedo_diffuse = calculate_albedo_ground(dz, density, grain_radius, cfs, mp)
    sw_net = shortwave_net_surface(albedo_broadband, albedo_diffuse, cfs)
    z0 = isempty(dz) ? mp.ground_roughness : surface_roughness(mp, density[1], water[1])
    table = _stability_table!(thermal_workspace, cfs.wind_observation_height,
        cfs.temperature_observation_height, z0, z0 * mp.surface_roughness_effective_ratio)

    # 4. Thermal step. `E_base` [J m-2] is the energy the reservoir supplies to the column;
    # `sw_base` [W m-2] the shortwave it absorbs.
    if !isempty(dz) && column_mass(cols) >= snow_layer_min_mass(cfs.dt)
        sw_snow = calculate_shortwave_radiation(dz, density, grain_radius, albedo_broadband,
            albedo_diffuse, cfs, mp)
        shortwave_flux = vcat(sw_snow, sw_net - sum(sw_snow))
        emissivity, melt_switch = _emissivity_initialize(grain_radius[1], mp)
        U = max(cfs.wind_speed, 0.01)
        sfc = _ThermalSurface(air_density(cfs.pressure_air, cfs.temperature_air), emissivity,
            melt_switch, VON_KARMAN^2 * U, GRAVITY * cfs.wind_observation_height / U^2, table)
        longwave_upward, heat_flux_sensible, heat_flux_latent, F_basal, EC, EC_liquid,
            refreeze_thermal = _thermal_solve_prescribed!(temperature, dz, density, water,
                T_base, shortwave_flux, sfc, emissivity, melt_switch, cfs, mp)
        refreeze += refreeze_thermal
        E_base = F_basal * cfs.dt
        sw_base = shortwave_flux[end]
    else
        if isempty(dz)
            longwave_upward, heat_flux_sensible, heat_flux_latent, EC, EC_liquid =
                _prescribed_bare_surface(T_base, mp.ground_emissivity,
                    mp.ground_evaporation_factor, cfs, mp, table)
            E_change = 0.0
        else
            # Trace snow at the base temperature, its surface balance that of snow there.
            E_change, melt_trace, runoff_trace, refreeze_trace =
                _trace_snow_on_base!(cols, T_base, mp)
            melt += melt_trace; runoff += runoff_trace; refreeze += refreeze_trace
            longwave_upward, heat_flux_sensible, heat_flux_latent, EC, EC_liquid =
                _trace_snow_surface(cols, min(T_base, CtoK), cfs, mp, table)
        end
        Q_surface = (sw_net + cfs.longwave_downward - longwave_upward + heat_flux_sensible +
                     heat_flux_latent) * cfs.dt
        # The column's enthalpy changes by `E_change` (the trace snow's, with its melt water
        # leaving at `h_water`), so the reservoir supplies what the surface does not.
        E_base = E_change - Q_surface
        sw_base = 0.0
    end

    snow_new, melt, runoff, refreeze, rain, E_vapor, EC_ground, E_ground_added,
        mass_blowing_snow, E_blowing_snow, melt_surface, percolation_depth, viscosity =
        _snow_after_thermal!(cols, EC, EC_liquid, cfs, mp, verbose, column_workspace,
            melt, runoff, refreeze)
    E_base -= E_ground_added
    heat_flux_basal = E_base / cfs.dt
    (; dz, density, temperature, water) = snow_new

    if verbose
        M_final = column_mass(snow_new)
        M_delta = M_final - M_initial + runoff - cfs.precipitation - (EC - EC_ground) -
                  mass_blowing_snow
        abs(M_delta) > M_TOLERANCE && error("gemb_core_prescribed: mass not conserved: " *
                                            "M_delta = $(M_delta)")
        E_final = column_enthalpy(mp, dz .* density, temperature, water)
        E_snowfall = (cfs.precipitation - rain) * specific_enthalpy(mp, cfs.temperature_air)
        E_rain = rain * specific_enthalpy_water(mp, cfs.temperature_air)
        E_supplied = (sw_net - sw_base + cfs.longwave_downward - longwave_upward +
                      heat_flux_sensible + heat_flux_latent + heat_flux_basal) * cfs.dt +
                     E_snowfall + E_rain + E_vapor + E_blowing_snow - runoff * h_water
        E_delta = (E_final - E_initial) - E_supplied
        abs(E_delta) > energy_tolerance(E_initial) && error("gemb_core_prescribed: energy " *
            "not conserved: E_delta = $(E_delta) J m-2")
    end

    ice_slab_thickness, ice_slab_depth = ice_slab_diagnostics(dz, density, mp)
    flux = (
        albedo_broadband=effective_albedo(albedo_broadband, albedo_diffuse, cfs),
        shortwave_net=sw_net,
        heat_flux_sensible=heat_flux_sensible,
        heat_flux_latent=heat_flux_latent,
        heat_flux_basal=heat_flux_basal,
        heat_flux_ground=heat_flux_basal,
        longwave_upward=longwave_upward,
        evaporation_condensation=EC,
        evaporation_condensation_ground=EC_ground,
        rain=rain,
        melt=melt,
        melt_surface=melt_surface,
        runoff=runoff,
        refreeze=refreeze,
        mass_blowing_snow=mass_blowing_snow,
        percolation_depth=percolation_depth,
        ice_slab_thickness=ice_slab_thickness,
        ice_slab_depth=ice_slab_depth,
        viscosity=viscosity,
        temperature_base=T_base,
    )
    return snow_new, flux
end
