# One time step of a `lower_boundary = :ground` column. Every `:ice` path is in `gemb_core`; this
# file is reached only under `:ground`.

"""
    _ground_vapor_exchange!(cols, EC, EC_liquid, mp) -> (E, EC_ground)

Exchange the step's vapor mass `EC` [kg m-2, positive for deposition], of which `EC_liquid` at
`LV`, with the snow of `cols`, as [`expose_vapor_exchange_layer!`](@ref) and
[`apply_vapor_exchange!`](@ref) do for a glacier column, and return the enthalpy `E` [J m-2] the
exchanged mass carried into the snow and the part `EC_ground` [kg m-2] exchanged with the ground.

- With no snow all of it is exchanged with the ground: `EC_ground = EC`.
- A loss the snow cannot supply, once merged into one cell, takes that whole cell (its mass leaves
  at its enthalpy), and the rest, `EC_ground`, is exchanged with the ground.

Vapor exchanged with the ground changes neither the column's mass nor its enthalpy beyond the
latent heat flux the surface balance already charged: the soil water is held fixed, so
evaporation from the ground stands for water the soil resupplies, as Crocus's `PEVAPCOR` takes a
thin snow cover's excess sublimation from the soil (`snowcro.F90`, `SNOWCROGONE`).
"""
function _ground_vapor_exchange!(cols::NamedTuple, EC::Float64, EC_liquid::Float64,
    mp::ModelParameters)
    isempty(cols.dz) && return 0.0, EC
    expose_vapor_exchange_layer!(cols, EC, EC_liquid, mp)
    (; temperature, dz, density, water) = cols
    ice_change = _vapor_ice_change(EC, EC_liquid, water[1])
    M1 = dz[1] * density[1]
    supported = M1 + ice_change > 0.0 && dz[1] + ice_change / density[1] > 0.0 &&
                _vapor_latent_draw_supported(temperature[1], M1, water[1], EC, EC_liquid, mp)
    if !supported
        # Only one cell is left (`expose_vapor_exchange_layer!` merges until it supports the
        # loss or one remains), and it sublimates or evaporates entirely.
        E = -(M1 * specific_enthalpy(mp, temperature[1]) + water[1] * specific_enthalpy_water(mp))
        EC_ground = EC + M1 + water[1]
        close_slot!(cols, 1)
        return E, EC_ground
    end
    S_surface = M1 + water[1]
    E = apply_vapor_exchange!(temperature, dz, density, water, EC, EC_liquid, mp)
    cols.age[1] = dilute_age(cols.age[1], S_surface, EC)
    return E, 0.0
end

# Merge every snow cell of `cols` into one, conserving mass and enthalpy.
function _merge_all_snow!(cols::NamedTuple, mp::ModelParameters)
    while length(cols.dz) > 1
        _merge_thin_cell!(cols, 1, 2, cols.dz[1] * cols.density[1],
            cols.dz[2] * cols.density[2], mp)
        close_slot!(cols, 1)
    end
    return nothing
end

"""
    _snow_after_thermal!(cols, EC, EC_liquid, cfs, mp, verbose, column_workspace,
                         melt, runoff, refreeze)
        -> (snow, melt, runoff, refreeze, rain, E_vapor, EC_ground, E_ground_added,
            mass_blowing_snow, E_blowing_snow, melt_surface, percolation_depth, viscosity)

Steps 5 to 8 of [`gemb_core_ground`](@ref), shared by every column whose snow rests on ground
(`lower_boundary = :ground` and `:prescribed`): vapor exchange, snowfall and rain, blowing snow,
melt and percolation, thickness bands and densification of the snow `cols`, after the thermal
solve has supplied the step's vapor mass `EC` (of which `EC_liquid` at `LV`).

`melt`, `runoff` and `refreeze` [kg m-2] are the step's totals so far, returned with this part
added. `E_ground_added` [J m-2] is the energy the ground below receives: the sensible heat above
the melting point of rain on snow-free ground, which runs off at 273.15 K, and the melt energy
beyond what melts the deepest snow cell. `EC_ground` [kg m-2] is the vapor exchanged with the
ground rather than the snow.
"""
function _snow_after_thermal!(cols::NamedTuple, EC::Float64, EC_liquid::Float64,
    cfs::ClimateForcingStep, mp::ModelParameters, verbose::Bool,
    column_workspace::ColumnWorkspace, melt::Float64, runoff::Float64, refreeze::Float64)

    (; temperature, dz, density, water, grain_radius, grain_dendricity,
       grain_sphericity, age) = cols
    h_water = specific_enthalpy_water(mp)
    rain = 0.0

    # 5. Vapor exchange
    E_vapor, EC_ground = _ground_vapor_exchange!(cols, EC, EC_liquid, mp)

    # 6. Snowfall and rain
    E_ground_added = 0.0
    if isempty(dz)
        if cfs.precipitation > 0
            if cfs.temperature_air <= (mp.rain_temperature_threshold + T_TOLERANCE)
                ρ_new, r_new, d_new, s_new = fresh_snow_properties(cfs, mp, mp.ground_roughness)
                push!(temperature, cfs.temperature_air); push!(dz, cfs.precipitation / ρ_new)
                push!(density, ρ_new); push!(water, 0.0); push!(grain_radius, r_new)
                push!(grain_dendricity, d_new); push!(grain_sphericity, s_new); push!(age, 0.0)
            else
                rain = cfs.precipitation
                runoff += rain
                E_ground_added = rain * (specific_enthalpy_water(mp, cfs.temperature_air) - h_water)
            end
        end
    else
        temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity,
            age, rain = calculate_accumulation(temperature, dz, density, water, grain_radius,
                grain_dendricity, grain_sphericity, age, cfs, mp, verbose)
    end
    cols = column_state(temperature, dz, density, water, grain_radius, grain_dendricity,
        grain_sphericity, age)

    # 7. Blowing snow, melt and percolation
    mass_blowing_snow = 0.0; E_blowing_snow = 0.0
    if !isempty(dz)
        m1, e1 = apply_blowing_snow!(cols, cfs, mp)
        mass_blowing_snow += m1; E_blowing_snow += e1
    end
    basal_surplus = Ref(0.0)
    melt_surface = 0.0
    percolation_depth = 0.0
    if !isempty(dz)
        temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity,
            age, melt_snow, melt_surface, runoff_snow, refreeze_snow, percolation_depth =
            calculate_melt(temperature, dz, density, water, grain_radius, grain_dendricity,
                grain_sphericity, age, mp, verbose; workspace=column_workspace, basal_surplus)
        melt += melt_snow; runoff += runoff_snow; refreeze += refreeze_snow
        runoff += apply_lateral_drainage!(column_state(temperature, dz, density, water,
            grain_radius, grain_dendricity, grain_sphericity, age), cfs.dt, mp)
    end
    E_ground_added += basal_surplus[]

    # 8. Thickness bands and densification of the snow
    viscosity = NO_VISCOSITY
    if !isempty(dz)
        temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity,
            age, _ = manage_layer_thickness(temperature, dz, density, water, grain_radius,
                grain_dendricity, grain_sphericity, age, mp, verbose; dt=cfs.dt,
                n_target=length(dz), workspace=column_workspace,
                n_max=ground_snow_slots(mp), dirichlet_bottom=false)
        viscosity = mp.output_viscosity ? Vector{Float64}(undef, length(dz)) : NO_VISCOSITY
        dz, density = calculate_density(temperature, dz, density, grain_radius, water, cfs, mp;
            grain_dendricity, grain_sphericity, viscosity)
    end
    snow_new = column_state(temperature, dz, density, water, grain_radius, grain_dendricity,
        grain_sphericity, age)

    return snow_new, melt, runoff, refreeze, rain, E_vapor, EC_ground, E_ground_added,
        mass_blowing_snow, E_blowing_snow, melt_surface, percolation_depth, viscosity
end

"""
    gemb_core_ground(snow, temperature_ground, water_ground, g, cfs, mp, verbose;
                     thermal_workspace=ThermalWorkspace(), column_workspace=ColumnWorkspace())
        -> (snow, temperature_ground, water_ground, flux)

One time step of the snow over the soil column `g` (`lower_boundary = :ground`). `snow` is the
[`column_state`](@ref) of the snow cells, possibly empty; `temperature_ground` [K] and
`water_ground` (unfrozen pore water [kg m-2]) the ground cells. Returns the new snow bundle, the
ground state and the step's `flux` terms.

The sequence is [`gemb_core`](@ref)'s, applied to the snow, with the ground in the thermal solve:

1. grain metamorphism of the snow;
2. albedo by [`calculate_albedo_ground`](@ref): snow-free ground at `mp.ground_albedo`, thin snow
   blended to it;
3. shortwave absorbed in the snow as on a glacier, and what the snow transmits, or with no snow
   all of the net shortwave, absorbed in the top ground cell;
4. the thermal solve of snow and ground together ([`_thermal_solve_ground!`](@ref)). Snow lighter
   than [`snow_layer_min_mass`](@ref)`(dt)` is merged into one cell and solved as frozen (or
   liquid) water of the top ground node, at that node's temperature, and separated again after
   the solve: the snow takes the node's temperature capped at 273.15 K and the node's frozen
   water first, so it stays snow while the node is on its freezing plateau. Melted snow runs off. This keeps the explicit sub-step at the scale of a snow cell
   of that mass;
5. vapor exchange ([`_ground_vapor_exchange!`](@ref));
6. snowfall and rain. Over snow, [`calculate_accumulation`](@ref). On snow-free ground, snowfall
   forms a new snow cell with [`fresh_snow_properties`](@ref), and rain runs off at the melting
   point, its sensible heat above it entering the top ground cell;
7. blowing snow, melt, percolation and refreeze in the snow. Water leaving the base of the snow
   runs off: the ground is impermeable to it. Energy beyond what melts the deepest snow cell enters
   the top ground cell;
8. the snow's thickness bands ([`manage_layer_thickness`](@ref) with a ceiling of
   `mp.ground_snow_cells_max` cells and no fixed-temperature bottom) and densification.

The ground cells are never merged, split, densified or melted: their thickness, solids and pore
water are fixed, and only their enthalpy changes, the pore water freezing and thawing at 273.15 K
on the plateau of [`ground_state`](@ref). The deepest ground cell is the Dirichlet reservoir at
its initial temperature. There is no basal mass flux and no horizontal strain.

With `verbose` the step's mass and energy budgets are checked as in `gemb_core`.
"""
function gemb_core_ground(snow::NamedTuple, temperature_ground::Vector{Float64},
    water_ground::Vector{Float64}, g::GroundColumn, cfs::ClimateForcingStep,
    mp::ModelParameters, verbose::Bool;
    thermal_workspace::ThermalWorkspace=ThermalWorkspace(),
    column_workspace::ColumnWorkspace=ColumnWorkspace())

    (; temperature, dz, density, water, grain_radius, grain_dendricity,
       grain_sphericity, age) = snow
    ng = length(g.dz)
    T_bottom = temperature_ground[ng]
    age .+= cfs.dt / 86400.0
    h_water = specific_enthalpy_water(mp)

    H_ground = [ground_enthalpy(mp, g.C_solid[j], g.water_total[j], temperature_ground[j],
                    water_ground[j]) for j in 1:ng]

    if verbose
        M_initial = column_mass(snow)
        E_initial = column_enthalpy(mp, dz .* density, temperature, water) + sum(H_ground)
    end

    melt = 0.0; runoff = 0.0; refreeze = 0.0; rain = 0.0

    # 1. Grain metamorphism
    if !isempty(dz)
        grain_radius, grain_dendricity, grain_sphericity =
            calculate_grain_size(temperature, dz, density, water, grain_radius,
                grain_dendricity, grain_sphericity, cfs, mp)
    end
    cols = column_state(temperature, dz, density, water, grain_radius, grain_dendricity,
        grain_sphericity, age)

    # 2. Albedo
    albedo_broadband, albedo_diffuse = calculate_albedo_ground(dz, density, grain_radius, cfs, mp)

    # 4 (set-up). Snow lighter than the least layer mass is solved inside the top ground node.
    lumped = false
    if !isempty(dz) && column_mass(cols) < snow_layer_min_mass(cfs.dt)
        _merge_all_snow!(cols, mp)
        lumped = true
    end
    ns = length(dz)
    ns_solve = lumped ? 0 : ns

    # 3. Shortwave
    sw_net = shortwave_net_surface(albedo_broadband, albedo_diffuse, cfs)
    shortwave_flux = zeros(ns_solve + ng)
    if ns_solve > 0
        sw_snow = calculate_shortwave_radiation(dz, density, grain_radius, albedo_broadband,
            albedo_diffuse, cfs, mp)
        shortwave_flux[1:ns] .= sw_snow
        shortwave_flux[ns+1] = sw_net - sum(sw_snow)
    else
        shortwave_flux[1] = sw_net
    end

    # 4. Thermal solve
    Cs = copy(g.C_solid)
    Wn = copy(g.water_total)
    K = [ground_conductivity(g, water_ground[j], g.water_total[j]) for j in 1:ng]
    if lumped
        S = dz[1] * density[1] + water[1]
        Wn[1] += S
        H_ground[1] += dz[1] * density[1] * specific_enthalpy(mp, temperature[1]) +
                       water[1] * h_water
    end
    gn = _GroundNodes(Cs, Wn, H_ground, K, g.dz, T_bottom)

    density_air = air_density(cfs.pressure_air, cfs.temperature_air)
    z0 = ns > 0 ? surface_roughness(mp, density[1], water[1]) : mp.ground_roughness
    zT = z0 * mp.surface_roughness_effective_ratio
    if ns_solve > 0
        emissivity, melt_switch = _emissivity_initialize(grain_radius[1], mp)
    else
        emissivity, melt_switch = (ns > 0 ? mp.emissivity : mp.ground_emissivity), false
    end
    U = max(cfs.wind_speed, 0.01)
    table = _stability_table!(thermal_workspace, cfs.wind_observation_height,
        cfs.temperature_observation_height, z0, zT)
    sfc = _ThermalSurface(density_air, emissivity, melt_switch, VON_KARMAN^2 * U,
        GRAVITY * cfs.wind_observation_height / U^2, table)
    β = ns > 0 ? 1.0 : mp.ground_evaporation_factor

    snow_solve = lumped ? column_state(Float64[], Float64[], Float64[], Float64[], Float64[],
        Float64[], Float64[], Float64[]) : cols
    longwave_upward, heat_flux_sensible, heat_flux_latent, heat_flux_basal, EC, EC_liquid,
        heat_flux_ground, refreeze_thermal = _thermal_solve_ground!(snow_solve.temperature,
            snow_solve.dz, snow_solve.density, snow_solve.water, gn, shortwave_flux, sfc,
            emissivity, melt_switch, β, cfs, mp)
    refreeze += refreeze_thermal

    if lumped
        # Separate the snow from the top ground node.
        W_total = Wn[1]
        S = W_total - g.water_total[1]
        M_s = dz[1] * density[1]
        W_s = water[1]
        T, L = ground_state(mp, Cs[1], W_total, H_ground[1])
        # The node's frozen water is the snow's first, so a thin snow over soil on its freezing
        # plateau stays snow; above the plateau it is all liquid and runs off.
        snow_liquid = S - min(S, W_total - L)
        T_snow = min(T, CtoK)
        H_snow = (S - snow_liquid) * specific_enthalpy(mp, T_snow) + snow_liquid * h_water
        H_ground[1] -= H_snow
        melt += max(0.0, snow_liquid - W_s)
        refreeze += max(0.0, W_s - snow_liquid)
        if S - snow_liquid <= WATER_TOLERANCE
            runoff += S
            close_slot!(cols, 1)
        else
            M_new = S - snow_liquid
            temperature[1] = T_snow
            water[1] = snow_liquid
            dz[1] = M_new / density[1]
        end
    end
    for j in 1:ng-1
        temperature_ground[j], water_ground[j] = ground_state(mp, g.C_solid[j],
            g.water_total[j], H_ground[j])
    end

    snow_new, melt, runoff, refreeze, rain, E_vapor, EC_ground, E_ground_added,
        mass_blowing_snow, E_blowing_snow, melt_surface, percolation_depth, viscosity =
        _snow_after_thermal!(cols, EC, EC_liquid, cfs, mp, verbose, column_workspace,
            melt, runoff, refreeze)
    (; temperature, dz, density, water) = snow_new
    if E_ground_added != 0.0
        H_ground[1] = ground_enthalpy(mp, g.C_solid[1], g.water_total[1], temperature_ground[1],
            water_ground[1]) + E_ground_added
        temperature_ground[1], water_ground[1] = ground_state(mp, g.C_solid[1],
            g.water_total[1], H_ground[1])
    end

    temperature_ground[ng] == T_bottom || error("gemb_core_ground: the bottom ground cell's " *
        "temperature changed: $(T_bottom) K to $(temperature_ground[ng]) K")
    length(dz) <= ground_snow_slots(mp) || error("gemb_core_ground: $(length(dz)) snow cells " *
        "exceed ground_snow_cells_max = $(ground_snow_slots(mp))")

    if verbose
        M_final = column_mass(snow_new)
        M_delta = M_final - M_initial + runoff - cfs.precipitation - (EC - EC_ground) -
                  mass_blowing_snow
        abs(M_delta) > M_TOLERANCE && error("gemb_core_ground: mass not conserved: " *
                                            "M_delta = $(M_delta)")
        H_g = sum(ground_enthalpy(mp, g.C_solid[j], g.water_total[j], temperature_ground[j],
                  water_ground[j]) for j in 1:ng)
        E_final = column_enthalpy(mp, dz .* density, temperature, water) + H_g
        E_snowfall = (cfs.precipitation - rain) * specific_enthalpy(mp, cfs.temperature_air)
        E_rain = rain * specific_enthalpy_water(mp, cfs.temperature_air)
        E_supplied = (sw_net + cfs.longwave_downward - longwave_upward + heat_flux_sensible +
                      heat_flux_latent + heat_flux_basal) * cfs.dt + E_snowfall + E_rain +
                     E_vapor + E_blowing_snow - runoff * h_water
        E_delta = (E_final - E_initial) - E_supplied
        abs(E_delta) > energy_tolerance(E_initial) && error("gemb_core_ground: energy not " *
            "conserved: E_delta = $(E_delta) J m-2")
    end

    ice_slab_thickness, ice_slab_depth = ice_slab_diagnostics(dz, density, mp)
    flux = (
        albedo_broadband=effective_albedo(albedo_broadband, albedo_diffuse, cfs),
        shortwave_net=sw_net,
        heat_flux_sensible=heat_flux_sensible,
        heat_flux_latent=heat_flux_latent,
        heat_flux_basal=heat_flux_basal,
        heat_flux_ground=heat_flux_ground,
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
    )
    return snow_new, temperature_ground, water_ground, flux
end
