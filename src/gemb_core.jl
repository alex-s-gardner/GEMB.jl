"""
    apply_vapor_exchange!(temperature, dz, density, water, evaporation_condensation,
                          evaporation_condensation_liquid, mp) -> E [J m-2]

Add or remove the vapor mass of one timestep at the surface cell, in the phase whose latent
heat the surface flux charged, and return the enthalpy the mass carried into the column.
Positive masses are deposition or condensation, negative sublimation or evaporation.

- The part exchanged at `LS`, `evaporation_condensation − evaporation_condensation_liquid`, is
  ice: it changes `dz[1]` at constant density and carries `h_ice(T₁)` per kg.
- The part exchanged at `LV`, `evaporation_condensation_liquid`, is liquid water. Condensate is
  added to `water[1]` and evaporation takes `water[1]` first, each at the melt-point water
  enthalpy. The thermal solve charges `LV` only for evaporation its pore water can supply
  ([`_surface_liquid`](@ref)). Evaporation beyond `water[1]` arises only where the solve itself
  froze some of that water at the end ([`_wet_cell_handoff!`](@ref)); it comes from the ice of
  cell 1, the remaining `LS − LV = LF` per kg is drawn from the cell's enthalpy, and the mass
  leaves the column at `h_ice(T₁) + LF`.

Throws if cell 1 would lose all its ice, or if that draw would take its ice below
`WET_HANDOFF_TEMPERATURE_MIN`, the coldest temperature the heat-capacity fit is supported at.
"""
function apply_vapor_exchange!(temperature::Vector{Float64}, dz::Vector{Float64},
    density::Vector{Float64}, water::Vector{Float64}, evaporation_condensation::Float64,
    evaporation_condensation_liquid::Float64, mp::ModelParameters)

    h_ice = specific_enthalpy(mp, temperature[1])
    h_water = specific_enthalpy_water(mp)
    ice_change = evaporation_condensation - evaporation_condensation_liquid
    E = ice_change * h_ice

    if evaporation_condensation_liquid >= 0.0
        water[1] += evaporation_condensation_liquid
        E += evaporation_condensation_liquid * h_water
    else
        evaporated = -evaporation_condensation_liquid
        from_water = min(evaporated, water[1])
        water[1] -= from_water
        E -= from_water * h_water
        from_ice = evaporated - from_water
        if from_ice > 0.0
            M1 = dz[1] * density[1]
            M1_new = M1 + ice_change - from_ice
            # The cell keeps its enthalpy less the departing ice and the extra LF it needs.
            H1_new = M1 * h_ice + ice_change * h_ice - from_ice * (h_ice + LF)
            M1_new > 0.0 || error("apply_vapor_exchange!: evaporation of $(evaporated) kg m-2 " *
                                  "exceeds the surface cell's water and ice")
            H1_new >= M1_new * specific_enthalpy(mp, WET_HANDOFF_TEMPERATURE_MIN) ||
                error("apply_vapor_exchange!: evaporating $(from_ice) kg m-2 of the surface " *
                      "cell's ice at LV leaves LF·$(from_ice) = $(from_ice * LF) J m-2 to be drawn " *
                      "from a cell of $(M1) kg m-2 (dz $(dz[1]) m, density $(density[1]) kg m-3, " *
                      "T $(temperature[1]) K, water $(water[1] + from_water) kg m-2), which would " *
                      "take it below $(WET_HANDOFF_TEMPERATURE_MIN) K")
            temperature[1] = temperature_from_specific_enthalpy(mp, H1_new / M1_new)
            ice_change -= from_ice
            E -= from_ice * (h_ice + LF)
        end
    end

    dz[1] += ice_change / density[1]
    dz[1] > 0.0 || error("apply_vapor_exchange!: sublimation of $(-ice_change) kg m-2 " *
                         "removes the whole surface cell")
    return E
end

# Net ice mass change of the surface cell under `apply_vapor_exchange!` [kg m-2]: the part
# exchanged at LS, less the ice evaporation beyond `water[1]` takes.
@inline function _vapor_ice_change(evaporation_condensation::Float64,
    evaporation_condensation_liquid::Float64, water_surface::Float64)
    ice_change = evaporation_condensation - evaporation_condensation_liquid
    evaporation_condensation_liquid < 0.0 &&
        (ice_change -= max(-evaporation_condensation_liquid - water_surface, 0.0))
    return ice_change
end

# Whether the surface cell's enthalpy can pay the LF that `apply_vapor_exchange!` draws from it
# for evaporation beyond `water[1]` and stay at or above `WET_HANDOFF_TEMPERATURE_MIN`; the same
# test as that function's guard. True when no ice is evaporated.
@inline function _vapor_latent_draw_supported(temperature1::Float64, M1::Float64,
    water_surface::Float64, evaporation_condensation::Float64,
    evaporation_condensation_liquid::Float64, mp::ModelParameters)
    from_ice = max(-evaporation_condensation_liquid - water_surface, 0.0)
    from_ice > 0.0 || return true
    h_ice = specific_enthalpy(mp, temperature1)
    sublimated = evaporation_condensation - evaporation_condensation_liquid
    M1_new = M1 + sublimated - from_ice
    H1_new = M1 * h_ice + sublimated * h_ice - from_ice * (h_ice + LF)
    return M1_new > 0.0 && H1_new >= M1_new * specific_enthalpy(mp, WET_HANDOFF_TEMPERATURE_MIN)
end

"""
    expose_vapor_exchange_layer!(cols, evaporation_condensation,
                                 evaporation_condensation_liquid, mp) -> n_merged

Merge the surface cell into the cell below, as [`manage_layer_thickness`](@ref) merges a thin
cell, while the surface cell cannot supply the timestep's vapor loss; return the number of
cells merged. Mass, enthalpy and pore water are conserved, so [`apply_vapor_exchange!`](@ref)
then takes the loss from the merged cell. The surface cell cannot supply the loss when

- the loss would remove all of its ice, or
- evaporation beyond its pore water takes LF per kg of its ice from its enthalpy, and that
  would take it below `WET_HANDOFF_TEMPERATURE_MIN`.

The surface cell can be as light as [`snow_layer_min_mass`](@ref)`(dt)`, which covers a
sublimation flux of 262 W m-2 sustained over the step. Larger fluxes occur over thin snow on
ablating ice under strong wind, and the step's loss is then supplied by the cell below as well.
A column of one cell is left unchanged, and `apply_vapor_exchange!` throws.
"""
function expose_vapor_exchange_layer!(cols::NamedTuple, evaporation_condensation::Float64,
    evaporation_condensation_liquid::Float64, mp::ModelParameters)
    (; temperature, dz, density, water) = cols
    n_merged = 0
    while length(dz) > 1
        ice_change = _vapor_ice_change(evaporation_condensation,
            evaporation_condensation_liquid, water[1])
        M1 = dz[1] * density[1]
        (M1 + ice_change > 0.0 && dz[1] + ice_change / density[1] > 0.0 &&
         _vapor_latent_draw_supported(temperature[1], M1, water[1], evaporation_condensation,
             evaporation_condensation_liquid, mp)) && break
        _merge_thin_cell!(cols, 1, 2, M1, dz[2] * density[2], mp)
        close_slot!(cols, 1)
        n_merged += 1
    end
    return n_merged
end

"""
    gemb_core(state, cfs::ClimateForcingStep, mp::ModelParameters, verbose::Bool;
              n_target=length(state.dz), z_target=sum(state.dz),
              thermal_workspace=ThermalWorkspace(), column_workspace=ColumnWorkspace())

Perform a single time-step of the GEMB model.

The column is returned with exactly `n_target` cells summing to exactly `z_target` metres,
so both are invariant across the whole run (see `grid_ops.jl` for the two controllers that
enforce them). Cell count is restored by [`manage_layer_thickness`](@ref) at step 9; column depth is
pinned at step 11, *after* [`calculate_density`](@ref), which is the last thing in the
timestep to change `dz`.

Returns `(state, flux)` where `state` is the updated column state and
`flux` contains energy/mass budget terms for output accumulation.

`thermal_workspace` is the reusable scratch for the thermal solve (see
[`ThermalWorkspace`](@ref)). The default allocates a fresh one per call; [`gemb`](@ref)
passes one down so the buffers are grown once for the whole run. One workspace belongs to
one column being stepped — give each concurrent thread its own.
"""
function gemb_core(state, cfs::ClimateForcingStep, mp::ModelParameters, verbose::Bool;
    n_target::Int=length(state.dz), z_target::Float64=sum(state.dz),
    thermal_workspace::ThermalWorkspace=ThermalWorkspace(),
    column_workspace::ColumnWorkspace=ColumnWorkspace())

    # Destructure state - arrays are mutated in-place by physics functions, and the names are
    # rebound to the new arrays the physics functions return. This is safe because gemb_driver
    # rebinds state to new_state after each call. `evaporation_condensation` is deliberately
    # absent: step 5 assigns it unconditionally before any read, so binding it here would be dead.
    (; temperature, dz, density, water, grain_radius, grain_dendricity,
       grain_sphericity, age, melt_surface) = state

    # 0. Advance the age clock, before any mass moves, so mass added later in this timestep
    #    reads as age 0 in this timestep's output. The column carries no absolute time, so
    #    this counter *is* the model clock (days from profile initialization).
    age .+= cfs.dt / 86400.0

    if verbose
        M = dz .* density
        M_total_initial = sum(M) + sum(water)
        E_total_initial = column_enthalpy(mp, M, temperature, water)
        T_bottom = temperature[end]
    end

    # 1. Snow grain metamorphism
    grain_radius, grain_dendricity, grain_sphericity =
        calculate_grain_size(temperature, dz, density, water, grain_radius,
            grain_dendricity, grain_sphericity, cfs, mp)

    # 2. Calculate snow, firn, and ice albedo. Both are scalars diagnosed from the current
    #    column, not carried state. The direct and diffuse albedos combine into the
    #    shortwave-weighted `effective_albedo` returned in `flux` for output.
    albedo_broadband, albedo_diffuse =
        calculate_albedo(dz, density, grain_radius, cfs, mp)

    # 3. Determine distribution of absorbed SW radiation with depth
    _resize_workspace!(column_workspace, length(dz))
    shortwave_flux = calculate_shortwave_radiation!(column_workspace.shortwave, dz, density,
        grain_radius, albedo_broadband, albedo_diffuse, cfs, mp)

    # 4. Calculate net shortwave [W m-2]
    shortwave_net = sum(shortwave_flux)

    # 5. Calculate new temperature-depth profile and turbulent heat flux
    temperature, longwave_upward, heat_flux_sensible, heat_flux_latent, heat_flux_basal,
        evaporation_condensation, evaporation_condensation_liquid, refreeze_thermal =
        calculate_temperature(temperature, dz, density, water[1], grain_radius,
            shortwave_flux, cfs, mp, verbose, thermal_workspace; water)
    # Pore water the solve froze itself (see `_wet_cell_handoff!`) is refreeze like the rest;
    # under the SNICAR aging it enters the grain radius as `calculate_melt`'s refreeze does.
    refreeze_thermal_count = 0
    if refreeze_thermal > 0.0
        refrozen_cells = view(thermal_workspace.refrozen, 1:length(dz))
        refreeze_thermal_count = count(>(0.0), refrozen_cells)
        if independent_grain_radius(mp)
            snicar_refreeze!(grain_radius, grain_dendricity, refrozen_cells, dz .* density,
                water, density)
        end
    end

    # 6. Exchange the vapor mass with the surface cell, in the phase its latent heat implied.
    #    Deposition and condensation (positive) are new mass at age 0, so they dilute the
    #    surface cell's age; removal takes a fraction of the cell and is age-neutral, which
    #    `dilute_age` handles by returning the age unchanged for a non-positive input.
    #    A loss larger than the surface cell's ice is supplied by the cells below it; the cell
    #    count this removes is restored at step 9.
    expose_vapor_exchange_layer!(column_state(temperature, dz, density, water, grain_radius,
            grain_dendricity, grain_sphericity, age),
        evaporation_condensation, evaporation_condensation_liquid, mp)
    S_surface = dz[1] * density[1] + water[1]
    E_evaporation_condensation = apply_vapor_exchange!(temperature, dz, density, water,
        evaporation_condensation, evaporation_condensation_liquid, mp)
    age[1] = dilute_age(age[1], S_surface, evaporation_condensation)

    # 7. Add snow/rain to top grid cell
    temperature, dz, density, water, grain_radius, grain_dendricity,
        grain_sphericity, age, rain =
        calculate_accumulation(temperature, dz, density, water, grain_radius,
            grain_dendricity, grain_sphericity, age, cfs, mp, verbose)

    # 7b. Wind rework of snow already on the ground (both schemes are in `blowing_snow.jl`).
    # Here because wind acts on snow once it has landed, and this timestep's snowfall is the most
    # driftable in the column; before `calculate_melt`, so a wind slab formed now can block this
    # timestep's percolation; before the grid controllers (step 9), so a cell driven out of its
    # thickness band is fixed up immediately. The two schemes are independent and may both be on.
    # A no-op at the `blowing_snow_method = :none`, `drift_rate = 0.0` defaults.
    cols_drift = column_state(temperature, dz, density, water, grain_radius,
        grain_dendricity, grain_sphericity, age)

    # Both signed as a flux *into* the column, like `mass_lateral`.
    mass_sublimated, E_sublimated = apply_blowing_snow!(cols_drift, cfs, mp)
    mass_drift, E_drift = apply_prescribed_drift!(cols_drift, cfs, mp)
    mass_blowing_snow = mass_sublimated + mass_drift
    E_blowing_snow = E_sublimated + E_drift

    # 8. Melt and wet compaction
    densification_from_melt = sum(dz)

    temperature, dz, density, water, grain_radius, grain_dendricity,
        grain_sphericity, age, melt, melt_surface, runoff, refreeze,
        percolation_depth =
        calculate_melt(temperature, dz, density, water, grain_radius,
            grain_dendricity, grain_sphericity, age, mp, verbose;
            workspace=column_workspace)
    refreeze += refreeze_thermal

    densification_from_melt = densification_from_melt - sum(dz)

    # 8b. Drain water standing above irreducible saturation out of the column laterally, over
    # the timescale `mp.runoff_method` prescribes. Immediately after `calculate_melt`, which is
    # what put the water there, and before the grid controllers, so a cell about to be merged
    # has already drained. Only `water` changes, so no grid invariant is at stake and the mass
    # joins `runoff` directly — the runoff enthalpy below is derived from that scalar, and pore
    # water carries the same melt-point enthalpy, so no new energy term is needed either.
    # A no-op at the `runoff_method = :instantaneous` default.
    runoff += apply_lateral_drainage!(
        column_state(temperature, dz, density, water, grain_radius,
            grain_dendricity, grain_sphericity, age), cfs.dt, mp)

    # 9. Return cells to their thickness bands and restore the cell count to n_target
    #    (exactly conservative)
    temperature, dz, density, water, grain_radius, grain_dendricity,
        grain_sphericity, age, E_added =
        manage_layer_thickness(temperature, dz, density, water, grain_radius,
            grain_dendricity, grain_sphericity, age, mp, verbose;
            dt=cfs.dt, n_target=n_target, workspace=column_workspace)

    # 10. Allow non-melt densification
    densification_from_compaction = sum(dz)

    # Opt-in viscosity diagnostic — a per-step diagnostic, not column state, so it is allocated
    # here rather than carried in `state`. Allocated after the grid controllers (step 9), and
    # nothing below merges or splits cells, so these indices match the `dz`/`density` profiles
    # this timestep reports.
    #
    # "Off" is empty `NO_VISCOSITY`, not `nothing`: `viscosity` lands in `flux`, and a
    # `Union{Nothing,Vector{Float64}}` field would make `flux`'s type depend on
    # `output_viscosity`, boxing its ~40 `Float64` fields every timestep (measured: 42 allocations
    # per step). Length, not type, is what says "not requested".
    viscosity = mp.output_viscosity ? Vector{Float64}(undef, length(dz)) : NO_VISCOSITY

    dz, density = calculate_density(temperature, dz, density, grain_radius, water, cfs, mp;
        grain_dendricity, grain_sphericity, viscosity)

    densification_from_compaction = densification_from_compaction - sum(dz)

    # 11. Thin the column by the horizontal (ice-dynamic) strain rate, at constant density.
    # After `calculate_density`, which rebuilds `dz` from the conserved cell mass and would
    # overwrite this; before `trim_bottom!`, which restores the fixed column depth. The mass it
    # exports leaves laterally, not basally, so it is reported separately. A no-op at the
    # `horizontal_strain_rate = 0.0` default.
    cols = column_state(temperature, dz, density, water, grain_radius,
        grain_dendricity, grain_sphericity, age)
    mass_lateral, strain_energy = apply_horizontal_strain!(cols, cfs.dt, mp)
    E_added += strain_energy

    # 12. Pin the total column depth to z_target. This runs last because
    # `calculate_density` rebuilds `dz` from the conserved cell mass and, with step 11, is
    # the final step to change column thickness. The adjustment is the model's only basal
    # mass and energy flux, and is signed: negative `mass_added` under accumulation (material
    # leaves through the base), positive under net ablation (basal accretion).
    mass_added, trim_energy = trim_bottom!(cols, z_target, mp)
    E_added += trim_energy

    # Restore the cell count, which `trim_bottom!` can lower: an accumulation step thicker than the
    # deepest cell exports whole cells through the base, and how thick that step is only becomes
    # known after `calculate_density` and the strain thinning, so step 9 cannot anticipate it.
    # Splitting is exact in depth — `dz[i] /= 2` and duplicating the slot sums back to `dz[i]` in
    # binary floating point — so the count is repaired without disturbing the depth pinned above.
    # Free when nothing was consumed: the controller returns on its first comparison.
    enforce_column_length!(cols, n_target, mp; workspace=column_workspace)

    # 13. Ice-slab diagnostics. Taken here, after every step that can change `dz` or
    # `density`, so they describe the same column the profile output records at this
    # timestamp — a user recomputing slab depth from the output `dz`/`density` gets this
    # answer. Run unconditionally: a slab is a property of the column, present whether or not
    # anything melted this timestep. Read-only, and absent from the budgets above.
    ice_slab_thickness, ice_slab_depth = ice_slab_diagnostics(dz, density, mp)
    aquifer_thickness, aquifer_depth = aquifer_diagnostics(dz, density, water, mp)

    if verbose
        dt = cfs.dt
        M = dz .* density
        M_total_final = sum(M) + sum(water)
        M_delta = M_total_final - M_total_initial + runoff - cfs.precipitation -
                  evaporation_condensation - mass_added - mass_lateral -
                  mass_blowing_snow

        if abs(M_delta) > M_TOLERANCE
            error("total system mass not conserved: M_delta = $(M_delta)")
        end

        longwave_net = cfs.longwave_downward - longwave_upward
        E_snow = (cfs.precipitation - rain) * specific_enthalpy(mp, cfs.temperature_air)
        E_rain = rain * specific_enthalpy_water(mp, cfs.temperature_air)
        E_runoff = runoff * specific_enthalpy_water(mp)
        E_thermal = column_enthalpy(mp, M, temperature)
        E_water = sum(water) * specific_enthalpy_water(mp)
        E_shortwave = shortwave_net * dt
        E_longwave = longwave_net * dt
        E_thf = (heat_flux_sensible + heat_flux_latent) * dt
        E_basal = heat_flux_basal * dt

        E_total_final = E_thermal + E_water + E_runoff
        E_used = E_total_final - E_total_initial
        E_supplied = E_shortwave + E_longwave + E_thf + E_snow + E_rain + E_basal +
                     E_evaporation_condensation + E_added + E_blowing_snow
        E_delta = E_used - E_supplied

        if abs(E_delta) > energy_tolerance(E_total_initial)
            error("total system energy not conserved: E_delta = $(E_delta)")
        end

        if abs(temperature[end] - T_bottom) > T_BOTTOM_TOLERANCE
            error("temperature of bottom grid cell changed")
        end

        for (name, a) in pairs(cols)
            if length(a) != n_target
                error("state field $(name) has length $(length(a)), expected $(n_target)")
            end
        end
    end

    # Grid invariants. Output arrays are sized exactly `n_target` rows on the strength of these
    # two equalities, so a violation is a wrong-shaped write. Checked every timestep, not just
    # on output steps: under `output_frequency=:last` (which `gemb_spinup` forces) the latter
    # would hide a violation for an entire spinup cycle.
    if length(dz) != n_target
        error("gemb_core: column length not conserved: $(length(dz)) != $(n_target). " *
              "The count controller (enforce_column_length!) did not restore the fixed " *
              "cell count.")
    end
    z_err = sum(dz) - z_target
    if abs(z_err) > Z_TOLERANCE
        error("gemb_core: column depth not conserved: Σdz - z_target = $(z_err) m " *
              "(target $(z_target) m). The mass controller (trim_bottom!) did not pin the " *
              "fixed column depth.")
    end

    new_state = (
        temperature=temperature,
        dz=dz,
        density=density,
        water=water,
        grain_radius=grain_radius,
        grain_dendricity=grain_dendricity,
        grain_sphericity=grain_sphericity,
        age=age,
        evaporation_condensation=evaporation_condensation,
        melt_surface=melt_surface,
    )

    flux = (
        albedo_broadband=effective_albedo(albedo_broadband, albedo_diffuse, cfs),
        shortwave_net=shortwave_net,
        heat_flux_sensible=heat_flux_sensible,
        heat_flux_latent=heat_flux_latent,
        heat_flux_basal=heat_flux_basal,
        longwave_upward=longwave_upward,
        rain=rain,
        melt=melt,
        runoff=runoff,
        refreeze=refreeze,
        mass_added=mass_added,
        mass_lateral=mass_lateral,
        mass_blowing_snow=mass_blowing_snow,
        E_added=E_added,
        densification_from_compaction=densification_from_compaction,
        densification_from_melt=densification_from_melt,
        percolation_depth=percolation_depth,
        ice_slab_thickness=ice_slab_thickness,
        ice_slab_depth=ice_slab_depth,
        aquifer_thickness=aquifer_thickness,
        aquifer_depth=aquifer_depth,
        # Part of `refreeze`: the water the thermal solve froze itself (`_wet_cell_handoff!`), and
        # the number of cells it froze in.
        refreeze_thermal=refreeze_thermal,
        refreeze_thermal_count=refreeze_thermal_count,
        viscosity=viscosity,
    )

    return new_state, flux
end
