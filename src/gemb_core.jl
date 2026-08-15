"""
    gemb_core(state, cfs::ClimateForcingStep, mp::ModelParameters, verbose::Bool;
              n_target=length(state.dz), z_target=sum(state.dz))

Perform a single time-step of the GEMB model.
Matches MATLAB's `gemb_core.m`.

The column is returned with exactly `n_target` cells summing to exactly `z_target` metres,
so both are invariant across the whole run (see `grid_ops.jl` for the two controllers that
enforce them). Cell count is restored by [`manage_layer_thickness`](@ref) at step 9; column depth is
pinned at step 11, *after* [`calculate_density`](@ref), which is the last thing in the
timestep to change `dz`.

Returns `(state, flux)` where `state` is the updated column state and
`flux` contains energy/mass budget terms for output accumulation.
"""
function gemb_core(state, cfs::ClimateForcingStep, mp::ModelParameters, verbose::Bool;
    n_target::Int=length(state.dz), z_target::Float64=sum(state.dz))
    # Destructure state - arrays are mutated in-place by physics functions.
    # This is safe because gemb_driver rebinds state to new_state after each call.
    temperature = state.temperature
    dz = state.dz
    density = state.density
    water = state.water
    grain_radius = state.grain_radius
    grain_dendricity = state.grain_dendricity
    grain_sphericity = state.grain_sphericity
    age = state.age
    evaporation_condensation = state.evaporation_condensation
    melt_surface = state.melt_surface

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
    #    column, not carried state: `albedo_broadband` is returned in `flux` for output and
    #    is not read again.
    albedo_broadband, albedo_diffuse =
        calculate_albedo(dz, density, water, grain_radius, melt_surface, cfs, mp)

    # 3. Determine distribution of absorbed SW radiation with depth
    shortwave_flux = calculate_shortwave_radiation(dz, density, grain_radius,
        albedo_broadband, albedo_diffuse, cfs, mp)

    # 4. Calculate net shortwave [W m-2]
    shortwave_net = sum(shortwave_flux)

    # 5. Calculate new temperature-depth profile and turbulent heat flux
    temperature, longwave_upward, heat_flux_sensible, heat_flux_latent, ghf, evaporation_condensation =
        calculate_temperature(temperature, dz, density, water[1], grain_radius,
            shortwave_flux, cfs, mp, verbose)

    # 6. Change in thickness of top cell due to evaporation/condensation.
    #    Deposition (positive) is new mass at age 0, so it dilutes the surface cell's age;
    #    sublimation (negative) removes a fraction of the cell and is age-neutral, which
    #    `dilute_age` handles by returning the age unchanged for a non-positive input.
    S_surface = dz[1] * density[1] + water[1]
    dz[1] = dz[1] + evaporation_condensation / density[1]
    age[1] = dilute_age(age[1], S_surface, evaporation_condensation)

    if verbose
        E_evaporation_condensation = evaporation_condensation * specific_enthalpy(mp, temperature[1])
    end

    # 7. Add snow/rain to top grid cell
    temperature, dz, density, water, grain_radius, grain_dendricity,
        grain_sphericity, age, rain =
        calculate_accumulation(temperature, dz, density, water, grain_radius,
            grain_dendricity, grain_sphericity, age, cfs, mp, verbose)

    # 8. Melt and wet compaction
    densification_from_melt = sum(dz)

    temperature, dz, density, water, grain_radius, grain_dendricity,
        grain_sphericity, age, melt, melt_surface, runoff, refreeze,
        percolation_depth =
        calculate_melt(temperature, dz, density, water, grain_radius,
            grain_dendricity, grain_sphericity, age, rain, mp, verbose)

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
            n_target=n_target)

    # 10. Allow non-melt densification
    densification_from_compaction = sum(dz)

    dz, density = calculate_density(temperature, dz, density, grain_radius, water, cfs, mp)

    densification_from_compaction = densification_from_compaction - sum(dz)

    # 11. Thin the column by the horizontal (ice-dynamic) strain rate, at constant density.
    # After `calculate_density`, which rebuilds `dz` from the conserved cell mass and would
    # overwrite this; before `trim_bottom!`, which restores the fixed column depth. The mass
    # this exports leaves laterally, not through the base, and is reported separately.
    # A no-op when `horizontal_strain_rate == 0` (the default).
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
                  evaporation_condensation - mass_added - mass_lateral

        if abs(M_delta) > 1e-3
            error("total system mass not conserved: M_delta = $(M_delta)")
        end

        longwave_net = cfs.longwave_downward - longwave_upward
        E_snow = (cfs.precipitation - rain) * specific_enthalpy(mp, cfs.temperature_air)
        E_rain = rain * (specific_enthalpy(mp, cfs.temperature_air) + LF)
        E_runoff = runoff * specific_enthalpy_water(mp)
        E_thermal = column_enthalpy(mp, M, temperature)
        E_water = sum(water) * specific_enthalpy_water(mp)
        E_shortwave = shortwave_net * dt
        E_longwave = longwave_net * dt
        E_thf = (heat_flux_sensible + heat_flux_latent) * dt
        E_ghf = ghf * dt

        E_total_final = E_thermal + E_water + E_runoff
        E_used = E_total_final - E_total_initial
        E_supplied = E_shortwave + E_longwave + E_thf + E_snow + E_rain + E_ghf + E_evaporation_condensation + E_added
        E_delta = E_used - E_supplied

        if abs(E_delta) > energy_tolerance(E_total_initial)
            error("total system energy not conserved: E_delta = $(E_delta)")
        end

        if abs(temperature[end] - T_bottom) > 1e-3
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
    if abs(z_err) > 1e-9
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
        albedo_broadband=albedo_broadband,
        shortwave_net=shortwave_net,
        heat_flux_sensible=heat_flux_sensible,
        heat_flux_latent=heat_flux_latent,
        longwave_upward=longwave_upward,
        rain=rain,
        melt=melt,
        runoff=runoff,
        refreeze=refreeze,
        mass_added=mass_added,
        mass_lateral=mass_lateral,
        E_added=E_added,
        densification_from_compaction=densification_from_compaction,
        densification_from_melt=densification_from_melt,
        # Diagnostics. Read-only: deliberately absent from the conservation checks above
        # because they move no mass and carry no energy. `percolation_depth` is necessarily
        # mid-timestep (it describes an event during percolation); the slab terms are taken
        # at step 13, on the final column.
        percolation_depth=percolation_depth,
        ice_slab_thickness=ice_slab_thickness,
        ice_slab_depth=ice_slab_depth,
        aquifer_thickness=aquifer_thickness,
        aquifer_depth=aquifer_depth,
    )

    return new_state, flux
end
