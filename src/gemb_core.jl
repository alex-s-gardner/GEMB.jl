"""
    gemb_core(temperature, dz, density, water, grain_radius, grain_dendricity,
              grain_sphericity, albedo, albedo_diffuse, evaporation_condensation,
              melt_surface, cfs::ClimateForcingStep, mp::ModelParameters, verbose::Bool)

Perform a single time-step of the GEMB model.
Matches MATLAB's `gemb_core.m`.

Returns a tuple of all updated state vectors and diagnostic scalars.
"""
function gemb_core(temperature, dz, density, water, grain_radius, grain_dendricity,
    grain_sphericity, albedo, albedo_diffuse, evaporation_condensation,
    melt_surface, cfs::ClimateForcingStep, mp::ModelParameters, verbose::Bool)

    if verbose
        M = dz .* density
        M_total_initial = sum(M) + sum(water)
        E_total_initial = sum(M .* temperature * C_ICE) +
                          sum(water .* (LF + CtoK * C_ICE))
        T_bottom = temperature[end]
    end

    # 1. Snow grain metamorphism
    grain_radius, grain_dendricity, grain_sphericity =
        calculate_grain_size(temperature, dz, density, water, grain_radius,
            grain_dendricity, grain_sphericity, cfs, mp)

    # 2. Calculate snow, firn, and ice albedo
    albedo, albedo_diffuse =
        calculate_albedo(temperature, dz, density, water, grain_radius,
            albedo, albedo_diffuse, evaporation_condensation, melt_surface, cfs, mp)

    # 3. Determine distribution of absorbed SW radiation with depth
    shortwave_flux = calculate_shortwave_radiation(dz, density, grain_radius,
        albedo[1], albedo_diffuse[1], cfs, mp)

    # 4. Calculate net shortwave [W m-2]
    shortwave_net = sum(shortwave_flux)

    # 5. Calculate new temperature-depth profile and turbulent heat fluxes
    temperature, longwave_upward, heat_flux_sensible, heat_flux_latent, ghf, evaporation_condensation =
        calculate_temperature(temperature, dz, density, water[1], grain_radius,
            shortwave_flux, cfs, mp, verbose)

    # 6. Change in thickness of top cell due to evaporation/condensation
    dz = copy(dz)
    dz[1] = dz[1] + evaporation_condensation / density[1]

    if verbose
        E_evaporation_condensation = evaporation_condensation * temperature[1] * C_ICE
    end

    # 7. Add snow/rain to top grid cell
    temperature, dz, density, water, grain_radius, grain_dendricity,
        grain_sphericity, albedo, albedo_diffuse, rain =
        calculate_accumulation(temperature, dz, density, water, grain_radius,
            grain_dendricity, grain_sphericity, albedo, albedo_diffuse, cfs, mp, verbose)

    # 8. Melt and wet compaction
    densification_from_melt = sum(dz)

    temperature, dz, density, water, grain_radius, grain_dendricity,
        grain_sphericity, albedo, albedo_diffuse, melt, melt_surface, runoff, refreeze =
        calculate_melt(temperature, dz, density, water, grain_radius,
            grain_dendricity, grain_sphericity, albedo, albedo_diffuse, rain, mp, verbose)

    densification_from_melt = densification_from_melt - sum(dz)

    # 9. Manage the layering
    temperature, dz, density, water, grain_radius, grain_dendricity,
        grain_sphericity, albedo, albedo_diffuse, mass_added, E_added =
        manage_layers(temperature, dz, density, water, grain_radius,
            grain_dendricity, grain_sphericity, albedo, albedo_diffuse, mp, verbose)

    # 10. Allow non-melt densification
    densification_from_compaction = sum(dz)

    dz, density = calculate_density(temperature, dz, density, grain_radius, cfs, mp)

    densification_from_compaction = densification_from_compaction - sum(dz)

    if verbose
        dt = cfs.dt
        M = dz .* density
        M_total_final = sum(M) + sum(water)
        M_delta = M_total_final - M_total_initial + runoff - cfs.precipitation - evaporation_condensation - mass_added

        if abs(M_delta) > 1e-3
            error("total system mass not conserved: M_delta = $(M_delta)")
        end

        longwave_net = cfs.longwave_downward - longwave_upward
        E_snow = (cfs.precipitation - rain) * cfs.temperature_air * C_ICE
        E_rain = rain * (cfs.temperature_air * C_ICE + LF)
        E_runoff = runoff * (LF + CtoK * C_ICE)
        E_thermal = sum((dz .* density) .* temperature * C_ICE)
        E_water = sum(water .* (LF + CtoK * C_ICE))
        E_shortwave = shortwave_net * dt
        E_longwave = longwave_net * dt
        E_thf = (heat_flux_sensible + heat_flux_latent) * dt
        E_ghf = ghf * dt

        E_total_final = E_thermal + E_water + E_runoff
        E_used = E_total_final - E_total_initial
        E_supplied = E_shortwave + E_longwave + E_thf + E_snow + E_rain + E_ghf + E_evaporation_condensation + E_added
        E_delta = E_used - E_supplied

        if abs(E_delta) > 1e-3
            error("total system energy not conserved: E_delta = $(E_delta)")
        end

        if abs(temperature[end] - T_bottom) > 1e-3
            error("temperature of bottom grid cell changed")
        end
    end

    return (temperature, dz, density, water, grain_radius, grain_dendricity,
        grain_sphericity, albedo, albedo_diffuse, evaporation_condensation,
        melt_surface, shortwave_net, heat_flux_sensible, heat_flux_latent,
        longwave_upward, rain, melt, runoff, refreeze, mass_added, E_added,
        densification_from_compaction, densification_from_melt)
end
