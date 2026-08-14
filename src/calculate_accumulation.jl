"""
    fresh_snow_density(mp::ModelParameters, T_air_mean, precip_mean, wind_speed_mean)

Return the density of fresh snow [kg m-3] for the configured `mp.new_snow_method`.

Depends only on climatological means (temperature [K], accumulation [kg m-2 yr-1],
wind speed [m s-1]), so it is shared by `calculate_accumulation` (per timestep) and
`initialize_profile` (as the surface density ρ₀ of the steady-state firn profile).

Methods: `150kgm2` → 150; `350kgm2` → 350; `:Fausto` → 315; `:Kaspers` and
`:KuipersMunneke` → temperature/accumulation/wind-dependent fits.
"""
function fresh_snow_density(mp::ModelParameters, T_air_mean::Real,
    precip_mean::Real, wind_speed_mean::Real)
    T_tolerance = 1e-10
    if mp.new_snow_method == Symbol("150kgm2")
        return 150.0
    elseif mp.new_snow_method == Symbol("350kgm2")
        return 350.0
    elseif mp.new_snow_method == :Fausto
        return 315.0
    elseif mp.new_snow_method == :Kaspers
        return (7.36e-2 + 1.06e-3 * min(T_air_mean, CtoK - T_tolerance) +
                6.69e-2 * precip_mean / 1000.0 + 4.77e-3 * wind_speed_mean) * 1000.0
    elseif mp.new_snow_method == :KuipersMunneke
        return 481.0 + 4.834 * (T_air_mean - CtoK)
    end
    return 0.0
end

"""
    calculate_accumulation(temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, cfs::ClimateForcingStep, mp::ModelParameters, verbose::Bool)

Add precipitation and deposition to the model column.

Precipitation is classified as snow or rain based on `mp.rain_temperature_threshold`.
Snow is added as a new layer (if depth > dzmin) or merged into the top cell.
Rain is added by increasing the mass and temperature of the top grid cell,
with temperature adjusted to account for latent heat of fusion.

Returns `(temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, rain)`.
Arrays grow by one cell (via [`open_slot!`](@ref)) when snow depth > dzmin; the extra cell
is reclaimed later in the timestep by [`manage_layer_thickness`](@ref)'s count controller.

Albedo is absent: it is diagnosed from the column at the top of the *next* timestep (see
[`calculate_albedo`](@ref)), which reads the fresh-snow grain size and density this
function sets, so new snow brightens the surface without an albedo being stored here.
"""
function calculate_accumulation(temperature::Vector{Float64}, dz::Vector{Float64},
    density::Vector{Float64}, water::Vector{Float64},
    grain_radius::Vector{Float64}, grain_dendricity::Vector{Float64},
    grain_sphericity::Vector{Float64},
    cfs::ClimateForcingStep, mp::ModelParameters, verbose::Bool)

    # Note: arrays are modified in place, and grow by one cell via `open_slot!` when new
    # snow is deep enough to warrant its own layer.

    # Define tolerances
    T_tolerance = 1e-10
    d_tolerance = 1e-11
    gdn_tolerance = 1e-10

    # Specify constants (shared with the steady-state initial guess, so a freshly
    # initialized column and freshly fallen snow agree by construction)
    re_new_snow = RE_NEW_SNOW    # new snow grain size [mm]
    gdn_new_snow = GDN_NEW_SNOW  # new snow dendricity
    gsp_new_snow = GSP_NEW_SNOW  # new snow sphericity
    rain = 0.0               # rainfall [mm w.e. or kg m^-2]

    if verbose
        M = dz .* density
        M_total_initial = sum(M)
        E_total_initial = column_enthalpy(mp, M, temperature, water)
    end

    # Density of fresh snow [kg m-3] (shared with initialize_profile)
    density_new_snow = fresh_snow_density(mp, cfs.temperature_air_mean,
        cfs.precipitation_mean, cfs.wind_speed_mean)
    if mp.new_snow_method == :Fausto
        # From Vionnet et al., 2012 (Crocus): wind-dependent grain properties.
        gdn_new_snow = min(max(1.29 - 0.17 * cfs.wind_speed, 0.20), 1.0)
        gsp_new_snow = min(max(0.08 * cfs.wind_speed + 0.38, 0.5), 0.9)
        re_new_snow = max(1e-1 * (gdn_new_snow / 0.99 + (1.0 - 1.0 * gdn_new_snow / 0.99) * (gsp_new_snow / 0.99 * 3.0 + (1.0 - gsp_new_snow / 0.99) * 4.0)) / 2.0, gdn_tolerance)
    end

    M_surface = dz[1] * density[1]

    if cfs.precipitation > 0
        # if snow
        if cfs.temperature_air <= (mp.rain_temperature_threshold + T_tolerance)
            z_snow = cfs.precipitation / density_new_snow          # depth of snow
            dfall = gdn_new_snow
            sfall = gsp_new_snow
            refall = re_new_snow

            # if snow depth is greater than specified min dz, new cell created
            if z_snow > mp.column_dzmin + d_tolerance
                cols = column_state(temperature, dz, density, water, grain_radius,
                    grain_dendricity, grain_sphericity)
                open_slot!(cols, 1)
                @inbounds begin
                    temperature[1] = cfs.temperature_air
                    dz[1] = z_snow
                    density[1] = density_new_snow
                    water[1] = 0.0
                    grain_radius[1] = refall
                    grain_dendricity[1] = dfall
                    grain_sphericity[1] = sfall
                end
            else
                # if snow depth is less than specified minimum dz
                M_surface_new = M_surface + cfs.precipitation
                dz[1] = dz[1] + cfs.precipitation / density_new_snow
                density[1] = M_surface_new / dz[1]

                # adjust temperature (assume precipitation is same temp as air)
                temperature[1] = mix_temperature(mp, cfs.temperature_air, cfs.precipitation,
                    temperature[1], M_surface)

                grain_dendricity[1] = dfall
                grain_sphericity[1] = sfall
                grain_radius[1] = max(0.1 * (grain_dendricity[1] / 0.99 + (1.0 - 1.0 * grain_dendricity[1] / 0.99) * (grain_sphericity[1] / 0.99 * 3.0 + (1.0 - grain_sphericity[1] / 0.99) * 4.0)) / 2, gdn_tolerance)
            end

        else
            # rain
            # grid cell adjusted mass
            M_surface_new = M_surface + cfs.precipitation

            # adjust temperature (liquid: must account for latent heat of fusion)
            temperature[1] = mix_temperature_liquid(mp, cfs.temperature_air, cfs.precipitation,
                temperature[1], M_surface)

            # adjust grid cell density
            density[1] = M_surface_new / dz[1]

            # if density > the density of ice
            if density[1] > mp.density_ice - d_tolerance
                density[1] = mp.density_ice
                dz[1] = M_surface_new / density[1]
            end

            rain = cfs.precipitation
        end

        if verbose
            # Check for conservation of mass
            M = dz .* density
            M_total_final = sum(M)
            M_delta = M_total_final - M_total_initial - cfs.precipitation

            E_total_final = column_enthalpy(mp, M, temperature, water)

            E_snow = (cfs.precipitation - rain) * specific_enthalpy(mp, cfs.temperature_air)
            E_rain = rain * (specific_enthalpy(mp, cfs.temperature_air) + LF)

            E_delta = E_total_final - E_total_initial - E_snow - E_rain

            E_tol = energy_tolerance(E_total_initial)

            if (abs(M_delta) > 1e-3) || (abs(E_delta) > E_tol)
                error("Mass and/or energy are not conserved:\n M_delta: $(M_delta) E_delta: $(E_delta)\n")
            end
        end
    end

    return temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, rain
end
