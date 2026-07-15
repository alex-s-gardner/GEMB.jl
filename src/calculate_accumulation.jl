"""
    calculate_accumulation(temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, albedo, albedo_diffuse, cfs::ClimateForcingStep, mp::ModelParameters, verbose::Bool)

Add precipitation and deposition to the model column.

Precipitation is classified as snow or rain based on `mp.rain_temperature_threshold`.
Snow is added as a new layer (if depth > dzmin) or merged into the top cell.
Rain is added by increasing the mass and temperature of the top grid cell,
with temperature adjusted to account for latent heat of fusion.

Returns `(temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, albedo, albedo_diffuse, rain)`.
Arrays may grow (prepend) when snow depth > dzmin.
"""
function calculate_accumulation(temperature::Vector{Float64}, dz::Vector{Float64},
    density::Vector{Float64}, water::Vector{Float64},
    grain_radius::Vector{Float64}, grain_dendricity::Vector{Float64},
    grain_sphericity::Vector{Float64}, albedo::Vector{Float64},
    albedo_diffuse::Vector{Float64},
    cfs::ClimateForcingStep, mp::ModelParameters, verbose::Bool)

    # Note: arrays are modified in-place or extended via pushfirst! when new snow is added.

    # Define tolerances
    T_tolerance = 1e-10
    d_tolerance = 1e-11
    gdn_tolerance = 1e-10

    # Specify constants
    re_new_snow = 0.05       # new snow grain size [mm]
    gdn_new_snow = 1.0       # new snow dendricity
    gsp_new_snow = 0.5       # new snow sphericity
    rain = 0.0               # rainfall [mm w.e. or kg m^-2]

    if verbose
        M = dz .* density
        M_total_initial = sum(M)
        E_total_initial = sum(M .* temperature .* C_ICE) +
            sum(water .* (LF + CtoK * C_ICE))
    end

    # Density of fresh snow [kg m-3]
    density_new_snow = 0.0
    if mp.new_snow_method == Symbol("150kgm2")
        density_new_snow = 150.0
    elseif mp.new_snow_method == Symbol("350kgm2")
        density_new_snow = 350.0
    elseif mp.new_snow_method == :Fausto
        density_new_snow = 315.0
        # From Vionnet et al., 2012 (Crocus)
        gdn_new_snow = min(max(1.29 - 0.17 * cfs.wind_speed, 0.20), 1.0)
        gsp_new_snow = min(max(0.08 * cfs.wind_speed + 0.38, 0.5), 0.9)
        re_new_snow = max(1e-1 * (gdn_new_snow / 0.99 + (1.0 - 1.0 * gdn_new_snow / 0.99) * (gsp_new_snow / 0.99 * 3.0 + (1.0 - gsp_new_snow / 0.99) * 4.0)) / 2.0, gdn_tolerance)
    elseif mp.new_snow_method == :Kaspers
        density_new_snow = (7.36e-2 + 1.06e-3 * min(cfs.temperature_air_mean, CtoK - T_tolerance) + 6.69e-2 * cfs.precipitation_mean / 1000.0 + 4.77e-3 * cfs.wind_speed_mean) * 1000.0
    elseif mp.new_snow_method == :KuipersMunneke
        density_new_snow = 481.0 + 4.834 * (cfs.temperature_air_mean - CtoK)
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
                temperature = vcat([cfs.temperature_air], temperature)
                dz = vcat([z_snow], dz)
                density = vcat([density_new_snow], density)
                water = vcat([0.0], water)
                albedo = vcat([mp.albedo_snow], albedo)
                albedo_diffuse = vcat([mp.albedo_snow], albedo_diffuse)
                grain_radius = vcat([refall], grain_radius)
                grain_dendricity = vcat([dfall], grain_dendricity)
                grain_sphericity = vcat([sfall], grain_sphericity)
            else
                # if snow depth is less than specified minimum dz
                M_surface_new = M_surface + cfs.precipitation
                dz[1] = dz[1] + cfs.precipitation / density_new_snow
                density[1] = M_surface_new / dz[1]

                # adjust temperature (assume precipitation is same temp as air)
                temperature[1] = ((cfs.temperature_air * cfs.precipitation) + (temperature[1] * M_surface)) / M_surface_new

                # adjust albedo
                if mp.albedo_method != "150kgm2"
                    albedo[1] = (mp.albedo_snow * cfs.precipitation + albedo[1] * M_surface) / M_surface_new
                end

                grain_dendricity[1] = dfall
                grain_sphericity[1] = sfall
                grain_radius[1] = max(0.1 * (grain_dendricity[1] / 0.99 + (1.0 - 1.0 * grain_dendricity[1] / 0.99) * (grain_sphericity[1] / 0.99 * 3.0 + (1.0 - grain_sphericity[1] / 0.99) * 4.0)) / 2, gdn_tolerance)
            end

        else
            # rain
            # grid cell adjusted mass
            M_surface_new = M_surface + cfs.precipitation

            # adjust temperature (liquid: must account for latent heat of fusion)
            temperature[1] = ((cfs.precipitation * (cfs.temperature_air + LF / C_ICE)) +
                (temperature[1] * M_surface)) / M_surface_new

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

            E_total_final = sum(M .* temperature .* C_ICE) +
                sum(water .* (LF + CtoK * C_ICE))

            E_snow = ((cfs.precipitation - rain) * (cfs.temperature_air) * C_ICE)
            E_rain = (rain * (cfs.temperature_air * C_ICE + LF))

            E_delta = E_total_final - E_total_initial - E_snow - E_rain

            if (abs(M_delta) > 1e-3) || (abs(E_delta) > 1e-3)
                error("Mass and/or energy are not conserved:\n M_delta: $(M_delta) E_delta: $(E_delta)\n")
            end
        end
    end

    return temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, albedo, albedo_diffuse, rain
end
