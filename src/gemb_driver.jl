"""
    gemb(profile::DimStack, cf::ClimateForcing, mp::ModelParameters; verbose=false)

Run the Glacier Energy and Mass Balance (GEMB) model.
Matches MATLAB's `gemb.m`.

Returns a DimStack containing time series of surface fluxes (monolevel)
and vertical profiles at the specified output frequency.
"""
function gemb(profile::DimStack, cf::ClimateForcing, mp::ModelParameters; verbose::Bool=false)
    # Extract column state from DimStack
    temperature = parent(profile[:temperature])
    dz = parent(profile[:dz])
    density = parent(profile[:density])
    water = parent(profile[:water])
    grain_radius = parent(profile[:grain_radius])
    grain_dendricity = parent(profile[:grain_dendricity])
    grain_sphericity = parent(profile[:grain_sphericity])
    albedo = parent(profile[:albedo])
    albedo_diffuse = parent(profile[:albedo_diffuse])

    # Get time information
    time_dim = dims(cf.temperature_air, Ti)
    times = Vector{DateTime}(time_dim.val)
    n_steps = length(times)
    dt = Dates.value(times[2] - times[1]) / 1000.0  # milliseconds to seconds

    # Round dt if not exact integer
    if rem(dt, 1) != 0
        dt = round(dt)
    end
    dt_int = round(Int, dt)

    # Pre-compute dt_divisors for thermal sub-stepping
    mp_with_divisors = ModelParameters(;
        (field => getfield(mp, field) for field in fieldnames(ModelParameters) if field != :dt_divisors)...,
        dt_divisors=fast_divisors(dt_int * 10000) ./ 10000
    )

    # Initialize output structure
    output_index = _compute_output_index(times, mp.output_frequency)
    n_outputs = sum(output_index)
    column_length = length(dz)
    profile_size = column_length + mp.output_padding

    # Pre-allocate output arrays
    out_time = times[output_index]

    # Monolevel outputs
    out_melt = fill(NaN, n_outputs)
    out_runoff = fill(NaN, n_outputs)
    out_refreeze = fill(NaN, n_outputs)
    out_evaporation_condensation = fill(NaN, n_outputs)
    out_shortwave_net = fill(NaN, n_outputs)
    out_longwave_net = fill(NaN, n_outputs)
    out_heat_flux_sensible = fill(NaN, n_outputs)
    out_heat_flux_latent = fill(NaN, n_outputs)
    out_albedo_surface = fill(NaN, n_outputs)
    out_densification_from_compaction = fill(NaN, n_outputs)
    out_densification_from_melt = fill(NaN, n_outputs)
    out_thickness_cumulative = fill(NaN, n_outputs)
    out_firn_air_content = fill(NaN, n_outputs)
    out_valid_profile_length = fill(0, n_outputs)

    # Forcing summary outputs (pre-computed from forcing data)
    out_temperature_air = fill(NaN, n_outputs)
    out_precipitation = fill(NaN, n_outputs)

    # Profile outputs
    out_temperature = fill(NaN, profile_size, n_outputs)
    out_dz = fill(NaN, profile_size, n_outputs)
    out_density = fill(NaN, profile_size, n_outputs)
    out_water = fill(NaN, profile_size, n_outputs)
    out_grain_radius = fill(NaN, profile_size, n_outputs)
    out_grain_dendricity = fill(NaN, profile_size, n_outputs)
    out_grain_sphericity = fill(NaN, profile_size, n_outputs)
    out_albedo = fill(NaN, profile_size, n_outputs)
    out_albedo_diffuse = fill(NaN, profile_size, n_outputs)

    # Initialize cumulative trackers
    cum_melt = 0.0
    cum_runoff = 0.0
    cum_refreeze = 0.0
    cum_ec = 0.0
    cum_rain = 0.0
    cum_mass_added = 0.0
    cum_shortwave_net = 0.0
    cum_longwave_net = 0.0
    cum_shf = 0.0
    cum_lhf = 0.0
    cum_albedo_surface = 0.0
    cum_densification_compaction = 0.0
    cum_densification_melt = 0.0
    cum_firn_air_content = 0.0
    cum_thickness = 0.0
    cum_temperature_air = 0.0
    cum_precipitation = 0.0
    cum_count = 0

    # Initialize state variables
    evaporation_condensation = 0.0
    melt_surface = 0.0
    thickness_added_total = 0.0

    # Output counter
    out_idx = 0

    # Main time loop
    for date_ind in 1:n_steps
        # Extract single timestep forcing
        cfs = _extract_forcing_step(date_ind, Float64(dt_int), cf, mp_with_divisors)

        # Run physics for single timestep
        temperature, dz, density, water, grain_radius, grain_dendricity,
            grain_sphericity, albedo, albedo_diffuse, evaporation_condensation,
            melt_surface, shortwave_net, heat_flux_sensible, heat_flux_latent,
            longwave_upward, rain, melt, runoff, refreeze, mass_added, _E_added,
            densification_from_compaction, densification_from_melt =
            gemb_core(temperature, dz, density, water, grain_radius, grain_dendricity,
                grain_sphericity, albedo, albedo_diffuse, evaporation_condensation,
                melt_surface, cfs, mp_with_divisors, verbose)

        # Calculate net longwave
        longwave_net = cfs.longwave_downward - longwave_upward

        # Sum total thickness
        thickness_added_total += mass_added / mp.density_ice

        # Accumulate outputs
        cum_melt += melt
        cum_runoff += runoff
        cum_refreeze += refreeze
        cum_ec += evaporation_condensation
        cum_rain += rain
        cum_mass_added += mass_added
        cum_shortwave_net += shortwave_net
        cum_longwave_net += longwave_net
        cum_shf += heat_flux_sensible
        cum_lhf += heat_flux_latent
        cum_albedo_surface += albedo[1]
        cum_densification_compaction += densification_from_compaction
        cum_densification_melt += densification_from_melt
        cum_firn_air_content += sum(dz .* (mp.density_ice .- min.(density, mp.density_ice))) / 1000
        cum_thickness += thickness_added_total
        cum_temperature_air += cfs.temperature_air
        cum_precipitation += cfs.precipitation
        cum_count += 1

        # Store output at designated intervals
        if output_index[date_ind]
            out_idx += 1

            # Cumulative variables
            out_melt[out_idx] = cum_melt
            out_runoff[out_idx] = cum_runoff
            out_refreeze[out_idx] = cum_refreeze
            out_evaporation_condensation[out_idx] = cum_ec
            out_densification_from_compaction[out_idx] = cum_densification_compaction
            out_densification_from_melt[out_idx] = cum_densification_melt

            # Averaged variables
            out_shortwave_net[out_idx] = cum_shortwave_net / cum_count
            out_longwave_net[out_idx] = cum_longwave_net / cum_count
            out_heat_flux_sensible[out_idx] = cum_shf / cum_count
            out_heat_flux_latent[out_idx] = cum_lhf / cum_count
            out_albedo_surface[out_idx] = cum_albedo_surface / cum_count
            out_firn_air_content[out_idx] = cum_firn_air_content / cum_count
            out_thickness_cumulative[out_idx] = cum_thickness / cum_count

            # Forcing summary
            out_temperature_air[out_idx] = cum_temperature_air / cum_count
            out_precipitation[out_idx] = cum_precipitation

            # Profile data (stored from bottom up, matching MATLAB convention)
            m = length(dz)
            out_valid_profile_length[out_idx] = m

            if m > profile_size
                error("Column length ($m) exceeds output array size ($profile_size). Increase output_padding.")
            end

            offset = profile_size - m
            out_temperature[(offset+1):profile_size, out_idx] = temperature
            out_dz[(offset+1):profile_size, out_idx] = dz
            out_density[(offset+1):profile_size, out_idx] = density
            out_water[(offset+1):profile_size, out_idx] = water
            out_grain_radius[(offset+1):profile_size, out_idx] = grain_radius
            out_grain_dendricity[(offset+1):profile_size, out_idx] = grain_dendricity
            out_grain_sphericity[(offset+1):profile_size, out_idx] = grain_sphericity
            out_albedo[(offset+1):profile_size, out_idx] = albedo
            out_albedo_diffuse[(offset+1):profile_size, out_idx] = albedo_diffuse

            # Reset accumulators
            cum_melt = 0.0
            cum_runoff = 0.0
            cum_refreeze = 0.0
            cum_ec = 0.0
            cum_rain = 0.0
            cum_mass_added = 0.0
            cum_shortwave_net = 0.0
            cum_longwave_net = 0.0
            cum_shf = 0.0
            cum_lhf = 0.0
            cum_albedo_surface = 0.0
            cum_densification_compaction = 0.0
            cum_densification_melt = 0.0
            cum_firn_air_content = 0.0
            cum_thickness = 0.0
            cum_temperature_air = 0.0
            cum_precipitation = 0.0
            cum_count = 0
        end
    end

    # Build output DimStack
    ti_dim = Ti(out_time)
    z_dim = Z(1:profile_size)

    return DimStack((
        # Monolevel
        melt=DimArray(out_melt, (ti_dim,)),
        runoff=DimArray(out_runoff, (ti_dim,)),
        refreeze=DimArray(out_refreeze, (ti_dim,)),
        evaporation_condensation=DimArray(out_evaporation_condensation, (ti_dim,)),
        shortwave_net=DimArray(out_shortwave_net, (ti_dim,)),
        longwave_net=DimArray(out_longwave_net, (ti_dim,)),
        heat_flux_sensible=DimArray(out_heat_flux_sensible, (ti_dim,)),
        heat_flux_latent=DimArray(out_heat_flux_latent, (ti_dim,)),
        albedo_surface=DimArray(out_albedo_surface, (ti_dim,)),
        densification_from_compaction=DimArray(out_densification_from_compaction, (ti_dim,)),
        densification_from_melt=DimArray(out_densification_from_melt, (ti_dim,)),
        thickness_cumulative=DimArray(out_thickness_cumulative, (ti_dim,)),
        firn_air_content=DimArray(out_firn_air_content, (ti_dim,)),
        valid_profile_length=DimArray(out_valid_profile_length, (ti_dim,)),
        temperature_air=DimArray(out_temperature_air, (ti_dim,)),
        precipitation=DimArray(out_precipitation, (ti_dim,)),
        # Profile
        temperature=DimArray(out_temperature, (z_dim, ti_dim)),
        dz=DimArray(out_dz, (z_dim, ti_dim)),
        density=DimArray(out_density, (z_dim, ti_dim)),
        water=DimArray(out_water, (z_dim, ti_dim)),
        grain_radius=DimArray(out_grain_radius, (z_dim, ti_dim)),
        grain_dendricity=DimArray(out_grain_dendricity, (z_dim, ti_dim)),
        grain_sphericity=DimArray(out_grain_sphericity, (z_dim, ti_dim)),
        albedo=DimArray(out_albedo, (z_dim, ti_dim)),
        albedo_diffuse=DimArray(out_albedo_diffuse, (z_dim, ti_dim)),
    ))
end

"""
Compute output indices based on output frequency.
"""
function _compute_output_index(times::Vector{DateTime}, frequency::String)
    n = length(times)
    if frequency == "all"
        return trues(n)
    elseif frequency == "last"
        idx = falses(n)
        idx[end] = true
        return idx
    elseif frequency == "daily"
        idx = falses(n)
        for i in 1:(n-1)
            if Dates.day(times[i]) != Dates.day(times[i+1]) ||
               Dates.month(times[i]) != Dates.month(times[i+1]) ||
               Dates.year(times[i]) != Dates.year(times[i+1])
                idx[i] = true
            end
        end
        idx[end] = true
        return idx
    elseif frequency == "monthly"
        idx = falses(n)
        for i in 1:(n-1)
            if Dates.month(times[i]) != Dates.month(times[i+1]) ||
               Dates.year(times[i]) != Dates.year(times[i+1])
                idx[i] = true
            end
        end
        idx[end] = true
        return idx
    else
        error("output_frequency must be one of: all, daily, monthly, last")
    end
end

"""
Extract a single timestep of forcing from ClimateForcing.
"""
function _extract_forcing_step(index::Int, dt::Float64, cf::ClimateForcing, mp::ModelParameters)
    return ClimateForcingStep(
        dt,
        parent(cf.temperature_air)[index],
        parent(cf.pressure_air)[index],
        parent(cf.precipitation)[index],
        parent(cf.wind_speed)[index],
        parent(cf.shortwave_downward)[index],
        parent(cf.longwave_downward)[index],
        parent(cf.vapor_pressure)[index],
        cf.temperature_air_mean,
        cf.wind_speed_mean,
        cf.precipitation_mean,
        cf.temperature_observation_height,
        cf.wind_observation_height,
        mp.black_carbon_snow,
        mp.black_carbon_ice,
        mp.cloud_optical_thickness,
        mp.solar_zenith_angle,
        mp.shortwave_downward_diffuse,
        mp.cloud_fraction,
    )
end
