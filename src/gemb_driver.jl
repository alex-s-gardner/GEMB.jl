"""
    gemb(profile::DimStack, climate_forcing::ClimateForcing, mp::ModelParameters; verbose=false)

Run the Glacier Energy and Mass Balance (GEMB) model.
Matches MATLAB's `gemb.m`.

Returns a DimStack containing time series of surface flux (monolevel)
and vertical profiles at the specified output frequency.
"""
function gemb(profile::DimStack, climate_forcing::ClimateForcing, mp::ModelParameters; verbose::Bool=false)
    # Get time information
    time_dim = dims(climate_forcing.temperature_air, Ti)
    times = Vector{DateTime}(time_dim.val)
    dt_int = climate_forcing.time_step

    # Merge model parameters into climate_forcing as Fill arrays
    n = length(times)
    tdim = Ti(times)

    # TODO: time evolving parameters should all be handled in ClimateForcing, not ModelParameters, as inputs to gemb.
    climate_forcing = ClimateForcing(
        climate_forcing.temperature_air,
        climate_forcing.pressure_air,
        climate_forcing.precipitation,
        climate_forcing.wind_speed,
        climate_forcing.shortwave_downward,
        climate_forcing.longwave_downward,
        climate_forcing.vapor_pressure,
        DimArray(Fill(mp.black_carbon_snow, n), (tdim,)),
        DimArray(Fill(mp.black_carbon_ice, n), (tdim,)),
        DimArray(Fill(mp.cloud_optical_thickness, n), (tdim,)),
        DimArray(Fill(mp.solar_zenith_angle, n), (tdim,)),
        DimArray(Fill(mp.shortwave_downward_diffuse, n), (tdim,)),
        DimArray(Fill(mp.cloud_fraction, n), (tdim,)),
        climate_forcing.time_step,
        climate_forcing.temperature_air_mean,
        climate_forcing.wind_speed_mean,
        climate_forcing.precipitation_mean,
        climate_forcing.temperature_observation_height,
        climate_forcing.wind_observation_height,
    )

    # Pre-compute dt_divisors for thermal sub-stepping
    model_parameters = ModelParameters(;
        (field => getfield(mp, field) for field in fieldnames(ModelParameters) if field != :dt_divisors)...,
        dt_divisors=fast_divisors(dt_int * 10000) ./ 10000
    )

    # Initialize column state from profile
    state = (
        temperature = Vector{Float64}(profile[:temperature]),
        dz = Vector{Float64}(profile[:dz]),
        density = Vector{Float64}(profile[:density]),
        water = Vector{Float64}(profile[:water]),
        grain_radius = Vector{Float64}(profile[:grain_radius]),
        grain_dendricity = Vector{Float64}(profile[:grain_dendricity]),
        grain_sphericity = Vector{Float64}(profile[:grain_sphericity]),
        albedo = Vector{Float64}(profile[:albedo]),
        albedo_diffuse = Vector{Float64}(profile[:albedo_diffuse]),
        evaporation_condensation = 0.0,
        melt_surface = 0.0,
    )

    # Compute output times
    output_times = Set(_compute_output_times(times, mp.output_frequency))
    out_time = sort(collect(output_times))
    n_outputs = length(out_time)
    column_length = length(state.dz)
    profile_size = column_length + mp.output_padding

    # Create output time coordinate and dimensions
    ti_dim = Ti(out_time)
    z_dim = Z(1:profile_size)

    # Initialize output as DimStack with NaN-filled arrays
    output = DimStack((
        # Monolevel outputs
        melt=DimArray(fill(NaN, n_outputs), (ti_dim,)),
        runoff=DimArray(fill(NaN, n_outputs), (ti_dim,)),
        refreeze=DimArray(fill(NaN, n_outputs), (ti_dim,)),
        evaporation_condensation=DimArray(fill(NaN, n_outputs), (ti_dim,)),
        shortwave_net=DimArray(fill(NaN, n_outputs), (ti_dim,)),
        longwave_net=DimArray(fill(NaN, n_outputs), (ti_dim,)),
        heat_flux_sensible=DimArray(fill(NaN, n_outputs), (ti_dim,)),
        heat_flux_latent=DimArray(fill(NaN, n_outputs), (ti_dim,)),
        albedo_surface=DimArray(fill(NaN, n_outputs), (ti_dim,)),
        densification_from_compaction=DimArray(fill(NaN, n_outputs), (ti_dim,)),
        densification_from_melt=DimArray(fill(NaN, n_outputs), (ti_dim,)),
        thickness_cumulative=DimArray(fill(NaN, n_outputs), (ti_dim,)),
        firn_air_content=DimArray(fill(NaN, n_outputs), (ti_dim,)),
        valid_profile_length=DimArray(fill(0, n_outputs), (ti_dim,)),

        # Forcing summary outputs
        temperature_air=DimArray(fill(NaN, n_outputs), (ti_dim,)),
        precipitation=DimArray(fill(NaN, n_outputs), (ti_dim,)),

        # Profile outputs (2D: vertical × time)
        temperature=DimArray(fill(NaN, profile_size, n_outputs), (z_dim, ti_dim)),
        dz=DimArray(fill(NaN, profile_size, n_outputs), (z_dim, ti_dim)),
        density=DimArray(fill(NaN, profile_size, n_outputs), (z_dim, ti_dim)),
        water=DimArray(fill(NaN, profile_size, n_outputs), (z_dim, ti_dim)),
        grain_radius=DimArray(fill(NaN, profile_size, n_outputs), (z_dim, ti_dim)),
        grain_dendricity=DimArray(fill(NaN, profile_size, n_outputs), (z_dim, ti_dim)),
        grain_sphericity=DimArray(fill(NaN, profile_size, n_outputs), (z_dim, ti_dim)),
        albedo=DimArray(fill(NaN, profile_size, n_outputs), (z_dim, ti_dim)),
        albedo_diffuse=DimArray(fill(NaN, profile_size, n_outputs), (z_dim, ti_dim)),
    ))

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
    thickness_added_total = 0.0

    # Main time loop
    for t in times
        # Extract single timestep forcing
        forcing_step = climate_forcing[Ti=At(t)]

        # Run physics for single timestep
        state, flux = gemb_core(state, forcing_step, model_parameters, verbose)

        # Sum total thickness
        thickness_added_total += flux.mass_added / mp.density_ice

        # Accumulate outputs
        cum_melt += flux.melt
        cum_runoff += flux.runoff
        cum_refreeze += flux.refreeze
        cum_ec += state.evaporation_condensation
        cum_rain += flux.rain
        cum_mass_added += flux.mass_added
        cum_shortwave_net += flux.shortwave_net
        cum_longwave_net += forcing_step.longwave_downward - flux.longwave_upward
        cum_shf += flux.heat_flux_sensible
        cum_lhf += flux.heat_flux_latent
        cum_albedo_surface += state.albedo[1]
        cum_densification_compaction += flux.densification_from_compaction
        cum_densification_melt += flux.densification_from_melt
        cum_firn_air_content += sum(state.dz .* (mp.density_ice .- min.(state.density, mp.density_ice))) / 1000
        cum_thickness += thickness_added_total
        cum_temperature_air += forcing_step.temperature_air
        cum_precipitation += forcing_step.precipitation
        cum_count += 1

        # Store output at designated intervals
        if t in output_times
            # Cumulative variables (1D arrays indexed by time)
            output[:melt][Ti=At(t)] = cum_melt
            output[:runoff][Ti=At(t)] = cum_runoff
            output[:refreeze][Ti=At(t)] = cum_refreeze
            output[:evaporation_condensation][Ti=At(t)] = cum_ec
            output[:densification_from_compaction][Ti=At(t)] = cum_densification_compaction
            output[:densification_from_melt][Ti=At(t)] = cum_densification_melt

            # Averaged variables
            output[:shortwave_net][Ti=At(t)] = cum_shortwave_net / cum_count
            output[:longwave_net][Ti=At(t)] = cum_longwave_net / cum_count
            output[:heat_flux_sensible][Ti=At(t)] = cum_shf / cum_count
            output[:heat_flux_latent][Ti=At(t)] = cum_lhf / cum_count
            output[:albedo_surface][Ti=At(t)] = cum_albedo_surface / cum_count
            output[:firn_air_content][Ti=At(t)] = cum_firn_air_content / cum_count
            output[:thickness_cumulative][Ti=At(t)] = cum_thickness / cum_count

            # Forcing summary
            output[:temperature_air][Ti=At(t)] = cum_temperature_air / cum_count
            output[:precipitation][Ti=At(t)] = cum_precipitation

            # Profile data (stored from bottom up, matching MATLAB convention)
            m = length(state.dz)
            output[:valid_profile_length][Ti=At(t)] = m

            if m > profile_size
                error("Column length ($m) exceeds output array size ($profile_size). Increase output_padding.")
            end

            offset = profile_size - m
            output[:temperature][Z=(offset+1):profile_size, Ti=At(t)] = state.temperature
            output[:dz][Z=(offset+1):profile_size, Ti=At(t)] = state.dz
            output[:density][Z=(offset+1):profile_size, Ti=At(t)] = state.density
            output[:water][Z=(offset+1):profile_size, Ti=At(t)] = state.water
            output[:grain_radius][Z=(offset+1):profile_size, Ti=At(t)] = state.grain_radius
            output[:grain_dendricity][Z=(offset+1):profile_size, Ti=At(t)] = state.grain_dendricity
            output[:grain_sphericity][Z=(offset+1):profile_size, Ti=At(t)] = state.grain_sphericity
            output[:albedo][Z=(offset+1):profile_size, Ti=At(t)] = state.albedo
            output[:albedo_diffuse][Z=(offset+1):profile_size, Ti=At(t)] = state.albedo_diffuse

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

    # Return the DimStack (already populated during time loop)
    return output
end

"""
Compute output times based on output frequency.
Returns the last timestep of each day/month, all timesteps, or just the final one.
"""
function _compute_output_times(times::Vector{DateTime}, frequency::Symbol)
    frequency == :all && return times
    frequency == :last && return [times[end]]
    groupfn = frequency == :daily ? Date :
              frequency == :monthly ? (t -> (year(t), month(t))) :
              error("output_frequency must be one of: :all, :daily, :monthly, :last")
    out = [times[i] for i in 1:(length(times)-1) if groupfn(times[i]) != groupfn(times[i+1])]
    push!(out, times[end])
    return out
end

