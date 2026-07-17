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

    # Run the time loop through a function barrier. ClimateForcing's fields are
    # typed `::DimArray` (a UnionAll, not a concrete type), so indexing them
    # directly infers to `Any` and dispatches at runtime on every timestep. By
    # passing the underlying concrete arrays (via `parent`) as arguments, the
    # barrier specializes on their real element/dimension types and the inner
    # loop becomes fully type-stable.
    _gemb_time_loop!(output, state, model_parameters, mp, verbose,
        times, output_times, profile_size, Float64(dt_int),
        parent(climate_forcing.temperature_air),
        parent(climate_forcing.pressure_air),
        parent(climate_forcing.precipitation),
        parent(climate_forcing.wind_speed),
        parent(climate_forcing.shortwave_downward),
        parent(climate_forcing.longwave_downward),
        parent(climate_forcing.vapor_pressure),
        parent(climate_forcing.black_carbon_snow),
        parent(climate_forcing.black_carbon_ice),
        parent(climate_forcing.cloud_optical_thickness),
        parent(climate_forcing.solar_zenith_angle),
        parent(climate_forcing.shortwave_downward_diffuse),
        parent(climate_forcing.cloud_fraction),
        climate_forcing.temperature_air_mean,
        climate_forcing.wind_speed_mean,
        climate_forcing.precipitation_mean,
        climate_forcing.temperature_observation_height,
        climate_forcing.wind_observation_height)

    # Return the DimStack (already populated during time loop)
    return output
end

"""
    _gemb_time_loop!(output, state, model_parameters, mp, verbose, times,
                     output_times, profile_size, dt_f, <forcing arrays/scalars>)

Function-barrier inner loop for [`gemb`](@ref). Receives the forcing series as
concrete arrays (already unwrapped with `parent`) so the compiler specializes on
their true types, eliminating the per-timestep runtime dispatch that indexing
the `::DimArray`-typed `ClimateForcing` fields would otherwise incur. Mutates
`output` in place; numerically identical to the previous inline loop.
"""
function _gemb_time_loop!(output, state, model_parameters, mp, verbose::Bool,
    times::Vector{DateTime}, output_times, profile_size::Int, dt_f::Float64,
    f_temperature_air::AbstractVector, f_pressure_air::AbstractVector,
    f_precipitation::AbstractVector, f_wind_speed::AbstractVector,
    f_shortwave_downward::AbstractVector, f_longwave_downward::AbstractVector,
    f_vapor_pressure::AbstractVector, f_black_carbon_snow::AbstractVector,
    f_black_carbon_ice::AbstractVector, f_cloud_optical_thickness::AbstractVector,
    f_solar_zenith_angle::AbstractVector, f_shortwave_downward_diffuse::AbstractVector,
    f_cloud_fraction::AbstractVector,
    temperature_air_mean::Float64, wind_speed_mean::Float64,
    precipitation_mean::Float64, temperature_observation_height::Float64,
    wind_observation_height::Float64)

    # Extract concrete output arrays once (avoids per-write DimStack/At dispatch).
    out_melt = parent(output[:melt])
    out_runoff = parent(output[:runoff])
    out_refreeze = parent(output[:refreeze])
    out_ec = parent(output[:evaporation_condensation])
    out_dcomp = parent(output[:densification_from_compaction])
    out_dmelt = parent(output[:densification_from_melt])
    out_swnet = parent(output[:shortwave_net])
    out_lwnet = parent(output[:longwave_net])
    out_shf = parent(output[:heat_flux_sensible])
    out_lhf = parent(output[:heat_flux_latent])
    out_albsurf = parent(output[:albedo_surface])
    out_fac = parent(output[:firn_air_content])
    out_thick = parent(output[:thickness_cumulative])
    out_ta = parent(output[:temperature_air])
    out_precip = parent(output[:precipitation])
    out_vpl = parent(output[:valid_profile_length])
    out_temperature = parent(output[:temperature])
    out_dz = parent(output[:dz])
    out_density = parent(output[:density])
    out_water = parent(output[:water])
    out_grain_radius = parent(output[:grain_radius])
    out_grain_dendricity = parent(output[:grain_dendricity])
    out_grain_sphericity = parent(output[:grain_sphericity])
    out_albedo = parent(output[:albedo])
    out_albedo_diffuse = parent(output[:albedo_diffuse])

    density_ice = mp.density_ice

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

    # Output writes occur in chronological order, matching the sorted output
    # time axis, so a single advancing index tracks the output column.
    oi = 0

    for i in eachindex(times)
        # Construct ClimateForcingStep from concrete forcing arrays
        @inbounds forcing_step = ClimateForcingStep(
            dt_f,
            f_temperature_air[i],
            f_pressure_air[i],
            f_precipitation[i],
            f_wind_speed[i],
            f_shortwave_downward[i],
            f_longwave_downward[i],
            f_vapor_pressure[i],
            temperature_air_mean,
            wind_speed_mean,
            precipitation_mean,
            temperature_observation_height,
            wind_observation_height,
            f_black_carbon_snow[i],
            f_black_carbon_ice[i],
            f_cloud_optical_thickness[i],
            f_solar_zenith_angle[i],
            f_shortwave_downward_diffuse[i],
            f_cloud_fraction[i])

        # Run physics for single timestep
        state, flux = gemb_core(state, forcing_step, model_parameters, verbose)

        # Sum total thickness
        thickness_added_total += flux.mass_added / density_ice

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
        # Compute firn air content without temporary array allocation
        _fac = 0.0
        @inbounds for j in eachindex(state.dz)
            _fac += state.dz[j] * (density_ice - min(state.density[j], density_ice))
        end
        cum_firn_air_content += _fac / 1000
        cum_thickness += thickness_added_total
        cum_temperature_air += forcing_step.temperature_air
        cum_precipitation += forcing_step.precipitation
        cum_count += 1

        # Store output at designated intervals
        t = times[i]
        if t in output_times
            oi += 1

            @inbounds begin
                # Cumulative variables (1D arrays indexed by time)
                out_melt[oi] = cum_melt
                out_runoff[oi] = cum_runoff
                out_refreeze[oi] = cum_refreeze
                out_ec[oi] = cum_ec
                out_dcomp[oi] = cum_densification_compaction
                out_dmelt[oi] = cum_densification_melt

                # Averaged variables (division preserved for bit-identical results)
                out_swnet[oi] = cum_shortwave_net / cum_count
                out_lwnet[oi] = cum_longwave_net / cum_count
                out_shf[oi] = cum_shf / cum_count
                out_lhf[oi] = cum_lhf / cum_count
                out_albsurf[oi] = cum_albedo_surface / cum_count
                out_fac[oi] = cum_firn_air_content / cum_count
                out_thick[oi] = cum_thickness / cum_count

                # Forcing summary
                out_ta[oi] = cum_temperature_air / cum_count
                out_precip[oi] = cum_precipitation
            end

            # Profile data (stored from bottom up, matching MATLAB convention)
            m = length(state.dz)
            @inbounds out_vpl[oi] = m

            if m > profile_size
                error("Column length ($m) exceeds output array size ($profile_size). Increase output_padding.")
            end

            offset = profile_size - m
            @inbounds for k in 1:m
                r = offset + k
                out_temperature[r, oi] = state.temperature[k]
                out_dz[r, oi] = state.dz[k]
                out_density[r, oi] = state.density[k]
                out_water[r, oi] = state.water[k]
                out_grain_radius[r, oi] = state.grain_radius[k]
                out_grain_dendricity[r, oi] = state.grain_dendricity[k]
                out_grain_sphericity[r, oi] = state.grain_sphericity[k]
                out_albedo[r, oi] = state.albedo[k]
                out_albedo_diffuse[r, oi] = state.albedo_diffuse[k]
            end

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

