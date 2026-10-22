# The `gemb` time loop of a `lower_boundary = :ground` column. `gemb` hands over to `_gemb_ground`
# before it builds anything of the glacier path.

# Monolevel layers of a ground run in output order: the glacier run's, then the ground-only ones.
const GROUND_MONOLEVEL_LAYERS = (:melt, :runoff, :refreeze, :evaporation_condensation,
    :shortwave_net, :longwave_net, :heat_flux_sensible, :heat_flux_latent, :heat_flux_basal,
    :albedo_broadband, :blowing_snow, :percolation_depth, :ice_slab_thickness, :ice_slab_depth,
    :temperature_air, :precipitation, :rain,
    :heat_flux_ground, :snow_depth, :snow_water_equivalent,
    :soil_temperature_005, :soil_temperature_010, :soil_temperature_020,
    :soil_temperature_050, :soil_temperature_100)
const GROUND_PROFILE_LAYERS = (RESTART_LAYERS..., :material)

"""
    _gemb_ground(profile, climate_forcing, mp, model_parameters; verbose, thermal_workspace,
                 column_workspace) -> DimStack

[`gemb`](@ref) for `lower_boundary = :ground`, stepping [`gemb_core_ground`](@ref).

The profile layers have the layout of `ground_column.jl`: snow cells, ground cells, empty slots,
with `material` saying which is which. Output layers:

- the glacier run's surface fluxes and mass terms, with the same meaning and interval reductions;
  `heat_flux_basal` is the flux the bottom ground cell supplies;
- `heat_flux_ground` [W m-2, positive upward, interval mean], the conductive flux from the ground
  into the snow (see [`_thermal_solve_ground!`](@ref));
- `snow_depth` [m] and `snow_water_equivalent` [kg m-2], instantaneous, of the snow alone;
- `soil_temperature_005` … `soil_temperature_100` [K], instantaneous, at 0.05, 0.1, 0.2, 0.5 and
  1.0 m below the ground surface ([`ground_soil_temperature`](@ref));
- the profile layers, `ground_snow_cells_max + ng` rows.

`evaporation_condensation` is the vapor exchanged with the snow and the ground together; the
whole-run mass budget, checked under `verbose`, counts only the part exchanged with the snow,
since the ground's water is held fixed.
"""
function _gemb_ground(profile::DimStack, climate_forcing::ClimateForcing, mp::ModelParameters,
    model_parameters::ModelParameters; verbose::Bool, thermal_workspace::ThermalWorkspace,
    column_workspace::ColumnWorkspace)

    haskey(profile, :material) || error("gemb: lower_boundary = :ground needs a profile with a " *
        "`material` layer; build it with `initialize_profile(mp, cf)` under that setting")
    g = ground_column(mp, ground_grid(mp.ground_depth))
    full = NamedTuple(l => Vector{Float64}(parent(profile[l])) for l in GROUND_PROFILE_LAYERS)
    length(full.dz) == ground_snow_slots(mp) + length(g.dz) || error("gemb: the profile has " *
        "$(length(full.dz)) rows; a ground column of these parameters has " *
        "$(ground_snow_slots(mp) + length(g.dz))")
    snow, temperature_ground, water_ground = _ground_state_split(full, g)
    abs(temperature_ground[end] - mp.ground_temperature_deep) <= T_TOLERANCE || error("gemb: " *
        "the profile's bottom ground cell is at $(temperature_ground[end]) K, but " *
        "ground_temperature_deep is $(mp.ground_temperature_deep) K")

    times = Vector{DateTime}(dims(climate_forcing.temperature_air, Ti).val)
    dt_f = Float64(climate_forcing.time_step)
    out_time = _compute_output_times(times, mp.output_frequency)
    output_times = Set(out_time)
    n_out = length(out_time)
    N = length(full.dz)
    ti_dim = Ti(out_time; metadata=cf_time_attributes())
    z_dim = Z(1:N; metadata=cf_layer_index_attributes())
    layers = merge(
        NamedTuple(k => DimArray(fill(NaN, n_out), (ti_dim,)) for k in GROUND_MONOLEVEL_LAYERS),
        NamedTuple(k => DimArray(fill(NaN, N, n_out), (z_dim, ti_dim))
                   for k in GROUND_PROFILE_LAYERS))
    output = DimStack(layers; layermetadata=cf_layermetadata(layers),
        metadata=merge(GEMB_CF_GLOBAL_ATTRIBUTES,
            Dict{String,Any}("dataset" => climate_forcing.dataset,
                "latitude" => climate_forcing.latitude, "longitude" => climate_forcing.longitude,
                "elevation" => climate_forcing.elevation,
                "lower_boundary" => "ground", "density_ice" => mp.density_ice),
            _profile_provenance(profile)))

    cf = climate_forcing
    f = (cf.temperature_air, cf.pressure_air, cf.precipitation, cf.wind_speed,
        cf.shortwave_downward, cf.longwave_downward, cf.vapor_pressure, cf.black_carbon_snow,
        cf.black_carbon_ice, cf.cloud_optical_thickness, cf.solar_zenith_angle,
        cf.shortwave_downward_diffuse, cf.cloud_fraction, cf.snow_drift)
    sums = Dict(k => 0.0 for k in (:melt, :runoff, :refreeze, :evaporation_condensation,
        :shortwave_net, :longwave_net, :heat_flux_sensible, :heat_flux_latent, :heat_flux_basal,
        :albedo_broadband, :blowing_snow, :temperature_air, :precipitation, :rain,
        :heat_flux_ground))
    mean_layers = (:shortwave_net, :longwave_net, :heat_flux_sensible, :heat_flux_latent,
        :heat_flux_basal, :albedo_broadband, :temperature_air, :heat_flux_ground)
    max_percolation = 0.0
    count = 0
    oi = 0
    run_mass_initial = column_mass(snow)
    run_in = 0.0

    for i in eachindex(times)
        cfs = ClimateForcingStep(dt_f, f[1][i], f[2][i], f[3][i], f[4][i], f[5][i], f[6][i],
            f[7][i], cf.temperature_air_mean, cf.wind_speed_mean, cf.precipitation_mean,
            cf.temperature_observation_height, cf.wind_observation_height, f[8][i], f[9][i],
            f[10][i], f[11][i], f[12][i], f[13][i], cf.accumulation_mean,
            cf.temperature_air_effective, f[14][i])
        snow, temperature_ground, water_ground, flux = gemb_core_ground(snow, temperature_ground,
            water_ground, g, cfs, model_parameters, verbose; thermal_workspace, column_workspace)

        for k in (:melt, :runoff, :refreeze, :evaporation_condensation, :shortwave_net,
                  :heat_flux_sensible, :heat_flux_latent, :heat_flux_basal, :albedo_broadband,
                  :rain, :heat_flux_ground)
            sums[k] += getfield(flux, k)
        end
        sums[:longwave_net] += cfs.longwave_downward - flux.longwave_upward
        sums[:blowing_snow] += flux.mass_blowing_snow
        sums[:temperature_air] += cfs.temperature_air
        sums[:precipitation] += cfs.precipitation
        max_percolation = max(max_percolation, flux.percolation_depth)
        count += 1
        run_in += cfs.precipitation + (flux.evaporation_condensation -
                  flux.evaporation_condensation_ground) + flux.mass_blowing_snow - flux.runoff

        if times[i] in output_times
            oi += 1
            for k in keys(sums)
                parent(output[k])[oi] = k in mean_layers ? sums[k] / count : sums[k]
                sums[k] = 0.0
            end
            parent(output[:percolation_depth])[oi] = max_percolation
            parent(output[:ice_slab_thickness])[oi] = flux.ice_slab_thickness
            parent(output[:ice_slab_depth])[oi] = flux.ice_slab_depth
            parent(output[:snow_depth])[oi] = sum(snow.dz)
            parent(output[:snow_water_equivalent])[oi] = column_mass(snow)
            for (k, d) in zip((:soil_temperature_005, :soil_temperature_010,
                    :soil_temperature_020, :soil_temperature_050, :soil_temperature_100),
                    GROUND_OUTPUT_DEPTHS)
                parent(output[k])[oi] = ground_soil_temperature(temperature_ground, g, d)
            end
            col = _ground_state_vectors(mp, g, snow.temperature, snow.dz, snow.density,
                snow.water, snow.grain_radius, snow.grain_dendricity, snow.grain_sphericity,
                snow.age, temperature_ground, water_ground)
            for l in GROUND_PROFILE_LAYERS
                parent(output[l])[:, oi] .= getfield(col, l)
            end
            max_percolation = 0.0
            count = 0
        end
    end

    if verbose
        residual = column_mass(snow) - run_mass_initial - run_in
        tol = max(1e-6, 1e-9 * (abs(run_in) + run_mass_initial))
        abs(residual) > tol && error("gemb: whole-run mass budget of the snow over the ground " *
                                     "does not close: residual = $(residual) kg m-2")
    end
    return output
end
