# The `gemb` time loop of a `lower_boundary = :prescribed` column. `gemb` hands over to
# `_gemb_prescribed` before it builds anything of the glacier path.
#
# The profile holds `ground_snow_cells_max` rows, surface first: the snow cells
# (`material = MATERIAL_SNOW`, possibly none), then empty slots (`MATERIAL_EMPTY`). It is the
# layout of `ground_column.jl` with no ground cells, and is built and split by the same helpers
# with a ground column of zero cells.

# Monolevel layers of a prescribed-base run in output order.
const PRESCRIBED_MONOLEVEL_LAYERS = (:melt, :runoff, :refreeze, :evaporation_condensation,
    :shortwave_net, :longwave_net, :heat_flux_sensible, :heat_flux_latent, :heat_flux_basal,
    :albedo_broadband, :blowing_snow, :percolation_depth, :ice_slab_thickness, :ice_slab_depth,
    :temperature_air, :precipitation, :rain, :temperature_base, :snow_depth,
    :snow_water_equivalent)

# The ground column of zero cells that gives the prescribed-base profile its layout.
_no_ground(mp::ModelParameters) = ground_column(mp, Float64[])

"""
    initialize_profile_prescribed(mp) -> profile::DimStack

The initial column of `lower_boundary = :prescribed`: no snow, `mp.ground_snow_cells_max` empty
slots. The profile carries the `RESTART_LAYERS` and `material`.
"""
function initialize_profile_prescribed(mp::ModelParameters)
    validate_parameters(mp)
    e = Float64[]
    state = _ground_state_vectors(mp, _no_ground(mp), e, e, e, e, e, e, e, e, e, e)
    zdim = Z(1:length(state.dz); metadata=cf_layer_index_attributes())
    layers = NamedTuple(l => DimArray(getfield(state, l), (zdim,))
                        for l in (RESTART_LAYERS..., :material))
    return DimStack(layers; layermetadata=cf_layermetadata(layers; time_axis=false))
end

"""
    _gemb_prescribed(profile, climate_forcing, mp, model_parameters; verbose, thermal_workspace,
                     column_workspace) -> DimStack

[`gemb`](@ref) for `lower_boundary = :prescribed`, stepping [`gemb_core_prescribed`](@ref) with
the base temperature of each step from the `temperature_ground` forcing layer, which must be
finite at every step. Output layers:

- the glacier run's surface fluxes and mass terms, with the same meaning and interval reductions;
  `heat_flux_basal` is the energy the prescribed base supplies (see `gemb_core_prescribed`);
- `temperature_base` [K, interval mean], the prescribed base temperature;
- `snow_depth` [m] and `snow_water_equivalent` [kg m-2], instantaneous;
- the profile layers, `ground_snow_cells_max` rows.

Under `verbose` the whole-run snow mass budget is checked as under `:ground`.
"""
function _gemb_prescribed(profile::DimStack, climate_forcing::ClimateForcing,
    mp::ModelParameters, model_parameters::ModelParameters; verbose::Bool,
    thermal_workspace::ThermalWorkspace, column_workspace::ColumnWorkspace)

    haskey(profile, :material) || error("gemb: lower_boundary = :prescribed needs a profile " *
        "with a `material` layer; build it with `initialize_profile(mp, cf)` under that setting")
    g = _no_ground(mp)
    full = NamedTuple(l => Vector{Float64}(parent(profile[l])) for l in GROUND_PROFILE_LAYERS)
    length(full.dz) == ground_snow_slots(mp) || error("gemb: the profile has " *
        "$(length(full.dz)) rows; a prescribed-base column of these parameters has " *
        "$(ground_snow_slots(mp))")
    snow, _, _ = _ground_state_split(full, g)

    T_base = parent(climate_forcing.temperature_ground)
    i_bad = findfirst(!isfinite, T_base)
    i_bad === nothing || error("gemb: lower_boundary = :prescribed needs the forcing layer " *
        "temperature_ground at every step; it is $(T_base[i_bad]) at step $(i_bad)")

    times = Vector{DateTime}(dims(climate_forcing.temperature_air, Ti).val)
    dt_f = Float64(climate_forcing.time_step)
    out_time = _compute_output_times(times, mp.output_frequency)
    output_times = Set(out_time)
    n_out = length(out_time)
    N = length(full.dz)
    ti_dim = Ti(out_time; metadata=cf_time_attributes())
    z_dim = Z(1:N; metadata=cf_layer_index_attributes())
    layers = merge(
        NamedTuple(k => DimArray(fill(NaN, n_out), (ti_dim,)) for k in PRESCRIBED_MONOLEVEL_LAYERS),
        NamedTuple(k => DimArray(fill(NaN, N, n_out), (z_dim, ti_dim))
                   for k in GROUND_PROFILE_LAYERS))
    output = DimStack(layers; layermetadata=cf_layermetadata(layers),
        metadata=merge(GEMB_CF_GLOBAL_ATTRIBUTES,
            Dict{String,Any}("dataset" => climate_forcing.dataset,
                "latitude" => climate_forcing.latitude, "longitude" => climate_forcing.longitude,
                "elevation" => climate_forcing.elevation,
                "lower_boundary" => "prescribed", "density_ice" => mp.density_ice),
            _profile_provenance(profile)))

    cf = climate_forcing
    f = (cf.temperature_air, cf.pressure_air, cf.precipitation, cf.wind_speed,
        cf.shortwave_downward, cf.longwave_downward, cf.vapor_pressure, cf.black_carbon_snow,
        cf.black_carbon_ice, cf.cloud_optical_thickness, cf.solar_zenith_angle,
        cf.shortwave_downward_diffuse, cf.cloud_fraction, cf.snow_drift)
    sum_keys = (:melt, :runoff, :refreeze, :evaporation_condensation, :shortwave_net,
        :longwave_net, :heat_flux_sensible, :heat_flux_latent, :heat_flux_basal,
        :albedo_broadband, :blowing_snow, :temperature_air, :precipitation, :rain,
        :temperature_base)
    sums = Dict(k => 0.0 for k in sum_keys)
    mean_layers = (:shortwave_net, :longwave_net, :heat_flux_sensible, :heat_flux_latent,
        :heat_flux_basal, :albedo_broadband, :temperature_air, :temperature_base)
    max_percolation = 0.0
    count = 0
    oi = 0
    run_mass_initial = column_mass(snow)
    run_in = 0.0
    e = Float64[]

    for i in eachindex(times)
        cfs = ClimateForcingStep(dt_f, f[1][i], f[2][i], f[3][i], f[4][i], f[5][i], f[6][i],
            f[7][i], cf.temperature_air_mean, cf.wind_speed_mean, cf.precipitation_mean,
            cf.temperature_observation_height, cf.wind_observation_height, f[8][i], f[9][i],
            f[10][i], f[11][i], f[12][i], f[13][i], cf.accumulation_mean,
            cf.temperature_air_effective, f[14][i])
        snow, flux = gemb_core_prescribed(snow, T_base[i], cfs, model_parameters, verbose;
            thermal_workspace, column_workspace)
        # The base is an argument of the step, not a cell: within a step it is this step's
        # prescribed value, and it changes only between steps.
        flux.temperature_base === T_base[i] || error("gemb: step $(i) ran on a base of " *
            "$(flux.temperature_base) K, not its prescribed $(T_base[i]) K")

        for k in (:melt, :runoff, :refreeze, :evaporation_condensation, :shortwave_net,
                  :heat_flux_sensible, :heat_flux_latent, :heat_flux_basal, :albedo_broadband,
                  :rain)
            sums[k] += getfield(flux, k)
        end
        sums[:longwave_net] += cfs.longwave_downward - flux.longwave_upward
        sums[:blowing_snow] += flux.mass_blowing_snow
        sums[:temperature_air] += cfs.temperature_air
        sums[:precipitation] += cfs.precipitation
        sums[:temperature_base] += T_base[i]
        max_percolation = max(max_percolation, flux.percolation_depth)
        count += 1
        run_in += cfs.precipitation + (flux.evaporation_condensation -
                  flux.evaporation_condensation_ground) + flux.mass_blowing_snow - flux.runoff

        if times[i] in output_times
            oi += 1
            for k in sum_keys
                parent(output[k])[oi] = k in mean_layers ? sums[k] / count : sums[k]
                sums[k] = 0.0
            end
            parent(output[:percolation_depth])[oi] = max_percolation
            parent(output[:ice_slab_thickness])[oi] = flux.ice_slab_thickness
            parent(output[:ice_slab_depth])[oi] = flux.ice_slab_depth
            parent(output[:snow_depth])[oi] = sum(snow.dz)
            parent(output[:snow_water_equivalent])[oi] = column_mass(snow)
            col = _ground_state_vectors(mp, g, snow.temperature, snow.dz, snow.density,
                snow.water, snow.grain_radius, snow.grain_dendricity, snow.grain_sphericity,
                snow.age, e, e)
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
        abs(residual) > tol && error("gemb: whole-run mass budget of the snow over the " *
                                     "prescribed base does not close: residual = $(residual) kg m-2")
    end
    return output
end
