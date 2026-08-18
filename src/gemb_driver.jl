"""
    gemb(profile::DimStack, climate_forcing::ClimateForcing, mp::ModelParameters;
         verbose=false, thermal_workspace=ThermalWorkspace())

Run the Glacier Energy and Mass Balance (GEMB) model.
Matches MATLAB's `gemb.m`.

Returns a DimStack containing time series of surface flux (monolevel)
and vertical profiles at the specified output frequency.

# Threading
`profile`, `mp`, and `thermal_workspace` divide cleanly into what a run reads and what it
writes. `mp` is an immutable description and holds no scratch, so a single `mp` can be read
by any number of threads stepping columns concurrently. `thermal_workspace` is the mutable
half: it is reused across timesteps within a run, so each concurrent run needs its own (the
default gives every call a fresh one). See [`ThermalWorkspace`](@ref).

# Provenance
The output metadata records where the forcing came from (`dataset`, `latitude`,
`longitude`, `elevation_offset`) and how the initial column was prepared. If
`profile` was produced by [`gemb_spinup`](@ref), its `spinup_*` / `climatology_*`
provenance is copied onto the output and `spinup_performed => true` is set. If
`profile` came straight from [`initialize_profile`](@ref) (no spinup), the output
records `spinup_performed => false`, so the run is unambiguous about having
started from an un-spun-up column.

# Metadata
Every output layer carries CF-style attributes — `units`, `long_name`, the interval
reduction as `cell_methods` (`time: sum` / `time: mean` / `time: point`), and a
`standard_name` where a CF standard name matches the quantity exactly — and the stack
carries `Conventions = "CF-1.11"` alongside the provenance above. See
[`GEMB_CF_ATTRIBUTES`](@ref) for the table and [`cf_attributes`](@ref) to read one
layer's attributes. Conformance is at the attribute level: GEMB.jl ships no NetCDF
writer, so the attributes are there for whichever writer or plotting code consumes them.
"""
function gemb(profile::DimStack, climate_forcing::ClimateForcing, mp::ModelParameters;
    verbose::Bool=false, thermal_workspace::ThermalWorkspace=ThermalWorkspace())
    _assert_restart_layers(profile)

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
        climate_forcing.wind_observation_height;
        accumulation_mean=climate_forcing.accumulation_mean,
        temperature_air_effective=climate_forcing.temperature_air_effective,
        dataset=climate_forcing.dataset,
        latitude=climate_forcing.latitude,
        longitude=climate_forcing.longitude,
        elevation=climate_forcing.elevation,
        elevation_native=climate_forcing.elevation_native,
        elevation_offset=climate_forcing.elevation_offset,
    )

    # Pre-compute dt_divisors for thermal sub-stepping. Excluded by way of `DERIVED_PARAMETERS`,
    # the exported statement of which fields this overwrites, so a caller comparing two runs'
    # settings reads the exclusion from here rather than restating it.
    model_parameters = ModelParameters(;
        (field => getfield(mp, field) for field in fieldnames(ModelParameters)
         if !(field in DERIVED_PARAMETERS))...,
        dt_divisors=fast_divisors(dt_int * 10000) ./ 10000
    )

    # Initialize column state from profile. Every layer read here is in `RESTART_LAYERS`, which
    # `_assert_restart_layers` has already checked for above — so a profile that reaches this point
    # cannot fail on a missing name.
    state = (
        temperature = Vector{Float64}(profile[:temperature]),
        dz = Vector{Float64}(profile[:dz]),
        density = Vector{Float64}(profile[:density]),
        water = Vector{Float64}(profile[:water]),
        grain_radius = Vector{Float64}(profile[:grain_radius]),
        grain_dendricity = Vector{Float64}(profile[:grain_dendricity]),
        grain_sphericity = Vector{Float64}(profile[:grain_sphericity]),
        age = Vector{Float64}(profile[:age]),
        evaporation_condensation = 0.0,
        melt_surface = 0.0,
    )

    # Compute output times. `_compute_output_times` already returns a sorted, unique
    # subset of `times`, so it becomes the output time axis directly; the `Set` is kept
    # only for O(1) membership testing inside the time loop below.
    out_time = _compute_output_times(times, mp.output_frequency)
    output_times = Set(out_time)
    n_outputs = length(out_time)

    # Period-based frequencies emit only *complete* days/weeks/months. A forcing record
    # shorter than one period therefore yields no output at all; warn rather than
    # silently returning empty arrays.
    if n_outputs == 0
        @warn "gemb: output_frequency=:$(mp.output_frequency) produced no output times — " *
              "the forcing record ($(length(times)) steps, $(times[1]) to $(times[end])) " *
              "does not span a complete $(mp.output_frequency == :monthly ? "month" :
                                          mp.output_frequency == :weekly ? "week" : "day"). " *
              "Returning empty output."
    end
    # The column now holds a fixed cell count and a fixed total depth for the whole run
    # (enforced each timestep by the two controllers in `grid_ops.jl`), so the profile
    # output is sized exactly to the column — no padding, no NaN rows, and no possibility
    # of overflow.
    profile_size = length(state.dz)
    z_target = sum(state.dz)
    _assert_grid_feasible(state.dz, z_target, mp)

    # Create output time coordinate and dimensions. `Z` is a bare 1-based cell index,
    # not a height, so its attributes say so rather than claiming a vertical coordinate.
    ti_dim = Ti(out_time; metadata=cf_time_attributes())
    z_dim = Z(1:profile_size; metadata=cf_layer_index_attributes())

    # Initialize output as NaN-filled arrays. Held as a NamedTuple first so the CF
    # attribute table can be indexed by the same keys in the same order, which is what
    # `cf_layermetadata` relies on.
    layers = (
        # Monolevel outputs
        melt=DimArray(fill(NaN, n_outputs), (ti_dim,)),
        runoff=DimArray(fill(NaN, n_outputs), (ti_dim,)),
        refreeze=DimArray(fill(NaN, n_outputs), (ti_dim,)),
        evaporation_condensation=DimArray(fill(NaN, n_outputs), (ti_dim,)),
        shortwave_net=DimArray(fill(NaN, n_outputs), (ti_dim,)),
        longwave_net=DimArray(fill(NaN, n_outputs), (ti_dim,)),
        heat_flux_sensible=DimArray(fill(NaN, n_outputs), (ti_dim,)),
        heat_flux_latent=DimArray(fill(NaN, n_outputs), (ti_dim,)),
        heat_flux_basal=DimArray(fill(NaN, n_outputs), (ti_dim,)),
        albedo_broadband=DimArray(fill(NaN, n_outputs), (ti_dim,)),
        densification_from_compaction=DimArray(fill(NaN, n_outputs), (ti_dim,)),
        densification_from_melt=DimArray(fill(NaN, n_outputs), (ti_dim,)),
        strain_thinning=DimArray(fill(NaN, n_outputs), (ti_dim,)),
        ice_flux=DimArray(fill(NaN, n_outputs), (ti_dim,)),
        firn_air_content=DimArray(fill(NaN, n_outputs), (ti_dim,)),
        firn_air_content_10m=DimArray(fill(NaN, n_outputs), (ti_dim,)),
        firn_air_content_20m=DimArray(fill(NaN, n_outputs), (ti_dim,)),
        close_off_age=DimArray(fill(NaN, n_outputs), (ti_dim,)),
        percolation_depth=DimArray(fill(NaN, n_outputs), (ti_dim,)),
        ice_slab_thickness=DimArray(fill(NaN, n_outputs), (ti_dim,)),
        ice_slab_depth=DimArray(fill(NaN, n_outputs), (ti_dim,)),
        aquifer_thickness=DimArray(fill(NaN, n_outputs), (ti_dim,)),
        aquifer_depth=DimArray(fill(NaN, n_outputs), (ti_dim,)),

        # Forcing summary outputs
        temperature_air=DimArray(fill(NaN, n_outputs), (ti_dim,)),
        precipitation=DimArray(fill(NaN, n_outputs), (ti_dim,)),
        rain=DimArray(fill(NaN, n_outputs), (ti_dim,)),

        # Profile outputs (2D: vertical × time)
        temperature=DimArray(fill(NaN, profile_size, n_outputs), (z_dim, ti_dim)),
        dz=DimArray(fill(NaN, profile_size, n_outputs), (z_dim, ti_dim)),
        density=DimArray(fill(NaN, profile_size, n_outputs), (z_dim, ti_dim)),
        water=DimArray(fill(NaN, profile_size, n_outputs), (z_dim, ti_dim)),
        grain_radius=DimArray(fill(NaN, profile_size, n_outputs), (z_dim, ti_dim)),
        grain_dendricity=DimArray(fill(NaN, profile_size, n_outputs), (z_dim, ti_dim)),
        grain_sphericity=DimArray(fill(NaN, profile_size, n_outputs), (z_dim, ti_dim)),
        age=DimArray(fill(NaN, profile_size, n_outputs), (z_dim, ti_dim)),
    )

    # Optional layers, appended last so the base stack's layer order — and therefore the
    # order `cf_layermetadata` walks — is unchanged whether or not they are present.
    if mp.output_viscosity
        layers = merge(layers, (;
            viscosity=DimArray(fill(NaN, profile_size, n_outputs), (z_dim, ti_dim)),))
    end

    output = DimStack(layers;
        # Per-layer CF attributes: units, long_name, the interval reduction
        # (`cell_methods`), and a standard_name only where a CF standard name matches
        # GEMB's quantity exactly. See `cf_metadata.jl`.
        layermetadata=cf_layermetadata(layers),
        # CF global attributes, then forcing provenance so downstream consumers
        # (e.g. `plot_output`) can report where the forcing came from, plus
        # spinup/climatology provenance from the initial profile (or a flag noting
        # the run started from an un-spun-up profile). Globals go first so a
        # run-specific key would win on any collision.
        metadata=merge(
            GEMB_CF_GLOBAL_ATTRIBUTES,
            Dict{String,Any}(
                "dataset" => climate_forcing.dataset,
                "latitude" => climate_forcing.latitude,
                "longitude" => climate_forcing.longitude,
                "elevation" => climate_forcing.elevation,
                "elevation_native" => climate_forcing.elevation_native,
                "elevation_offset" => climate_forcing.elevation_offset,
                # Ice density this run used [kg m-3]. Recorded so a consumer can convert the
                # mass-flux outputs to metres of ice with the same value the run did, rather
                # than assuming a default that may not match `mp` (`plot_output` does this for
                # its detrended height axis).
                "density_ice" => mp.density_ice,
            ),
            _profile_provenance(profile),
        ),
    )

    # Run the time loop through a function barrier so it specializes on the concrete
    # forcing types rather than re-dispatching each timestep. The layers are passed
    # as-is: `ClimateForcing` is a `DimStack`, so `climate_forcing.temperature_air`
    # already infers to a concrete `DimVector`, and integer-indexing one is free
    # (measured indistinguishable from the bare `Vector`).
    _gemb_time_loop!(output, state, model_parameters, mp, verbose, thermal_workspace,
        times, output_times, profile_size, z_target, Float64(dt_int),
        climate_forcing.temperature_air,
        climate_forcing.pressure_air,
        climate_forcing.precipitation,
        climate_forcing.wind_speed,
        climate_forcing.shortwave_downward,
        climate_forcing.longwave_downward,
        climate_forcing.vapor_pressure,
        climate_forcing.black_carbon_snow,
        climate_forcing.black_carbon_ice,
        climate_forcing.cloud_optical_thickness,
        climate_forcing.solar_zenith_angle,
        climate_forcing.shortwave_downward_diffuse,
        climate_forcing.cloud_fraction,
        climate_forcing.temperature_air_mean,
        climate_forcing.wind_speed_mean,
        climate_forcing.precipitation_mean,
        climate_forcing.temperature_observation_height,
        climate_forcing.wind_observation_height,
        climate_forcing.accumulation_mean,
        climate_forcing.temperature_air_effective)

    # Return the DimStack (already populated during time loop)
    return output
end

"""
    _assert_restart_layers(profile)

Check that `profile` carries every layer in [`RESTART_LAYERS`](@ref) before the run starts.

Without this the first missing layer surfaces as a bare `FieldError` from deep inside the time
loop, naming a field of an internal `state` NamedTuple rather than the profile that is actually
incomplete — and it surfaces only once the run is already underway, which for a resumed
long-record run can be well after the caller has committed real work. The common cause is a
column persisted to disk before a layer joined the state, so the message names the gap and the
remedy instead.
"""
function _assert_restart_layers(profile::DimStack)
    missing_layers = filter(l -> !haskey(profile, l), RESTART_LAYERS)
    isempty(missing_layers) && return nothing
    error("gemb: the profile is missing the state layer(s) " *
          join(missing_layers, ", ") * ". A complete column state is " *
          join(RESTART_LAYERS, ", ") * " (see GEMB.RESTART_LAYERS). A profile from " *
          "`gemb_profile` or `initialize_profile` always carries all of them; one short of a " *
          "layer has most likely been round-tripped through storage that predates it, and the " *
          "column must be rebuilt rather than continued.")
end

"""
    _assert_grid_feasible(dz, z_target, mp)

Check at startup that the column the two grid controllers are being handed is one they can
actually hold: at least two cells (the thermal solve indexes `temperature[2]`), all
thicknesses positive and finite, a positive target depth, and — the substantive check —
that `z_target` lies inside `[Σdzmin_i, Σdzmax_i]`, the depth range the `N` cells can span
at their band limits.

That last condition is what makes the two controllers compatible. The count controller
holds the cell count at `N` while keeping cells inside their bands; the mass controller
holds total depth at `z_target`. If `z_target` fell outside the band-limited range, no
admissible `N`-cell grid would have the required depth and the two controllers would fight
indefinitely — the count controller pushing cells toward their bands, the mass controller
pulling total depth away from what those bands permit.

Stable only because [`column_bands!`](@ref) references each cell's band to its **depth** (see
there for why). The range then holds steady over a run — measured on the default grid, `Σdzmin`
118–127 m and `Σdzmax` 355–382 m against a 254.6 m target — so the check can be trusted rather
than rejecting legitimate configurations after the first cycle.
"""
function _assert_grid_feasible(dz::Vector{Float64}, z_target::Float64, mp::ModelParameters)
    n = length(dz)
    if n < 2
        error("gemb: the column has $n cell(s); at least 2 are required (the thermal " *
              "solve differences adjacent cells). Check column_depth_max / column_dztop.")
    end
    if !all(d -> isfinite(d) && d > 0.0, dz)
        error("gemb: the initial column contains a non-positive or non-finite cell " *
              "thickness. Every dz must be finite and > 0.")
    end
    if !(isfinite(z_target) && z_target > 0.0)
        error("gemb: the fixed column depth must be finite and positive, got $z_target m.")
    end

    dzmin = Vector{Float64}(undef, n)
    dzmax = Vector{Float64}(undef, n)
    column_bands!(dzmin, dzmax, dz, mp)
    z_lo = sum(dzmin)
    z_hi = sum(dzmax)
    if !(z_lo - D_TOLERANCE <= z_target <= z_hi + D_TOLERANCE)
        error("gemb: the fixed column depth ($(z_target) m) is outside the range $n cells " *
              "can span at their band limits ([$(z_lo), $(z_hi)] m), so no admissible " *
              "grid has the required depth and the count and depth controllers cannot " *
              "both be satisfied. Adjust column_depth_max / column_dzmin / column_dzmax / " *
              "column_ztop / column_zy.")
    end
    return nothing
end

# Extract spinup / climatology provenance from the initial `profile` for stamping
# onto the `gemb` output. `gemb_spinup` attaches a NamedTuple of `spinup_*` /
# `climatology_*` keys; a profile straight from `initialize_profile` (no spinup)
# carries `NoMetadata`. In that case record `spinup_performed => false` so the
# output is unambiguous about having started from an un-spun-up column.
function _profile_provenance(profile::DimStack)
    md = DD.metadata(profile)
    if md isa DD.Dimensions.Lookups.NoMetadata || isempty(keys(md))
        return Dict{String,Any}("spinup_performed" => false)
    end
    prov = Dict{String,Any}("spinup_performed" => true)
    for k in keys(md)
        prov[String(k)] = md[k]
    end
    return prov
end

"""
    _gemb_time_loop!(output, state, model_parameters, mp, verbose, thermal_workspace,
                     times, output_times, profile_size, z_target, dt_f,
                     <forcing arrays/scalars>)

Function-barrier inner loop for [`gemb`](@ref). Receives the forcing series as
arguments so the compiler specializes on their concrete types instead of
re-dispatching every timestep. Mutates `output` in place.

`profile_size` doubles as the fixed cell count and the profile output row count;
`z_target` is the fixed total column depth. Both are passed to every `gemb_core` call and
hold at every timestep boundary.
"""
function _gemb_time_loop!(output, state, model_parameters, mp, verbose::Bool,
    thermal_workspace::ThermalWorkspace,
    times::Vector{DateTime}, output_times, profile_size::Int, z_target::Float64,
    dt_f::Float64,
    f_temperature_air::AbstractVector, f_pressure_air::AbstractVector,
    f_precipitation::AbstractVector, f_wind_speed::AbstractVector,
    f_shortwave_downward::AbstractVector, f_longwave_downward::AbstractVector,
    f_vapor_pressure::AbstractVector, f_black_carbon_snow::AbstractVector,
    f_black_carbon_ice::AbstractVector, f_cloud_optical_thickness::AbstractVector,
    f_solar_zenith_angle::AbstractVector, f_shortwave_downward_diffuse::AbstractVector,
    f_cloud_fraction::AbstractVector,
    temperature_air_mean::Float64, wind_speed_mean::Float64,
    precipitation_mean::Float64, temperature_observation_height::Float64,
    wind_observation_height::Float64,
    accumulation_mean::Float64, temperature_air_effective::Float64)

    # Extract concrete output arrays once (avoids per-write DimStack/At dispatch).
    out_melt = parent(output[:melt])
    out_runoff = parent(output[:runoff])
    out_refreeze = parent(output[:refreeze])
    out_ec = parent(output[:evaporation_condensation])
    out_dcomp = parent(output[:densification_from_compaction])
    out_dmelt = parent(output[:densification_from_melt])
    out_strain = parent(output[:strain_thinning])
    out_swnet = parent(output[:shortwave_net])
    out_lwnet = parent(output[:longwave_net])
    out_shf = parent(output[:heat_flux_sensible])
    out_lhf = parent(output[:heat_flux_latent])
    out_bhf = parent(output[:heat_flux_basal])
    out_albedo_broadband = parent(output[:albedo_broadband])
    out_fac = parent(output[:firn_air_content])
    out_fac10 = parent(output[:firn_air_content_10m])
    out_fac20 = parent(output[:firn_air_content_20m])
    out_close_off_age = parent(output[:close_off_age])
    out_percolation_depth = parent(output[:percolation_depth])
    out_slab_thickness = parent(output[:ice_slab_thickness])
    out_slab_depth = parent(output[:ice_slab_depth])
    out_aquifer_thickness = parent(output[:aquifer_thickness])
    out_aquifer_depth = parent(output[:aquifer_depth])
    out_ice_flux = parent(output[:ice_flux])
    out_ta = parent(output[:temperature_air])
    out_precip = parent(output[:precipitation])
    out_rain = parent(output[:rain])
    out_temperature = parent(output[:temperature])
    out_dz = parent(output[:dz])
    out_density = parent(output[:density])
    out_water = parent(output[:water])
    out_grain_radius = parent(output[:grain_radius])
    out_grain_dendricity = parent(output[:grain_dendricity])
    out_grain_sphericity = parent(output[:grain_sphericity])
    out_age = parent(output[:age])
    # Empty when the layer was not allocated, which is what gates the write below. Empty
    # rather than `nothing` for the same reason as `NO_VISCOSITY`: this is a local in the
    # time loop, so a `Union{Nothing,Matrix{Float64}}` would leave its element type
    # unknowable at the write site.
    out_viscosity = mp.output_viscosity ? parent(output[:viscosity]) : Matrix{Float64}(undef, 0, 0)

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
    cum_bhf = 0.0
    cum_albedo_broadband = 0.0
    cum_densification_compaction = 0.0
    cum_densification_melt = 0.0
    cum_strain_thinning = 0.0
    cum_firn_air_content = 0.0
    cum_firn_air_content_10m = 0.0
    cum_firn_air_content_20m = 0.0
    # Interval *maximum*, not a sum: the question a wetting-front depth answers is how deep
    # water got during the interval, which is what upward-looking radar measures.
    max_percolation_depth = 0.0
    cum_ice_flux = 0.0
    cum_temperature_air = 0.0
    cum_precipitation = 0.0
    cum_count = 0

    # Whole-run mass budget (verbose only). `gemb_core` checks the budget every timestep,
    # but only against that timestep's own terms — a bias far below the per-step tolerance
    # would pass every step and still accumulate over a multi-decade run. These run-level
    # accumulators are never reset by the output interval, so they close the budget across
    # the entire run. Reported once at the end rather than per timestep.
    run_mass_initial = verbose ? column_mass(state) : 0.0
    run_precipitation = 0.0
    run_runoff = 0.0
    run_ec = 0.0
    run_mass_added = 0.0
    run_mass_lateral = 0.0

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
            f_cloud_fraction[i],
            accumulation_mean,
            temperature_air_effective)

        # Run physics for single timestep, holding the column at its fixed cell count
        # and fixed total depth.
        state, flux = gemb_core(state, forcing_step, model_parameters, verbose;
            n_target=profile_size, z_target=z_target,
            thermal_workspace=thermal_workspace)

        # Basal ice flux over this interval, in metres of ice. An interval *sum*, like every
        # other mass flux in the output, so `cumsum` recovers the cumulative flux exactly;
        # storing the cumulative directly (as this did before) meant reporting the interval
        # mean of a running total, which lagged the true cumulative by half an interval.
        cum_ice_flux += flux.mass_added / density_ice

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
        cum_bhf += flux.heat_flux_basal
        cum_albedo_broadband += flux.albedo_broadband
        cum_densification_compaction += flux.densification_from_compaction
        cum_densification_melt += flux.densification_from_melt
        # Metres of thinning, sign-flipped so a divergent column reports a positive number.
        # `flux.mass_lateral` is negative when mass leaves laterally.
        cum_strain_thinning -= flux.mass_lateral / density_ice
        # Firn air content [m of air]: Σ dz (1 - ρ/ρ_ice). Normalizing by `density_ice` (not
        # 1000) is what makes this metres of air rather than metres of water equivalent — see
        # upstream GEMB issue #198. The depth-limited variants are what the firn-core
        # literature reports and are the only form comparable to it.
        cum_firn_air_content += firn_air_content(state.dz, state.density, density_ice)
        cum_firn_air_content_10m += firn_air_content(state.dz, state.density, density_ice, 10.0)
        cum_firn_air_content_20m += firn_air_content(state.dz, state.density, density_ice, 20.0)
        max_percolation_depth = max(max_percolation_depth, flux.percolation_depth)
        cum_temperature_air += forcing_step.temperature_air
        cum_precipitation += forcing_step.precipitation
        cum_count += 1

        if verbose
            run_precipitation += forcing_step.precipitation
            run_runoff += flux.runoff
            run_ec += state.evaporation_condensation
            run_mass_added += flux.mass_added
            run_mass_lateral += flux.mass_lateral
        end

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
                out_strain[oi] = cum_strain_thinning
                out_ice_flux[oi] = cum_ice_flux

                # Averaged variables (division preserved for bit-identical results)
                out_swnet[oi] = cum_shortwave_net / cum_count
                out_lwnet[oi] = cum_longwave_net / cum_count
                out_shf[oi] = cum_shf / cum_count
                out_lhf[oi] = cum_lhf / cum_count
                out_bhf[oi] = cum_bhf / cum_count
                out_albedo_broadband[oi] = cum_albedo_broadband / cum_count
                out_fac[oi] = cum_firn_air_content / cum_count
                out_fac10[oi] = cum_firn_air_content_10m / cum_count
                out_fac20[oi] = cum_firn_air_content_20m / cum_count

                # Interval maximum
                out_percolation_depth[oi] = max_percolation_depth

                # Instantaneous: slab geometry is slowly-varying column state, and
                # `ice_slab_depth` is NaN when no slab qualifies, which an interval mean would
                # poison for the whole interval. `close_off_age` is NaN for an open column for
                # the same reason, and is likewise read off the current column rather than
                # `flux` — it is pure state, not a budget term.
                out_close_off_age[oi] = close_off_age(state.density, state.age)
                out_slab_thickness[oi] = flux.ice_slab_thickness
                out_slab_depth[oi] = flux.ice_slab_depth

                # Instantaneous for the same reasons: the water table is column state, and
                # `aquifer_depth` is NaN for a dry column.
                out_aquifer_thickness[oi] = flux.aquifer_thickness
                out_aquifer_depth[oi] = flux.aquifer_depth

                # Forcing summary
                out_ta[oi] = cum_temperature_air / cum_count
                out_precip[oi] = cum_precipitation
                out_rain[oi] = cum_rain
            end

            # Profile data. The column is a fixed `profile_size` cells (guaranteed by the
            # grid controllers), so rows map one-to-one onto cells with the surface at row
            # 1 and no padding. Consumers can index rows directly instead of scanning for
            # the first non-NaN row.
            m = length(state.dz)

            if m != profile_size
                error("Column length ($m) does not match the fixed profile size " *
                      "($profile_size). The grid controllers in manage_layer_thickness / " *
                      "trim_bottom! should make this impossible.")
            end

            @inbounds for k in 1:m
                out_temperature[k, oi] = state.temperature[k]
                out_dz[k, oi] = state.dz[k]
                out_density[k, oi] = state.density[k]
                out_water[k, oi] = state.water[k]
                out_grain_radius[k, oi] = state.grain_radius[k]
                out_grain_dendricity[k, oi] = state.grain_dendricity[k]
                out_grain_sphericity[k, oi] = state.grain_sphericity[k]
                out_age[k, oi] = state.age[k]
            end

            # From `flux`, not `state`: viscosity is a per-step diagnostic of the
            # densification call, not column state carried forward. Written at the output
            # step only, so it is the instantaneous value at this timestamp rather than an
            # interval mean — averaging it would mix Pa s across cells whose density changed,
            # and NaN under any non-Crocus scheme would poison the whole interval anyway.
            if !isempty(out_viscosity)
                @inbounds for k in 1:m
                    out_viscosity[k, oi] = flux.viscosity[k]
                end
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
            cum_bhf = 0.0
            cum_albedo_broadband = 0.0
            cum_densification_compaction = 0.0
            cum_densification_melt = 0.0
            cum_strain_thinning = 0.0
            cum_firn_air_content = 0.0
            cum_firn_air_content_10m = 0.0
            cum_firn_air_content_20m = 0.0
            max_percolation_depth = 0.0
            cum_ice_flux = 0.0
            cum_temperature_air = 0.0
            cum_precipitation = 0.0
            cum_count = 0
        end
    end

    # Whole-run mass budget. Everything that entered the column (precipitation, basal flux,
    # condensation) minus everything that left (runoff, evaporation) must equal the change
    # in column mass. Tolerance scales with the total mass turned over: the per-step 1e-3
    # absolute check in `gemb_core` is the strict one, and this looks for accumulated bias
    # across a run of arbitrary length, where floating-point summation error grows with the
    # number of steps.
    if verbose
        run_mass_final = column_mass(state)
        supplied = run_precipitation + run_ec + run_mass_added + run_mass_lateral - run_runoff
        residual = (run_mass_final - run_mass_initial) - supplied
        turnover = abs(run_precipitation) + abs(run_runoff) + abs(run_ec) +
                   abs(run_mass_added) + abs(run_mass_lateral) + abs(run_mass_initial)
        tol = max(1e-6, 1e-9 * turnover)
        if abs(residual) > tol
            error("gemb: whole-run mass budget does not close: residual = $(residual) " *
                  "kg m-2 (tolerance $(tol)). Column mass changed by " *
                  "$(run_mass_final - run_mass_initial); sources supplied $(supplied) " *
                  "(precipitation $(run_precipitation), basal $(run_mass_added), " *
                  "evap/cond $(run_ec), runoff $(run_runoff), " *
                  "horizontal strain $(run_mass_lateral)).")
        end
        @info "GEMB whole-run mass budget closed" residual basal_flux=run_mass_added
    end

    return output
end

"""
Compute output times based on output frequency.
Returns the last timestep of each day/month, all timesteps, or just the final one.

For `:daily`/`:weekly`/`:monthly`, a timestep is emitted when the next timestep
falls in a different day/week/month. The final timestep is compared against a
synthetic successor (`times[end] + dt`) rather than emitted unconditionally, so a
trailing *partial* period (e.g. a lone midnight boundary step belonging to a new
day) is not saved — only complete days/weeks/months are written. Weeks are grouped
by Monday-anchored calendar week (`floor(Date, Week)`). This matches MATLAB's `gemb.m` output
indexing (which appends `dates(end) + (dates(end)-dates(end-1))` before diffing).
"""
function _compute_output_times(times::Vector{DateTime}, frequency::Symbol)
    frequency == :all && return times
    frequency == :last && return [times[end]]
    groupfn = frequency == :daily ? Date :
              frequency == :weekly ? (t -> floor(Date(t), Week)) :
              frequency == :monthly ? (t -> (year(t), month(t))) :
              error("output_frequency must be one of: :all, :daily, :weekly, :monthly, :last")
    n = length(times)
    n == 1 && return copy(times)
    # Synthetic successor for the final step (uniform time axis assumed).
    ext_last = times[n] + (times[n] - times[n-1])
    return [times[i] for i in 1:n
            if groupfn(times[i]) != groupfn(i < n ? times[i+1] : ext_last)]
end

