"""
    gemb_spinup(profile, cf, mp; max_iterations=100, convergence_depth=nothing,
                convergence_delta_density=nothing, verbose=false)

Run GEMB for multiple spinup cycles to reach quasi-steady state.

Forces `output_frequency=:last` internally to minimize memory usage during spinup.
Returns the spun-up profile DimStack.

# Provenance
The returned profile carries a metadata `NamedTuple` (accessible via
`DimensionalData.metadata(profile)`) recording how the spinup ran and which
climatology it used: `spinup_cycles`, `spinup_converged`,
`spinup_final_delta_density`, the convergence parameters (`spinup_max_iterations`,
`spinup_convergence_delta_density`, `spinup_convergence_depth`), and the
climatology fields copied forward from `cf` (`climatology_window_start`,
`climatology_window_stop`, `climatology_n_years`, `climatology_steps_per_year`).
This provenance is propagated onto the [`gemb`](@ref) output when the spun-up
profile is used to start a transient run.

# Keyword arguments
- `max_iterations`: maximum number of spinup cycles (default 100). The spinup always
  exits after this many cycles even if convergence has not been reached.
- `convergence_depth`: depth [m] over which to compute depth-averaged density for
  convergence testing. Defaults to `mp.column_zmin`. If the column depth ever drops
  below this value an error is thrown.
- `convergence_delta_density`: if provided, the spinup exits early when the absolute
  change in depth-averaged density between consecutive cycles is less than this
  value [kg/m³]. When `nothing` (default) no convergence check is performed and the
  spinup always runs `max_iterations` cycles.
- `verbose`: print a convergence message when early exit occurs.

Matches MATLAB's `gemb_spinup.m`.
"""
function gemb_spinup(profile::DimStack, cf::ClimateForcing, mp::ModelParameters;
                     max_iterations::Int=100,
                     convergence_depth=nothing,
                     convergence_delta_density=nothing,
                     verbose::Bool=false)

    # Force output_frequency to :last for spinup efficiency
    mp_spinup = ModelParameters(;
        (field => getfield(mp, field) for field in fieldnames(ModelParameters) if field != :output_frequency)...,
        output_frequency=:last
    )

    # Resolve convergence_depth (default = mp.column_zmin)
    initial_depth = sum(parent(profile[:dz]))
    cdepth = convergence_depth === nothing ? Float64(mp.column_zmin) : Float64(convergence_depth)

    if cdepth > initial_depth
        error("convergence_depth ($cdepth m) exceeds initial column depth ($initial_depth m)")
    end

    prev_avg_density = nothing
    current_profile  = profile

    # Convergence provenance, captured for the returned profile.
    cycles_run          = 0
    converged           = false
    final_delta_density = NaN

    for cycle in 1:max_iterations
        cycles_run = cycle
        out = gemb(current_profile, cf, mp_spinup)
        current_profile = gemb_profile(out)

        if convergence_delta_density !== nothing
            col_depth = sum(parent(current_profile[:dz]))
            if col_depth < cdepth
                error("Column depth ($col_depth m) dropped below convergence_depth ($cdepth m) at cycle $cycle")
            end

            avg_rho = _depth_avg_density(current_profile, cdepth)
            if prev_avg_density !== nothing
                final_delta_density = abs(avg_rho - prev_avg_density)
                if final_delta_density < Float64(convergence_delta_density)
                    verbose && @info "gemb_spinup converged at cycle $cycle " *
                                     "(Δρ = $(round(final_delta_density, digits=4)) kg/m³)"
                    converged = true
                    break
                end
            end
            prev_avg_density = avg_rho
        end
    end

    @info "GEMB Spinup" climatology_window=(cf.climatology_window_start, cf.climatology_window_stop) cycles=cycles_run converged=converged final_delta_density=final_delta_density

    return _attach_spinup_provenance(current_profile, cf;
        cycles=cycles_run, converged=converged,
        final_delta_density=final_delta_density,
        max_iterations=max_iterations,
        convergence_delta_density=convergence_delta_density,
        convergence_depth=cdepth)
end

# Attach spinup + climatology provenance to the spun-up profile so it is a
# self-describing artifact. Uses plain `DimStack` metadata (free-form NamedTuple);
# climatology fields are copied forward from the forcing `cf` when present.
function _attach_spinup_provenance(profile::DimStack, cf::ClimateForcing;
        cycles, converged, final_delta_density,
        max_iterations, convergence_delta_density, convergence_depth)
    cf_meta = DD.metadata(cf)
    _cf(key) = haskey(cf_meta, key) ? cf_meta[key] : nothing
    prov = (
        spinup_cycles = cycles,
        spinup_converged = converged,
        spinup_final_delta_density = final_delta_density,
        spinup_max_iterations = max_iterations,
        spinup_convergence_delta_density = convergence_delta_density,
        spinup_convergence_depth = convergence_depth,
        climatology_window_start = _cf(:climatology_window_start),
        climatology_window_stop = _cf(:climatology_window_stop),
        climatology_n_years = _cf(:climatology_n_years),
        climatology_steps_per_year = _cf(:climatology_steps_per_year),
    )
    return rebuild(profile; metadata=prov)
end

# Compute depth-averaged density over [0, depth] using a cubic spline so the result
# is independent of the Lagrangian grid spacing at each cycle.
function _depth_avg_density(profile::DimStack, depth::Float64)
    z_pos = -dz2z(profile[:dz])  # negative depths → positive, shallowest first
    rho   = profile[:density]
    # Constant left-extrapolation extends the surface layer density to z=0
    # (z_pos[1] = dz[1]/2, not 0.0, because z_center holds cell-center depths)
    itp = DataInterpolations.CubicSpline(rho, z_pos;
              extrapolation_left=DataInterpolations.ExtrapolationType.Constant)
    return DataInterpolations.integral(itp, 0.0, depth) / depth
end
