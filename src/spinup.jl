"""
    gemb_spinup(profile, cf, mp; max_iterations=100, convergence_depth=nothing,
                convergence_delta_density=nothing, verbose=false)

Run GEMB for multiple spinup cycles to reach quasi-steady state.

Forces `output_frequency=:last` internally to minimize memory usage during spinup.
Returns the spun-up profile DimStack.

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

    for cycle in 1:max_iterations
        out = gemb(current_profile, cf, mp_spinup)
        current_profile = gemb_profile(out)

        if convergence_delta_density !== nothing
            col_depth = sum(parent(current_profile[:dz]))
            if col_depth < cdepth
                error("Column depth ($col_depth m) dropped below convergence_depth ($cdepth m) at cycle $cycle")
            end

            avg_rho = _depth_avg_density(current_profile, cdepth)
            if prev_avg_density !== nothing &&
               abs(avg_rho - prev_avg_density) < Float64(convergence_delta_density)
                verbose && @info "gemb_spinup converged at cycle $cycle " *
                                 "(Δρ = $(round(abs(avg_rho - prev_avg_density), digits=4)) kg/m³)"
                break
            end
            prev_avg_density = avg_rho
        end
    end

    return current_profile
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
