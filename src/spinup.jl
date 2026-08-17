# Number of trailing cycles the drift criterion regresses over. Ten is enough for the slope
# to see through the seasonal limit cycle the column settles into without being so long that
# early transient cycles keep polluting the window.
const SPINUP_DRIFT_WINDOW = 10

"""
    gemb_spinup(profile, cf, mp; max_iterations=100,
                convergence_delta_density=nothing, convergence_drift_density=nothing,
                drift_window=$SPINUP_DRIFT_WINDOW, verbose=false)

Run GEMB for multiple spinup cycles to reach quasi-steady state.

Forces `output_frequency=:last` internally to minimize memory usage during spinup.
Returns the spun-up profile DimStack.

Convergence is judged on the **column-mean density** — averaged over the whole column,
since the column depth is fixed for the run (chosen by [`initialize_profile`](@ref),
held by [`trim_bottom!`](@ref)) and so the same domain is compared at every cycle.

# Provenance
The returned profile carries a metadata `NamedTuple` (accessible via
`DimensionalData.metadata(profile)`) recording how the spinup ran and which
climatology it used: `spinup_cycles`, `spinup_converged`,
`spinup_final_delta_density`, `spinup_final_drift_density`, the convergence
parameters (`spinup_max_iterations`, `spinup_convergence_delta_density`,
`spinup_convergence_drift_density`, `spinup_drift_window`), and the climatology
fields copied forward from `cf` (`climatology_window_start`,
`climatology_window_stop`, `climatology_n_years`, `climatology_steps_per_year`).
This provenance is propagated onto the [`gemb`](@ref) output when the spun-up
profile is used to start a transient run.

# Convergence criteria

Two independent tests, either or both of which may be requested. When **both** are
given, **both** must hold — a small per-cycle step and no systematic trend.

- `convergence_delta_density`: exit when the absolute change in column-mean density
  between consecutive cycles is below this value [kg/m³]. Cheap, but a *step* test: a
  column creeping steadily at just under the tolerance passes it while still drifting.
- `convergence_drift_density`: exit when the magnitude of the least-squares slope of
  column-mean density against cycle number, over the trailing `drift_window` cycles,
  is below this value [kg/m³ per cycle]. This is the trend test the step test cannot
  make: it distinguishes a column oscillating about a settled mean (slope ≈ 0) from
  one still densifying monotonically, and by using every point in the window it is
  not fooled by a single anomalous cycle. Inactive until `drift_window` cycles have
  run.

When neither is given the spinup always runs `max_iterations` cycles.

# Keyword arguments
- `max_iterations`: maximum number of spinup cycles (default 100). The spinup always
  exits after this many cycles even if convergence has not been reached.
- `convergence_delta_density`, `convergence_drift_density`: see above.
- `drift_window`: trailing cycles the drift slope is fitted over (default
  $SPINUP_DRIFT_WINDOW). Must be ≥ 2.
- `verbose`: print a convergence message when early exit occurs.
- `thermal_workspace`: reusable thermal-solve scratch, shared across every cycle of this
  spinup. One belongs to one column; give each concurrent thread its own
  (see [`ThermalWorkspace`](@ref)).

Derived from MATLAB's `gemb_spinup.m`; the drift criterion and the whole-column mean
have no MATLAB counterpart.
"""
function gemb_spinup(profile::DimStack, cf::ClimateForcing, mp::ModelParameters;
                     max_iterations::Int=100,
                     convergence_delta_density=nothing,
                     convergence_drift_density=nothing,
                     drift_window::Int=SPINUP_DRIFT_WINDOW,
                     verbose::Bool=false,
                     thermal_workspace::ThermalWorkspace=ThermalWorkspace())

    if convergence_drift_density !== nothing && drift_window < 2
        error("drift_window must be at least 2 to fit a slope, got $drift_window")
    end

    # Force output_frequency to :last for spinup efficiency
    mp_spinup = ModelParameters(;
        (field => getfield(mp, field) for field in fieldnames(ModelParameters) if field != :output_frequency)...,
        output_frequency=:last
    )

    prev_avg_density = nothing
    current_profile  = profile
    # Column-mean density per cycle, for the drift regression. Only the trailing
    # `drift_window` entries are ever read, but a spinup is at most a few hundred
    # cycles of one Float64 each, so there is nothing to gain by ringbuffering it.
    history = Float64[]

    # Convergence provenance, captured for the returned profile.
    cycles_run          = 0
    converged           = false
    final_delta_density = NaN
    final_drift_density = NaN

    checking = convergence_delta_density !== nothing || convergence_drift_density !== nothing

    for cycle in 1:max_iterations
        cycles_run = cycle
        out = gemb(current_profile, cf, mp_spinup; thermal_workspace=thermal_workspace)
        current_profile = gemb_profile(out)

        checking || continue

        avg_rho = _column_mean_density(current_profile)
        push!(history, avg_rho)

        # Each criterion reports `nothing` while it lacks the cycles to judge, which is
        # distinct from reporting "not converged": an unmet criterion must not be able to
        # pass by abstention, so a `nothing` blocks convergence.
        if prev_avg_density !== nothing
            final_delta_density = abs(avg_rho - prev_avg_density)
        end
        prev_avg_density = avg_rho

        if convergence_drift_density !== nothing && length(history) >= drift_window
            final_drift_density = abs(_density_drift(history, drift_window))
        end

        delta_ok = convergence_delta_density === nothing ||
                   (!isnan(final_delta_density) &&
                    final_delta_density < Float64(convergence_delta_density))
        drift_ok = convergence_drift_density === nothing ||
                   (!isnan(final_drift_density) &&
                    final_drift_density < Float64(convergence_drift_density))

        if delta_ok && drift_ok
            verbose && @info "gemb_spinup converged at cycle $cycle " *
                             "(Δρ = $(round(final_delta_density, digits=4)) kg/m³, " *
                             "drift = $(round(final_drift_density, digits=6)) kg/m³/cycle)"
            converged = true
            break
        end
    end

    @info "GEMB Spinup" climatology_window=(cf.climatology_window_start, cf.climatology_window_stop) cycles=cycles_run converged=converged final_delta_density=final_delta_density final_drift_density=final_drift_density

    return _attach_spinup_provenance(current_profile, cf;
        cycles=cycles_run, converged=converged,
        final_delta_density=final_delta_density,
        final_drift_density=final_drift_density,
        max_iterations=max_iterations,
        convergence_delta_density=convergence_delta_density,
        convergence_drift_density=convergence_drift_density,
        drift_window=drift_window)
end

# Attach spinup + climatology provenance to the spun-up profile so it is a
# self-describing artifact. Uses plain `DimStack` metadata (free-form NamedTuple);
# climatology fields are copied forward from the forcing `cf` when present.
function _attach_spinup_provenance(profile::DimStack, cf::ClimateForcing;
        cycles, converged, final_delta_density, final_drift_density,
        max_iterations, convergence_delta_density, convergence_drift_density,
        drift_window)
    cf_meta = DD.metadata(cf)
    _cf(key) = haskey(cf_meta, key) ? cf_meta[key] : nothing
    prov = (
        spinup_cycles = cycles,
        spinup_converged = converged,
        spinup_final_delta_density = final_delta_density,
        spinup_final_drift_density = final_drift_density,
        spinup_max_iterations = max_iterations,
        spinup_convergence_delta_density = convergence_delta_density,
        spinup_convergence_drift_density = convergence_drift_density,
        spinup_drift_window = drift_window,
        climatology_window_start = _cf(:climatology_window_start),
        climatology_window_stop = _cf(:climatology_window_stop),
        climatology_n_years = _cf(:climatology_n_years),
        climatology_steps_per_year = _cf(:climatology_steps_per_year),
    )
    return rebuild(profile; metadata=prov)
end

# Mass-weighted mean density of the whole column [kg/m³] — total column mass over total
# column depth. Exact, and grid-independent without interpolation: the Lagrangian grid
# redistributes the same mass over the same (fixed) depth, so no spline is needed here.
# A partial-depth average would need one, because the averaging depth would then fall
# mid-cell; that is a large part of why the convergence depth is gone.
function _column_mean_density(profile::DimStack)
    dz  = parent(profile[:dz])
    rho = parent(profile[:density])
    return sum(rho[i] * dz[i] for i in eachindex(dz)) / sum(dz)
end

# Least-squares slope of the last `window` column-mean densities against cycle number
# [kg/m³ per cycle]. Cycles are equally spaced, so the design matrix is known a priori and
# the slope reduces to a closed form: with x centered on the window, x̄ = 0 and the slope is
# Σxᵢyᵢ / Σxᵢ². No allocation, no least-squares solve.
function _density_drift(history::Vector{Float64}, window::Int)
    n = min(window, length(history))
    n < 2 && return NaN
    first_i = length(history) - n + 1
    x_mean = (n - 1) / 2                    # mean of 0:(n-1)
    sxy = 0.0
    sxx = 0.0
    y_mean = sum(@view history[first_i:end]) / n
    for k in 0:(n - 1)
        dx = k - x_mean
        sxy += dx * (history[first_i + k] - y_mean)
        sxx += dx * dx
    end
    return sxy / sxx
end
