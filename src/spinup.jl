# Number of trailing cycles the drift criterion regresses over. Ten is enough for the slope
# to see through the seasonal limit cycle the column settles into without being so long that
# early transient cycles keep polluting the window.
const SPINUP_DRIFT_WINDOW = 10

"""
    _SpinupTracker(delta_tolerance, drift_tolerance)

Step-and-trend convergence bookkeeping for one scalar diagnostic of the column.

One of these per convergence *quantity* (column-mean density, firn air content), so the two
share this logic rather than each carrying its own copy of the history, the previous value, the
`NaN` sentinels and the abstention rule. Adding a third quantity is a third instance and one
`_update!` call, with no new branching in [`gemb_spinup`](@ref)'s loop.

Either tolerance may be `nothing`, meaning "not requested"; a tracker with both unset is inert
and [`_is_satisfied`](@ref) returns `true` for it, which is what makes the criteria conjunctive
without special-casing the unused ones.
"""
mutable struct _SpinupTracker
    delta_tolerance::Union{Nothing,Float64}
    drift_tolerance::Union{Nothing,Float64}
    history::Vector{Float64}
    previous::Union{Nothing,Float64}
    final_delta::Float64
    final_drift::Float64
end

_SpinupTracker(delta_tolerance, drift_tolerance) = _SpinupTracker(
    delta_tolerance === nothing ? nothing : Float64(delta_tolerance),
    drift_tolerance === nothing ? nothing : Float64(drift_tolerance),
    Float64[], nothing, NaN, NaN)

# Whether this quantity is being tested at all.
_is_checking(t::_SpinupTracker) =
    t.delta_tolerance !== nothing || t.drift_tolerance !== nothing

"""
    _update!(t::_SpinupTracker, value, drift_window)

Record this cycle's `value` and refresh the tracker's step and drift measures.

The step is computed unconditionally once there is a previous cycle to compare against, but the
drift regression is only run when its tolerance was requested *and* the window is full — the
slope of a partial window is not the slope the criterion is about.

A no-op for an inert tracker, so an unrequested quantity costs neither the `firn_air_content`
integral that produced it nor a growing history.
"""
function _update!(t::_SpinupTracker, value::Float64, drift_window::Int)
    _is_checking(t) || return t
    push!(t.history, value)
    if t.previous !== nothing
        t.final_delta = abs(value - t.previous)
    end
    t.previous = value
    if t.drift_tolerance !== nothing && length(t.history) >= drift_window
        t.final_drift = abs(_least_squares_slope(t.history, drift_window))
    end
    return t
end

"""
    _is_satisfied(t::_SpinupTracker) -> Bool

Whether every criterion this tracker was given currently holds. `true` for an inert tracker
(nothing was asked of it), and `false` while a requested measure is still `NaN` — an unmet
criterion must not pass by abstention.
"""
function _is_satisfied(t::_SpinupTracker)
    delta_ok = t.delta_tolerance === nothing ||
               (!isnan(t.final_delta) && t.final_delta < t.delta_tolerance)
    drift_ok = t.drift_tolerance === nothing ||
               (!isnan(t.final_drift) && t.final_drift < t.drift_tolerance)
    return delta_ok && drift_ok
end

"""
    gemb_spinup(profile, cf, mp; max_iterations=100,
                convergence_delta_density=nothing, convergence_drift_density=nothing,
                convergence_delta_fac=nothing, convergence_drift_fac=nothing,
                drift_window=$SPINUP_DRIFT_WINDOW, verbose=false)

Run GEMB for multiple spinup cycles to reach quasi-steady state.

Forces `output_frequency=:last` internally to minimize memory usage during spinup.
Returns the spun-up profile DimStack.

Convergence is judged on **column-mean density**, on whole-column **firn air content**, or on
both — see the criteria below. Either quantity is averaged or integrated over the whole column,
which is a fair comparison at every cycle because the column depth is fixed for the run (chosen
by [`initialize_profile`](@ref), held by [`trim_bottom!`](@ref)).

**Which quantity to converge is an application question, not a detail.** Column-mean density and
firn air content are not interchangeable tolerances on the same thing. They are related by
`FAC = Σ dz·(ρᵢ − ρ)/ρᵢ`, so at fixed total depth `Z`, `∂FAC/∂ρ̄ = −Z/ρᵢ`: **the FAC residual a
given density tolerance admits is proportional to the depth of the column.** Measured across a
21-site synthetic fleet, `convergence_delta_density = 1e-3` admits 0.02 mm of FAC drift in a
14 m column and 0.23 mm in a 212 m one — a 12x spread, loosest exactly at the deep, cold sites.
Work that reads FAC directly, such as altimetry requiring no residual trend in firn air, should
therefore converge on FAC and not rely on a density tolerance to imply it.

# Provenance
The returned profile carries a metadata `NamedTuple` (accessible via
`DimensionalData.metadata(profile)`) recording how the spinup ran and which
climatology it used: `spinup_cycles`, `spinup_converged`,
`spinup_final_delta_density`, `spinup_final_drift_density`,
`spinup_final_delta_fac`, `spinup_final_drift_fac`, the convergence
parameters (`spinup_max_iterations`, `spinup_convergence_delta_density`,
`spinup_convergence_drift_density`, `spinup_convergence_delta_fac`,
`spinup_convergence_drift_fac`, `spinup_drift_window`),
`spinup_smb_rate` (mean surface mass balance over the final cycle [m of ice per
year] from the surface fluxes, positive under accumulation — see
[`_cycle_smb_rate`](@ref); `NaN` if undeterminable), and the climatology
fields copied forward from `cf` (`climatology_window_start`,
`climatology_window_stop`, `climatology_n_years`, `climatology_steps_per_year`).
This provenance is propagated onto the [`gemb`](@ref) output when the spun-up
profile is used to start a transient run.

# Convergence criteria

Four independent tests, any subset of which may be requested. Every criterion given must
hold — they are conjunctive, so adding one can only make convergence stricter, never looser.
Two quantities (density, FAC) times two forms:

- **Step** (`convergence_delta_*`): exit when the absolute change between consecutive cycles
  is below this value. Cheap, but a column creeping steadily at just under the tolerance
  passes it while still drifting.
- **Drift** (`convergence_drift_*`): exit when the magnitude of the least-squares slope
  against cycle number, over the trailing `drift_window` cycles, is below this value. This is
  the trend test the step test cannot make: it distinguishes a column oscillating about a
  settled mean (slope ≈ 0) from one still evolving monotonically, and by using every point in
  the window it is not fooled by a single anomalous cycle. Inactive until `drift_window`
  cycles have run.

Units follow the quantity: density criteria are [kg/m³] and [kg/m³ per cycle], FAC criteria are
[m of air] and [m of air per cycle]. A useful FAC tolerance is therefore numerically much
smaller than a density one — millimetres of air, not thousandths of a kg/m³ (see the depth
scaling above).

When none is given the spinup always runs `max_iterations` cycles.

# Keyword arguments
- `max_iterations`: maximum number of spinup cycles (default 100). The spinup always
  exits after this many cycles even if convergence has not been reached.
- `convergence_delta_density`, `convergence_drift_density`: density criteria, see above.
- `convergence_delta_fac`, `convergence_drift_fac`: whole-column firn-air-content criteria
  ([`firn_air_content`](@ref)), see above.
- `drift_window`: trailing cycles the drift slope is fitted over (default
  $SPINUP_DRIFT_WINDOW). Must be ≥ 2. Shared by both drift criteria.
- `verbose`: print a convergence message when early exit occurs.
- `thermal_workspace`: reusable thermal-solve scratch, shared across every cycle of this
  spinup. One belongs to one column; give each concurrent thread its own
  (see [`ThermalWorkspace`](@ref)).

"""
function gemb_spinup(profile::DimStack, cf::ClimateForcing, mp::ModelParameters;
                     max_iterations::Int=100,
                     convergence_delta_density=nothing,
                     convergence_drift_density=nothing,
                     convergence_delta_fac=nothing,
                     convergence_drift_fac=nothing,
                     drift_window::Int=SPINUP_DRIFT_WINDOW,
                     verbose::Bool=false,
                     thermal_workspace::ThermalWorkspace=ThermalWorkspace())

    if (convergence_drift_density !== nothing || convergence_drift_fac !== nothing) &&
       drift_window < 2
        error("drift_window must be at least 2 to fit a slope, got $drift_window")
    end

    # Force output_frequency to :last for spinup efficiency
    mp_spinup = ModelParameters(;
        (field => getfield(mp, field) for field in fieldnames(ModelParameters) if field != :output_frequency)...,
        output_frequency=:last
    )

    current_profile = profile

    # One tracker per convergence quantity, so density and FAC share the step/drift
    # bookkeeping instead of duplicating it. Only the trailing `drift_window` entries of each
    # history are ever read, but a spinup is at most a few hundred cycles of one Float64 each,
    # so there is nothing to gain by ringbuffering them.
    rho_track = _SpinupTracker(convergence_delta_density, convergence_drift_density)
    fac_track = _SpinupTracker(convergence_delta_fac, convergence_drift_fac)

    # Convergence provenance, captured for the returned profile.
    cycles_run = 0
    converged  = false

    checking = _is_checking(rho_track) || _is_checking(fac_track)

    # Mean SMB of the most recent cycle, refreshed each pass so the equilibrated rate is to
    # hand after the loop. Recorded rather than holding the cycle's whole output stack: under
    # the forced `output_frequency=:last` each cycle is a single step, so the recompute is a
    # handful of scalar sums.
    smb_rate = NaN

    for cycle in 1:max_iterations
        cycles_run = cycle
        out = gemb(current_profile, cf, mp_spinup; thermal_workspace=thermal_workspace)
        current_profile = gemb_profile(out)
        smb_rate = _cycle_smb_rate(out, cf, mp)

        checking || continue

        _update!(rho_track, _column_mean_density(current_profile), drift_window)
        _update!(fac_track, _column_fac(current_profile, mp), drift_window)

        # Conjunctive across every requested criterion. Each reports `nothing` while it lacks
        # the cycles to judge, which is distinct from reporting "not converged": an unmet
        # criterion must not be able to pass by abstention, so a `nothing` blocks convergence.
        if _is_satisfied(rho_track) && _is_satisfied(fac_track)
            verbose && @info "gemb_spinup converged at cycle $cycle " *
                             "(Δρ = $(round(rho_track.final_delta, digits=4)) kg/m³, " *
                             "ρ drift = $(round(rho_track.final_drift, digits=6)) kg/m³/cycle, " *
                             "ΔFAC = $(round(fac_track.final_delta, digits=8)) m, " *
                             "FAC drift = $(round(fac_track.final_drift, digits=10)) m/cycle)"
            converged = true
            break
        end
    end

    @info "GEMB Spinup" climatology_window=(cf.climatology_window_start, cf.climatology_window_stop) cycles=cycles_run converged=converged final_delta_density=rho_track.final_delta final_drift_density=rho_track.final_drift final_delta_fac=fac_track.final_delta final_drift_fac=fac_track.final_drift

    return _attach_spinup_provenance(current_profile, cf;
        cycles=cycles_run, converged=converged,
        final_delta_density=rho_track.final_delta,
        final_drift_density=rho_track.final_drift,
        final_delta_fac=fac_track.final_delta,
        final_drift_fac=fac_track.final_drift,
        max_iterations=max_iterations,
        convergence_delta_density=convergence_delta_density,
        convergence_drift_density=convergence_drift_density,
        convergence_delta_fac=convergence_delta_fac,
        convergence_drift_fac=convergence_drift_fac,
        drift_window=drift_window,
        smb_rate=smb_rate)
end

"""
    _cycle_smb_rate(out, cf, mp) -> m of ice per year

Mean surface mass balance over the final spinup cycle, as metres of ice per year (positive
under net accumulation).

Derived from the **surface** mass fluxes alone:

    SMB = precipitation + evaporation_condensation - runoff

in kg m-2 over the cycle, divided by `mp.density_ice` and by the cycle length. Snowfall and
rain enter through `precipitation`, deposition and sublimation through the signed
`evaporation_condensation`, and meltwater leaving the column through `runoff`; refreezing is
internal and so does not appear. This is SMB by definition, and it is a property of the
climate rather than of the column's state — spinup changes it only through the small
feedback of albedo and surface temperature on melt.

It is deliberately **not** read from `ice_flux`. That output is the basal ice flux
the fixed-depth column performs, and its negation is surface *elevation* change, which is SMB
plus a firn-air-content term (see [`trim_bottom!`](@ref) for the identity): densification
drives it as surely as accumulation does, so it equals SMB only where compaction has
equilibrated, and in general it can even carry the opposite sign — at an accumulating
synthetic site with SMB of +0.135 m ice/yr the elevation trend was -0.157 m/yr, compaction
lowering the surface faster than accumulation raised it.

A spinup cycle is one climatological year by construction ([`forcing_climatology`](@ref)),
but the length is taken from the cycle's own step count and step length rather than assumed,
so a non-annual cycle gives a correct rate rather than a silently mis-scaled one.

Returns `NaN` when the rate cannot be determined (an empty or missing flux layer, a
degenerate cycle length, or a non-finite total), which callers treat as "unavailable" rather
than as a rate of zero.
"""
function _cycle_smb_rate(out, cf::ClimateForcing, mp::ModelParameters)
    isempty(dims(out, Ti)) && return NaN

    # `output_frequency=:last` leaves a single output step whose value is the whole cycle's
    # total, so the cycle length must come from the forcing rather than from the output times.
    # It is the number of steps times the step length, not the span from the first forcing time
    # to the last: the last step integrates a full `time_step` of its own, which a first-to-last
    # span omits, and dropping it would leave the rate biased high by one step in `n`.
    years = length(dims(cf, Ti)) * cf.time_step / SECONDS_PER_YEAR
    return _smb_rate(out, years, mp.density_ice)
end

"""
    _smb_rate(out, years, density_ice) -> m of ice per year

Mean surface mass balance of the `gemb` output `out` over a record of `years` years, as metres
of ice per year (positive under net accumulation). The flux algebra of
[`_cycle_smb_rate`](@ref), factored out so the spinup cycle and a transient run share one
definition of SMB — they differ only in which record they measure and how its length is known.

Returns `NaN` for a missing flux layer, a non-positive length, or a non-finite total.
"""
function _smb_rate(out, years::Real, density_ice::Real)
    for v in (:precipitation, :evaporation_condensation, :runoff)
        haskey(out, v) || return NaN
    end
    years > 0 || return NaN

    # Sum the whole record rather than taking its last element: these are per-output-step
    # accumulations, not running totals, so under `:last` (one step holding the cycle total)
    # and under any finer frequency alike the sum is the record total.
    tot(v) = sum(x -> isfinite(x) ? x : 0.0, parent(out[v]))
    smb_mass = tot(:precipitation) + tot(:evaporation_condensation) - tot(:runoff)
    isfinite(smb_mass) || return NaN

    return smb_mass / density_ice / years
end

# Attach spinup + climatology provenance to the spun-up profile so it is a
# self-describing artifact. Uses plain `DimStack` metadata (free-form NamedTuple);
# climatology fields are copied forward from the forcing `cf` when present.
function _attach_spinup_provenance(profile::DimStack, cf::ClimateForcing;
        cycles, converged, final_delta_density, final_drift_density,
        max_iterations, convergence_delta_density, convergence_drift_density,
        drift_window, smb_rate=NaN,
        final_delta_fac=NaN, final_drift_fac=NaN,
        convergence_delta_fac=nothing, convergence_drift_fac=nothing)
    cf_meta = DD.metadata(cf)
    _cf(key) = haskey(cf_meta, key) ? cf_meta[key] : nothing
    prov = (
        spinup_cycles = cycles,
        spinup_converged = converged,
        spinup_final_delta_density = final_delta_density,
        spinup_final_drift_density = final_drift_density,
        # FAC measures [m of air] and [m of air per cycle]. `NaN` when the FAC criteria were
        # not requested, which is distinct from a measured zero.
        spinup_final_delta_fac = final_delta_fac,
        spinup_final_drift_fac = final_drift_fac,
        spinup_max_iterations = max_iterations,
        spinup_convergence_delta_density = convergence_delta_density,
        spinup_convergence_drift_density = convergence_drift_density,
        spinup_convergence_delta_fac = convergence_delta_fac,
        spinup_convergence_drift_fac = convergence_drift_fac,
        spinup_drift_window = drift_window,
        # Equilibrated mean surface mass balance [m of ice per year], positive under accumulation.
        # Read by `plot_output`'s detrended height axis; `NaN` when undeterminable.
        spinup_smb_rate = smb_rate,
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

# Whole-column firn air content [m of air], via the same `firn_air_content` the `gemb` output
# reports, so the convergence criterion and the diagnostic cannot disagree about what FAC means.
# No depth cutoff: the column depth is fixed for the run, so the whole column is already a
# consistent domain across cycles, and a cutoff would answer a different question (see the
# depth-scaling note in `gemb_spinup`).
function _column_fac(profile::DimStack, mp::ModelParameters)
    return firn_air_content(parent(profile[:dz]), parent(profile[:density]), mp.density_ice)
end

# Least-squares slope of the last `window` history entries against cycle number, in the
# history's own units per cycle. Quantity-agnostic: density [kg/m³ per cycle] and firn air
# content [m per cycle] use the same regression.
#
# Cycles are equally spaced, so the design matrix is known a priori and the slope reduces to a
# closed form: with x centered on the window, x̄ = 0 and the slope is Σxᵢyᵢ / Σxᵢ². No
# allocation, no least-squares solve.
function _least_squares_slope(history::Vector{Float64}, window::Int)
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
