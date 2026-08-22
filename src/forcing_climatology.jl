"""
    forcing_climatology(cf::ClimateForcing; method=:average, kwargs...)
    forcing_climatology(cf, datetime_range::Tuple{DateTime,DateTime}; kwargs...)
    forcing_climatology(cf, selector; kwargs...)

Build a one-year repeating forcing cycle from a [`ClimateForcing`](@ref), for
[`gemb_spinup`](@ref).

Both methods first subset to the requested window (if any), drop leap days (day 366), and
discard partial years. They differ in what they do with the complete years that remain:

- `method = :average` (default) — average every field across the years, slot by slot. The
  historical behaviour, and bit-identical to it.
- `method = :representative` — select a block of `n_years` consecutive **real** years whose mean
  melt matches the record's, and repeat that block. Requires `model_parameters`; see below.

Returns a new `ClimateForcing` holding one year under `:average`, or `n_years` under
`:representative`. The time step and scalar metadata
(`temperature_air_mean`, `wind_speed_mean`, `precipitation_mean`,
`temperature_observation_height`, `wind_observation_height`) are carried forward unchanged
under both methods.

The window can be given as a `(start, stop)` `DateTime` tuple, or — since `ClimateForcing` is
an `AbstractDimStack` — as any DimensionalData `Ti` selector, e.g. a closed interval
`DateTime(1950,1,1) .. DateTime(1980,12,31)`, `At`, `Near`, or an index range. Both forms
subset via `cf[Ti(...)]` first.

# Which method to use: averaging destroys melt

`:average` preserves the mean of every field but shrinks its variance by roughly `1/n_years`,
and melt is not a linear function of the forcing. It is **rectified** — identically zero until
the surface energy balance drives the skin to the melt point, and growing only above it — so
averaging pairs each year's warm excursions against other years' cold ones in the same slot,
and the peaks that carried the melt stop crossing the threshold. This is Jensen's inequality on
a convex, one-sided response: the mean of the melt is not the melt of the mean, and for a
threshold process the loss is *total*, not a percentage.

Measured on **real ERA5-Land forcing** at 12 glacierized Greenland sites, 2000-2019
(`bench/greenland_spinup_forcing.jl`), as melt on the cycle over melt on the full record. Report
the per-site distribution, not the aggregate: the aggregate is dominated by the highest-melt
sites and hides the ones where averaging fails completely.

| method | melt retained, median | sites below 10% | accumulation error, median |
|---|---|---|---|
| `:average` | **12%** | **5 of 10** | 6.2% |
| `:representative`, `n_years=1` | 104% | 0 of 10 | 18.5% |
| `:representative`, `n_years=3` | 102% | 1 of 10 | 6.9% |

Averaging's failure tracks **elevation and absolute melt magnitude, not the melt/accumulation
ratio**: retention is 0% at Saddle (2456 m, 23.5 kg m-2 yr-1), 0% at DYE-2 (2094 m, 81.6), 24% at
KAN-M (1300 m, 857.5) and 86-115% at the lowest three sites (271-577 m). High sites melt from the
*tail* of the distribution, which averaging removes; low sites melt from the mean, which it
preserves. A high melting site spun up under `:average` therefore equilibrates as a melt-free
site, with its firn air content and refreeze terms wrong by construction — which matters most for
work reading those terms directly, such as altimetry.

!!! note "An earlier version of this docstring claimed `:average` recovers 0% of melt"
    That was measured on a synthetic fleet derived from a single parent site, whose one seasonal
    cycle plus noise averages away almost completely. Real forcing has a coherent seasonal cycle
    that survives averaging at low elevations, so the true figure is a median of 12% and not zero.
    The synthetic 21-site "critical band" of melt/accumulation 0.3-1.0 did not reproduce either:
    only one real site landed in it, and no real site lost its firn there.

# Which method to use

- **Negligible melt: use `:average`.** It preserves accumulation essentially exactly (0.1-0.5% at
  the three dry Greenland sites, against 22% and 15% errors for a melt-ranked block there), and it
  integrates less than half the model-years. This is why it remains the default. Averaging is not
  *exactly* neutral even here — the densification rate carries a convex Arrhenius factor — but the
  measured temperature offset of both methods was 0.00 K at every site, so that channel is inert.
- **Significant melt: use `:representative`.** Melt is recovered to a median 102% (`n_years=3`)
  against 12% for averaging.
- **Ablation (melt above accumulation): use `:representative`.** The column equilibrates to bare
  ice, which is correct there — accumulation resets each year and what matters is the mean melt
  driving latent-heat warming, which `:average` can set to zero.

A regime-specific composite (average / SMB-ranked / melt-ranked by zone) was measured and is
**worse than simply using `n_years=3` everywhere** (median joint error 7.9% against 4.8%), so it
is not offered.

# What this does not tell you

Every number above is **forcing fidelity** — how well the cycle reproduces the record's melt and
accumulation. The equilibrated-column comparison that would close the loop was attempted and
abandoned: under strict convergence criteria 6 of these 12 sites do not converge within 800
cycles, and a multi-variant sweep ran for over an hour and 4000+ model-years without finishing
(see `bench/greenland_block_length.jl`). So it is not established that better forcing fidelity
translates proportionally into a better spun-up column.

`:representative` also costs about **2.2x the model-years** of `:average` (13116 against 6096 over
these 12 sites). Cycle counts are not comparable between the methods — `:average` integrates one
year per cycle and `:representative` integrates `n_years`, so a longer cycle lowers the count
without lowering the work.

# Keyword arguments

- `method`: `:average` or `:representative` (see above).
- `model_parameters`: a [`ModelParameters`](@ref). **Required** for `:representative`, which
  needs the physics to rank candidate years; ignored by `:average`.
- `n_years`: for `:representative`, the block length (default 3). Clamped to the number of
  complete years available, with a warning. 3 is the measured optimum on the Greenland sites: a
  single year matches melt marginally better (median 104% against 102%) but is a noisy sample of
  accumulation (median error 18.5% against 6.9%, worst 37.8%), and on a joint melt-and-accumulation
  score 3 years beat both 1 and 5 (median 4.8% against 10.9% and 5.9%). **Cycling the whole record
  is deliberately not offered**: it would give one convergence comparison per record length, so a
  10-cycle drift window would span 10 records.
- `statistic`: for `:representative`, which year to call representative — `:mean` (default) or
  `:median` of the candidate years' score. `:mean` is the default because the column equilibrates
  to a long-run mean and the annual melt distribution is right-skewed (measured at one synthetic
  site: mean 150.8 against median 127.3 kg m-2 yr-1), so selecting on the median biases the spinup
  low. Measured on the Greenland sites, `:median` was worse than `:mean` under SMB ranking
  (median joint error 6.1% against 4.7%).
- `rank_by`: for `:representative`, what each candidate year is scored on.
  - `:model` (default) — annual melt, integrated by running one year of [`gemb`](@ref) per
    candidate. The default because melt is the term averaging destroys, and because it is cheap:
    32 candidate years cost 2.4 s (0.074 s/yr) against the hundreds of cycles the spinup it feeds
    will run.
  - `:smb` — annual surface mass balance instead, via the same `_smb_rate` [`gemb_spinup`](@ref)
    reports. Combines accumulation and melt with the weighting the mass balance gives them, so it
    matches accumulation far better (median 2.5% against 6.9% for `:model` over 12 Greenland
    sites) at some cost in melt fidelity (median 95% against 102%). The cost concentrates where
    accumulation dominates SMB: at a percolation site with 485 kg m-2 yr-1 of accumulation against
    55 of melt, an SMB-matched year says little about melt and recovered only 33% of it. Prefer it
    when accumulation and firn air content matter more than melt.
  - `:estimate` — the surface-energy-balance melt from [`initialize_climate_summary`](@ref), which
    needs no model integration and so is much the cheapest. Its melt is biased ~2.6x high, and
    while a uniform bias cancels when *ordering* years, it does not cancel when picking the year
    nearest a target value, so prefer `:model` unless the ranking cost matters.
- `fallback_accumulation_tolerance`: guard on the selected block's accumulation, as a fraction
  (default `0.05`); `nothing` disables it. Melt-ranking optimizes melt and leaves accumulation to
  chance — measured, that chance is sometimes bad, missing the record's accumulation by 30.2% at
  QAS-L and 18.0% at KAN-U. When the block exceeds this tolerance the selection is re-run with
  `rank_by = :smb`, and the result is kept only if it actually improves the accumulation match.
  Ignored under `rank_by = :smb`, which is already the fallback. Measured at the default tolerance
  over 9 melting Greenland sites: median joint (melt, accumulation) error 4.0% -> 3.5%, firing at
  4 sites — a small average gain whose value is bounding the 30% tail rather than shifting the
  median.
- `profile`: for `:representative`, the column the candidate years are integrated on. Defaults to
  [`initialize_profile`](@ref)`(model_parameters, cf)`. Only the *ranking* depends on it, never the
  forcing returned. Built even under `rank_by = :estimate`, which does not use it, because the
  accumulation fallback may re-select on `:smb`, which does.
- `verbose`: for `:representative`, report the selected year and how typical it is (default
  `true`). See the caveat below.

# The short-block caveat

A block is a handful of real weather sequences repeated for the whole spinup, and scoring it on
one variable cannot make it typical in every other. The measured failure mode is **accumulation**:
melt-ranking left the block 30.2% off the record's accumulation at QAS-L and 18.0% at KAN-U, which
`fallback_accumulation_tolerance` exists to bound. Temperature, by contrast, was *not* a problem —
the mean offset of both methods was 0.00 K at all 12 Greenland sites.

`verbose=true` reports the chosen block, its score against the target, and its accumulation and
temperature offsets, and warns when either is large. The averaged year is equally a single
sequence, with the difference that it is not one that ever occurred.

# Provenance

Both methods set `climatology_window_start`/`_stop`, `climatology_n_years` and
`climatology_steps_per_year`. `:representative` additionally sets
`climatology_representative_year` to the calendar year chosen; it is absent for `:average`, so
the two are distinguishable after the fact.

# Examples
```julia
cf = initialize_forcing(ds)
mp = initialize_parameters()

cycle = forcing_climatology(cf)                                   # averaged (default)
cycle = forcing_climatology(cf; method=:representative,           # one real year
                            model_parameters=mp)
```
"""
function forcing_climatology(cf::ClimateForcing, datetime_range::Tuple{DateTime,DateTime};
                             kwargs...)
    # Record the requested window as climatology provenance on the result.
    return forcing_climatology(cf[Ti(datetime_range[1] .. datetime_range[2])];
                               window=datetime_range, kwargs...)
end

# Idiomatic DimensionalData form: forward any `Ti` selector (interval, `At`, `Near`,
# index range, ...) through the stack's own indexing, then average. The tuple method
# above is strictly more specific, so it always wins for `(start, stop)` calls.
function forcing_climatology(cf::ClimateForcing, selector; kwargs...)
    return forcing_climatology(cf[Ti(selector)]; kwargs...)
end

"""
    _complete_years(cf::ClimateForcing) -> NamedTuple

Identify the complete, leap-day-free years of `cf` — the common first step of both
[`forcing_climatology`](@ref) methods.

Returns `non_leap` (a mask over the original series), `forcing_index` (a mask over the
leap-day-filtered series selecting the complete years), `times_complete`, `complete_years`,
`n_complete_years` and `steps_per_year`.

A "complete" year is one holding the maximum number of steps any year holds, which identifies
partial years at either end of a record without assuming a step length. Leap day (day 366) is
dropped first so that every year has the same step count, which is what lets `:average` reshape
on `steps_per_year` and `:representative` slice a single year of exactly that length.
"""
function _complete_years(cf::ClimateForcing)
    times = collect(lookup(dims(cf.temperature_air, Ti)))

    # Eliminate leap days (day 366 of the year)
    non_leap = [Dates.dayofyear(t) != 366 for t in times]
    times_noleap = times[non_leap]

    # Count timesteps per year
    years_all = Dates.year.(times_noleap)
    unique_years = sort(unique(years_all))
    counts_per_year = [count(==(yr), years_all) for yr in unique_years]

    # Find years with the maximum number of entries (complete years)
    max_count = maximum(counts_per_year)
    complete_years = unique_years[counts_per_year .== max_count]

    # Indices (within the leap-day-filtered series) for complete years only
    forcing_index = [yr in complete_years for yr in years_all]

    # `complete_idx` composes both masks into one index vector over the *original* series, so a
    # slice is a single gather rather than two record-sized intermediate copies (see
    # `_years_slice`). `years_complete` is the calendar year per retained step, precomputed
    # because every slice would otherwise re-broadcast `Dates.year` over the whole series.
    complete_idx = findall(non_leap)[forcing_index]

    return (non_leap=non_leap, forcing_index=forcing_index,
            complete_idx=complete_idx,
            times_complete=times_noleap[forcing_index],
            years_complete=Dates.year.(times_noleap[forcing_index]),
            complete_years=complete_years,
            n_complete_years=length(complete_years),
            steps_per_year=max_count)
end

# `window` is the requested `(start, stop)` averaging window recorded as
# provenance; when `nothing` (e.g. called on an already-subset stack) it falls
# back to the extent of the complete years actually averaged.
function forcing_climatology(cf::ClimateForcing; window=nothing,
                             method::Symbol=:average,
                             model_parameters=nothing,
                             n_years::Integer=3,
                             statistic::Symbol=:mean,
                             rank_by::Symbol=:model,
                             fallback_accumulation_tolerance=0.05,
                             profile=nothing,
                             verbose::Bool=true)

    if method === :representative
        model_parameters === nothing && throw(ArgumentError(
            "forcing_climatology(method=:representative) requires `model_parameters`: " *
            "ranking candidate years needs the physics that turns forcing into melt."))
        return _representative_block(cf, model_parameters; n_years=n_years, window=window,
                                     statistic=statistic, rank_by=rank_by,
                                     fallback_accumulation_tolerance=fallback_accumulation_tolerance,
                                     profile=profile, verbose=verbose)
    elseif method !== :average
        throw(ArgumentError(
            "method must be :average or :representative, got $(repr(method))"))
    end

    cy = _complete_years(cf)
    non_leap = cy.non_leap
    forcing_index = cy.forcing_index
    n_complete_years = cy.n_complete_years
    steps_per_year = cy.steps_per_year
    times_complete = cy.times_complete

    # Reshape into (steps_per_year × n_years) and average across years
    reshape_avg(arr) = vec(Statistics.mean(
        reshape(arr[non_leap][forcing_index], steps_per_year, n_complete_years), dims=2))

    # Use times from the first complete year
    tdim = Ti(times_complete[1:steps_per_year])
    _avg(a) = DimArray(reshape_avg(parent(a)), (tdim,))

    # Climatology provenance: the requested window (or, when not given, the extent
    # of the complete years actually averaged) plus the year / step counts.
    window_start = window === nothing ? first(times_complete) : window[1]
    window_stop  = window === nothing ? last(times_complete)  : window[2]

    return ClimateForcing(
        _avg(cf.temperature_air),
        _avg(cf.pressure_air),
        _avg(cf.precipitation),
        _avg(cf.wind_speed),
        _avg(cf.shortwave_downward),
        _avg(cf.longwave_downward),
        _avg(cf.vapor_pressure),
        _avg(cf.black_carbon_snow),
        _avg(cf.black_carbon_ice),
        _avg(cf.cloud_optical_thickness),
        _avg(cf.solar_zenith_angle),
        _avg(cf.shortwave_downward_diffuse),
        _avg(cf.cloud_fraction),
        cf.time_step,
        cf.temperature_air_mean,
        cf.wind_speed_mean,
        cf.precipitation_mean,
        cf.temperature_observation_height,
        cf.wind_observation_height;
        snow_drift=_avg(cf.snow_drift),
        accumulation_mean=cf.accumulation_mean,
        temperature_air_effective=cf.temperature_air_effective,
        dataset=cf.dataset,
        latitude=cf.latitude,
        longitude=cf.longitude,
        elevation=cf.elevation,
        elevation_native=cf.elevation_native,
        elevation_offset=cf.elevation_offset,
        climatology_window_start=window_start,
        climatology_window_stop=window_stop,
        climatology_n_years=n_complete_years,
        climatology_steps_per_year=steps_per_year,
    )
end

"""
    _representative_block(cf, mp; n_years, window, statistic, rank_by, profile, verbose)

The `method = :representative` path of [`forcing_climatology`](@ref): return a **block of
`n_years` consecutive real years** whose mean melt best matches the record, unmodified, so every
field's variance, cross-field covariance and day-to-day sequencing are intact.

## Why a block rather than one year — and what it does not fix

A block carries some interannual variability while staying short enough that `gemb_spinup`'s
per-cycle convergence criteria are still judged often. That is the only claim made for it.

It was introduced on the hypothesis that alternating warm and cool years would stop a repeated
cycle from melting a firn column away, and **the measurement refuted that**: on the fleet, 3- and
5-year blocks lost the firn at as many sites as a single year (4 of 21, equilibrated FAC below
1 m, against 0 for `:average`), and the two marginal percolation sites went to ice at every block
length tried. Where mean melt approaches mean accumulation, a short cycle of real years cannot
sustain the column, and lengthening it within the range that keeps convergence judgeable does not
help. See [`_warn_cycle_melt_ratio`](@ref).

What a block **does** buy, measured on real Greenland forcing, is accumulation fidelity. A single
year matches melt marginally better (median 104% against 102%) but samples accumulation badly
(median error 18.5%, worst 37.8%, against 6.9% and 30.2% for three years). On a joint
melt-and-accumulation score, 3 years beat both 1 and 5 (median 4.8% against 10.9% and 5.9%). That
is the case for the default, and it is about accumulation rather than about sustaining the column.

**Consecutive** years, not the `n` individually-closest ones: a hand-assembled set would break
the transition between successive years — the carry-over of a warm summer into the following
winter — which is exactly the multi-year memory the block exists to preserve.

## Selection

Each candidate window of `n_years` consecutive years is scored by the mean of its years' annual
melt (see `rank_by` in [`forcing_climatology`](@ref)), and the window closest to the `statistic`
of *all* candidate years' melt is returned. Melt is the ranking variable because it is the term
averaging destroys outright; the other fields survive averaging far better (see the regime table
in [`forcing_climatology`](@ref)).
"""
function _representative_block(cf::ClimateForcing, mp::ModelParameters;
                               n_years::Integer=3, window=nothing,
                               statistic::Symbol=:mean, rank_by::Symbol=:model,
                               fallback_accumulation_tolerance=0.05,
                               profile=nothing, verbose::Bool=true)

    statistic in (:mean, :median) ||
        throw(ArgumentError("statistic must be :mean or :median, got $(repr(statistic))"))
    rank_by in (:smb, :model, :estimate) ||
        throw(ArgumentError("rank_by must be :smb, :model or :estimate, got $(repr(rank_by))"))
    n_years >= 1 ||
        throw(ArgumentError("n_years must be at least 1, got $n_years"))

    cy = _complete_years(cf)
    isempty(cy.complete_years) &&
        error("forcing_climatology(method=:representative): forcing has no complete years")

    years = cy.complete_years
    n_block = min(Int(n_years), length(years))
    if n_block < n_years
        @warn "forcing_climatology: only $(length(years)) complete years available; " *
              "using a $(n_block)-year block instead of the requested $n_years."
    end

    # The ranking column. Built unless the caller supplied one, and needed even under
    # `rank_by = :estimate` — which does not use it itself — because the accumulation guard below
    # may re-select on `:smb`, which does. Building it is cheap next to the spinup this feeds.
    # The record's climate summary, computed once. The guard needs its accumulation, and
    # `_report_representative_block` its melt — and building the ranking column needs the same
    # summary internally, so it is built here and handed to `initialize_profile` rather than
    # letting that recompute it. Measured on a 32-year 3-hourly record: 346 ms per summary, of a
    # 3020 ms total, so the duplicate was ~11% of the whole call.
    #
    # Deliberately *not* short-cut to `cf.accumulation_mean`: that is derived from total
    # precipitation and is 3x larger here (1177 against 386 kg m-2 yr-1), and a sliced block
    # carries the record's value through unchanged, so comparing block metadata to record metadata
    # would compare a number to itself.
    record_summary = initialize_climate_summary(cf, mp)
    col = profile === nothing ? initialize_profile(mp, cf; climate_summary=record_summary) :
          profile

    # `(melt, smb)` per candidate year, from one `gemb` pass each, computed once and shared by
    # both the initial selection and the fallback's re-selection. One thermal workspace is reused
    # across all years, as `gemb_spinup` does across cycles. Skipped entirely under `:estimate`,
    # which integrates nothing — but the fallback may still need SMB, so it is computed lazily
    # there rather than never.
    mp_last = ModelParameters(;
        (f => getfield(mp, f) for f in fieldnames(ModelParameters)
         if f != :output_frequency)..., output_frequency=:last)
    ws = ThermalWorkspace()
    year_scores = rank_by === :estimate ? nothing :
        [_year_scores(col, _year_slice(cf, cy, y), mp_last, ws) for y in years]

    sel = _select_block(cf, cy, years, n_block, mp, statistic, rank_by, col, year_scores)

    # Accumulation guard. Melt-ranking optimizes melt and leaves accumulation to chance, and on
    # real Greenland forcing that chance is sometimes bad: the melt-ranked 3-year block missed
    # the record's accumulation by 30.2% at QAS-L and 18.0% at KAN-U. Since accumulation drives
    # the firn column directly, a block that far off is a poor cycle whatever its melt fidelity.
    #
    # The re-selection ranks on SMB instead, which matched accumulation to ~2.5% median across
    # the same sites because SMB carries accumulation with the weighting the mass balance gives
    # it. Measured effect of the guard at a 5% tolerance over 9 melting sites: median joint
    # (melt, accumulation) error 4.0% -> 3.5%, median accumulation error 4.1% -> 3.6%, firing at
    # 4 sites. A small average gain whose real value is bounding the 30% tail.
    #
    # **The trigger is the accumulation error itself, not the block's SMB bias.** SMB bias was
    # tried first and is the wrong signal: SMB = accumulation - melt, so at an ablation site it is
    # dominated by melt and says little about accumulation. Measured, the two are effectively
    # uncorrelated — SwissCamp had the largest SMB bias of any site (-52.2%) with a 0.7%
    # accumulation error, while QAS-L had the worst accumulation error (30.2%) at -16.3% SMB bias.
    # An SMB-triggered guard fired on the wrong sites and made the median *worse* at every
    # threshold tried.
    # The cycle's own summary, for the guard, the diagnostics and the warning. The record's was
    # computed above.
    cycle_summary = initialize_climate_summary(sel.cycle, mp)
    acc_record = record_summary.accumulation
    effective_rank = rank_by
    accumulation_error = acc_record > 0 ?
        abs(cycle_summary.accumulation / acc_record - 1) : NaN

    if fallback_accumulation_tolerance !== nothing && rank_by !== :smb &&
       isfinite(accumulation_error) &&
       accumulation_error > Float64(fallback_accumulation_tolerance)

        # Reuses `year_scores`, so the re-selection costs no further integrations — only under
        # `:estimate`, which computed none, is the scoring pass paid here.
        alt_scores = year_scores === nothing ?
            [_year_scores(col, _year_slice(cf, cy, y), mp_last, ws) for y in years] : year_scores
        alt = _select_block(cf, cy, years, n_block, mp, statistic, :smb, col, alt_scores)
        alt_summary = initialize_climate_summary(alt.cycle, mp)
        alt_error = abs(alt_summary.accumulation / acc_record - 1)
        # Only accept the fallback if it actually improves the quantity that triggered it. SMB
        # ranking is better on accumulation *on average*, not at every site, and a guard that can
        # make its own trigger worse is not a guard.
        if alt_error < accumulation_error
            verbose && @info "forcing_climatology: melt-ranked block accumulation is " *
                "$(round(100 * accumulation_error, digits=1))% off the record (tolerance " *
                "$(round(100 * Float64(fallback_accumulation_tolerance), digits=1))%); " *
                "re-selected on SMB (now $(round(100 * alt_error, digits=1))%). " *
                "Set `fallback_accumulation_tolerance=nothing` to disable."
            sel, cycle_summary, accumulation_error = alt, alt_summary, alt_error
            effective_rank = :smb
        end
    end

    if verbose
        _report_representative_block(sel, record_summary, cycle_summary, accumulation_error,
                                     statistic, effective_rank)
        _warn_cycle_melt_ratio(cycle_summary, sel.block_years)
    end

    return _stamp_representative(sel.cycle, sel.block_years, cy.steps_per_year, window;
                                 rank_by=effective_rank,
                                 accumulation_error=accumulation_error,
                                 melt_ratio=cycle_summary.accumulation > 0 ?
                                     cycle_summary.melt / cycle_summary.accumulation : NaN)
end

"""
    _select_block(cf, cy, years, n_block, mp, statistic, rank_by, col) -> NamedTuple

Score every candidate year, then return the `n_block`-long consecutive run whose mean score is
closest to the `statistic` of the per-year scores.

Factored out of [`_representative_block`](@ref) because the accumulation guard re-runs the whole
selection under a different `rank_by`, and a guard that duplicated the scoring logic could drift
from it. Returns the chosen `cycle`, its `block_years`, the per-year `scores`, the chosen block's
own mean `block_score`, and the `target` it was matched against.
"""
function _select_block(cf::ClimateForcing, cy, years, n_block::Int, mp::ModelParameters,
                       statistic::Symbol, rank_by::Symbol, col, year_scores)
    # The ranking score per candidate year. `:model` and `:smb` both read `year_scores`, which
    # holds `(melt, smb)` per year from a single `gemb` pass each — the two are different
    # reductions of the same integration, so scoring one after the other would double the
    # model-years for nothing. That matters because the accumulation guard re-selects under
    # `:smb`. `:estimate` instead reads the summary's surface-energy-balance melt, which needs no
    # integration at all but is biased ~2.6x high.
    scores = if rank_by === :estimate
        [initialize_climate_summary(_year_slice(cf, cy, y), mp).melt for y in years]
    else
        take = rank_by === :smb ? last : first
        [take(year_scores[i]) for i in eachindex(years)]
    end

    # Target is the statistic over every candidate *year*, so it describes the record rather
    # than the blocks; each block is then scored by its own mean against it.
    target = statistic === :mean ? Statistics.mean(scores) : Statistics.median(scores)
    n_windows = length(years) - n_block + 1
    block_scores = [Statistics.mean(@view scores[i:(i + n_block - 1)]) for i in 1:n_windows]
    # Non-allocating: the broadcast form would build two temporary vectors to pick one index.
    start = argmin(i -> abs(block_scores[i] - target), eachindex(block_scores))
    block_years = years[start:(start + n_block - 1)]

    return (cycle = _years_slice(cf, cy, block_years), block_years = block_years,
            scores = scores, block_score = block_scores[start], target = target)
end

# Melt-to-accumulation ratio of a constructed cycle above which `_warn_cycle_melt_ratio` speaks.
# Deliberately below 1.0: the useful signal is "approaching the threshold", because a cycle whose
# *mean* ratio is 0.8 can contain individual years above 1.
const CYCLE_MELT_RATIO_WARN = 0.6

"""
    _warn_cycle_melt_ratio(cycle_summary::ClimateSummary, block_years)

Warn when the selected cycle's own melt-to-accumulation ratio puts it at or past the point where
repeating it drives the column to solid ice.

The test is on the **constructed cycle**, not on the record, and that distinction is the whole
point. `melt/accumulation > 1` over a record means the site is ablating, and a cycle that
equilibrates to ice is then the *correct* answer. But the record's ratio is a mean, and
interannual variability puts individual years above 1 while the mean sits below it — so a site
whose record ratio is 0.9 can still yield a cycle that ablates.

Not an error, and not a refusal to return the cycle: at a genuinely ablating site the ice column
*is* the answer, and what matters there is that the mean melt — and so the latent-heat warming
it drives — is right, which is exactly what this path gets right and `:average` does not.

The ratio that motivated the sub-1.0 threshold was measured on the **synthetic** fleet, where two
sites with record ratios of 0.83 and 0.93 lost their firn to a repeated block. That firn loss did
**not** reproduce on real Greenland forcing (see [`forcing_climatology`](@ref)), so the threshold
is a conservative guard rather than a calibrated boundary.
"""

# `cycle_summary` is a `ClimateSummary`, left unannotated because
# `initialize_climate_summary.jl` (which defines that type) is included after this file.
function _warn_cycle_melt_ratio(cycle_summary, block_years)
    acc = cycle_summary.accumulation
    acc > 0.0 || return nothing

    # The ratio is formed from the cycle's *own* summary melt, deliberately not from the block's
    # selection score. Those are the same quantity only under `rank_by = :model`: under `:smb`
    # the score is a surface mass balance in m of ice per year, and under `:estimate` it is a
    # melt biased ~2.6x high. Dividing either by an accumulation in kg m-2 and comparing to
    # `CYCLE_MELT_RATIO_WARN` is meaningless — and it would fail in the worst possible way,
    # because the accumulation guard switches the effective ranking to `:smb` at precisely the
    # melting sites this warning exists for, so the units break would silence it exactly where
    # it is needed. Reading melt from the summary makes the test correct under every ranking.
    ratio = cycle_summary.melt / acc
    ratio < CYCLE_MELT_RATIO_WARN && return nothing

    if ratio > 1.0
        @warn "Representative cycle $(block_years) melts more than it accumulates " *
              "(melt/accumulation = $(round(ratio, digits=2))): this is the ablation regime, " *
              "and a spinup on it will equilibrate to bare ice. That is the physically " *
              "correct outcome — accumulation resets annually there and what matters is the " *
              "mean melt driving latent-heat warming, which this cycle reproduces and an " *
              "averaged climatology sets to zero. Reported so the ice column is not a surprise."
    else
        @warn "Representative cycle $(block_years) has melt/accumulation = " *
              "$(round(ratio, digits=2)), close enough to 1 that repeating it may drive the " *
              "column to bare ice even though the site accumulates on average. Interannual " *
              "variability puts individual years above 1 while the mean sits below it. " *
              "Compare the equilibrated firn air content against a `method=:average` spinup " *
              "before trusting either." maxlog=2
    end
    return nothing
end

"""
    _year_scores(col, year_forcing, mp_last) -> (melt, smb)

Annual melt [kg m-2 yr-1] and surface mass balance [m of ice per year] of one candidate year,
integrated by the model itself on `col`.

Both from **one** `gemb` pass. They are different reductions of the same integration — melt sums
the `melt` flux, SMB is `precipitation + evaporation_condensation - runoff` via
[`_smb_rate`](@ref) — so computing them separately would double the model-years spent ranking,
which the accumulation fallback would then double again when it re-selects on SMB.

`output_frequency=:last` leaves a single step holding the year's total, so the sum is that total,
and the year length comes from the forcing rather than from the output times.

SMB shares `_smb_rate` with [`gemb_spinup`](@ref)'s own diagnostic, so the ranking and the
reported `spinup_smb_rate` cannot disagree about what SMB means.
"""
function _year_scores(col, year_forcing::ClimateForcing, mp_last::ModelParameters,
                      ws::ThermalWorkspace)::Tuple{Float64,Float64}
    # Return type annotated rather than left to inference. `gemb`'s output stack and the `col`
    # profile are untyped here, so both reductions infer as `Any` and the tuple with them; the
    # annotation makes what callers see concrete without narrowing either argument. Both
    # quantities are `Float64` by construction (`_smb_rate` returns one, and the melt sum is over
    # a `Float64` layer), so this converts nothing.
    out = gemb(col, year_forcing, mp_last; thermal_workspace=ws)
    years = length(dims(year_forcing, Ti)) * year_forcing.time_step / SECONDS_PER_YEAR
    years > 0 || return (NaN, NaN)
    return (sum(parent(out[:melt])) / years, _smb_rate(out, years, mp_last.density_ice))
end

# One or more complete years of `cf` as a standalone `ClimateForcing`, carrying every layer and
# all scalar metadata forward. Sliced from the leap-day-filtered, complete-year series so the
# result is exactly `n * steps_per_year` long and can serve as a spinup cycle.
#
# The `years` method takes a *consecutive* run, so the slice is contiguous in time and the
# transition from one year into the next is the real one — the multi-year memory a block exists
# to preserve (see `_representative_block`).
_year_slice(cf::ClimateForcing, cy, year::Integer) = _years_slice(cf, cy, (year,))

function _years_slice(cf::ClimateForcing, cy, years)
    # One composed index vector, then one gather per layer. Indexing with the two masks in
    # succession (`parent(a)[non_leap][forcing_index][keep]`) would allocate two *record-sized*
    # temporaries per layer before producing the small result — 14 layers x 3 allocations per
    # slice, and a slice happens once per candidate year plus once per selected block.
    keep = findall(in(Set(Int.(years))), cy.years_complete)
    idx = cy.complete_idx[keep]
    tdim = Ti(cy.times_complete[keep])
    _slice(a) = DimArray(parent(a)[idx], (tdim,))

    return ClimateForcing(
        _slice(cf.temperature_air), _slice(cf.pressure_air), _slice(cf.precipitation),
        _slice(cf.wind_speed), _slice(cf.shortwave_downward), _slice(cf.longwave_downward),
        _slice(cf.vapor_pressure), _slice(cf.black_carbon_snow), _slice(cf.black_carbon_ice),
        _slice(cf.cloud_optical_thickness), _slice(cf.solar_zenith_angle),
        _slice(cf.shortwave_downward_diffuse), _slice(cf.cloud_fraction),
        cf.time_step, cf.temperature_air_mean, cf.wind_speed_mean, cf.precipitation_mean,
        cf.temperature_observation_height, cf.wind_observation_height;
        snow_drift=_slice(cf.snow_drift),
        accumulation_mean=cf.accumulation_mean,
        temperature_air_effective=cf.temperature_air_effective,
        dataset=cf.dataset, latitude=cf.latitude, longitude=cf.longitude,
        elevation=cf.elevation, elevation_native=cf.elevation_native,
        elevation_offset=cf.elevation_offset,
    )
end

# Stamp the selected block with climatology provenance. `climatology_n_years` is the number of
# years the *cycle contains*, which for this path is the block length rather than the number of
# years averaged; `climatology_representative_year` records the years chosen (a single `Int` for
# a one-year block, a `Vector{Int}` otherwise), and its presence is what distinguishes this path
# from `:average` after the fact. The requested `window` is preserved when given, so the
# provenance still says what was asked for and not only what was selected from it.
function _stamp_representative(cf::ClimateForcing, years, steps_per_year::Integer, window;
                               rank_by::Symbol, accumulation_error, melt_ratio)
    times = collect(lookup(dims(cf.temperature_air, Ti)))
    ys = Int.(collect(years))
    meta = merge(DD.metadata(cf), (
        climatology_window_start = window === nothing ? first(times) : window[1],
        climatology_window_stop = window === nothing ? last(times) : window[2],
        climatology_n_years = length(ys),
        climatology_steps_per_year = Int(steps_per_year),
        climatology_representative_year = length(ys) == 1 ? ys[1] : ys,
        # Diagnostics recorded rather than only logged. `verbose` controls *printing*; these are
        # always measured, so a batch caller running with `verbose=false` can still act on them —
        # in particular on `climatology_melt_ratio`, which is the "this cycle will equilibrate to
        # bare ice" signal that would otherwise exist only as warning text.
        # `climatology_rank_by` is the ranking actually used, which differs from the requested one
        # when the accumulation fallback fired.
        climatology_rank_by = rank_by,
        climatology_accumulation_error = Float64(accumulation_error),
        climatology_melt_ratio = Float64(melt_ratio),
    ))
    return rebuild(cf; metadata=meta)
end

# Report the selection, and warn when the block chosen on melt is atypical in a variable it was
# *not* chosen on. Selecting on one statistic cannot make a block typical in every term, so the
# honest thing is to say how far off the others are rather than to imply it is representative of
# everything.
#
# The accumulation and temperature checks compare the block's mean against the record's, as a
# fraction and as a difference respectively, rather than as percentile ranks: a block of `n`
# consecutive years has only `n_years - n + 1` candidates, so a rank over them is too coarse to
# mean anything once `n` is a meaningful fraction of the record.
"""
    _report_representative_block(sel, record_summary, cycle_summary, accumulation_error,
                                 statistic, rank_by)

Report which block was selected, how closely its score matched the target, and how far its
accumulation is from the record's.

Takes the two summaries already computed by [`_representative_block`](@ref) rather than deriving
its own. An earlier version built one `ClimateSummary` **per candidate year** — 20 extra
surface-energy-balance sweeps on a 20-year record — to produce two numbers, and computed a
second, differently-defined accumulation offset (mean-of-year-means rather than the block cycle's
own) that could disagree with the guard's in the same log line.

No temperature diagnostic is reported. `_years_slice` carries `cf.temperature_air_mean` through
unchanged and `initialize_climate_summary` reads that scalar rather than recomputing it, so a
per-block temperature offset is **identically zero by construction** — not a physical result. The
earlier version reported it and warned above 0.5 K, a branch that could never fire.
"""
# The two summaries are `ClimateSummary`s, left unannotated for the include-order reason above.
function _report_representative_block(sel, record_summary, cycle_summary, accumulation_error,
                                      statistic::Symbol, rank_by::Symbol)
    @info "Representative spinup block selected" years=sel.block_years n_years=length(sel.block_years) statistic=statistic rank_by=rank_by n_candidate_years=length(sel.scores) block_score=sel.block_score score_target=sel.target accumulation_error=accumulation_error record_melt=record_summary.melt cycle_melt=cycle_summary.melt

    # 10%: large enough that a spinup would equilibrate a visibly different column, loose enough
    # not to fire on the ordinary sampling spread of a few years drawn from a few decades. This
    # fires only when the accumulation fallback did not already fix it — either because it was
    # disabled, or because SMB ranking did not improve matters at this site.
    if isfinite(accumulation_error) && accumulation_error > 0.1
        @warn "Representative block $(sel.block_years) was chosen on `$(rank_by)` but its mean " *
              "accumulation is $(round(100 * accumulation_error, digits=1))% off the record's. " *
              "The spinup will equilibrate that term to an atypical value; consider a longer " *
              "block, a different window, or `rank_by=:smb`." maxlog=2
    end
end
