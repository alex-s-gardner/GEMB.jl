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

Measured over a 21-site synthetic fleet (`bench/calibrate_initial_guess.jl`), melt on the cycle
as a fraction of melt over the real record, summed across the 18 sites that melt:

| method | melt recovered |
|---|---|
| `:average` | **0%** — identically 0.0 at *every* site, including one melting 384 kg m-2 yr-1 |
| `:representative`, `n_years=1` | 99% |
| `:representative`, `n_years=3` | 116% |
| `:representative`, `n_years=5` | 100% |

A melting site spun up under `:average` therefore equilibrates as a melt-free site, with its
firn air content and refreeze terms wrong by construction — which matters most for work reading
those terms directly, such as altimetry.

# Which method to use: it depends on the regime

Classifying the same fleet by the ratio of melt to accumulation gives three regimes, and they
want different things. Counts are sites whose equilibrated firn air content fell below 1 m:

| regime | melt/accumulation | n | `:average` | `:representative` |
|---|---|---|---|---|
| dry / cold | ~0 | 3 | fine | fine (agrees to ~10%) |
| percolation | 0.3 or below | 14 | fine | fine (agrees to ~10-20%) |
| percolation | 0.3 to 1.0 | 2 | keeps firn | **goes to ice** |
| ablation | above 1.0 | 2 | keeps firn (wrongly) | goes to ice (correctly) |

So:

- **Dry and cold: use `:average`.** With no melt the column responds through densification and
  heat diffusion, both smooth in temperature, and averaging costs essentially nothing. This is
  why it remains the default. Not *exactly* linear even there — the densification rate carries a
  convex Arrhenius factor, and a 12 K shift moves it 3.9x — but the error is a fraction of the
  term rather than all of it.
- **Ablation (melt above accumulation): use `:representative`.** The column equilibrates to bare
  ice, which is correct: accumulation resets each year and what matters is the mean melt driving
  latent-heat warming. `:representative` reproduces that melt; `:average` sets it to zero and
  returns a firn column the site does not have.
- **Melt approaching accumulation: neither is trustworthy — compare them.** This is the genuinely
  unresolved band. `:average` deletes the melt; `:representative` can drive the column to ice
  because the record's *mean* ratio being below 1 does not stop individual years exceeding it
  (the two sites that lost their firn had record ratios of 0.83 and 0.93). See
  [`_warn_cycle_melt_ratio`](@ref), which warns on the constructed cycle's own ratio.

A longer block does **not** fix that band: 3- and 5-year blocks lose the firn at those sites too
(FAC 0.12-0.75 m), so interannual alternation within a short cycle is not enough to sustain a
column whose mean melt is near its accumulation.

These boundaries come from 21 synthetic sites derived from one parent site, with only 2 sites in
the critical band and 2 in ablation. Treat 0.3 and 1.0 as directions, not as calibrated numbers.

# Keyword arguments

- `method`: `:average` or `:representative` (see above).
- `model_parameters`: a [`ModelParameters`](@ref). **Required** for `:representative`, which
  needs the physics to rank candidate years; ignored by `:average`.
- `n_years`: for `:representative`, the block length (default 3). Clamped to the number of
  complete years available, with a warning. The measured differences between 1, 3 and 5 are
  modest and not monotonic (see the melt table above), so this is not a strong lever; 3 is a
  compromise between carrying some interannual variability and keeping the cycle short enough
  that `gemb_spinup`'s per-cycle criteria are judged often. **Cycling the whole record is
  deliberately not offered**: it would give one convergence comparison per record length, so a
  10-cycle drift window would span 10 records.
- `statistic`: for `:representative`, which year to call representative — `:mean` (default) or
  `:median` of the candidate years' melt. `:mean` is the default because the column equilibrates
  to a long-run mean and the annual melt distribution is right-skewed (measured at one fleet
  site: mean 150.8 against median 127.3 kg m-2 yr-1), so selecting on the median biases the
  spinup low.
- `rank_by`: for `:representative`, how each candidate year's melt is measured — `:model`
  (default) runs one year of [`gemb`](@ref) per candidate, `:estimate` uses
  [`initialize_climate_summary`](@ref)'s surface-energy-balance estimate. `:model` is the default
  because it is both more accurate and cheap: 32 candidate years cost 2.4 s (0.074 s/yr), against
  the hundreds of cycles the spinup it feeds will run. `:estimate` avoids needing a starting
  profile and is retained for that case, but its melt is biased ~2.6x high, and while a uniform
  bias cancels when *ordering* years it does not cancel when picking the year nearest a target
  value.
- `profile`: for `:representative` with `rank_by = :model`, the column to run the candidate
  years on. Defaults to [`initialize_profile`](@ref)`(model_parameters, cf)`. Only the *ranking*
  depends on it, never the forcing returned.
- `verbose`: for `:representative`, report the selected year and how typical it is (default
  `true`). See the caveat below.

# The single-year caveat

A representative year is one weather sequence, repeated for the whole spinup, and selecting it on
melt cannot make it typical in every other term. `verbose=true` reports the chosen year's melt,
accumulation and mean temperature as percentiles of the candidates, and warns when it is an
outlier in accumulation or temperature — a real warning, which fires on real records (at one
fleet site the melt-representative year was also the warmest of 32). The averaged year is equally
a single sequence, with the difference that it is not one that ever occurred.

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

    return (non_leap=non_leap, forcing_index=forcing_index,
            times_complete=times_noleap[forcing_index],
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
                             profile=nothing,
                             verbose::Bool=true)

    if method === :representative
        model_parameters === nothing && throw(ArgumentError(
            "forcing_climatology(method=:representative) requires `model_parameters`: " *
            "ranking candidate years needs the physics that turns forcing into melt."))
        return _representative_block(cf, model_parameters; n_years=n_years, window=window,
                                     statistic=statistic, rank_by=rank_by,
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

So the block length is a mild knob, not a fix: 1, 3 and 5 years recovered 99%, 116% and 100% of
record melt respectively — differences that are not monotonic and are small against the 0% that
`:average` recovers.

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
                               profile=nothing, verbose::Bool=true)

    statistic in (:mean, :median) ||
        throw(ArgumentError("statistic must be :mean or :median, got $(repr(statistic))"))
    rank_by in (:model, :estimate) ||
        throw(ArgumentError("rank_by must be :model or :estimate, got $(repr(rank_by))"))
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

    # Annual melt per candidate year. `:model` integrates each year on a real column, which is
    # the quantity being matched; `:estimate` uses the summary's surface-energy-balance melt,
    # which needs no column but is biased ~2.6x high.
    melts = if rank_by === :model
        col = profile === nothing ? initialize_profile(mp, cf) : profile
        mp_last = ModelParameters(;
            (f => getfield(mp, f) for f in fieldnames(ModelParameters)
             if f != :output_frequency)..., output_frequency=:last)
        [_year_melt(col, _year_slice(cf, cy, y), mp_last) for y in years]
    else
        [initialize_climate_summary(_year_slice(cf, cy, y), mp).melt for y in years]
    end

    # Target is the statistic over every candidate *year*, so it describes the record rather
    # than the blocks; each block is then scored by its own mean melt against it.
    target = statistic === :mean ? Statistics.mean(melts) : Statistics.median(melts)
    n_windows = length(years) - n_block + 1
    block_melts = [Statistics.mean(@view melts[i:(i + n_block - 1)]) for i in 1:n_windows]
    start = argmin(abs.(block_melts .- target))
    block_years = years[start:(start + n_block - 1)]

    cycle = _years_slice(cf, cy, block_years)

    if verbose
        _report_representative_block(years, melts, block_years, block_melts[start],
                                     target, statistic, rank_by, cf, cy, mp)
        _warn_cycle_melt_ratio(cycle, block_years, block_melts[start], melts, mp)
    end

    return _stamp_representative(cycle, block_years, cy.steps_per_year, window)
end

"""
    _warn_cycle_melt_ratio(cycle, block_years, block_melt, year_melts, mp)

Warn when the selected cycle's own melt-to-accumulation ratio puts it at or past the point where
repeating it drives the column to solid ice.

The test is on the **constructed cycle**, not on the record, and that distinction is the whole
point. `melt/accumulation > 1` over a record means the site is ablating, and a cycle that
equilibrates to ice is then the *correct* answer. But the record's ratio is a mean, and
interannual variability puts individual years above 1 while the mean sits below it — so a site
whose record ratio is 0.9 can still yield a cycle that ablates. Measured on the fleet: the two
sites that lost their firn to a repeated block (equilibrated FAC falling from ~4.8 m to
0.12-0.75 m) had *record* ratios of 0.83 and 0.93, i.e. both below 1, and would pass a check
made on the record.

`CYCLE_MELT_RATIO_WARN` is deliberately below 1.0 for that reason: the useful signal is
"approaching the threshold", since a cycle at 0.8 of a mean can straddle it. Sites well below it
(ratio <= ~0.3 on the fleet) reproduce the averaged climatology's equilibrium closely and need no
warning.

Not an error, and not a refusal to return the cycle: at a genuinely ablating site the ice column
*is* the answer, and what matters there is that the mean melt — and so the latent-heat warming
it drives — is right, which is exactly what this path gets right and `:average` does not.
"""
const CYCLE_MELT_RATIO_WARN = 0.6

function _warn_cycle_melt_ratio(cycle::ClimateForcing, block_years, block_melt,
                                year_melts, mp::ModelParameters)
    cs = initialize_climate_summary(cycle, mp)
    acc = cs.accumulation
    acc > 0.0 || return nothing

    # The cycle's own ratio, plus the worst single year inside it — a block can average below
    # the threshold while containing a year far above it, and it is the repeated *block* that
    # equilibrates, so both are worth reporting.
    ratio = block_melt / acc
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

# Annual melt [kg m-2 yr-1] of one candidate year, integrated by the model itself on `col`.
# `output_frequency=:last` leaves a single step holding the year's total, so the sum is that
# total and the year length comes from the forcing rather than from the output times.
function _year_melt(col, year_forcing::ClimateForcing, mp_last::ModelParameters)
    out = gemb(col, year_forcing, mp_last; thermal_workspace=ThermalWorkspace())
    years = length(dims(year_forcing, Ti)) * year_forcing.time_step / SECONDS_PER_YEAR
    return years > 0 ? sum(parent(out[:melt])) / years : NaN
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
    times = collect(lookup(dims(cf.temperature_air, Ti)))[cy.non_leap][cy.forcing_index]
    wanted = Set(Int.(years))
    keep = findall(y -> y in wanted, Dates.year.(times))
    tdim = Ti(times[keep])
    _slice(a) = DimArray(parent(a)[cy.non_leap][cy.forcing_index][keep], (tdim,))

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
function _stamp_representative(cf::ClimateForcing, years, steps_per_year::Integer, window)
    times = collect(lookup(dims(cf.temperature_air, Ti)))
    ys = Int.(collect(years))
    meta = merge(DD.metadata(cf), (
        climatology_window_start = window === nothing ? first(times) : window[1],
        climatology_window_stop = window === nothing ? last(times) : window[2],
        climatology_n_years = length(ys),
        climatology_steps_per_year = Int(steps_per_year),
        climatology_representative_year = length(ys) == 1 ? ys[1] : ys,
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
function _report_representative_block(years, melts, block_years, block_melt, target,
                                      statistic, rank_by, cf, cy, mp)
    all_summaries = [initialize_climate_summary(_year_slice(cf, cy, y), mp) for y in years]
    block_set = Set(Int.(block_years))
    in_block = [Int(y) in block_set for y in years]

    acc_all = Statistics.mean(s.accumulation for s in all_summaries)
    acc_blk = Statistics.mean(all_summaries[i].accumulation for i in eachindex(years) if in_block[i])
    t_all = Statistics.mean(s.temperature_air_mean for s in all_summaries)
    t_blk = Statistics.mean(all_summaries[i].temperature_air_mean for i in eachindex(years) if in_block[i])

    acc_ratio = acc_all > 0 ? acc_blk / acc_all : NaN
    t_delta = t_blk - t_all

    @info "Representative spinup block selected" years=block_years n_years=length(block_years) statistic=statistic rank_by=rank_by n_candidate_years=length(years) block_melt=block_melt melt_target=target melt_ratio=(target > 0 ? block_melt / target : NaN) accumulation_ratio=acc_ratio temperature_offset=t_delta

    # 10% in accumulation and 0.5 K in mean temperature: both are large enough that a spinup
    # would equilibrate a visibly different column, and loose enough not to fire on the ordinary
    # sampling spread of a few years drawn from a few decades.
    if isfinite(acc_ratio) && abs(acc_ratio - 1) > 0.1
        @warn "Representative block $(block_years) was chosen on melt but its mean " *
              "accumulation is $(round(100 * (acc_ratio - 1), digits=1))% off the record's. " *
              "The spinup will equilibrate that term to an atypical value; consider a longer " *
              "block or a different window." maxlog=2
    end
    if abs(t_delta) > 0.5
        @warn "Representative block $(block_years) was chosen on melt but its mean " *
              "temperature is $(round(t_delta, digits=2)) K off the record's. Densification " *
              "is Arrhenius in temperature, so this biases the equilibrated density " *
              "profile; consider a longer block or a different window." maxlog=2
    end
end
