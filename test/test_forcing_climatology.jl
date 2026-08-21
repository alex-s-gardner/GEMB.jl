# Tests for `forcing_climatology`, both methods.
#
# The `:average` method is covered incidentally by several other files (provenance in
# `test_gemb_spinup.jl`, metadata carry-over in `test_calculate_density.jl`, the drift layer in
# `test_blowing_snow.jl`). What is tested here is the contract shared by both methods and the
# `:representative` path in full.
#
# The motivating measurement, from `bench/calibrate_initial_guess.jl` over a 21-site fleet: melt
# integrated over an *averaged* one-year cycle is identically 0.0 at every site, including sites
# melting up to 384 kg m-2 yr-1 over their real record. Melt is rectified — zero until the skin
# reaches the melt point — so averaging years slot-by-slot cancels the warm excursions that carry
# it. `:representative` exists for that, and the tests below pin the properties that make it a
# real alternative: the forcing it returns is unmodified real data, and its selection is on melt.

# `Statistics` via GEMB rather than as a direct dependency, matching
# `test_synthetic_regression.jl` — `runtests.jl` does not import it into `Main`.
using GEMB: DimensionalData, Statistics

# Forcing whose years differ, so averaging is distinguishable from selecting and the melt
# ranking has something to rank. `warm_years` get a temperature offset (and the radiation to
# make the surface actually melt), so melt is concentrated in them.
function _multiyear_forcing(; years=1990:1994, warm_years=(1992,), warm_offset=14.0,
                            base_temperature=258.0, precip=1.5)
    days = Date(first(years), 1, 1):Day(1):Date(last(years), 12, 31)
    days = filter(d -> dayofyear(d) != 366, days)      # keep every year 365 steps
    time = DateTime.(days)

    seasonal(d) = 8.0 * cos(2π * dayofyear(d) / 365.0 - π)
    warm(d) = year(d) in warm_years ? warm_offset : 0.0

    temperature = [base_temperature + seasonal(d) + warm(d) for d in days]
    # Radiation follows the warmth, as it must for the surface energy balance to reach melting
    # in the warm years without being pushed there in the cold ones.
    shortwave = [90.0 + 6.0 * warm(d) for d in days]
    longwave = [200.0 + 4.0 * warm(d) for d in days]

    return initialize_forcing(
        time, temperature, fill(80000.0, length(time)), fill(precip, length(time)),
        fill(3.0, length(time)), shortwave, longwave, fill(90.0, length(time));
        temperature_air_mean=Statistics.mean(temperature), wind_speed_mean=3.0,
        precipitation_mean=precip * 365.25,
        temperature_observation_height=2.0, wind_observation_height=10.0)
end

@testset "forcing_climatology: shared contract" begin
    cf = _multiyear_forcing()
    mp = GEMB.ModelParameters()

    # `:average` is the default, and remains the historical behaviour.
    @test all(collect(forcing_climatology(cf)[:temperature_air]) .==
              collect(forcing_climatology(cf; method=:average)[:temperature_air]))

    # An unknown method is rejected rather than silently averaging.
    @test_throws ArgumentError forcing_climatology(cf; method=:nonsense)

    # `:representative` needs the physics to rank years, so omitting it is an error rather than
    # a silent fallback to averaging.
    @test_throws ArgumentError forcing_climatology(cf; method=:representative)

    # Both methods return a forcing usable as a spinup cycle: same time step, same scalar
    # metadata, and a whole number of years.
    for cycle in (forcing_climatology(cf),
                  forcing_climatology(cf; method=:representative, model_parameters=mp,
                                      n_years=2, verbose=false))
        @test cycle.time_step == cf.time_step
        @test cycle.temperature_air_mean == cf.temperature_air_mean
        @test cycle.precipitation_mean == cf.precipitation_mean
        @test length(DimensionalData.dims(cycle, Ti)) % 365 == 0
    end
end

@testset "forcing_climatology: :representative returns real, unmodified forcing" begin
    cf = _multiyear_forcing()
    mp = GEMB.ModelParameters()

    cycle = forcing_climatology(cf; method=:representative, model_parameters=mp,
                                n_years=1, verbose=false)
    md = DimensionalData.metadata(cycle)
    chosen = md[:climatology_representative_year]

    # The point of the method: every value is the value that actually occurred in that year, not
    # an average of anything. This is what preserves the variance and cross-field covariance
    # that melt depends on, so it is asserted exactly rather than approximately.
    times = collect(DimensionalData.lookup(DimensionalData.dims(cf, Ti)))
    src = findall(t -> year(t) == chosen, times)
    for layer in (:temperature_air, :precipitation, :shortwave_downward, :longwave_downward,
                  :wind_speed, :vapor_pressure, :pressure_air)
        @test collect(cycle[layer]) == collect(cf[layer])[src]
    end

    # ...and the averaged cycle is *not* the raw data, which is what makes the above a real
    # distinction rather than a tautology on this forcing.
    @test collect(forcing_climatology(cf)[:temperature_air]) != collect(cf[:temperature_air])[src]
end

@testset "forcing_climatology: :representative selects on melt" begin
    # One year is much warmer than the rest, so the melt-ranked choice is knowable a priori: the
    # mean-melt target sits below the warm year and above the cold ones, and with four cold
    # years and one warm one the closest single year to the mean is a cold one.
    cf = _multiyear_forcing(years=1990:1994, warm_years=(1992,))
    mp = GEMB.ModelParameters()

    one_year = forcing_climatology(cf; method=:representative, model_parameters=mp,
                                   n_years=1, verbose=false)
    @test DimensionalData.metadata(one_year)[:climatology_representative_year] != 1992

    # The two ranking paths are independent implementations of "how much does this year melt",
    # so they are not required to agree on the year — only to both return a year that exists.
    for rank_by in (:model, :estimate)
        y = DimensionalData.metadata(
            forcing_climatology(cf; method=:representative, model_parameters=mp,
                                n_years=1, rank_by=rank_by, verbose=false)
        )[:climatology_representative_year]
        @test y in 1990:1994
    end

    # Selecting on the median rather than the mean is a different question and may pick a
    # different year; both must still be real years.
    for statistic in (:mean, :median)
        y = DimensionalData.metadata(
            forcing_climatology(cf; method=:representative, model_parameters=mp,
                                n_years=1, statistic=statistic, verbose=false)
        )[:climatology_representative_year]
        @test y in 1990:1994
    end

    @test_throws ArgumentError forcing_climatology(cf; method=:representative,
        model_parameters=mp, statistic=:mode)
    @test_throws ArgumentError forcing_climatology(cf; method=:representative,
        model_parameters=mp, rank_by=:vibes)
    @test_throws ArgumentError forcing_climatology(cf; method=:representative,
        model_parameters=mp, n_years=0)
end

@testset "forcing_climatology: block length and consecutiveness" begin
    cf = _multiyear_forcing(years=1990:1994)
    mp = GEMB.ModelParameters()

    for n in 1:3
        cycle = forcing_climatology(cf; method=:representative, model_parameters=mp,
                                    n_years=n, verbose=false)
        md = DimensionalData.metadata(cycle)

        # The cycle is exactly `n` years long, in both its length and its provenance.
        @test length(DimensionalData.dims(cycle, Ti)) == n * 365
        @test md[:climatology_n_years] == n

        ys = md[:climatology_representative_year]
        if n == 1
            @test ys isa Int                       # scalar for a one-year block
        else
            @test length(ys) == n
            # Consecutive, which is the property that preserves the real transition from one
            # year into the next — the multi-year memory a block exists to carry.
            @test ys == collect(first(ys):(first(ys) + n - 1))
        end

        # The selected years are the years present in the cycle, so provenance is not merely
        # decorative.
        present = sort(unique(year.(collect(DimensionalData.lookup(
            DimensionalData.dims(cycle, Ti))))))
        # `ys` is a scalar Int for a one-year block and a Vector otherwise, so normalize before
        # comparing rather than assuming either shape.
        @test present == sort(Int[ys...])
    end

    # Asking for more years than exist clamps, with a warning, rather than erroring or
    # silently returning a shorter cycle without saying so.
    cycle = @test_logs (:warn, r"complete years available") match_mode = :any (
        forcing_climatology(cf; method=:representative, model_parameters=mp,
                            n_years=99, verbose=false))
    @test DimensionalData.metadata(cycle)[:climatology_n_years] == 5
end

@testset "forcing_climatology: provenance distinguishes the methods" begin
    cf = _multiyear_forcing(years=1990:1994)
    mp = GEMB.ModelParameters()

    avg = DimensionalData.metadata(forcing_climatology(cf))
    rep = DimensionalData.metadata(forcing_climatology(cf; method=:representative,
                                                      model_parameters=mp, n_years=2,
                                                      verbose=false))

    # `climatology_representative_year` is the discriminator: present for the selected block,
    # absent for the average, so a spun-up profile records which cycle produced it.
    @test !haskey(avg, :climatology_representative_year)
    @test haskey(rep, :climatology_representative_year)

    # `climatology_n_years` means "years in the cycle" for `:representative` and "years
    # averaged" for `:average`, which is why the discriminator above is needed to read it.
    @test avg[:climatology_n_years] == 5
    @test rep[:climatology_n_years] == 2

    # A requested window is preserved verbatim by both, so the provenance says what was asked
    # for and not only what was selected from it.
    window = (DateTime(1990, 1, 1), DateTime(1993, 12, 31))
    for m in (DimensionalData.metadata(forcing_climatology(cf, window)),
              DimensionalData.metadata(forcing_climatology(cf, window;
                  method=:representative, model_parameters=mp, n_years=1, verbose=false)))
        @test m[:climatology_window_start] == window[1]
        @test m[:climatology_window_stop] == window[2]
    end
end

@testset "forcing_climatology: :representative preserves melt that :average destroys" begin
    # The measurement this whole method exists for, in miniature. A site that melts in one year
    # of five: the averaged cycle dilutes that year's warmth across all five and can lose the
    # melt entirely, while a cycle containing the warm year retains it.
    cf = _multiyear_forcing(years=1990:1994, warm_years=(1992,), warm_offset=16.0)
    mp = GEMB.ModelParameters(output_frequency=:last)
    profile = initialize_profile(mp, cf)

    melt_rate(f) = begin
        out = gemb(profile, f, mp)
        years = length(DimensionalData.dims(f, Ti)) * f.time_step / GEMB.SECONDS_PER_YEAR
        sum(parent(out[:melt])) / years
    end

    melt_record = melt_rate(cf)
    @test melt_record > 0.0                        # the premise: this site melts

    # The warm year on its own melts; the average of all five melts strictly less, since
    # averaging is what removes the excursion.
    warm_cycle = forcing_climatology(cf, (DateTime(1992, 1, 1), DateTime(1992, 12, 31));
                                     method=:representative, model_parameters=mp,
                                     n_years=1, verbose=false)
    @test melt_rate(warm_cycle) > melt_rate(forcing_climatology(cf))
end

@testset "forcing_climatology: cycle melt/accumulation guard" begin
    # The guard reports on the *constructed cycle*, not the record, because the record's mean
    # ratio being below 1 does not stop individual years exceeding it: on the fleet, the two
    # sites whose firn was lost to a repeated cycle had record ratios of 0.83 and 0.93.
    mp = GEMB.ModelParameters()

    # A site that melts far more than it accumulates: the ablation regime, where the warning
    # says the ice column is expected rather than a failure.
    #
    # Built directly rather than via `_multiyear_forcing`, because the regime is narrow and this
    # fixture has to sit inside it: it needs `melt > accumulation > 0`, and simply making the
    # site hotter overshoots into `accumulation == 0` (every step falls as rain), which would
    # skip the ratio test entirely rather than exercise it. Measured here: melt 1661,
    # accumulation 293, ratio 5.7.
    days = filter(d -> dayofyear(d) != 366, Date(1990, 1, 1):Day(1):Date(1993, 12, 31))
    hot_time = DateTime.(days)
    hot_temperature = [265.0 + 10.0 * cos(2π * dayofyear(d) / 365.0 - π) for d in days]
    hot = initialize_forcing(
        hot_time, hot_temperature, fill(80000.0, length(days)), fill(1.0, length(days)),
        fill(3.0, length(days)), fill(240.0, length(days)), fill(320.0, length(days)),
        fill(90.0, length(days));
        temperature_air_mean=Statistics.mean(hot_temperature), wind_speed_mean=3.0,
        precipitation_mean=365.25,
        temperature_observation_height=2.0, wind_observation_height=10.0)

    # Guard the fixture's premise, so a change that moves it out of the ablation regime fails
    # here rather than quietly making the warning assertion below unreachable.
    cs = GEMB.initialize_climate_summary(hot, mp)
    @test cs.melt > cs.accumulation > 0.0

    @test_logs (:warn, r"melts more than it accumulates") match_mode = :any (
        forcing_climatology(hot; method=:representative, model_parameters=mp,
                            n_years=1, verbose=true))

    # A cold, dry site is far below the threshold, so no ratio warning is emitted.
    cold = _multiyear_forcing(years=1990:1993, warm_years=(), base_temperature=245.0)
    @test GEMB.initialize_climate_summary(cold, mp).melt == 0.0
    cycle = forcing_climatology(cold; method=:representative, model_parameters=mp,
                                n_years=1, verbose=true)
    @test DimensionalData.metadata(cycle)[:climatology_representative_year] in 1990:1993

    # `verbose=false` silences the reporting entirely, including the guard.
    @test_logs (forcing_climatology(cold; method=:representative, model_parameters=mp,
                                    n_years=1, verbose=false))
end
