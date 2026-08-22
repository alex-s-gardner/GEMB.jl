# Tests for initialize_profile
# DimArray and Ti come from GEMB (re-exported from DimensionalData)
using GEMB: DimArray, Ti

@testset "Default parameters" begin
    mp = GEMB.ModelParameters()

    # Create minimal ClimateForcing for initialization
    times = [DateTime(2000, 1, 1), DateTime(2000, 1, 1, 3)]
    cf = initialize_forcing(times,
        [253.15, 253.15], [80000.0, 80000.0], [0.001, 0.001],
        [5.0, 5.0], [100.0, 100.0], [250.0, 250.0], [300.0, 300.0];
        temperature_air_mean=253.15, wind_speed_mean=5.0,
        precipitation_mean=200.0, temperature_observation_height=2.0,
        wind_observation_height=10.0)

    # Both escape-hatch flags → the MATLAB pure-ice initialization, which also
    # keeps the configured column limits for the grid validation below.
    profile = GEMB.initialize_profile(mp, cf; constant_density=true, constant_temperature=true)

    # Check that profile has expected fields
    @test haskey(profile, :temperature)
    @test haskey(profile, :dz)
    @test haskey(profile, :density)
    @test haskey(profile, :water)
    @test haskey(profile, :grain_radius)

    # Check values
    dz = collect(profile[:dz])
    temperature = collect(profile[:temperature])
    density = collect(profile[:density])

    # All temperatures should be the mean
    @test all(temperature .≈ 253.15)

    # All densities should be density_ice
    @test all(density .≈ mp.density_ice)

    # Top layers should have constant spacing
    n_top = round(Int, mp.column_ztop / mp.column_dztop)
    @test all(dz[1:n_top] .≈ mp.column_dztop)

    # Total depth should be at or slightly above column_depth_max
    # (last cell extends past boundary in MATLAB implementation)
    @test sum(dz) >= mp.column_depth_max
    @test sum(dz) < mp.column_depth_max + dz[end] + 1e-10
end

@testset "Grid stretching" begin
    mp = GEMB.ModelParameters(column_ztop=5.0, column_dztop=0.05, column_depth_max=50.0, column_zy=1.10)

    times = [DateTime(2000, 1, 1), DateTime(2000, 1, 1, 3)]
    cf = initialize_forcing(times,
        [260.0, 260.0], [80000.0, 80000.0], [0.001, 0.001],
        [5.0, 5.0], [100.0, 100.0], [250.0, 250.0], [300.0, 300.0];
        temperature_air_mean=260.0, wind_speed_mean=5.0,
        precipitation_mean=200.0, temperature_observation_height=2.0,
        wind_observation_height=10.0)

    profile = GEMB.initialize_profile(mp, cf)
    dz = collect(profile[:dz])

    n_top = round(Int, mp.column_ztop / mp.column_dztop)

    # Below top zone, layers should increase geometrically
    for i in (n_top+2):length(dz)
        ratio = dz[i] / dz[i-1]
        @test ratio ≈ mp.column_zy atol = 1e-10
    end
end

@testset "z_center calculation" begin
    mp = GEMB.ModelParameters()

    times = [DateTime(2000, 1, 1), DateTime(2000, 1, 1, 3)]
    cf = initialize_forcing(times,
        [253.15, 253.15], [80000.0, 80000.0], [0.001, 0.001],
        [5.0, 5.0], [100.0, 100.0], [250.0, 250.0], [300.0, 300.0];
        temperature_air_mean=253.15, wind_speed_mean=5.0,
        precipitation_mean=200.0, temperature_observation_height=2.0,
        wind_observation_height=10.0)

    profile = GEMB.initialize_profile(mp, cf)
    # z_center is no longer stored on the profile; recompute it from dz.
    dz = collect(profile[:dz])
    z_center = GEMB.dz2z(dz)

    # First center should be at -dz[1]/2
    @test z_center[1] ≈ -dz[1] / 2 atol = 1e-12

    # All centers should be negative (below surface)
    @test all(z_center .< 0)

    # Centers should be monotonically decreasing
    @test all(diff(z_center) .< 0)
end

# ---------------------------------------------------------------------------
# The universal, threshold-free steady-state initialization.
#
# `initialize_profile` derives every profile field from the climate forcing via
# `initialize_climate_summary` -> `steady_state_profile`, with the sign of the net
# annual mass balance as the only regime discriminator. These tests cover the
# claims that motivate that design; the escape-hatch (MATLAB-fidelity) path is
# covered by the pure-ice testsets above and by `test_ablation_regime.jl`.
# ---------------------------------------------------------------------------

# Synthetic single-site forcing with a real seasonal cycle, so the annual harmonic
# fit has something to fit. Defaults describe a cold, dry, accumulating site
# (dry-snow zone); the kwargs walk it toward the melt-dominated end.
#
# `vapor_fraction` is the vapor pressure as a fraction of saturation over ice at
# `T_mean`. It has to scale with temperature: a fixed value that is reasonable at
# 250 K is desert-dry at 271 K, and the resulting latent-heat loss would hold the
# skin below freezing at any radiation.
function _dry_snow_forcing(; n=365*2, T_mean=250.0, T_amp=10.0, precip=0.8,
    shortwave=80.0, longwave=190.0, vapor_fraction=0.7)
    time = DateTime(2000, 1, 1) .+ Day.(0:n-1)
    doy = collect(1:n)
    temp = T_mean .+ T_amp .* cos.(2π .* doy ./ 365.25 .- π)

    # Saturation vapor pressure over ice [Pa], Clausius-Clapeyron about the triple
    # point — enough for a plausible test forcing.
    e_sat = 611.0 * exp((GEMB.LS / 461.5) * (1.0 / GEMB.CtoK - 1.0 / T_mean))

    return initialize_forcing(
        time, temp, fill(85000.0, n), fill(precip, n), fill(4.0, n),
        fill(shortwave, n), fill(longwave, n), fill(vapor_fraction * e_sat, n);
        temperature_air_mean=T_mean, wind_speed_mean=4.0,
        precipitation_mean=precip * 365.25,
        temperature_observation_height=2.0, wind_observation_height=10.0)
end

@testset "steady_state_density: monotonic and method-dependent" begin
    # Fine near the surface so `d[1]` is genuinely the surface value rather than a
    # cell whose center has already densified.
    zc = vcat(-collect(0.025:0.05:1.0), -collect(1.25:0.5:200.0))
    ρ0, ρi = 350.0, 910.0

    profiles = Dict{Symbol,Vector{Float64}}()
    for method in (:HerronLangway, :Arthern, :Ligtenberg)
        mp = GEMB.initialize_parameters(densification_method=method)
        d = GEMB.steady_state_density(zc, 253.0, 300.0, ρ0, ρi, mp)
        profiles[method] = d

        @test isapprox(d[1], ρ0; atol=5.0)          # surface is fresh snow
        @test all(diff(d) .>= -1e-9)                # monotonic with depth
        @test all(ρ0 - 1e-9 .<= d .<= ρi + 1e-9)    # bounded by snow and ice
        @test d[end] > 0.98 * ρi                    # asymptotes to ice
    end

    # Each law is genuinely plumbed through rather than all resolving to H–L.
    @test !isapprox(profiles[:Arthern], profiles[:HerronLangway]; rtol=1e-3)
    @test !isapprox(profiles[:Ligtenberg], profiles[:HerronLangway]; rtol=1e-3)
end

@testset "Steady-state init: firn column from climate" begin
    cf = _dry_snow_forcing()
    mp = GEMB.ModelParameters()
    profile = GEMB.initialize_profile(mp, cf)

    cs = GEMB.initialize_climate_summary(cf, mp)
    @test cs.balance > 0.0                          # accumulating site

    density = collect(profile[:density])
    temperature = collect(profile[:temperature])
    water = collect(profile[:water])
    re = collect(profile[:grain_radius])
    gdn = collect(profile[:grain_dendricity])
    gsp = collect(profile[:grain_sphericity])

    # A graded firn column, not a block of ice.
    @test !all(density .== mp.density_ice)
    @test density[1] < 400.0
    @test all(diff(density) .>= -1e-9)
    @test density[end] > 0.98 * mp.density_ice

    # Surface density is the fresh-snow density for the configured method.
    ρ0 = GEMB.fresh_snow_density(mp, cs.temperature_air_mean,
        cs.accumulation_effective, cs.wind_speed_mean)
    @test isapprox(density[1], ρ0; rtol=0.05)

    # Grain size starts near fresh snow (0.05 mm) rather than the 2.5 mm ice cap,
    # and coarsens with depth.
    @test re[1] < 0.5
    @test re[end] > re[1]
    @test all(GEMB.RE_NEW_SNOW - 1e-9 .<= re .<= GEMB.GRAIN_RADIUS_ICE + 1e-9)
    @test all(0.0 - 1e-9 .<= gdn .<= 1.0 + 1e-9)
    @test all(0.0 - 1e-9 .<= gsp .<= 1.0 + 1e-9)

    # Nothing above the melt point; a cold column holds no water.
    @test all(temperature .<= GEMB.CtoK + 1e-9)
    @test all(water .>= 0.0)
    @test all(water[temperature .< GEMB.CtoK - GEMB.T_TOLERANCE] .== 0.0)

    # A cold site needs the depth it was given: the derived column is not truncated
    # above the depth where firn reaches ice.
    @test sum(profile[:dz]) > 20.0
end

@testset "Damped annual thermal wave" begin
    cf = _dry_snow_forcing(T_mean=250.0, T_amp=12.0)
    mp = GEMB.ModelParameters()
    cs = GEMB.initialize_climate_summary(cf, mp)

    # The harmonic fit recovers the amplitude that was imposed.
    @test isapprox(cs.temperature_amplitude, 12.0; rtol=0.05)

    profile = GEMB.initialize_profile(mp, cf)
    T = collect(profile[:temperature])
    depth = -GEMB.dz2z(collect(profile[:dz]))
    # The deep asymptote is the mean *surface* temperature, not the mean air temperature: the
    # column is coupled to the atmosphere only through the surface energy balance. Latent
    # warming does not appear because it decays over the annual accumulation layer, which is
    # metres from the base of a column this deep. See `_steady_state_temperature`.
    T_deep = cs.temperature_surface_mean

    # The wave decays with depth: the deviation from the deep mean shrinks, and the
    # deepest cell sits at the mean.
    dev = abs.(T .- T_deep)
    @test dev[1] > dev[end]
    @test isapprox(T[end], min(T_deep, GEMB.CtoK); atol=0.1)

    # Decay is set by the damping depth, so it is essentially complete well before
    # the bottom of a column this deep.
    i10 = findfirst(>=(10.0), depth)
    @test dev[i10] < 0.1 * max(dev[1], 1e-12)
end

@testset "Irreducible water matches the runtime formula" begin
    # A *temperate accumulating* site, which is the only regime that initializes with pore
    # water — and a narrow target. Two conditions must hold at once:
    #
    #   * the surface energy balance must reach the melt point, since `_irreducible_water`
    #     gates on the marched temperature (not the air temperature), and
    #   * `balance` must stay positive, or `steady_state_profile` takes its ablation branch and
    #     returns a dry block of ice.
    #
    # Radiation is what warms the surface, but it also drives the melt that pushes `balance`
    # negative, so the two conditions pull against each other: at 2 kg m-2 d-1 of precipitation
    # there is no radiation that satisfies both (150/280 leaves the surface 2.9 K short of
    # melting; 200/310 melts 3466 kg m-2 yr-1 against 731 of snowfall and the site ablates).
    # The heavy precipitation here is what buys the margin — 7300 kg m-2 yr-1 keeps `balance` at
    # +3804 while the surface sits at the melt point, giving water in 64 of 248 cells.
    cf = _dry_snow_forcing(T_mean=272.0, T_amp=3.0, precip=20.0,
                           shortwave=200.0, longwave=310.0)
    mp = GEMB.ModelParameters()
    profile = GEMB.initialize_profile(mp, cf)

    # Guard the premise: if a change ever makes this site ablate, the assertions below would
    # pass vacuously on a dry ice column rather than testing the water formula.
    @test GEMB.initialize_climate_summary(cf, mp).balance > 0.0

    dz = collect(profile[:dz])
    density = collect(profile[:density])
    temperature = collect(profile[:temperature])
    water = collect(profile[:water])

    @test any(water .> 0.0)

    # Where water is present it equals `calculate_melt`'s irreducible content,
    # (ρi − ρ)·S·(M/ρ) with M = ρ·dz, so init and runtime agree by construction.
    # The saturation comes from `irreducible_saturation`, not from
    # `mp.water_irreducible_saturation` directly, so this follows whichever
    # `water_irreducible_method` is configured.
    for i in eachindex(water)
        water[i] == 0.0 && continue
        M = density[i] * dz[i]
        expected = (mp.density_ice - density[i]) *
            GEMB.irreducible_saturation(mp, density[i]) * (M / density[i])
        @test isapprox(water[i], expected; rtol=1e-12)
        @test temperature[i] >= GEMB.CtoK - GEMB.T_TOLERANCE
    end
end

@testset "SEB melt responds to radiation" begin
    mp = GEMB.ModelParameters()

    # Cold polar site: the surface never reaches melting.
    cf_cold = _dry_snow_forcing(T_mean=240.0, T_amp=8.0, shortwave=40.0, longwave=150.0)
    @test GEMB.initialize_climate_summary(cf_cold, mp).melt == 0.0

    # A near-melting site with weak incoming longwave *also* gets no melt, even at
    # high shortwave: at snow albedo, εσT⁴ at the melt point (≈306 W m-2) exceeds
    # what the surface receives. A degree-day scheme, blind to radiation and
    # albedo, would predict substantial melt at these air temperatures — this is
    # the case that motivates the SEB estimate.
    cf_dim = _dry_snow_forcing(T_mean=271.0, T_amp=4.0, shortwave=300.0, longwave=280.0)
    @test GEMB.initialize_climate_summary(cf_dim, mp).melt == 0.0

    # With realistic incoming longwave the same site melts, increasingly with
    # shortwave.
    melts = [GEMB.initialize_climate_summary(
                 _dry_snow_forcing(T_mean=271.0, T_amp=4.0, shortwave=sw, longwave=320.0),
                 mp).melt
             for sw in (150.0, 250.0, 350.0)]
    @test all(melts .> 0.0)
    @test issorted(melts)
    @test melts[end] > 2 * melts[1]
end

@testset "Continuity across the balance sign change" begin
    # Walk the net annual balance through zero by varying snowfall at fixed melt
    # (melt is radiation-driven here, so it stays put while accumulation sweeps).
    # The old threshold-based scheme jumped from a 250 m firn column to a 25 m ice
    # block across this crossing; the continuous scheme must not.
    mp = GEMB.ModelParameters()
    make(p) = _dry_snow_forcing(T_mean=270.5, T_amp=4.0, precip=p,
                                shortwave=100.0, longwave=320.0)
    precips = collect(range(0.8, 2.4; length=17))

    summaries = [GEMB.initialize_climate_summary(make(p), mp) for p in precips]
    balances = [cs.balance for cs in summaries]

    # The sweep genuinely crosses the regime boundary, and does so monotonically.
    @test balances[1] < 0.0
    @test balances[end] > 0.0
    @test issorted(balances)

    # Melt is *not* constant across the sweep: more snowfall means more snow cover,
    # a higher albedo and so less melt (here 395 -> 134 kg m-2 yr-1). That coupling
    # is the albedo/melt feedback `initialize_climate_summary` iterates to close, so
    # assert it runs the right way and stays smooth rather than assuming melt is an
    # independent knob.
    melts = [cs.melt for cs in summaries]
    @test issorted(melts; rev=true)
    @test melts[1] > 2 * melts[end]
    # Smooth: no single 0.1 mm/day step in snowfall moves melt by more than a
    # quarter of the total range it covers.
    @test maximum(abs.(diff(melts))) < 0.25 * (melts[1] - melts[end])

    results = [GEMB.initialize_profile(mp, make(p)) for p in precips]
    surface_density = [collect(prof[:density])[1] for prof in results]
    # The derived depth is read off the grid it produced.
    depths = [sum(prof[:dz]) for prof in results]

    i_cross = findlast(<=(0.0), balances)     # last ablating case
    @test i_cross !== nothing && i_cross < length(precips)

    # Below the crossing the column is ice, and it stays ice all the way down.
    @test all(surface_density[1:i_cross] .== mp.density_ice)

    # Above the crossing the surface keeps lightening as accumulation grows, rather
    # than being pinned at one of two branch values.
    @test issorted(surface_density[i_cross+1:end]; rev=true)

    # Derived depth grows with accumulation rather than switching between two
    # constants, and never exceeds the configured ceiling.
    @test depths[1] < depths[end]
    @test issorted(depths)
    @test all(depths .<= mp.column_depth_max)

    # --- the continuity claim itself -------------------------------------------
    #
    # Surface density approaches ice *from above* as `b -> 0+`, and this is the
    # whole point: the age of the surface cell's mid-depth is `z·ρ/b`, which
    # diverges as burial stops, so the parcel there has had unbounded time to
    # densify. The ablation column (exactly ice) is therefore the continuous limit
    # of the firn column, not a separate case — which is what lets the scheme drop
    # the old threshold.
    #
    # A fixed sweep always shows *some* jump at the crossing simply because `b`
    # is sampled discretely, so the meaningful test is that the jump shrinks
    # without bound as the sampling refines. A genuine discontinuity would hold its
    # size instead.
    # Sample the first accumulating case at successively finer resolution and record
    # both the balance reached there and how far its surface sits below ice.
    function crossing_jump(n)
        ps = collect(range(1.1, 1.4; length=n))
        bs = [GEMB.initialize_climate_summary(make(p), mp).balance for p in ps]
        i = findlast(<=(0.0), bs)
        (i === nothing || i == length(ps)) && return nothing
        ρ_above = collect(GEMB.initialize_profile(mp, make(ps[i+1]))[:density])[1]
        return (balance=bs[i+1], gap=mp.density_ice - ρ_above)
    end

    steps = filter(!isnothing, [crossing_jump(n) for n in (17, 33, 65, 129, 257)])
    @test length(steps) == 5

    # Refining drives the sampled balance toward zero...
    @test issorted([s.balance for s in steps]; rev=true)
    @test steps[end].balance < 0.1 * steps[1].balance
    # ...and the surface density closes onto ice as it does, with no floor. Measured:
    # b = 9.8 -> 4.8 -> 2.6 -> 1.3 -> 0.61 kg m-2 yr-1 gives a gap below ice of
    # 520 -> 479 -> 408 -> 300 -> 162 kg m-3. Every refinement shrinks it, so the
    # ablation column (exactly ice) is the `b -> 0+` limit of the firn column rather
    # than a separate branch — the claim that lets the scheme drop the threshold.
    gaps = [s.gap for s in steps]
    @test issorted(gaps; rev=true)
    @test gaps[end] < 0.4 * gaps[1]
end

@testset "No initialized cell starts above the melt point" begin
    # An initialized column carries no enthalpy above `CtoK` — that energy would be
    # melt, and there is no water in the profile to hold it — so every path must
    # return T <= 273.15 K. The climate-derived march clamps in
    # `_steady_state_temperature`; the two fidelity flags fill `temperature_air_mean`
    # verbatim, so they are clamped in `initialize_profile` itself. A temperate site
    # (mean annual air temperature above freezing) exercises all of it.
    mp = GEMB.ModelParameters()

    for T_mean in (274.0, 280.0, 295.0)
        cf = _dry_snow_forcing(T_mean=T_mean, T_amp=6.0, precip=2.0,
                               shortwave=120.0, longwave=300.0)
        @test cf.temperature_air_mean > GEMB.CtoK   # the forcing really is warm

        # Climate-derived path: no warning, clamped by the march.
        prof = GEMB.initialize_profile(mp, cf)
        @test maximum(collect(prof[:temperature])) <= GEMB.CtoK

        # Each fidelity path warns as it clamps, and the clamp holds.
        for kw in ((; constant_temperature=true),
                   (; constant_density=true, constant_temperature=true))
            prof_f = @test_logs (:warn, r"above the melt point") match_mode = :any (
                GEMB.initialize_profile(mp, cf; kw...))
            T = collect(prof_f[:temperature])
            @test all(T .<= GEMB.CtoK)
            @test all(T .≈ GEMB.CtoK)   # clamped from the (warmer) mean, not below it
        end

        # `constant_density` alone leaves temperature to the march, so it needs no
        # clamp of its own and must not emit the warning.
        prof_d = GEMB.initialize_profile(mp, cf; constant_density=true)
        @test maximum(collect(prof_d[:temperature])) <= GEMB.CtoK
    end

    # A cold site is untouched: the clamp is a ceiling, not a rewrite.
    cf_cold = _dry_snow_forcing(T_mean=250.0)
    prof_cold = GEMB.initialize_profile(mp, cf_cold;
        constant_density=true, constant_temperature=true)
    @test all(collect(prof_cold[:temperature]) .≈ 250.0)
end

@testset "Climate summary cost is negligible" begin
    # The scheme's premise is that a best guess is rapidly computable from the
    # climate data alone. The meaningful comparison is against the integration it
    # seeds, not an absolute wall-clock bound: a spinup is many passes over this
    # same forcing, so anything well under one pass is free in context. Measured as
    # a ratio so the test does not depend on the machine.
    ds = GEMB_ClimateForcing.simulate_climate_forcing("test_1", 3)
    cf = initialize_forcing(ds)
    mp = GEMB.initialize_parameters()

    profile = GEMB.initialize_profile(mp, cf)
    gemb(profile, cf, mp; verbose=false)             # warm up
    t_pass = @elapsed gemb(profile, cf, mp; verbose=false)

    GEMB.initialize_profile(mp, cf)                  # warm up
    t_init = @elapsed GEMB.initialize_profile(mp, cf)

    # The whole initialization costs a fraction of a single forcing pass. (Measured
    # ~0.27x: the cost is dominated by the per-timestep SEB melt estimate, which is
    # a secant solve over the series.)
    @test t_init < t_pass

    # The march itself is negligible even against the summary: it is a few thousand
    # steps over a column, not a sweep of the series.
    dz = GEMB.initialize_grid(mp)
    cs = GEMB.initialize_climate_summary(cf, mp)
    GEMB.steady_state_profile(dz, cs, mp)            # warm up
    t_march = @elapsed GEMB.steady_state_profile(dz, cs, mp)
    @test t_march < 0.1 * t_pass
end


@testset "initialize_age" begin
    cf = _dry_snow_forcing()
    mp_ss = GEMB.ModelParameters(initialize_age=:steady_state)
    mp_zero = GEMB.ModelParameters(initialize_age=:zero)

    @test GEMB.ModelParameters().initialize_age === :steady_state

    age_ss = collect(GEMB.initialize_profile(mp_ss, cf)[:age])
    age_zero = collect(GEMB.initialize_profile(mp_zero, cf)[:age])

    # :zero is the old behaviour on every cell.
    @test all(age_zero .== 0.0)

    # The marched age is the residence time of a parcel buried at cs.balance: increasing
    # downward and finite everywhere. The surface cell is *not* age 0 — the profile is sampled
    # at cell centers, so cell 1 has already been buried by half its own thickness, which is
    # the right analogue of the transient run's mass-weighted cell age.
    @test all(isfinite, age_ss)
    @test issorted(age_ss)
    @test age_ss[end] > age_ss[1] > 0.0

    # Cell 1's age is exactly the burial time of half its own mass, which pins the units
    # (days, not years) as well as the cell-center convention.
    let cs = GEMB.initialize_climate_summary(cf, mp_ss),
        p = GEMB.initialize_profile(mp_ss, cf)

        half_mass = 0.5 * collect(p[:dz])[1] * collect(p[:density])[1]
        @test isapprox(age_ss[1], half_mass / cs.balance * (GEMB.SECONDS_PER_YEAR / 86400.0);
            rtol=0.05)
    end

    # Consistency with the burial rate rather than a hardcoded number: the age at the base is
    # the column's mass above it divided by the annual balance, to march resolution.
    profile = GEMB.initialize_profile(mp_ss, cf)
    cs = GEMB.initialize_climate_summary(cf, mp_ss)
    dz = collect(profile[:dz])
    density = collect(profile[:density])
    # Mass above the *center* of the deepest cell, matching the sampling convention.
    mass_above = sum(dz[1:end-1] .* density[1:end-1]) + 0.5 * dz[end] * density[end]
    expected_years = mass_above / cs.balance
    @test isapprox(age_ss[end] / (GEMB.SECONDS_PER_YEAR / 86400.0), expected_years; rtol=0.1)

    # The fidelity flags replace the marched column, so its age must not be inherited.
    @test all(collect(GEMB.initialize_profile(mp_ss, cf; constant_density=true)[:age]) .== 0.0)
    @test all(collect(GEMB.initialize_profile(mp_ss, cf;
        constant_density=true, constant_temperature=true)[:age]) .== 0.0)

    # An ablation column buries nothing, so the march never advances and there is no
    # residence time to report on either setting. Same forcing as the ablation end of the
    # balance sweep below.
    cf_ablation = _dry_snow_forcing(T_mean=270.5, T_amp=4.0, precip=0.8,
        shortwave=100.0, longwave=320.0)
    @test GEMB.initialize_climate_summary(cf_ablation, mp_ss).balance < 0.0
    @test all(collect(GEMB.initialize_profile(mp_ss, cf_ablation)[:age]) .== 0.0)
end

@testset "close_off_age" begin
    # NaN for an open column, the age of the shallowest closed-off cell otherwise. Same
    # NaN-means-absent convention as ice_slab_depth.
    @test isnan(GEMB.close_off_age([300.0, 500.0, 700.0], [1.0, 2.0, 3.0]))
    @test GEMB.close_off_age([300.0, 850.0, 900.0], [1.0, 2.0, 3.0]) == 2.0
    # Exactly at the threshold counts as closed off.
    @test GEMB.close_off_age([GEMB.DENSITY_PORE_CLOSEOFF], [7.0]) == 7.0
    @test isnan(GEMB.close_off_age(Float64[], Float64[]))
    # The threshold is a parameter, so a study can move it.
    @test GEMB.close_off_age([700.0, 800.0], [1.0, 2.0], 750.0) == 2.0

    # On a real initialized firn column it lies strictly between the surface and basal ages.
    cf = _dry_snow_forcing()
    profile = GEMB.initialize_profile(GEMB.ModelParameters(initialize_age=:steady_state), cf)
    density = collect(profile[:density])
    age = collect(profile[:age])
    coa = GEMB.close_off_age(density, age)
    @test maximum(density) >= GEMB.DENSITY_PORE_CLOSEOFF   # this column does close off
    @test 0.0 < coa < age[end]
end
