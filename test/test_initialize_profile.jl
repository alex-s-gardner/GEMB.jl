# Tests for initialize_profile - matches MATLAB test_model_initialize_profile.m
# DimArray and Ti come from GEMB (re-exported from DimensionalData)
using GEMB: DimArray, Ti

@testset "Default parameters" begin
    mp = GEMB.ModelParameters()

    # Create minimal ClimateForcing for initialization
    times = [DateTime(2000, 1, 1), DateTime(2000, 1, 1, 3)]
    cf = GEMB.initialize_forcing(times,
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
    cf = GEMB.initialize_forcing(times,
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
    cf = GEMB.initialize_forcing(times,
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

    return GEMB.initialize_forcing(
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
    T_deep = cs.temperature_air_mean + cs.latent_warming

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
    # Warm site: the surface reaches the melt point, so water is nonzero there.
    cf = _dry_snow_forcing(T_mean=272.0, T_amp=3.0, precip=2.0,
                           shortwave=150.0, longwave=280.0)
    mp = GEMB.ModelParameters()
    profile = GEMB.initialize_profile(mp, cf)

    dz = collect(profile[:dz])
    density = collect(profile[:density])
    temperature = collect(profile[:temperature])
    water = collect(profile[:water])

    @test any(water .> 0.0)

    # Where water is present it equals `calculate_melt`'s irreducible content,
    # (ρi − ρ)·S·(M/ρ) with M = ρ·dz, so init and runtime agree by construction.
    for i in eachindex(water)
        water[i] == 0.0 && continue
        M = density[i] * dz[i]
        expected = (mp.density_ice - density[i]) *
            mp.water_irreducible_saturation * (M / density[i])
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

@testset "Climate summary cost is negligible" begin
    # The scheme's premise is that a best guess is rapidly computable from the
    # climate data alone. The meaningful comparison is against the integration it
    # seeds, not an absolute wall-clock bound: a spinup is many passes over this
    # same forcing, so anything well under one pass is free in context. Measured as
    # a ratio so the test does not depend on the machine.
    ds = GEMB_ClimateForcing.simulate_climate_forcing("test_1", 3)
    cf = GEMB.initialize_forcing(ds)
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

# MATLAB validation test
matlab_validation_testset("initialize_profile", "initialize_profile.mat") do ref
    # This will validate the grid geometry matches MATLAB
    # Note: Full profile validation requires matching all parameters exactly
    # Here we just validate grid structure
    
    params = GEMB.initialize_parameters()
    ds = GEMB_ClimateForcing.simulate_climate_forcing("test_1", 3)
    forcing = GEMB.initialize_forcing(ds)
    # Both escape-hatch flags keep the configured column limits, which is what the
    # MATLAB reference grid was built from — the climate-derived depth is a
    # deliberate departure and would not be a like-for-like geometry comparison.
    profile = GEMB.initialize_profile(params, forcing;
        constant_density=true, constant_temperature=true)

    # Validate number of layers
    n_layers_julia = length(profile.dz)
    n_layers_matlab = Int(ref["n_layers_init"][1])
    
    @test n_layers_julia == n_layers_matlab
    
    # Validate grid structure (dz and z_center patterns)
    # Note: Exact values may differ slightly due to forcing differences
    # but structure should match
    @test length(profile.dz) > 0
    @test all(profile.dz .> 0)  # All layers have positive thickness
end
