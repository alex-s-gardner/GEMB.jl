# Tests for calculate_grain_size

# Helper to create a ClimateForcingStep for grain size tests
function _make_grain_cfs(; dt=86400.0)
    return GEMB.ClimateForcingStep(
        dt,               # dt [s]
        265.0,            # temperature_air
        100000.0,         # pressure_air
        0.0,              # precipitation
        5.0,              # wind_speed
        200.0,            # shortwave_downward
        300.0,            # longwave_downward
        400.0,            # vapor_pressure
        260.0,            # temperature_air_mean
        5.0,              # wind_speed_mean
        200.0,            # precipitation_mean
        2.0,              # temperature_observation_height
        10.0,             # wind_observation_height
        0.0,              # black_carbon_snow
        0.0,              # black_carbon_ice
        0.0,              # cloud_optical_thickness
        0.0,              # solar_zenith_angle
        0.0,              # shortwave_downward_diffuse
        0.1,              # cloud_fraction
    )
end

@testset "grain size evolves regardless of configuration" begin
    # Metamorphism is not contingent on model configuration. MATLAB skips
    # `calculate_grain_size` unless `albedo_method` is :GardnerSharp or :BrunLefebre, which
    # covers only two of the four schemes that read `grain_radius` — :ArthernB/:Crocus/
    # :CrocusPure densification (which goes as 1/r²) and the grain-radius emissivity
    # methods read it too, and select independently of `albedo_method`. Those configurations
    # silently ran with grain size frozen at the initial profile.
    n = 5
    temperature = 260.0 * ones(n)
    dz = 0.1 * ones(n)
    density = 300.0 * ones(n)
    water = zeros(n)
    grain_radius = 0.5 * ones(n)
    grain_dendricity = 0.5 * ones(n)
    grain_sphericity = 0.5 * ones(n)
    cfs = _make_grain_cfs(dt=86400.0 * 30)

    evolves(mp) = GEMB.calculate_grain_size(temperature, dz, density, water,
        copy(grain_radius), copy(grain_dendricity), copy(grain_sphericity),
        cfs, mp)[1] != grain_radius

    # Every albedo scheme, including the ones MATLAB skips on.
    for am in (:GardnerSharp, :BrunLefebre, :GreuellKonzelmann, :None)
        @test evolves(GEMB.ModelParameters(albedo_method=am))
    end

    # The two consumers that select independently of `albedo_method`, paired with an albedo
    # scheme that reads no grain size — the configurations the old gate broke.
    for dm in (:ArthernB, :Crocus, :CrocusPure)
        @test evolves(GEMB.ModelParameters(albedo_method=:None, densification_method=dm))
    end
    for em in (:grain_radius_threshold, :grain_radius_w_threshold)
        @test evolves(GEMB.ModelParameters(albedo_method=:None, emissivity_method=em))
    end

    # ...and a configuration where nothing reads grain size at all still evolves it, since
    # the work is now unconditional rather than gated on who happens to consume the result.
    @test evolves(GEMB.ModelParameters(albedo_method=:None, densification_method=:Arthern,
        emissivity_method=:uniform))
end

@testset "all three state arrays are updated in place" begin
    # `grain_radius` used to be returned as a fresh array while dendricity and sphericity
    # were mutated. It now shares their convention, which removes the largest allocation in
    # the function. Pinned here because callers rebinding the result would not notice the
    # difference, so nothing else would catch a regression to copy-on-return.
    mp = GEMB.ModelParameters(albedo_method=:GardnerSharp)
    cfs = _make_grain_cfs(dt=86400.0 * 10)
    n = 4
    temperature = [250.0, 252.0, 254.0, 256.0]
    dz = fill(0.1, n)
    density = fill(300.0, n)
    water = zeros(n)
    grain_radius = fill(0.5, n)
    grain_dendricity = fill(0.5, n)
    grain_sphericity = fill(0.5, n)

    (gs_out, gdn_out, gsp_out) = GEMB.calculate_grain_size(
        temperature, dz, density, water, grain_radius,
        grain_dendricity, grain_sphericity, cfs, mp)

    @test gs_out === grain_radius
    @test gdn_out === grain_dendricity
    @test gsp_out === grain_sphericity
    # ...and the aliased array really did change, so `===` is not passing on a no-op.
    @test !all(grain_radius .== 0.5)
end

@testset "Dendritic dry low gradient" begin
    mp = GEMB.ModelParameters(albedo_method=:GardnerSharp)
    cfs = _make_grain_cfs(dt=86400.0)

    temperature = [260.0, 260.0, 260.0]
    dz = [0.1, 0.1, 0.1]
    density = [200.0, 200.0, 200.0]
    water = [0.0, 0.0, 0.0]
    grain_radius = [0.2, 0.2, 0.2]
    grain_dendricity = [0.8, 0.8, 0.8]
    grain_sphericity = [0.2, 0.2, 0.2]

    gdn_before = copy(grain_dendricity)
    gsp_before = copy(grain_sphericity)
    gr_before = copy(grain_radius)

    (gs_out, gdn_out, gsp_out) = GEMB.calculate_grain_size(
        temperature, dz, density, water, grain_radius,
        grain_dendricity, grain_sphericity, cfs, mp)

    # Dendricity should decrease (decay)
    @test all(gdn_out .< gdn_before)
    # Sphericity should increase
    @test all(gsp_out .> gsp_before)
    # Grain radius should change
    @test gs_out != gr_before
end

@testset "Dendritic dry high gradient" begin
    mp = GEMB.ModelParameters(albedo_method=:GardnerSharp)
    cfs = _make_grain_cfs(dt=86400.0)

    # High temperature gradient (> 5 K/m)
    temperature = [260.0, 270.0, 280.0]
    dz = [0.1, 0.1, 0.1]
    density = [200.0, 200.0, 200.0]
    water = [0.0, 0.0, 0.0]
    grain_radius = [0.2, 0.2, 0.2]
    grain_dendricity = [0.8, 0.8, 0.8]
    grain_sphericity = [0.2, 0.2, 0.2]

    gdn_before = copy(grain_dendricity)
    gsp_before = copy(grain_sphericity)

    (_, gdn_out, gsp_out) = GEMB.calculate_grain_size(
        temperature, dz, density, water, grain_radius,
        grain_dendricity, grain_sphericity, cfs, mp)

    # Under high gradient: dendricity decreases, sphericity decreases
    @test all(gdn_out .< gdn_before)
    @test all(gsp_out .< gsp_before)
end

@testset "Dendritic wet snow" begin
    mp = GEMB.ModelParameters(albedo_method=:GardnerSharp)
    cfs = _make_grain_cfs(dt=86400.0)

    temperature = [GEMB.CtoK, GEMB.CtoK, GEMB.CtoK]
    dz = [0.1, 0.1, 0.1]
    density = [250.0, 250.0, 250.0]
    water = [1.0, 1.0, 1.0]
    grain_radius = [0.2, 0.2, 0.2]
    grain_dendricity = [0.8, 0.8, 0.8]
    grain_sphericity = [0.2, 0.2, 0.2]

    gdn_before = copy(grain_dendricity)
    gsp_before = copy(grain_sphericity)

    (_, gdn_out, gsp_out) = GEMB.calculate_grain_size(
        temperature, dz, density, water, grain_radius,
        grain_dendricity, grain_sphericity, cfs, mp)

    # Wet snow causes rapid rounding
    @test all(gdn_out .< gdn_before)
    @test all(gsp_out .> gsp_before)
end

@testset "Nondendritic dry (Marbouty)" begin
    mp = GEMB.ModelParameters(albedo_method=:GardnerSharp)
    cfs = _make_grain_cfs(dt=86400.0)

    # Moderate gradient (~20 K/m)
    temperature = [250.0, 252.0, 254.0]
    dz = [0.1, 0.1, 0.1]
    density = [300.0, 300.0, 300.0]  # must be < 400 for growth
    water = [0.0, 0.0, 0.0]
    grain_radius = [0.5, 0.5, 0.5]
    grain_dendricity = [0.0, 0.0, 0.0]
    grain_sphericity = [0.5, 0.5, 0.5]
    gr_before = copy(grain_radius)      # mutated in place

    (gs_out, gdn_out, _) = GEMB.calculate_grain_size(
        temperature, dz, density, water, grain_radius,
        grain_dendricity, grain_sphericity, cfs, mp)

    # Expect growth in grain size
    @test all(gs_out .> gr_before)
    # Dendricity stays at 0
    @test all(gdn_out .== 0.0)
end

@testset "Marbouty density limit" begin
    # `:Marbouty` explicitly — this ceiling is that method's behaviour, and the
    # `:Arthern` default has no density ceiling.
    mp = GEMB.ModelParameters(albedo_method=:GardnerSharp, grain_growth_method=:Marbouty)
    cfs = _make_grain_cfs(dt=86400.0)

    temperature = [250.0, 252.0, 254.0]
    dz = [0.1, 0.1, 0.1]
    density = [450.0, 450.0, 450.0]  # > 400 threshold
    water = [0.0, 0.0, 0.0]
    grain_radius = [0.5, 0.5, 0.5]
    grain_dendricity = [0.0, 0.0, 0.0]
    grain_sphericity = [0.5, 0.5, 0.5]
    gr_before = copy(grain_radius)      # mutated in place

    (gs_out, _, _) = GEMB.calculate_grain_size(
        temperature, dz, density, water, grain_radius,
        grain_dendricity, grain_sphericity, cfs, mp)

    # No growth expected above density threshold
    @test gs_out ≈ gr_before atol = 1e-10
end

@testset "Nondendritic wet snow (Brun)" begin
    mp = GEMB.ModelParameters(albedo_method=:GardnerSharp)
    cfs = _make_grain_cfs(dt=86400.0)

    temperature = [GEMB.CtoK, GEMB.CtoK, GEMB.CtoK]
    dz = [0.1, 0.1, 0.1]
    density = [300.0, 300.0, 300.0]
    water = [1.5, 1.5, 1.5]
    grain_radius = [0.5, 0.5, 0.5]
    grain_dendricity = [0.0, 0.0, 0.0]
    grain_sphericity = [0.5, 0.5, 0.5]
    gr_before = copy(grain_radius)      # mutated in place

    (gs_out, gdn_out, _) = GEMB.calculate_grain_size(
        temperature, dz, density, water, grain_radius,
        grain_dendricity, grain_sphericity, cfs, mp)

    # Expect growth via wet snow mechanism
    @test all(gs_out .> gr_before)
    # Dendricity stays at 0
    @test all(gdn_out .== 0.0)
end

@testset "Clamping limits" begin
    mp = GEMB.ModelParameters(albedo_method=:GardnerSharp)
    # Large dt to force dendricity toward 0
    cfs = _make_grain_cfs(dt=86400.0 * 100)

    temperature = [260.0, 260.0, 260.0]
    dz = [0.1, 0.1, 0.1]
    density = [200.0, 200.0, 200.0]
    water = [0.0, 0.0, 0.0]
    grain_radius = [0.2, 0.2, 0.2]
    grain_dendricity = [0.1, 0.1, 0.1]  # close to 0
    grain_sphericity = [0.2, 0.2, 0.2]

    (_, gdn_out, gsp_out) = GEMB.calculate_grain_size(
        temperature, dz, density, water, grain_radius,
        grain_dendricity, grain_sphericity, cfs, mp)

    @test all(gdn_out .>= 0.0)
    @test all(gsp_out .<= 1.0)
end

@testset "Grain size cap" begin
    mp = GEMB.ModelParameters(albedo_method=:GardnerSharp)
    cfs = _make_grain_cfs(dt=86400.0 * 50)

    temperature = [GEMB.CtoK, GEMB.CtoK, GEMB.CtoK]
    dz = [0.1, 0.1, 0.1]
    density = [300.0, 300.0, 300.0]
    water = [2.0, 2.0, 2.0]
    grain_radius = [1.9, 1.9, 1.9]  # radius 1.9mm
    grain_dendricity = [0.0, 0.0, 0.0]
    grain_sphericity = [1.0, 1.0, 1.0]

    (gs_out, _, _) = GEMB.calculate_grain_size(
        temperature, dz, density, water, grain_radius,
        grain_dendricity, grain_sphericity, cfs, mp)

    # Radius should be capped at 1.0mm for spherical grains
    @test all(gs_out .<= 1.0 + 1e-10)
end

@testset "grain_growth_method" begin
    # The dry non-dendritic branch is the only one this option touches. A column spanning the
    # Marbouty density ceiling (400) pins the handover: cells 1-2 are below it, cells 3-5 above.
    T = [263.0, 260.0, 255.0, 250.0, 248.0]
    dz = fill(0.05, 5)
    density = [200.0, 350.0, 450.0, 600.0, 800.0]
    water = zeros(5)
    gdn = zeros(5)   # non-dendritic
    gsp = zeros(5)
    cfs = _make_grain_cfs(dt=10800.0)

    _grow(method) = GEMB.calculate_grain_size(copy(T), copy(dz), copy(density), copy(water),
        fill(0.5, 5), copy(gdn), copy(gsp), cfs, GEMB.ModelParameters(grain_growth_method=method))[1]

    r_marbouty = _grow(:Marbouty)
    r_arthern = _grow(:Arthern)
    r_hybrid = _grow(:hybrid)

    # `:Arthern` is the default (it is what the Community Firn Model ships); `:Marbouty`
    # stops dead at DENSITY_MARBOUTY_MAX, so it shows no growth in cells 3-5.
    @test GEMB.ModelParameters().grain_growth_method === :Arthern
    @test all(r_marbouty[1:2] .> 0.5)
    @test r_marbouty[3:5] == fill(0.5, 3)

    # Arthern grows at every density, and monotonically faster the warmer the cell.
    @test all(r_arthern .> 0.5)
    @test issorted(r_arthern, rev=true)   # T decreases with depth, so growth does too

    # :hybrid is exactly Marbouty below the ceiling and exactly Arthern at or above it.
    @test r_hybrid[1:2] == r_marbouty[1:2]
    @test r_hybrid[3:5] == r_arthern[3:5]

    # The Arthern law integrates dr²/dt = kgr·exp(-Eg/RT) exactly over the step. Checked in SI
    # against the module's mm working units, which is where a unit slip would show.
    dt = 365 * 86400.0
    r0 = 0.5e-3   # m
    dr2 = GEMB.GRAIN_GROWTH_KGR * exp(-GEMB.GRAIN_GROWTH_EG / (GEMB.R_GAS * 250.0)) * dt
    expected = sqrt(r0^2 + dr2) * 1e3   # mm
    r_yr = GEMB.calculate_grain_size(fill(250.0, 1), [0.05], [600.0], [0.0], [0.5],
        [0.0], [0.0], _make_grain_cfs(dt=dt), GEMB.ModelParameters(grain_growth_method=:Arthern))[1]
    @test r_yr[1] ≈ expected rtol = 1e-12

    # Wet cells take the Brun branch on every setting, so the option cannot touch them.
    wet = fill(2.0, 5)
    r_wet = [GEMB.calculate_grain_size(copy(T), copy(dz), copy(density), copy(wet),
        fill(0.5, 5), copy(gdn), fill(1.0, 5), cfs,
        GEMB.ModelParameters(grain_growth_method=m))[1] for m in (:Marbouty, :Arthern, :hybrid)]
    @test r_wet[1] == r_wet[2] == r_wet[3]
end
