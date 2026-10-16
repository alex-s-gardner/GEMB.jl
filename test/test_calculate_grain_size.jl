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
    # `:Marbouty` explicitly — this ceiling is that method's behavior; the `:hybrid` default
    # hands over to Arthern growth above it.
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

    # `:hybrid` is the default; `:Marbouty` stops dead at DENSITY_MARBOUTY_MAX, so it shows no
    # growth in cells 3-5.
    @test GEMB.ModelParameters().grain_growth_method === :hybrid
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

@testset "Grain growth closed forms" begin
    dt = 86400.0 * 2
    dt_d = 2.0
    cfs = _make_grain_cfs(dt=dt)
    mp = GEMB.ModelParameters()
    ex(T) = exp(-6e3 / T)
    d_of(gd, gs) = 0.1 * (gd / 0.99 + (1 - gd / 0.99) * (gs / 0.99 * 3 + (1 - gs / 0.99) * 4))
    dz = [0.1, 0.1, 0.1]

    # Dendritic, dry, gradient below 5 K/m.
    T = fill(260.0, 3)
    r, gd, gs = GEMB.calculate_grain_size(T, dz, fill(200.0, 3), zeros(3), fill(0.2, 3),
        fill(0.8, 3), fill(0.2, 3), cfs, mp)
    @test gd[2] ≈ 0.8 - 2e8 * ex(260.0) * dt_d rtol = 1e-12
    @test gs[2] ≈ 0.2 + 1e9 * ex(260.0) * dt_d rtol = 1e-12
    @test 2 * r[2] ≈ d_of(gd[2], gs[2]) rtol = 1e-12

    # Dendritic, dry, gradient above 5 K/m: both fall by 2e8·e^(−6000/T)·|∇T|^0.4·dt.
    T = [258.0, 260.0, 262.0]        # 20 K/m: 4 K across the 0.2 m between cells 1 and 3
    r, gd, gs = GEMB.calculate_grain_size(T, dz, fill(200.0, 3), zeros(3), fill(0.2, 3),
        fill(0.8, 3), fill(0.5, 3), cfs, mp)
    Δ = 2e8 * ex(260.0) * 20.0^0.4 * dt_d
    @test gd[2] ≈ 0.8 - Δ rtol = 1e-12
    @test gs[2] ≈ 0.5 - Δ rtol = 1e-12

    # Dendritic, wet: θ = 100·W/(ρ·dz) %, capped at 9; dendricity falls by θ³/16·dt.
    W = 0.2
    r, gd, gs = GEMB.calculate_grain_size(fill(GEMB.CtoK, 3), dz, fill(250.0, 3), fill(W, 3),
        fill(0.2, 3), fill(0.8, 3), fill(0.2, 3), cfs, mp)
    θ = 100 * W / (250.0 * 0.1)
    @test gd[2] ≈ 0.8 - θ^3 / 16 * dt_d rtol = 1e-12
    @test gs[2] ≈ 0.2 + θ^3 / 16 * dt_d rtol = 1e-12

    # Brun (1989) wet growth: d³ grows by (6/π)(1.28e-8 + 4.22e-10·θ³)·dt [mm³ s⁻¹].
    W = 1.5
    d0 = 1.0
    r, _, _ = GEMB.calculate_grain_size(fill(GEMB.CtoK, 3), dz, fill(300.0, 3), fill(W, 3),
        fill(d0 / 2, 3), zeros(3), fill(1.0, 3), cfs, mp)
    θ = 100 * W / (300.0 * 0.1)
    @test 2 * r[2] ≈ (d0^3 + (6 / π) * (1.28e-8 + 4.22e-10 * θ^3) * dt)^(1 / 3) rtol = 1e-12

    # Marbouty at −10 °C, ρ = 200, ∇T = 30 K/m: F = 0.8, H = 0.8, G = 0.29, Q = 0.016704 mm/d.
    @test GEMB._marbouty_Q(GEMB.CtoK - 10.0, 200.0, 30.0) ≈ 0.8 * 0.8 * 0.29 * 0.09 rtol = 1e-12
    T = [GEMB.CtoK - 13.0, GEMB.CtoK - 10.0, GEMB.CtoK - 7.0]   # 30 K/m
    r, _, _ = GEMB.calculate_grain_size(T, dz, fill(200.0, 3), zeros(3), fill(0.5, 3),
        zeros(3), fill(0.5, 3), cfs, GEMB.ModelParameters(grain_growth_method=:hybrid))
    @test 2 * r[2] ≈ 1.0 + 0.016704 * dt_d rtol = 1e-9
end

@testset "Temperature gradient on a non-uniform grid" begin
    # Linear T(z) through the cell centers: the stencil returns the gradient exactly at every
    # cell, the two edge cells included.
    dz = [0.02, 0.05, 0.11, 0.07, 0.3, 0.2]
    zc = cumsum(dz) .- dz ./ 2
    g = 13.7
    T = 250.0 .+ g .* zc
    for i in eachindex(dz)
        @test GEMB._grain_gradient(T, dz, i, length(dz)) ≈ g rtol = 1e-12
    end
end

@testset "Grain size caps are exact" begin
    # Spherical grains cap at 2 mm diameter, non-spherical at 5 mm. Wet isothermal snow, so
    # Brun growth drives both past their caps and the sphericity is left unchanged.
    cfs = _make_grain_cfs(dt=86400.0 * 400)
    for (gs, cap) in ((1.0, 2.0), (0.0, 5.0))
        r, _, _ = GEMB.calculate_grain_size(fill(GEMB.CtoK, 3), fill(0.1, 3), fill(300.0, 3),
            fill(3.0, 3), fill(cap / 2 - 0.01, 3), zeros(3), fill(gs, 3), cfs,
            GEMB.ModelParameters())
        @test all(2 .* r .== cap)
    end
end

# Taillandier et al. (2007) eqs. 6-13 and Flanner and Zender (2006) eq. 16, written out in BigFloat:
# SSA [cm2 g-1] after t hours at constant temperature Tc [°C] from the initial SSA S0.
function _taillandier_closed(S0, Tc, t; tg=false)
    S0, Tc, t = big.((S0, Tc, t))
    A = tg ? big"0.0961" * S0 - big"3.44" * (Tc + big"1.90") : big"0.0760" * S0 - big"1.76" * (Tc - big"2.96")
    B = tg ? big"0.659" * S0 - big"27.2" * (Tc - big"2.03") : big"0.629" * S0 - big"15.0" * (Tc - big"11.2")
    return B - A * log(t + exp((B - S0) / A))
end

@testset "Taillandier (2007) dry-snow SSA decay" begin
    S0 = GEMB.SSA_NEW_SNOW
    @test S0 ≈ 3 / (0.091 * 0.05) rtol = 1e-15            # 659.3 cm2 g-1, r = 0.05 mm
    for (Tc, G, tg) in ((-15.0, 0.0, false), (-30.0, 10.0, false), (-5.0, 15.0, false),
                        (-15.0, 30.0, true), (-25.0, 60.0, true))
        T = Tc + 273.15
        for dt in (1.0, 3.0, 24.0)                         # 30 days in steps of dt hours
            s = S0
            for _ in 1:round(Int, 720 / dt)
                s = GEMB.taillandier_ssa(s, T, G, dt)
            end
            @test s ≈ Float64(_taillandier_closed(S0, Tc, 720; tg)) rtol = 1e-11
        end
    end
    # At −15 °C under ET, 30 and 90 days give SSA 269.2 and 179.9 cm2 g-1 (r = 0.122, 0.183 mm).
    @test Float64(_taillandier_closed(S0, -15, 720)) ≈ 269.2 rtol = 1e-3
    @test 3 / (0.091 * GEMB.taillandier_ssa(S0, 258.15, 0.0, 2160.0)) ≈ 0.183 rtol = 2e-3
    # t = 0 reproduces the initial SSA; the SSA only decreases; the TG law decays faster.
    @test GEMB.taillandier_ssa(S0, 258.15, 0.0, 0.0) ≈ S0 rtol = 1e-12
    @test GEMB.taillandier_ssa(300.0, 258.15, 0.0, 3.0) < 300.0
    @test GEMB.taillandier_ssa(S0, 258.15, 30.0, 720.0) < GEMB.taillandier_ssa(S0, 258.15, 0.0, 720.0)
    # Floor: the depth-hoar asymptote, 80 cm2 g-1, and nothing coarser is made finer.
    @test GEMB.taillandier_ssa(85.0, 268.15, 50.0, 1e7) == GEMB.SSA_DECAY_MIN
    @test GEMB.taillandier_ssa(60.0, 268.15, 0.0, 3.0) == 60.0
end

# The SNICAR dry-snow law integrated exactly in BigFloat: from r0 at time 0,
# (r − r0 + τ)^(n+1) = τ^(n+1) + (n+1)·(dr/dt)0·τ^n·t, n = 1/κ, t in hours.
function _snicar_closed(τ, κ, drdt0, t; r0=GEMB.SNICAR_RADIUS_FRESH)
    τ, κ, a, t = big.((τ, κ, drdt0, t))
    n = 1 / κ
    return big(r0) + (τ^(n + 1) + (n + 1) * a * τ^n * t)^(1 / (n + 1)) - τ
end

@testset "SNICAR aging table" begin
    τ, κ, a = GEMB.SNICAR_AGING_TAU, GEMB.SNICAR_AGING_KAPPA, GEMB.SNICAR_AGING_DRDT0
    @test size(τ) == size(κ) == size(a) == (8, 31, 11)
    @test !any(isnan, τ) && !any(isnan, κ) && !any(isnan, a)
    @test all(>(0), τ) && all(>(0), κ) && all(>(0), a)
    # Values of snicar_drdt_bst_fit_60_c070416.nc read by `ncdump -p 17,17`, at
    # [density, gradient, temperature] indices; (3, 11, 9) is the node CLM's SnowAge_init logs
    # (263 K, 100 K m-1, 150 kg m-3).
    for (ijk, v) in (((3, 11, 9), (6.320580029826786, 1.5669689418953008, 5.4872755027422642)),
                     ((1, 1, 1), (1.693772564487068, 0.8350176590222218, 0.15639772826553155)),
                     ((8, 31, 11), (0.94772251545821073, 1.9003584126169899, 7.9223126092813345)),
                     ((7, 3, 7), (4.4125164868207811, 0.55824171392066491, 2.7914618871761547)))
        @test (τ[ijk...], κ[ijk...], a[ijk...]) == v
    end
    @test GEMB.snicar_aging_parameters(263.15, 100.0, 150.0) == (τ[3, 11, 9], κ[3, 11, 9], a[3, 11, 9])
    # Nearest node with CLM's indexing, nint((T − 223)/5) + 1 etc., halves away from zero.
    p(T, G, ρ) = GEMB.snicar_aging_parameters(T, G, ρ)
    node(i, j, k) = (τ[i, j, k], κ[i, j, k], a[i, j, k])
    @test p(265.4, 0.0, 350.0) == node(7, 1, 9)        # (265.4 − 223)/5 = 8.48
    @test p(265.6, 0.0, 350.0) == node(7, 1, 10)       # 8.52
    @test p(225.5, 0.0, 350.0) == node(7, 1, 2)        # 0.5 rounds up
    @test p(258.15, 15.0, 350.0) == node(7, 3, 8)      # 1.5 rounds up, to 20 K m-1
    @test p(258.15, 14.9, 75.0) == node(2, 2, 8)       # ρ: 0.5 rounds up, to 100 kg m-3
    # Outside the table: edge nodes; density below 50 is read at 50; NaN fails.
    @test p(200.0, 0.0, 350.0) == node(7, 1, 1)
    @test p(280.0, 0.0, 350.0) == node(7, 1, 11)
    @test p(258.15, 1e6, 350.0) == node(7, 31, 8)
    @test p(258.15, 0.0, 10.0) == node(1, 1, 8)
    @test p(258.15, 0.0, 900.0) == node(8, 1, 8)
    @test_throws InexactError p(NaN, 0.0, 350.0)
end

@testset "SNICAR aging step" begin
    r0 = GEMB.SNICAR_RADIUS_FRESH
    # One dry step: (dr/dt)0·(τ/(r − r0 + τ))^(1/κ)·dt/3600, in BigFloat.
    for (T, G, ρ, r, dt) in ((258.15, 0.0, 350.0, 80.0, 3600.0), (243.15, 40.0, 150.0, 120.0, 10800.0),
                             (271.0, 210.0, 60.0, 300.0, 1800.0))
        τ, κ, a = GEMB.snicar_aging_parameters(T, G, ρ)
        expected = big(r) + big(a) * (big(τ) / (big(r) - big(r0) + big(τ)))^(1 / big(κ)) * big(dt) / 3600
        @test GEMB.snicar_aging_radius(r, T, G, ρ, 0.0, dt) ≈ Float64(expected) rtol = 1e-14
    end
    # Snow finer than r0 grows at the fresh-snow rate; a zero step changes nothing.
    τ, κ, a = GEMB.snicar_aging_parameters(258.15, 0.0, 350.0)
    @test GEMB.snicar_aging_radius(50.0, 258.15, 0.0, 350.0, 0.0, 3600.0) ≈ 50.0 + a rtol = 1e-15
    @test GEMB.snicar_aging_radius(80.0, 258.15, 0.0, 350.0, 0.0, 0.0) == 80.0
    # Wet growth 1e18·C2·f³·dt/(4π r²) on top of the dry step, f capped at 0.1; radius capped.
    dry = GEMB.snicar_aging_radius(200.0, 273.15, 0.0, 400.0, 0.0, 3600.0)
    for f in (0.03, 0.1, 0.4)
        wet = big"1e18" * big"4.22e-13" * big(min(f, 0.1))^3 * 3600 / (4 * big(pi) * big"200.0"^2)
        @test GEMB.snicar_aging_radius(200.0, 273.15, 0.0, 400.0, f, 3600.0) ≈ dry + Float64(wet) rtol = 1e-14
    end
    @test GEMB.snicar_aging_radius(1499.9, 273.15, 300.0, 50.0, 0.1, 86400.0) == GEMB.SNICAR_RADIUS_MAX
    # Stepping integrates the law: at 3-minute steps, 14 days from r0 match the exact solution.
    for (T, G, ρ) in ((253.15, 0.0, 350.0), (263.15, 50.0, 150.0), (273.15, 200.0, 50.0))
        r = r0
        for _ in 1:(14 * 24 * 20)
            r = GEMB.snicar_aging_radius(r, T, G, ρ, 0.0, 180.0)
        end
        @test r ≈ Float64(_snicar_closed(GEMB.snicar_aging_parameters(T, G, ρ)..., 336)) rtol = 2e-3
    end
    # Against the source: Flanner and Zender (2006) eq. 16 with the 16 parameter pairs of their
    # Table 4 (p. 10) gives the 14-day SSA ratio S/S0 = (τ/(336 + τ))^(1/κ). The table in GEMB is
    # a separate fit of the same microphysical model, in radius and from S0 = 60 m2 kg-1, so the
    # two agree to within 12 % (largest at 0 °C, 200 K m-1, 50 kg m-3), not exactly. Putting the
    # table's τ and κ into eq. 16 instead, as hours, gives ratios 5e-4 to 3 times the paper's.
    table4 = ((0, 150, ((43.6, 11.4), (7.1, 6.7), (4.5, 6.1), (3.2, 5.8))),
              (50, 150, ((27.5, 15.3), (47.1, 1.7), (21.0, 1.8), (11.9, 1.9))),
              (200, 50, ((370.6, 0.9), (5.2, 1.9), (2.5, 1.9), (1.5, 1.9))),
              (200, 300, ((47.2, 11.7), (35.0, 1.8), (15.5, 1.9), (8.8, 1.9))))
    for (G, ρ, rows) in table4, (Tc, (τ4, κ4)) in zip((-50, -20, -10, 0), rows)
        paper = (big(τ4) / (336 + big(τ4)))^(1 / big(κ4))
        model = r0 / _snicar_closed(GEMB.snicar_aging_parameters(273.15 + Tc, Float64(G), Float64(ρ))..., 336)
        @test Float64(model) ≈ Float64(paper) rtol = 0.12
    end
    # Refrozen water: mass-weighted with SNICAR_RADIUS_REFROZEN, the fraction capped at 1.
    @test GEMB.snicar_refrozen_radius(100.0, 2.0, 40.0) ≈ 100.0 + 0.05 * 900.0 rtol = 1e-15
    @test GEMB.snicar_refrozen_radius(100.0, 50.0, 40.0) == GEMB.SNICAR_RADIUS_REFROZEN
    @test GEMB.snicar_refrozen_radius(100.0, 0.0, 40.0) == 100.0
end

@testset "Dry-snow growth options in calculate_grain_size" begin
    # A uniform −15 °C column (no gradient) of non-dendritic snow at 350 kg m-3 over firn at 500:
    # under :hybrid Marbouty's G is 0, so dry snow does not grow; Taillandier coarsens it by its
    # closed form and SNICAR by one step of its law; the firn takes Arthern under every method.
    n = 10
    T = fill(258.15, n); dz = fill(0.05, n); ρ = [fill(350.0, 5); fill(500.0, 5)]
    cfs = GEMB.ClimateForcingStep(; dt=86400.0)
    run(method, r0=0.07, w=zeros(n), d=zeros(n)) = GEMB.calculate_grain_size(T, dz, ρ, w,
        fill(r0, n), copy(d), fill(0.5, n), cfs, GEMB.ModelParameters(grain_growth_method=method))
    r_h, = run(:hybrid)
    @test all(r_h[1:5] .== 0.07)
    r_t, = run(:Taillandier2007)
    s0 = 3 / (0.091 * 0.07)
    @test r_t[1] ≈ 3 / (0.091 * GEMB.taillandier_ssa(s0, 258.15, 0.0, 24.0)) rtol = 1e-14
    @test r_t[1] > 0.07
    r_f, = run(:FlannerZender2006)
    @test r_f[1] ≈ GEMB.snicar_aging_radius(70.0, 258.15, 0.0, 350.0, 0.0, 86400.0) / 1000 rtol = 1e-14
    for r in (r_t, r_f)
        @test r[6:10] == r_h[6:10]                       # Arthern at and above 400 kg m-3
    end
    # Taillandier leaves wet and dendritic snow to the common branches.
    w = [0.5; zeros(n - 1)]
    @test run(:Taillandier2007, 0.3, w)[1][1] == run(:hybrid, 0.3, w)[1][1]
    d = [0.8; zeros(n - 1)]
    @test run(:Taillandier2007, 0.05, zeros(n), d)[1][1] == run(:hybrid, 0.05, zeros(n), d)[1][1]
    # SNICAR ages wet snow by its own dry + wet law, at the liquid fraction W/(M + W) and the
    # density (M + W)/dz.
    M = 350.0 * 0.05
    @test run(:FlannerZender2006, 0.3, w)[1][1] ≈
          GEMB.snicar_aging_radius(300.0, 258.15, 0.0, (M + 0.5) / 0.05, 0.5 / (M + 0.5), 86400.0) / 1000 rtol = 1e-14
    # SNICAR ages dendritic snow from its own radius, while dendricity and sphericity evolve as
    # under :hybrid.
    rd, dd, sd = run(:FlannerZender2006, 0.05, zeros(n), d)
    rdh, ddh, sdh = run(:hybrid, 0.05, zeros(n), d)
    @test rd[1] ≈ GEMB.snicar_aging_radius(50.0, 258.15, 0.0, 350.0, 0.0, 86400.0) / 1000 rtol = 1e-14
    @test rd[1] != rdh[1]
    @test dd == ddh && sd == sdh
    @test_nowarn GEMB.validate_parameters(GEMB.ModelParameters(grain_growth_method=:Taillandier2007,
        dt_divisors=[1.0]))
    @test_nowarn GEMB.validate_parameters(GEMB.ModelParameters(grain_growth_method=:FlannerZender2006,
        dt_divisors=[1.0]))
end

@testset "Independent radius in mixing" begin
    # Two dendritic states: dendricity and sphericity mix by mass either way; the radius is
    # re-derived from them by default and mixed by the albedo-conserving rule when independent.
    d, s, r = GEMB.mix_grains(0.9, 0.5, 0.12, 30.0, 1.0, 0.5, 0.05, 10.0)
    di, si, ri = GEMB.mix_grains(0.9, 0.5, 0.12, 30.0, 1.0, 0.5, 0.05, 10.0, true)
    @test (di, si) == (d, s)
    @test d ≈ 0.925 rtol = 1e-15
    @test r == GEMB.dendritic_grain_radius(d, s)
    q = GEMB.ALBEDO_GRAIN_RADIUS_EXPONENT
    @test ri ≈ (0.75 * 0.12^q + 0.25 * 0.05^q)^(1 / q) rtol = 1e-14
    # A non-dendritic state is mixed identically either way.
    @test GEMB.mix_grains(0.0, 0.3, 0.2, 30.0, 1.0, 0.5, 0.05, 10.0, true) ==
          GEMB.mix_grains(0.0, 0.3, 0.2, 30.0, 1.0, 0.5, 0.05, 10.0)
    @test GEMB.independent_grain_radius(GEMB.ModelParameters(grain_growth_method=:FlannerZender2006))
    @test !GEMB.independent_grain_radius(GEMB.ModelParameters())
end
