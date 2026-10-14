# Tests for calculate_albedo
# Translated from MATLAB test_calculate_albedo.m

@testset "Method None" begin
    n = 10
    dz = fill(0.1, n)
    density = fill(400.0, n)
    grain_radius = fill(0.5, n)

    cfs = GEMB.ClimateForcingStep(
        3600.0, 260.0, 100000.0, 0.0, 5.0, 0.0, 0.0,
        100.0, 260.0, 5.0, 0.0, 2.0, 2.0,
        0.1, 0.1, 1.0, 60.0, 0.0, 0.5
    )

    mp = GEMB.ModelParameters(
        albedo_method=:None,
        albedo_fixed=0.75,
        albedo_density_threshold=1023.0,
        density_ice=917.0,
        albedo_snow=0.85,
        albedo_ice=0.48
    )

    a_out, _ = GEMB.calculate_albedo(
        dz, density, grain_radius, cfs, mp)

    @test a_out ≈ 0.75 atol=1e-10
end

@testset "Density threshold override" begin
    n = 10
    dz = fill(0.1, n)
    density = fill(400.0, n)
    density[1] = 350.0  # Higher than threshold
    grain_radius = fill(0.5, n)

    cfs = GEMB.ClimateForcingStep(
        3600.0, 260.0, 100000.0, 0.0, 5.0, 0.0, 0.0,
        100.0, 260.0, 5.0, 0.0, 2.0, 2.0,
        0.1, 0.1, 1.0, 60.0, 0.0, 0.5
    )

    mp = GEMB.ModelParameters(
        albedo_method=:GardnerSharp,
        albedo_fixed=0.4,
        albedo_density_threshold=300.0,  # Low threshold
        density_ice=917.0,
        albedo_snow=0.85,
        albedo_ice=0.48
    )

    a_out, _ = GEMB.calculate_albedo(
        dz, density, grain_radius, cfs, mp)

    @test a_out ≈ 0.4 atol=1e-10
end

@testset "GardnerSharp valid range" begin
    n = 10
    dz = fill(0.1, n)
    density = fill(400.0, n)
    grain_radius = fill(0.5, n)

    cfs = GEMB.ClimateForcingStep(
        3600.0, 260.0, 100000.0, 0.0, 5.0, 0.0, 0.0,
        100.0, 260.0, 5.0, 0.0, 2.0, 2.0,
        0.1, 0.1, 1.0, 60.0, 0.0, 0.5
    )

    mp = GEMB.ModelParameters(
        albedo_method=:GardnerSharp,
        albedo_fixed=0.8,
        albedo_density_threshold=Inf,
        density_ice=917.0,
        albedo_snow=0.85,
        albedo_ice=0.48
    )

    a_out, a_diff_out = GEMB.calculate_albedo(
        dz, density, grain_radius, cfs, mp)

    @test 0.0 <= a_out <= 1.0
    @test 0.0 <= a_diff_out <= 1.0
end

@testset "GreuellKonzelmann" begin
    n = 10
    dz = fill(0.1, n)
    density = fill(400.0, n)
    density[1] = 450.0
    grain_radius = fill(0.5, n)

    cfs = GEMB.ClimateForcingStep(
        3600.0, 260.0, 100000.0, 0.0, 5.0, 0.0, 0.0,
        100.0, 260.0, 5.0, 0.0, 2.0, 2.0,
        0.1, 0.1, 1.0, 60.0, 0.0, 0.5
    )

    mp = GEMB.ModelParameters(
        albedo_method=:GreuellKonzelmann,
        albedo_fixed=0.8,
        albedo_density_threshold=Inf,
        density_ice=900.0,
        albedo_snow=0.85,
        albedo_ice=0.5
    )

    a_out, _ = GEMB.calculate_albedo(
        dz, density, grain_radius, cfs, mp)

    # Expected: a_ice + (density - density_ice) * (a_snow - a_ice) / (300 - density_ice) + 0.05*(cloud - 0.5)
    expected = 0.5 + (450 - 900) * (0.85 - 0.5) / (300 - 900)
    @test a_out ≈ expected atol=1e-5
end

@testset "Exposed ice uses albedo_ice exactly" begin
    # An all-ice column returns `mp.albedo_ice` for both streams, whatever the surface water,
    # the timestep, or the grain-size method.
    n = 20
    dz = fill(0.05, n)
    grain_radius = fill(GEMB.GRAIN_RADIUS_ICE, n)
    for method in (:GardnerSharp, :BrunLefebre), dt in (3600.0, 86400.0)
        mp = GEMB.ModelParameters(albedo_method=method, albedo_ice=0.30)
        density = fill(mp.density_ice, n)
        cfs = GEMB.ClimateForcingStep(; dt, solar_zenith_angle=50.0,
            cloud_optical_thickness=2.0)
        α, α_d = GEMB.calculate_albedo(dz, density, grain_radius, cfs, mp)
        @test α === 0.30
        @test α_d === 0.30
    end
end

@testset "Exposed firn is linear between albedo_ice and 0.65" begin
    mp = GEMB.ModelParameters(albedo_ice=0.35)
    ρ_i, ρ_c = mp.density_ice, GEMB.DENSITY_PORE_CLOSEOFF
    @test GEMB.firn_ice_albedo(ρ_c, mp) ≈ 0.65 rtol = 1e-12
    @test GEMB.firn_ice_albedo(ρ_i, mp) === 0.35
    ρ_mid = 0.5 * (ρ_i + ρ_c)
    @test GEMB.firn_ice_albedo(ρ_mid, mp) ≈ 0.5 * (0.35 + 0.65) rtol = 1e-12
    @test_throws "below pore close-off" GEMB.firn_ice_albedo(600.0, mp)

    n = 20
    density = fill(ρ_mid, n)
    cfs = GEMB.ClimateForcingStep(dt=3600.0)
    α, α_d = GEMB.calculate_albedo(fill(0.05, n), density, fill(1.0, n), cfs, mp)
    @test α ≈ 0.5 rtol = 1e-12
    @test α_d ≈ 0.5 rtol = 1e-12
end

@testset "Thin snow over ice blends toward albedo_ice" begin
    # Snow of depth d over ice: albedo_ice at d → 0, the semi-infinite snow albedo at and above
    # 0.1 m, linear in between (so monotonic in d).
    mp = GEMB.ModelParameters(albedo_ice=0.30)
    cfs = GEMB.ClimateForcingStep(dt=3600.0)
    n = 30
    r_snow = 0.5
    α_snow = GEMB._albedo_gardner(r_snow, 0.0, 0.0, 0.0)
    α_snow_d = GEMB._albedo_gardner(r_snow, 0.0, 50.0, 0.0)
    depths = [1e-9, 0.01, 0.025, 0.05, 0.075, 0.0999, 0.1, 0.2, 0.5]
    αs = Float64[]
    for d in depths
        dz = [d; fill(0.05, n - 1)]
        density = [350.0; fill(mp.density_ice, n - 1)]
        grain_radius = [r_snow; fill(GEMB.GRAIN_RADIUS_ICE, n - 1)]
        α, α_d = GEMB.calculate_albedo(dz, density, grain_radius, cfs, mp)
        f = min(d / 0.1, 1.0)
        @test α ≈ 0.30 + (α_snow - 0.30) * f rtol = 1e-12
        @test α_d ≈ 0.30 + (α_snow_d - 0.30) * f rtol = 1e-12
        push!(αs, α)
    end
    @test αs[1] ≈ 0.30 atol = 1e-6
    @test αs[end] == α_snow
    @test issorted(αs)

    # Over firn, the d → 0 limit is the firn albedo instead.
    ρ_firn = 870.0
    dz = [1e-9; fill(0.05, n - 1)]
    density = [350.0; fill(ρ_firn, n - 1)]
    α, _ = GEMB.calculate_albedo(dz, density, fill(r_snow, n), cfs, mp)
    @test α ≈ GEMB.firn_ice_albedo(ρ_firn, mp) atol = 1e-6
end

@testset "Absorbed shortwave on bare ice scales with 1 - albedo_ice" begin
    # Twenty days of bare ice with no snowfall: the surface stays ice throughout, so the mean
    # absorbed shortwave is (1 − albedo_ice)·SW↓ and the ratio between two runs is
    # (1 − 0.30)/(1 − 0.48).
    n = 20
    time = DateTime(2000, 7, 1) .+ Day.(0:n-1)
    cf = initialize_forcing(time, fill(274.0, n), fill(85000.0, n), zeros(n), fill(4.0, n),
        fill(300.0, n), fill(300.0, n), fill(500.0, n);
        temperature_air_mean=268.0, wind_speed_mean=4.0, precipitation_mean=0.0,
        temperature_observation_height=2.0, wind_observation_height=10.0)
    sw_net = map((0.30, 0.48)) do albedo_ice
        mp = ModelParameters(output_frequency=:daily, albedo_ice=albedo_ice)
        profile = initialize_profile(mp, cf; constant_density=true, constant_temperature=true)
        out = gemb(profile, cf, mp)
        @test all(a -> isapprox(a, albedo_ice; rtol=1e-12), parent(out[:albedo_broadband]))
        sum(parent(out[:shortwave_net]))
    end
    @test sw_net[1] / sw_net[2] ≈ 0.70 / 0.52 rtol = 0.02
end

@testset "Gardner-Sharp closed forms (deep snow, r_e = 0.5 mm)" begin
    r = 0.5
    S = 3.0 / (0.091 * r)
    α_s = 1.48 - S^(-0.07)
    @test α_s ≈ 0.73413 atol = 5e-6
    # Overhead sun, clear sky, clean snow: the pure-snow albedo exactly.
    @test GEMB._albedo_gardner(r, 0.0, 0.0, 0.0) === α_s
    # Zenith term at 60°, clear sky: u = cos 60° = 0.5.
    @test GEMB._albedo_gardner(r, 0.0, 60.0, 0.0) ≈
          α_s + 0.53 * α_s * (1 - α_s) * 0.5^1.2 rtol = 1e-12
    # Cloud term at τ = 5, overhead sun: x = min(√(5/3), 1) = 1, so u = 0.64.
    τ = 5.0
    @test GEMB._albedo_gardner(r, 0.0, 0.0, τ) ≈
          α_s + 0.53 * α_s * (1 - α_s) * (1 - 0.64)^1.2 +
          0.1 * τ * α_s^1.3 / (1 + 1.5 * τ)^α_s rtol = 1e-12
    # Soot at c = 1 ppmw, overhead sun, clear sky.
    dα_c = max(0.04 - α_s, -1.0 / (0.16 + 0.6 * S^0.5 + 1.8 * S^(-0.25)))
    @test GEMB._albedo_gardner(r, 1.0, 0.0, 0.0) ≈ α_s + dα_c rtol = 1e-12
end

@testset "Gardner-Sharp limits and monotonicity" begin
    for r in (0.05, 0.2, 0.5, 1.0, 2.5)
        α_s = 1.48 - (3.0 / (0.091 * r))^(-0.07)
        @test GEMB._albedo_gardner(r, 0.0, 0.0, 0.0) === α_s
        @test issorted([GEMB._albedo_gardner(r, 0.0, sza, 0.0) for sza in 0.0:1.0:89.0])
        @test issorted([GEMB._albedo_gardner(r, 0.0, 30.0, τ) for τ in 0.0:0.5:30.0])
    end
    @test_throws DomainError GEMB._albedo_gardner(0.5, 0.0, 90.0, 1.0)
    @test_throws DomainError GEMB._albedo_gardner(0.5, 0.0, 30.0, -1.0)
end

@testset "Night steps return the diffuse albedo" begin
    n = 10
    dz = fill(0.1, n)
    density = fill(350.0, n)
    grain_radius = fill(0.5, n)
    mp = GEMB.ModelParameters(albedo_method=:GardnerSharp)
    for sza in (90.0, 95.0, 180.0)
        cfs = GEMB.ClimateForcingStep(; dt=3600.0, solar_zenith_angle=sza,
            cloud_optical_thickness=3.0)
        α, α_d = GEMB.calculate_albedo(dz, density, grain_radius, cfs, mp)
        @test α === α_d
        @test α_d == GEMB._albedo_gardner(0.5, 0.0, 50.0, 3.0)
    end
end
