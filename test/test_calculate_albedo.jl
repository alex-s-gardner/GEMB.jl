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

@testset "Brun-Lefebre: thin snow over ice blends toward albedo_ice" begin
    # Snow of depth d over ice: albedo_ice at d → 0, the semi-infinite snow albedo at and above
    # 0.1 m, linear in between (so monotonic in d). Lefebre et al. (2003) p. 3.
    mp = GEMB.ModelParameters(albedo_method=:BrunLefebre, albedo_ice=0.30)
    cfs = GEMB.ClimateForcingStep(dt=3600.0)
    n = 30
    r_snow = 0.5
    α_snow = GEMB._albedo_brun_lefebre(r_snow)
    α_snow_d = α_snow
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

# Gardner & Sharp (2010) eqs. 7–13 (pp. 10–11), written from the paper in BigFloat: a snow layer
# z [m w.e.] deep, grain radius r [mm] and soot c [ppmw], over a semi-infinite layer whose
# clear-sky overhead albedo is α_btm. Eq. 13 replaces dα_c in eqs. 8–12 with dα′_c.
function _gs2010_two_layer(r, c, sza, τ, α_btm, z)
    r, c, sza, τ, α_btm, z = big.((r, c, sza, τ, α_btm, z))
    S = 3 / (big"0.091" * r)                                    # cm² g⁻¹, ρ_ice 0.91 g cm⁻³
    α_S = big"1.48" - S^big"-0.07"                              # eq. 7
    dα_c = max(big"0.04" - α_S,                                 # eq. 8
        -c^big"0.55" / (big"0.16" + big"0.6" * sqrt(S) + big"1.8" * c^big"0.6" * S^big"-0.25"))
    A = min(big"1", big"2.1" * z^(big"1.35" * (1 - α_S) - big"0.1" * c - big"0.13"))
    dα′_c = (α_btm - α_S) + A * ((α_S + dα_c) - α_btm)          # eq. 13
    α′_c = α_S + dα′_c
    μ = cospi(sza / 180)
    x = min(sqrt(τ / (3μ)), big"1")
    u′ = big"0.64" * x + (1 - x) * μ                            # eq. 10
    dα_u = big"0.53" * α_S * (1 - α′_c) * (1 - u′)^big"1.2"     # eq. 9
    dα_τ = big"0.1" * τ * α′_c^big"1.3" / (1 + big"1.5" * τ)^α_S  # eq. 11
    return (α=α_S + dα′_c + dα_u + dα_τ, A=A, α_S=α_S)          # eq. 12
end

@testset "Gardner-Sharp eq. 13: two-layer closed forms" begin
    # A at three (r_e, z) points, read back from the clear-sky overhead albedo, where eqs. 9 and
    # 11 vanish: α = α_btm + A·(α_S − α_btm) for clean snow.
    α_btm = 0.30
    for (r, z) in ((1.0, 0.003), (0.5, 0.003), (0.2, 0.002))
        ref = _gs2010_two_layer(r, 0.0, 0.0, 0.0, α_btm, z)
        @test 0 < ref.A < 1                      # every point lies inside the thin-snow range
        α = GEMB._albedo_gardner(r, 0.0, 0.0, 0.0, α_btm, z)
        @test α ≈ Float64(ref.α) rtol = 1e-12
        @test (α - α_btm) / (GEMB._albedo_gardner(r, 0.0, 0.0, 0.0) - α_btm) ≈ Float64(ref.A) rtol = 1e-10
    end
    # r_e = 1 mm, z = 3 mm w.e.: A = 2.1·0.003^(1.35(1 − α_S) − 0.13), evaluated here.
    α_S = 1.48 - (3 / (0.091 * 1.0))^(-0.07)
    @test Float64(_gs2010_two_layer(1.0, 0.0, 0.0, 0.0, α_btm, 0.003).A) ≈
          2.1 * 0.003^(1.35 * (1 - α_S) - 0.13) rtol = 1e-14

    # Soot, zenith and cloud terms combine with the replaced dα_c: x < 1 and x = 1 regimes, soot
    # inside and on the edge of the fitted range, firn- and ice-like lower layers.
    for (r, c, sza, τ, α_b, z) in ((0.5, 0.0, 60.0, 0.0, 0.30, 0.003),
                                    (0.2, 0.3, 35.0, 1.0, 0.55, 0.001),
                                    (1.0, 0.1, 75.0, 0.25, 0.42, 0.05),
                                    (0.1, 0.0, 0.0, 10.0, 0.30, 0.0005))
        ref = _gs2010_two_layer(r, c, sza, τ, α_b, z)
        @test ref.A < 1
        @test GEMB._albedo_gardner(r, c, sza, τ, α_b, z) ≈ Float64(ref.α) rtol = 1e-12
    end
end

@testset "Gardner-Sharp eq. 13: limits and monotonicity" begin
    α_btm = 0.30
    for (r, c, sza, τ) in ((0.5, 0.0, 0.0, 0.0), (0.5, 0.0, 60.0, 0.0), (1.0, 0.2, 50.0, 10.0),
                           (0.1, 0.0, 30.0, 1.0))
        α_inf = GEMB._albedo_gardner(r, c, sza, τ)
        # z → 0: eqs. 9–12 with α_c = α_btm; at SZA 0 and τ 0 that is α_btm itself.
        lim = _gs2010_two_layer(r, c, sza, τ, α_btm, 0.0).α
        @test GEMB._albedo_gardner(r, c, sza, τ, α_btm, 1e-300) ≈ Float64(lim) rtol = 1e-14
        if sza == 0 && τ == 0
            @test Float64(lim) == α_btm
        end
        # z past the point where A reaches 1: the semi-infinite value exactly.
        p = 1.35 * (1 - _gs2010_two_layer(r, c, 0.0, 0.0, α_btm, 1.0).α_S) - 0.1 * c - 0.13
        z_full = 2.1^(-1 / p)
        for z in (z_full * (1 + 1e-9), 0.5, 10.0, Inf)
            @test GEMB._albedo_gardner(r, c, sza, τ, α_btm, Float64(z)) === α_inf
        end
        # Monotonic in z between the two limits, since α_c^top > α_btm.
        zs = [0.0; 10.0 .^ (-8:0.25:0)]
        αs = [GEMB._albedo_gardner(r, c, sza, τ, α_btm, z) for z in zs]
        @test issorted(αs)
        @test αs[end] === α_inf
    end
    @test_throws "non-negative" GEMB._albedo_gardner(0.5, 0.0, 0.0, 0.0, 0.3, -1e-3)
    @test_throws "must be in [0, 1]" GEMB._albedo_gardner(0.5, 0.0, 0.0, 0.0, 1.3, 0.01)
end

@testset "Gardner-Sharp eq. 13 over firn and over ice in calculate_albedo" begin
    # Two snow cells of different density over a firn or ice cell: z is the summed snow mass
    # in m w.e., and the lower layer is the exposed-surface albedo of that cell.
    mp = GEMB.ModelParameters(albedo_ice=0.30)
    n = 30
    r = 0.4
    for ρ_under in (870.0, mp.density_ice), (sza, τ) in ((0.0, 0.0), (55.0, 0.5), (40.0, 8.0))
        dz = [0.012; 0.008; fill(0.05, n - 2)]
        density = [320.0; 410.0; fill(ρ_under, n - 2)]
        grain_radius = [r; 0.7; fill(GEMB.GRAIN_RADIUS_ICE, n - 2)]
        z = (0.012 * 320.0 + 0.008 * 410.0) / 1000
        # Firn ramp from albedo_ice at 917 to 0.65 at 830 kg m-3, written out.
        α_under = ρ_under >= 917.0 ? 0.30 : 0.30 + (0.65 - 0.30) * (917.0 - ρ_under) / (917.0 - 830.0)
        cfs = GEMB.ClimateForcingStep(; dt=3600.0, solar_zenith_angle=sza,
            cloud_optical_thickness=τ)
        α, α_d = GEMB.calculate_albedo(dz, density, grain_radius, cfs, mp)
        @test _gs2010_two_layer(r, 0.0, sza, τ, α_under, z).A < 1
        @test α ≈ Float64(_gs2010_two_layer(r, 0.0, sza, τ, α_under, z).α) rtol = 1e-12
        @test α_d ≈ Float64(_gs2010_two_layer(r, 0.0, 50.0, τ, α_under, z).α) rtol = 1e-12
        # The term is active: the thin layer is darker than semi-infinite snow.
        @test α < GEMB._albedo_gardner(r, 0.0, sza, τ)
    end

    # Night: the direct albedo is the diffuse one, which carries eq. 13 at 50°.
    dz = [0.01; fill(0.05, n - 1)]
    density = [350.0; fill(mp.density_ice, n - 1)]
    grain_radius = [0.5; fill(GEMB.GRAIN_RADIUS_ICE, n - 1)]
    cfs = GEMB.ClimateForcingStep(; dt=3600.0, solar_zenith_angle=95.0, cloud_optical_thickness=3.0)
    α, α_d = GEMB.calculate_albedo(dz, density, grain_radius, cfs, mp)
    @test α === α_d
    @test α_d ≈ Float64(_gs2010_two_layer(0.5, 0.0, 50.0, 3.0, 0.30, 0.0035).α) rtol = 1e-12
end

@testset "Gardner-Sharp: deep snow is semi-infinite" begin
    # Snow over ice deep enough that A = 1, and a column with no firn or ice at all, both give
    # the semi-infinite albedo bit for bit, matching the paper's eqs. 7–12.
    mp = GEMB.ModelParameters(albedo_ice=0.30)
    n = 30
    cfs = GEMB.ClimateForcingStep(; dt=3600.0, solar_zenith_angle=48.0, cloud_optical_thickness=2.0)
    for density in ([fill(350.0, 4); fill(mp.density_ice, n - 4)], fill(350.0, n))
        dz = fill(0.05, n)        # 4 snow cells: z = 0.07 m w.e., past A = 1 at r 0.5 mm
        α, α_d = GEMB.calculate_albedo(dz, density, fill(0.5, n), cfs, mp)
        @test α === GEMB._albedo_gardner(0.5, 0.0, 48.0, 2.0)
        @test α_d === GEMB._albedo_gardner(0.5, 0.0, 50.0, 2.0)
        @test α ≈ Float64(_gs2010_two_layer(0.5, 0.0, 48.0, 2.0, 0.3, Inf).α) rtol = 1e-12
    end
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
