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
# z [m w.e.] deep, grain radius r [mm] and soot c [ppmw], over a semi-infinite layer of albedo
# α_btm. Eq. 13 replaces dα_c in eqs. 8–12 with dα′_c. `α_paper` is eq. 12 as printed; `α` is
# GEMB's form, which multiplies the zenith (eq. 9) and cloud (eq. 11) increments by A.
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
    k_u = big"0.53" * α_S * (1 - u′)^big"1.2"
    k_τ = big"0.1" * τ / (1 + big"1.5" * τ)^α_S
    dα_u = k_u * (1 - α′_c)                                     # eq. 9
    dα_τ = k_τ * α′_c^big"1.3"                                  # eq. 11
    return (α=α′_c + A * (dα_u + dα_τ), α_paper=α_S + dα′_c + dα_u + dα_τ,  # eq. 12
        A=A, α_S=α_S, α_c=α_S + dα_c, F=dα_u + dα_τ, k_u=k_u, k_τ=k_τ)
end

@testset "Gardner-Sharp eq. 13: two-layer closed forms" begin
    # A at three (r_e, z) points, read back from the clear-sky overhead albedo, where eqs. 9 and
    # 11 vanish and GEMB's form is the paper's: α = α_btm + A·(α_S − α_btm) for clean snow.
    α_btm = 0.30
    for (r, z) in ((1.0, 0.003), (0.5, 0.003), (0.2, 0.002))
        ref = _gs2010_two_layer(r, 0.0, 0.0, 0.0, α_btm, z)
        @test 0 < ref.A < 1                      # every point lies inside the thin-snow range
        @test ref.α == ref.α_paper
        α = GEMB._albedo_gardner(r, 0.0, 0.0, 0.0, α_btm, z)
        @test α ≈ Float64(ref.α) rtol = 1e-12
        @test (α - α_btm) / (GEMB._albedo_gardner(r, 0.0, 0.0, 0.0) - α_btm) ≈ Float64(ref.A) rtol = 1e-10
    end
    # r_e = 1 mm, z = 3 mm w.e.: A = 2.1·0.003^(1.35(1 − α_S) − 0.13), evaluated here.
    α_S = 1.48 - (3 / (0.091 * 1.0))^(-0.07)
    @test Float64(_gs2010_two_layer(1.0, 0.0, 0.0, 0.0, α_btm, 0.003).A) ≈
          2.1 * 0.003^(1.35 * (1 - α_S) - 0.13) rtol = 1e-14

    # Soot, zenith and cloud terms combine with the replaced dα_c: x < 1 and x = 1 regimes, soot
    # inside and on the edge of the fitted range, firn- and ice-like lower layers. GEMB's albedo
    # is below the paper's by exactly (1 − A)·(dα_u + dα_τ).
    for (r, c, sza, τ, α_b, z) in ((0.5, 0.0, 60.0, 0.0, 0.30, 0.003),
                                    (0.2, 0.3, 35.0, 1.0, 0.55, 0.001),
                                    (1.0, 0.1, 75.0, 0.25, 0.42, 0.05),
                                    (0.1, 0.0, 0.0, 10.0, 0.30, 0.0005),
                                    (1.0, 0.0, 60.0, 0.0, 0.48, 0.01))
        ref = _gs2010_two_layer(r, c, sza, τ, α_b, z)
        @test ref.A < 1
        α = GEMB._albedo_gardner(r, c, sza, τ, α_b, z)
        @test α ≈ Float64(ref.α) rtol = 1e-12
        @test Float64(ref.α_paper) - α ≈ Float64((1 - ref.A) * ref.F) rtol = 1e-9
    end
end

@testset "Gardner-Sharp eq. 13: exact at exposure, continuous, never throws" begin
    # z = 0 returns the exposed albedo exactly at every geometry, including dark ice at low sun
    # under clear sky with fresh snow on top, where eqs. 9–12 alone exceed the exposed albedo.
    for r in (0.05, 0.5, 1.5), sza in (0.0, 60.0, 75.0, 85.0), τ in (0.0, 10.0), c in (0.0, 0.2),
        α_x in (0.20, 0.30, 0.48, 0.62)

        ref0 = _gs2010_two_layer(r, c, sza, τ, α_x, 0)
        @test ref0.A == 0 && ref0.α == big(α_x)
        @test GEMB._albedo_gardner(r, c, sza, τ, α_x, 0.0) === α_x

        # 1e-300 m w.e. and small z: within A·(|α_c^top − α_x| + k_u + k_τ) of the exposed albedo
        # (attained at SZA 0 and τ 0, hence the allowance of a few ulp for rounding in α), and
        # equal to the BigFloat reference.
        for z in (1e-300, 1e-100, 1e-30, 1e-12, 1e-6)
            ref = _gs2010_two_layer(r, c, sza, τ, α_x, z)
            α = GEMB._albedo_gardner(r, c, sza, τ, α_x, z)
            @test α ≈ Float64(ref.α) rtol = 1e-12
            @test abs(α - α_x) <= Float64(ref.A * (abs(ref.α_c - α_x) + ref.k_u + ref.k_τ)) +
                                  4 * eps(α_x)
        end
        @test abs(GEMB._albedo_gardner(r, c, sza, τ, α_x, 1e-300) - α_x) <= 1e-12

        # Where the snow's α_c is at least the exposed albedo, the albedo rises with z from the
        # exposed value to the semi-infinite value, which it reaches exactly at A = 1.
        zs = [0.0; 10.0 .^ (-300:10:-10); 10.0 .^ (-9:0.25:0)]
        αs = [GEMB._albedo_gardner(r, c, sza, τ, α_x, z) for z in zs]
        α_inf = GEMB._albedo_gardner(r, c, sza, τ)
        if ref0.α_c >= α_x
            @test issorted(αs)
        end
        @test αs[end] === α_inf
        p = Float64(1.35 * (1 - ref0.α_S) - 0.1 * big(c) - 0.13)
        for z in (2.1^(-1 / p) * (1 + 1e-9), 0.5, 10.0, Inf)
            @test GEMB._albedo_gardner(r, c, sza, τ, α_x, z) === α_inf
        end
    end

    # Finite z against the BigFloat reference across the fitted depths.
    for r in (0.05, 0.5, 1.5), (sza, τ) in ((0.0, 0.0), (60.0, 0.0), (75.0, 10.0), (85.0, 0.0)),
        c in (0.0, 0.2), α_x in (0.20, 0.48), z in (1e-4, 0.003, 0.01, 0.05)
        ref = _gs2010_two_layer(r, c, sza, τ, α_x, z)
        ref.A < 1 || continue
        @test GEMB._albedo_gardner(r, c, sza, τ, α_x, z) ≈ Float64(ref.α) rtol = 1e-12
    end

    # Defined over the whole domain: finite, in [0, 1] for clean snow under clear sky, no throw.
    grid = Iterators.product((0.05, 0.1, 0.5, 1.0, 2.5), 0.0:5.0:89.0,
        (0.0, 0.5, 3.0, 30.0, 100.0), (0.0, 0.2, 2.0), (0.0, 0.2, 0.48, 0.65, 1.0),
        (0.0, 1e-300, 1e-9, 1e-3, 0.02, 1.0, Inf))
    αs = [(c == 0.0 && τ == 0.0, GEMB._albedo_gardner(r, c, sza, τ, α_x, z))
          for (r, sza, τ, c, α_x, z) in grid]
    @test all(isfinite ∘ last, αs)
    @test all(((clean, α),) -> !clean || 0.0 <= α <= 1.0, αs)
    @test_throws "non-negative" GEMB._albedo_gardner(0.5, 0.0, 0.0, 0.0, 0.3, -1e-3)
    @test_throws "must be in [0, 1]" GEMB._albedo_gardner(0.5, 0.0, 0.0, 0.0, 1.3, 0.01)
end

@testset "Gardner-Sharp eq. 13 over firn and over ice in calculate_albedo" begin
    # Two snow cells of different density over a firn or ice cell: z is the summed snow mass
    # in m w.e., and the lower layer's albedo is the exposed albedo of that cell.
    n = 30
    r = 0.4
    for α_ice in (0.20, 0.30), ρ_under in (870.0, 917.0),
        (sza, τ) in ((0.0, 0.0), (55.0, 0.5), (40.0, 8.0), (85.0, 0.0))
        mp = GEMB.ModelParameters(albedo_ice=α_ice)
        dz = [0.012; 0.008; fill(0.05, n - 2)]
        density = [320.0; 410.0; fill(ρ_under, n - 2)]
        grain_radius = [r; 0.7; fill(GEMB.GRAIN_RADIUS_ICE, n - 2)]
        z = (0.012 * 320.0 + 0.008 * 410.0) / 1000
        # Firn ramp from albedo_ice at 917 to 0.65 at 830 kg m-3, written out.
        α_under = ρ_under >= 917.0 ? α_ice :
                  α_ice + (0.65 - α_ice) * (917.0 - ρ_under) / (917.0 - 830.0)
        cfs = GEMB.ClimateForcingStep(; dt=3600.0, solar_zenith_angle=sza,
            cloud_optical_thickness=τ)
        α, α_d = GEMB.calculate_albedo(dz, density, grain_radius, cfs, mp)
        ref = _gs2010_two_layer(r, 0.0, sza, τ, α_under, z)
        @test ref.A < 1
        @test α ≈ Float64(ref.α) rtol = 1e-12
        @test α_d ≈ Float64(_gs2010_two_layer(r, 0.0, 50.0, τ, α_under, z).α) rtol = 1e-12
        # The term is active: the thin layer is darker than semi-infinite snow.
        @test α < GEMB._albedo_gardner(r, 0.0, sza, τ)

        # Snow of 1e-300 m w.e.: the exposed surface's albedo, both streams.
        exposed = GEMB.calculate_albedo(fill(0.05, n), fill(ρ_under, n),
            fill(GEMB.GRAIN_RADIUS_ICE, n), cfs, mp)
        @test exposed[1] === exposed[2] === GEMB.firn_ice_albedo(ρ_under, mp)
        dz = [1e-300 * 1000 / 350.0; fill(0.05, n - 1)]
        density = [350.0; fill(ρ_under, n - 1)]
        for r_top in (0.05, 0.5)
            α, α_d = GEMB.calculate_albedo(dz, density, [r_top; fill(GEMB.GRAIN_RADIUS_ICE, n - 1)],
                cfs, mp)
            @test abs(α - exposed[1]) <= 1e-12
            @test abs(α_d - exposed[2]) <= 1e-12
        end
    end

    # Night: the direct albedo is the diffuse one, which carries eq. 13 at 50°.
    mp = GEMB.ModelParameters(albedo_ice=0.30)
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
