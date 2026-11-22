using Test
using GEMB

@testset "albedo_update: minimum weighted change to the observed albedo" begin
    mp = GEMB.ModelParameters()
    alb(r, c) = GEMB._intrinsic_albedo(r, c, 0.0, Inf)
    for r0 in (0.1, 0.3, 0.8), target in (0.80, 0.70, 0.60, 0.45)
        r, c, ok = GEMB.albedo_update(target, r0, 0.02, 0.0, Inf, mp)
        @test ok
        @test abs(alb(r, c) - target) < 1e-6
    end
    # Minimum of the weighted cost along the albedo contour: moving along it costs more.
    σ1, σ2 = mp.albedo_update_sigma_log_ssa, mp.albedo_update_sigma_log_impurity
    y(r, c) = (log(GEMB._ssa_of_radius(r)), log(c + GEMB.BLACK_CARBON_UPDATE_OFFSET))
    r0, c0 = 0.3, 0.02
    cost(r, c) = ((y(r, c)[1] - y(r0, c0)[1]) / σ1)^2 + ((y(r, c)[2] - y(r0, c0)[2]) / σ2)^2
    r, c, _ = GEMB.albedo_update(0.6, r0, c0, 0.0, Inf, mp)
    for f in (0.9, 0.95, 1.05, 1.1)
        rp = r * f
        # the loading that keeps the albedo at 0.6 with radius rp
        lo, hi = 0.0, GEMB.BLACK_CARBON_EQUIVALENT_MAX
        for _ in 1:80
            m = (lo + hi) / 2
            alb(rp, m) > 0.6 ? (lo = m) : (hi = m)
        end
        @test cost(rp, (lo + hi) / 2) >= cost(r, c) - 1e-9
    end
    # No SSA weight: only the loading changes.
    mp0 = GEMB.ModelParameters(; albedo_update_sigma_log_ssa = 0.0)
    r, c, ok = GEMB.albedo_update(0.6, 0.3, 0.02, 0.0, Inf, mp0)
    @test ok && r ≈ 0.3 && abs(alb(r, c) - 0.6) < 1e-6
    # Out of reach: stops at the bound closest to the observation.
    r, c, ok = GEMB.albedo_update(0.99, 0.3, 0.02, 0.0, Inf, mp)
    @test !ok && c == 0.0 && r ≈ GEMB.ALBEDO_UPDATE_RADIUS_RANGE[1]
    r, c, ok = GEMB.albedo_update(0.05, 0.3, 0.02, 0.0, Inf, mp)
    @test !ok && c ≈ GEMB.BLACK_CARBON_EQUIVALENT_MAX
    @test_throws "albedo_update_sigma_log_impurity > 0" GEMB.albedo_update(0.6, 0.3, 0.02, 0.0, Inf,
        GEMB.ModelParameters(; albedo_update_sigma_log_impurity = 0.0))
end

@testset "assimilate_albedo!: sets the near-surface layers, persists under burial" begin
    mp = GEMB.ModelParameters(; albedo_impurity_depth = 0.02)
    n = 6
    dz = fill(0.01, n); density = fill(350.0, n); r = fill(0.3, n); bc = fill(0.02, n)
    cfs(obs) = GEMB.ClimateForcingStep(; dt = 3600.0, albedo_observed = obs)
    # No observation: nothing changes.
    @test GEMB.assimilate_albedo!(dz, density, r, bc, cfs(NaN), mp) == (0.0, 0.0, false)
    @test all(==(0.02), bc)
    # An observation updates the top 0.02 m (two 0.01 m layers) only.
    d1, d2, ok = GEMB.assimilate_albedo!(dz, density, r, bc, cfs(0.55), mp)
    @test ok && d2 > 0
    @test bc[1] == bc[2] > 0.02 && all(==(0.02), bc[3:end])
    @test abs(GEMB._intrinsic_albedo(r[1], GEMB.surface_black_carbon(dz, density, bc, mp), 0.0, Inf) - 0.55) < 1e-6
    dark = bc[1]
    # Fresh clean snow on top: the surface reads the new snow; the dark layer is kept below.
    pushfirst!(dz, 0.05); pushfirst!(density, 120.0); pushfirst!(r, 0.05); pushfirst!(bc, 0.0)
    @test GEMB.surface_black_carbon(dz, density, bc, mp) == 0.0
    # Melting the fresh snow away re-exposes the dark layer with no new observation.
    popfirst!(dz); popfirst!(density); popfirst!(r); popfirst!(bc)
    @test GEMB.surface_black_carbon(dz, density, bc, mp) == dark
    # Firn or ice at the surface: no update.
    density[1] = 850.0
    @test GEMB.assimilate_albedo!(dz, density, r, bc, cfs(0.55), mp)[3] == false
end

@testset "albedo with uniform background loading equals the scalar form" begin
    n = 10
    dz = fill(0.05, n); density = fill(300.0, n); grain_radius = fill(0.5, n)
    cfs = GEMB.ClimateForcingStep(; dt = 3600.0, cloud_optical_thickness = 1.0,
        solar_zenith_angle = 60.0, black_carbon_snow = 0.05)
    mp = GEMB.ModelParameters()
    @test GEMB.calculate_albedo(dz, density, grain_radius, fill(0.05, n), cfs, mp) ==
          GEMB.calculate_albedo(dz, density, grain_radius, cfs, mp)
end
