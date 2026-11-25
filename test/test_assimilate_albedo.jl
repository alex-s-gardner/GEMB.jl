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
    y(r, c) = (log(GEMB._ssa_of_radius(r)), log(c + GEMB.IMPURITY_UPDATE_OFFSET))
    r0, c0 = 0.3, 0.02
    cost(r, c) = ((y(r, c)[1] - y(r0, c0)[1]) / σ1)^2 + ((y(r, c)[2] - y(r0, c0)[2]) / σ2)^2
    r, c, _ = GEMB.albedo_update(0.6, r0, c0, 0.0, Inf, mp)
    for f in (0.9, 0.95, 1.05, 1.1)
        rp = r * f
        # the loading that keeps the albedo at 0.6 with radius rp
        lo, hi = 0.0, GEMB.IMPURITY_EFFECTIVE_MAX
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
    @test !ok && c ≈ GEMB.IMPURITY_EFFECTIVE_MAX
    @test_throws "albedo_update_sigma_log_impurity > 0" GEMB.albedo_update(0.6, 0.3, 0.02, 0.0, Inf,
        GEMB.ModelParameters(; albedo_update_sigma_log_impurity = 0.0))
end

@testset "assimilate_albedo!: sets the retained surface film, persists under burial" begin
    mp = GEMB.ModelParameters(; albedo_impurity_depth = 0.02)
    n = 6
    dz = fill(0.01, n); density = fill(350.0, n); r = fill(0.3, n); age = fill(10.0, n)
    bc = fill(0.02, n); bcr = zeros(n)
    total(dz, density, bc, bcr) = GEMB.surface_impurity_effective(dz, density, bc, bcr, mp)
    cfs(obs) = GEMB.ClimateForcingStep(; dt = 3600.0, albedo_observed = obs)
    # No observation: nothing changes.
    @test GEMB.assimilate_albedo!(dz, density, r, age, bc, bcr, cfs(NaN), mp) == (0.0, 0.0, false)
    @test all(==(0.02), bc) && all(==(0.0), bcr)
    # An observation sets the film in the surface layer only, so that the albedo's loading matches;
    # the mobile loading from snowfall is unchanged.
    d1, d2, ok = GEMB.assimilate_albedo!(dz, density, r, age, bc, bcr, cfs(0.55), mp)
    @test ok && d2 > 0
    @test all(==(0.02), bc)
    @test bcr[1] > 0 && all(==(0.0), bcr[2:end])
    @test abs(GEMB._intrinsic_albedo(r[1], total(dz, density, bc, bcr), 0.0, Inf) - 0.55) < 1e-6
    dark = total(dz, density, bc, bcr)
    # The film is an areal amount, read as that amount over albedo_impurity_depth.
    @test dark ≈ 0.02 + bcr[1] / mp.albedo_impurity_depth
    # Fresh clean snow on top: the surface reads the new snow; the film is kept below.
    pushfirst!(dz, 0.05); pushfirst!(density, 120.0); pushfirst!(r, 0.05); pushfirst!(age, 0.5)
    pushfirst!(bc, 0.0); pushfirst!(bcr, 0.0)
    @test total(dz, density, bc, bcr) == 0.0
    # The fresh snow is younger than albedo_update_snow_age_min: an observation is not assimilated.
    @test GEMB.assimilate_albedo!(dz, density, r, age, bc, bcr, cfs(0.55), mp)[3] == false
    # Melting the fresh snow away re-exposes the film with no new observation.
    popfirst!(dz); popfirst!(density); popfirst!(r); popfirst!(age); popfirst!(bc); popfirst!(bcr)
    @test total(dz, density, bc, bcr) == dark
    # An observation brighter than the mobile loading alone allows: the film is 0.
    bc2 = fill(5.0, n); bcr2 = fill(1.0, n)
    GEMB.assimilate_albedo!(copy(dz), copy(density), copy(r), copy(age), bc2, bcr2, cfs(0.8), mp)
    @test bc2 == fill(5.0, n) && bcr2[1] == 0.0 && bcr2[2:end] == fill(1.0, n - 1)
    # Firn or ice at the surface: no update.
    density[1] = 850.0
    @test GEMB.assimilate_albedo!(dz, density, r, age, bc, bcr, cfs(0.55), mp)[3] == false
end

@testset "the retained film keeps its strength through thinning, merging and splitting" begin
    mp = GEMB.ModelParameters(; albedo_impurity_depth = 0.02)
    film(dz, bcr) = GEMB.surface_impurity_effective(dz, fill(350.0, length(dz)), zero(dz), bcr, mp)
    # The film reads the surface layer's film over albedo_impurity_depth, whatever its thickness.
    @test film([0.04, 0.05], [0.04, 0.0]) ≈ 2.0
    @test film([0.01, 0.05], [0.04, 0.0]) ≈ 2.0
    # The film of buried layers is not read.
    @test film([0.04, 0.05], [0.0, 3.0]) == 0.0
    # Compaction, sublimation and erosion change the surface layer's thickness and density but move
    # no impurity, so the film reads the same: the reading takes neither.
    @test GEMB.surface_impurity_effective([0.004, 0.05], [600.0, 350.0], zeros(2), [0.04, 0.0], mp) ≈ 2.0
    # A thin surface layer merged into the cell below carries its film mass with it: read the same.
    cols = GEMB.column_state(fill(270.0, 3), [0.01, 0.05, 0.05], fill(350.0, 3), zeros(3), fill(0.5, 3),
        zeros(3), fill(0.5, 3), fill(10.0, 3), zeros(3), [4.0, 0.0, 0.0])
    before = film(cols.dz, cols.impurity_effective_retained)
    M1, M2 = cols.dz[1] * 350.0, cols.dz[2] * 350.0
    GEMB.merge_snow_pair!(cols, 1, 2, M1, M2, mp)
    GEMB.close_slot!(cols, [1])
    @test film(cols.dz, cols.impurity_effective_retained) ≈ before rtol = 1e-12
    # Splitting the surface layer keeps the film in the upper half; the total film is conserved.
    cols = GEMB.column_state(fill(270.0, 2), [0.10, 0.05], fill(350.0, 2), zeros(2), fill(0.5, 2),
        zeros(2), fill(0.5, 2), fill(10.0, 2), zeros(2), [1.0, 0.0])
    mass(c) = sum(c.impurity_effective_retained)
    m0 = mass(cols); f0 = film(cols.dz, cols.impurity_effective_retained)
    GEMB.split_cell!(cols, 1)
    @test mass(cols) ≈ m0 rtol = 1e-14
    @test cols.impurity_effective_retained[2] == 0.0
    @test film(cols.dz, cols.impurity_effective_retained) ≈ f0 rtol = 1e-14
    # A split below the surface keeps that cell's film in its upper half too.
    cols = GEMB.column_state(fill(270.0, 3), [0.05, 0.10, 0.05], fill(350.0, 3), zeros(3), fill(0.5, 3),
        zeros(3), fill(0.5, 3), fill(10.0, 3), zeros(3), [0.0, 1.0, 0.0])
    GEMB.split_cell!(cols, 2)
    @test cols.impurity_effective_retained == [0.0, 1.0, 0.0, 0.0]
end

@testset "albedo with uniform background loading equals the scalar form" begin
    n = 10
    dz = fill(0.05, n); density = fill(300.0, n); grain_radius = fill(0.5, n)
    cfs = GEMB.ClimateForcingStep(; dt = 3600.0, cloud_optical_thickness = 1.0,
        solar_zenith_angle = 60.0, black_carbon_snow = 0.05)
    mp = GEMB.ModelParameters()
    @test GEMB.calculate_albedo(dz, density, grain_radius, fill(0.05, n), zeros(n), cfs, mp) ==
          GEMB.calculate_albedo(dz, density, grain_radius, cfs, mp)
end

@testset "albedo reads the mobile loading plus the surface film" begin
    n = 10
    dz = fill(0.05, n); density = fill(300.0, n); grain_radius = fill(0.5, n)
    cfs = GEMB.ClimateForcingStep(; dt = 3600.0, cloud_optical_thickness = 1.0,
        solar_zenith_angle = 60.0)
    mp = GEMB.ModelParameters()
    # A surface film f reads as f/albedo_impurity_depth of mobile loading spread through the
    # near-surface snow.
    f = 0.001
    equivalent = 0.03 + f / mp.albedo_impurity_depth
    @test GEMB.calculate_albedo(dz, density, grain_radius, fill(0.03, n), [f; zeros(n - 1)], cfs, mp) ==
          GEMB.calculate_albedo(dz, density, grain_radius, fill(equivalent, n), zeros(n), cfs, mp)
    @test GEMB.calculate_albedo_ground(dz, density, grain_radius, fill(0.03, n), [f; zeros(n - 1)],
        cfs, mp) == GEMB.calculate_albedo_ground(dz, density, grain_radius, fill(equivalent, n),
        zeros(n), cfs, mp)
end

@testset "albedo update is type-stable" begin
    # The trimmed C library cannot compile a dynamic call, so these must infer concretely.
    mp = GEMB.ModelParameters()
    @test (@inferred GEMB.albedo_update(0.55, 0.5, 0.05, 0.5, Inf, mp)) isa Tuple{Float64,Float64,Bool}
    n = 4
    cfs = GEMB.ClimateForcingStep(; dt = 3600.0, albedo_observed = 0.55)
    @test (@inferred GEMB.assimilate_albedo!(fill(0.01, n), fill(350.0, n), fill(0.3, n),
        fill(10.0, n), zeros(n), zeros(n), cfs, mp)) isa Tuple{Float64,Float64,Bool}
end
