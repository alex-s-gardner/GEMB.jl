using Test
using GEMB
using Dates

# A diffuse, cloud-free albedo of semi-infinite snow, as a stand-in for the albedo of a step.
const _ALB = (r, c) -> GEMB._albedo_gardner(r, c, GEMB.ALBEDO_DIFFUSE_ZENITH_ANGLE, 0.0, 0.0, Inf)

@testset "albedo_update: minimum weighted change to the observed albedo" begin
    mp = GEMB.ModelParameters()
    alb = _ALB
    for r0 in (0.1, 0.3, 0.8), target in (0.80, 0.70, 0.60, 0.45)
        r, c, ok = GEMB.albedo_update(target, r0, 0.02, alb, mp)
        @test ok
        @test abs(alb(r, c) - target) < 1e-6
    end
    # Minimum of the weighted cost along the albedo contour: moving along it costs more.
    σ1, σ2 = mp.albedo_update_sigma_log_ssa, mp.albedo_update_sigma_log_impurity
    y(r, c) = (log(GEMB._ssa_of_radius(r)), log(c + GEMB.IMPURITY_UPDATE_OFFSET))
    r0, c0 = 0.3, 0.02
    cost(r, c) = ((y(r, c)[1] - y(r0, c0)[1]) / σ1)^2 + ((y(r, c)[2] - y(r0, c0)[2]) / σ2)^2
    r, c, _ = GEMB.albedo_update(0.6, r0, c0, alb, mp)
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
    r, c, ok = GEMB.albedo_update(0.6, 0.3, 0.02, alb, mp0)
    @test ok && r ≈ 0.3 && abs(alb(r, c) - 0.6) < 1e-6
    # Out of reach: stops at the bound closest to the observation.
    r, c, ok = GEMB.albedo_update(0.99, 0.3, 0.02, alb, mp)
    @test !ok && c == 0.0 && r ≈ GEMB.ALBEDO_UPDATE_RADIUS_RANGE[1]
    r, c, ok = GEMB.albedo_update(0.05, 0.3, 0.02, alb, mp)
    @test !ok && c ≈ GEMB.IMPURITY_EFFECTIVE_MAX
    @test_throws "albedo_update_sigma_log_impurity > 0" GEMB.albedo_update(0.6, 0.3, 0.02, alb,
        GEMB.ModelParameters(; albedo_update_sigma_log_impurity = 0.0))
end

@testset "an update makes the albedo GEMB applies to the step the observed one" begin
    mp = GEMB.ModelParameters()
    n = 6; dz = fill(0.01, n); ρ = fill(300.0, n)
    # A sunny step (SZA 60°, 20 % diffuse), an overcast one (τ = 5, all diffuse) and a low-sun one.
    for (sza, τ, f) in ((60.0, 0.0, 0.2), (60.0, 5.0, 1.0), (75.0, 1.0, 0.5))
        cfs = GEMB.ClimateForcingStep(; dt = 3600.0, albedo_observed = 0.6, solar_zenith_angle = sza,
            cloud_optical_thickness = τ, shortwave_downward = 500.0, shortwave_downward_diffuse = 500.0 * f)
        r = fill(0.3, n); bc = fill(0.02, n); bcr = zeros(n)
        @test GEMB.assimilate_albedo!(dz, ρ, r, fill(10.0, n), bc, bcr, cfs, mp)[3]
        a, ad = GEMB.calculate_albedo(dz, ρ, r, bc, bcr, cfs, mp)
        @test abs(GEMB.effective_albedo(a, ad, cfs) - 0.6) < 1e-6
    end
end

@testset "a white-sky observation matches the clear-sky diffuse albedo, whatever the step's sky" begin
    mp = GEMB.ModelParameters()
    n = 6; dz = fill(0.01, n); ρ = fill(300.0, n)
    states = map(((60.0, 0.0, 0.2), (60.0, 5.0, 1.0), (75.0, 1.0, 0.5))) do (sza, τ, f)
        cfs = GEMB.ClimateForcingStep(; dt = 3600.0, albedo_observed_white_sky = 0.6, solar_zenith_angle = sza,
            cloud_optical_thickness = τ, shortwave_downward = 500.0, shortwave_downward_diffuse = 500.0 * f)
        r = fill(0.3, n); bc = fill(0.02, n); bcr = zeros(n)
        @test GEMB.assimilate_albedo!(dz, ρ, r, fill(10.0, n), bc, bcr, cfs, mp)[3]
        c = GEMB.surface_impurity_effective(dz, ρ, bc, bcr, mp)
        @test abs(GEMB.white_sky_albedo(0.0, Inf, Inf, r[1], c, mp) - 0.6) < 1e-6
        # The albedo then applied to an overcast step is brighter than the clear-sky observation.
        a, ad = GEMB.calculate_albedo(dz, ρ, r, bc, bcr, cfs, mp)
        τ > 0 && @test GEMB.effective_albedo(a, ad, cfs) > 0.6
        (r[1], c)
    end
    # The same snow state is found under every sky.
    @test all(s -> isapprox(s[1], states[1][1]; rtol = 1e-6) && isapprox(s[2], states[1][2]; rtol = 1e-6), states)
    # Under all-sky matching, an overcast step finds darker snow than the white-sky one does.
    cfs = GEMB.ClimateForcingStep(; dt = 3600.0, albedo_observed = 0.6, solar_zenith_angle = 60.0,
        cloud_optical_thickness = 5.0, shortwave_downward = 500.0, shortwave_downward_diffuse = 500.0)
    r = fill(0.3, n); bc = fill(0.02, n); bcr = zeros(n)
    GEMB.assimilate_albedo!(dz, ρ, r, fill(10.0, n), bc, bcr, cfs, mp)
    @test GEMB.white_sky_albedo(0.0, Inf, Inf, r[1], GEMB.surface_impurity_effective(dz, ρ, bc, bcr, mp), mp) < 0.6
    # A step may carry only one kind of observation.
    both = GEMB.ClimateForcingStep(; dt = 3600.0, albedo_observed = 0.6, albedo_observed_white_sky = 0.6)
    @test_throws "a step carries both albedo_observed and albedo_observed_white_sky" GEMB.assimilate_albedo!(
        dz, ρ, fill(0.3, n), fill(10.0, n), fill(0.02, n), zeros(n), both, mp)
    # Over exposed ice the observation sets the ice albedo exactly as an all-sky one does.
    ice = fill(900.0, n)
    w = GEMB.ClimateForcingStep(; dt = 3600.0, albedo_observed_white_sky = 0.35)
    a = GEMB.ClimateForcingStep(; dt = 3600.0, albedo_observed = 0.35)
    @test GEMB.assimilate_albedo!(dz, ice, fill(1.0, n), fill(10.0, n), zeros(n), zeros(n), w, mp) ==
          GEMB.assimilate_albedo!(dz, ice, fill(1.0, n), fill(10.0, n), zeros(n), zeros(n), a, mp)
end

@testset "assimilate_albedo!: sets the retained surface film, persists under burial" begin
    mp = GEMB.ModelParameters(; albedo_impurity_depth = 0.02)
    n = 6
    dz = fill(0.01, n); density = fill(350.0, n); r = fill(0.3, n); age = fill(10.0, n)
    bc = fill(0.02, n); bcr = zeros(n)
    total(dz, density, bc, bcr) = GEMB.surface_impurity_effective(dz, density, bc, bcr, mp)
    cfs(obs) = GEMB.ClimateForcingStep(; dt = 3600.0, albedo_observed = obs)
    # No observation: nothing changes.
    @test GEMB.assimilate_albedo!(dz, density, r, age, bc, bcr, cfs(NaN), mp) == (0.0, 0.0, false, mp.albedo_ice)
    @test all(==(0.02), bc) && all(==(0.0), bcr)
    # An observation sets the film in the surface layer only, so that the albedo's loading matches;
    # the mobile loading from snowfall is unchanged.
    d1, d2, ok = GEMB.assimilate_albedo!(dz, density, r, age, bc, bcr, cfs(0.55), mp)
    @test ok && d2 > 0
    @test all(==(0.02), bc)
    @test bcr[1] > 0 && all(==(0.0), bcr[2:end])
    a, ad = GEMB.calculate_albedo(dz, density, r, bc, bcr, cfs(0.55), mp)
    @test abs(GEMB.effective_albedo(a, ad, cfs(0.55)) - 0.55) < 1e-6
    dark = total(dz, density, bc, bcr)
    # The film is an areal amount, read as that amount over albedo_impurity_depth.
    @test dark ≈ 0.02 + bcr[1] / mp.albedo_impurity_depth
    # Fresh clean snow on top: the surface reads the new snow; the film is kept below.
    pushfirst!(dz, 0.05); pushfirst!(density, 120.0); pushfirst!(r, 0.05); pushfirst!(age, 0.5)
    pushfirst!(bc, 0.0); pushfirst!(bcr, 0.0)
    @test total(dz, density, bc, bcr) == 0.0
    # With a minimum age set, fresh snow younger than it is not updated.
    @test GEMB.assimilate_albedo!(dz, density, r, age, bc, bcr, cfs(0.55),
        GEMB.ModelParameters(; albedo_impurity_depth = 0.02, albedo_update_snow_age_min = 3.0))[3] == false
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
    @test (@inferred GEMB.albedo_update(0.55, 0.5, 0.05, _ALB, mp)) isa Tuple{Float64,Float64,Bool}
    n = 4
    cfs = GEMB.ClimateForcingStep(; dt = 3600.0, albedo_observed = 0.55)
    @test (@inferred GEMB.assimilate_albedo!(fill(0.01, n), fill(350.0, n), fill(0.3, n),
        fill(10.0, n), zeros(n), zeros(n), cfs, mp)) isa Tuple{Float64,Float64,Bool,Float64}
    @test (@inferred GEMB.assimilate_albedo!(fill(0.01, n), fill(mp.density_ice, n), fill(0.3, n),
        fill(10.0, n), zeros(n), zeros(n), cfs, mp)) isa Tuple{Float64,Float64,Bool,Float64}
end

@testset "an observation of exposed ice or firn sets the column's ice albedo" begin
    mp = GEMB.ModelParameters()
    n = 4; dz = fill(0.05, n)
    cfs(obs) = GEMB.ClimateForcingStep(; dt = 3600.0, albedo_observed = obs, solar_zenith_angle = 60.0,
        shortwave_downward = 500.0, shortwave_downward_diffuse = 100.0)
    col(ρ) = (dz, fill(ρ, n), fill(GEMB.GRAIN_RADIUS_ICE, n), fill(10.0, n), zeros(n), zeros(n))
    ρ_mid = mp.density_ice - 0.25 * (mp.density_ice - GEMB.DENSITY_PORE_CLOSEOFF)
    for ρ in (mp.density_ice, ρ_mid), obs in (0.25, 0.4, 0.55)
        c = col(ρ)
        _, _, ok, albedo_ice = GEMB.assimilate_albedo!(c..., cfs(obs), mp)
        @test ok
        # The albedo GEMB then applies is the observation, and the snow state is untouched.
        a, ad = GEMB.calculate_albedo(c[1:3]..., c[5:6]..., cfs(obs), mp; albedo_ice)
        @test GEMB.effective_albedo(a, ad, cfs(obs)) ≈ obs rtol = 1e-12
        @test c[3] == fill(GEMB.GRAIN_RADIUS_ICE, n) && c[6] == zeros(n)
    end
    # Out of range: clamped.
    @test GEMB.assimilate_albedo!(col(mp.density_ice)..., cfs(0.05), mp)[4] == GEMB.ALBEDO_ICE_OBSERVED_RANGE[1]
    @test GEMB.assimilate_albedo!(col(mp.density_ice)..., cfs(0.9), mp)[4] == GEMB.ALBEDO_ICE_OBSERVED_RANGE[2]
    # Firn near pore close-off, whose albedo hardly depends on the ice albedo, is not used.
    ρ_light = GEMB.DENSITY_PORE_CLOSEOFF + 10.0
    @test GEMB.assimilate_albedo!(col(ρ_light)..., cfs(0.4), mp; albedo_ice = 0.33) == (0.0, 0.0, false, 0.33)
    # No observation: the ice albedo passed in is returned.
    @test GEMB.assimilate_albedo!(col(mp.density_ice)..., cfs(NaN), mp; albedo_ice = 0.33)[3:4] == (false, 0.33)
    # A snow surface keeps the column's ice albedo and reads it under the snow.
    ρ = [300.0; fill(mp.density_ice, n - 1)]
    c = (dz, ρ, fill(0.3, n), fill(10.0, n), zeros(n), zeros(n))
    _, _, ok, albedo_ice = GEMB.assimilate_albedo!(c..., cfs(0.6), mp; albedo_ice = 0.3)
    @test ok && albedo_ice == 0.3
    a, ad = GEMB.calculate_albedo(c[1:3]..., c[5:6]..., cfs(0.6), mp; albedo_ice)
    @test GEMB.effective_albedo(a, ad, cfs(0.6)) ≈ 0.6 atol = 1e-6
end

@testset "gemb: an observed ice albedo holds until the next observation of exposed ice" begin
    n = 20
    time = DateTime(2000, 7, 1) .+ Day.(0:n-1)
    obs = fill(NaN, n); obs[3] = 0.30; obs[12] = 0.40
    cf = initialize_forcing(time, fill(274.0, n), fill(85000.0, n), zeros(n), fill(4.0, n),
        fill(300.0, n), fill(300.0, n), fill(500.0, n);
        temperature_air_mean=268.0, wind_speed_mean=4.0, precipitation_mean=0.0,
        temperature_observation_height=2.0, wind_observation_height=10.0, albedo_observed=obs)
    mp = ModelParameters(output_frequency=:daily)
    profile = initialize_profile(mp, cf; constant_density=true, constant_temperature=true)
    out = gemb(profile, cf, mp)
    a = parent(out[:albedo_broadband]); ai = parent(out[:albedo_ice])
    @test all(≈(mp.albedo_ice; rtol=1e-12), a[1:2])
    @test all(≈(0.30; rtol=1e-12), a[3:11]) && all(≈(0.30; rtol=1e-12), ai[3:11])
    @test all(≈(0.40; rtol=1e-12), a[12:end])
    @test parent(out[:albedo_update_count]) == ifelse.(isnan.(obs), 0.0, 1.0)
end
