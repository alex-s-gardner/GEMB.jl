# Tests for calculate_density - matches MATLAB test_calculate_density.m

# Helper to create a ClimateForcingStep for density tests
function _make_density_cfs(; dt=86400.0 * 30, precipitation_mean=200.0, temperature_air_mean=250.0)
    return GEMB.ClimateForcingStep(
        dt,             # dt [s]
        260.0,          # temperature_air
        80000.0,        # pressure_air
        0.0,            # precipitation
        5.0,            # wind_speed
        100.0,          # shortwave_downward
        250.0,          # longwave_downward
        300.0,          # vapor_pressure
        temperature_air_mean,  # temperature_air_mean
        5.0,            # wind_speed_mean
        precipitation_mean,    # precipitation_mean
        2.0,            # temperature_observation_height
        10.0,           # wind_observation_height
        0.0,            # black_carbon_snow
        0.0,            # black_carbon_ice
        0.0,            # cloud_optical_thickness
        0.0,            # solar_zenith_angle
        0.0,            # shortwave_downward_diffuse
        0.1,            # cloud_fraction
    )
end

@testset "Mass conservation (HerronLangway)" begin
    n = 10
    t_vec = 260.0 * ones(n)
    dz = 0.5 * ones(n)
    density = collect(range(300.0, 700.0, length=n))
    grain_radius = 0.5 * ones(n)

    mp = GEMB.ModelParameters(
        density_ice=917.0,
        densification_method=:HerronLangway,
        densification_coeffs_M01=:Ant_RACMO_GS_SW0,
    )
    cfs = _make_density_cfs()

    mass_initial = density .* dz
    density_before = copy(density)
    dz_before = copy(dz)

    (dz_out, density_out) = GEMB.calculate_density(t_vec, copy(dz), copy(density), grain_radius, cfs, mp)

    mass_final = density_out .* dz_out

    # Mass must be conserved during densification
    @test mass_final ≈ mass_initial atol = 1e-10

    # Density should increase over time
    @test all(density_out .> density_before)

    # Thickness should decrease over time
    @test all(dz_out .< dz_before)
end

@testset "Density clamping at ice density" begin
    n = 10
    t_vec = 260.0 * ones(n)
    dz = 0.5 * ones(n)
    grain_radius = 0.5 * ones(n)

    # Set density very close to ice density
    density_near_ice = fill(917.0 - 0.1, n)

    # Long timestep to force overshoot
    cfs = _make_density_cfs(dt=86400.0 * 365 * 100)

    mp = GEMB.ModelParameters(
        density_ice=917.0,
        densification_method=:HerronLangway,
        densification_coeffs_M01=:Ant_RACMO_GS_SW0,
    )

    (_, density_out) = GEMB.calculate_density(t_vec, copy(dz), copy(density_near_ice), grain_radius, cfs, mp)

    # Density must be clamped at density_ice
    @test all(density_out .<= 917.0 + 1e-10)
end

@testset "HerronLangway densification" begin
    n = 10
    t_vec = 260.0 * ones(n)
    dz = 0.5 * ones(n)
    density = collect(range(300.0, 700.0, length=n))
    grain_radius = 0.5 * ones(n)

    mp = GEMB.ModelParameters(
        density_ice=917.0,
        densification_method=:HerronLangway,
        densification_coeffs_M01=:Ant_RACMO_GS_SW0,
    )
    cfs = _make_density_cfs()

    density_before = copy(density)
    (_, density_out) = GEMB.calculate_density(t_vec, copy(dz), copy(density), grain_radius, cfs, mp)

    @test all(density_out .> density_before)
end

@testset "Arthern densification" begin
    n = 10
    t_vec = 260.0 * ones(n)
    dz = 0.5 * ones(n)
    density = collect(range(300.0, 700.0, length=n))
    grain_radius = 0.5 * ones(n)

    mp = GEMB.ModelParameters(
        density_ice=917.0,
        densification_method=:Arthern,
        densification_coeffs_M01=:Ant_RACMO_GS_SW0,
    )
    cfs = _make_density_cfs()

    density_before = copy(density)
    (_, density_out) = GEMB.calculate_density(t_vec, copy(dz), copy(density), grain_radius, cfs, mp)

    @test all(density_out .> density_before)
end

@testset "ArthernB grain size sensitivity" begin
    n = 10
    t_vec = 260.0 * ones(n)
    dz = 0.5 * ones(n)
    density = collect(range(300.0, 700.0, length=n))
    grain_radius = 0.5 * ones(n)

    mp = GEMB.ModelParameters(
        density_ice=917.0,
        densification_method=:ArthernB,
        densification_coeffs_M01=:Ant_RACMO_GS_SW0,
    )
    cfs = _make_density_cfs()

    # Standard grain size
    (_, density_std) = GEMB.calculate_density(t_vec, copy(dz), copy(density), grain_radius, cfs, mp)

    # Larger grain size -> slower densification
    grain_radius_large = grain_radius * 2
    (_, density_large) = GEMB.calculate_density(t_vec, copy(dz), copy(density), grain_radius_large, cfs, mp)

    # ArthernB: rate is proportional to 1/r^2
    # Larger grains -> slower densification -> lower final density
    # Exclude top layer (overburden = 0) and values near ice density
    valid_mask = (density_std .< 917.0 - 1.0)
    check_indices = findall(valid_mask)
    filter!(i -> i != 1, check_indices)

    if !isempty(check_indices)
        @test all(density_std[check_indices] .> density_large[check_indices])
    end
end

@testset "Ligtenberg with different coefficients" begin
    n = 10
    t_vec = 260.0 * ones(n)
    dz = 0.5 * ones(n)
    density = collect(range(300.0, 700.0, length=n))
    grain_radius = 0.5 * ones(n)
    cfs = _make_density_cfs()

    density_before = copy(density)

    # Test Case 1: Standard RACMO
    mp1 = GEMB.ModelParameters(
        density_ice=917.0,
        densification_method=:Ligtenberg,
        densification_coeffs_M01=:Ant_RACMO_GS_SW0,
    )
    (_, d1) = GEMB.calculate_density(t_vec, copy(dz), copy(density), grain_radius, cfs, mp1)

    # Test Case 2: ERA5 variant (different M coefficients)
    mp2 = GEMB.ModelParameters(
        density_ice=917.0,
        densification_method=:Ligtenberg,
        densification_coeffs_M01=:Ant_ERA5_BF_SW1,
    )
    (_, d2) = GEMB.calculate_density(t_vec, copy(dz), copy(density), grain_radius, cfs, mp2)

    # Ensure densification occurred
    @test all(d1 .> density_before)
    @test all(d2 .> density_before)

    # Different coefficients should produce different results
    @test d1 != d2
end

@testset "KuipersMunneke coefficients" begin
    n = 10
    t_vec = 260.0 * ones(n)
    dz = 0.5 * ones(n)
    density = collect(range(300.0, 700.0, length=n))
    grain_radius = 0.5 * ones(n)
    cfs = _make_density_cfs()

    density_before = copy(density)
    mp = GEMB.ModelParameters(
        density_ice=917.0,
        densification_method=:Ligtenberg,
        densification_coeffs_M01=:Gre_KuipersMunneke,
    )

    (_, density_out) = GEMB.calculate_density(t_vec, copy(dz), copy(density), grain_radius, cfs, mp)

    @test all(density_out .> density_before)
end

@testset "Zero time step (no change)" begin
    n = 10
    t_vec = 260.0 * ones(n)
    dz = 0.5 * ones(n)
    density = collect(range(300.0, 700.0, length=n))
    grain_radius = 0.5 * ones(n)

    cfs = _make_density_cfs(dt=0.0)

    mp = GEMB.ModelParameters(
        density_ice=917.0,
        densification_method=:HerronLangway,
        densification_coeffs_M01=:Ant_RACMO_GS_SW0,
    )

    density_before = copy(density)
    dz_before = copy(dz)
    (dz_out, density_out) = GEMB.calculate_density(t_vec, copy(dz), copy(density), grain_radius, cfs, mp)

    @test density_out ≈ density_before
    @test dz_out ≈ dz_before
end

@testset "Ligtenberg bare ice logic (density_ice = 820 vs 917)" begin
    n = 10
    t_vec = 260.0 * ones(n)
    dz = 0.5 * ones(n)
    density = collect(range(300.0, 700.0, length=n))
    grain_radius = 0.5 * ones(n)
    cfs = _make_density_cfs()

    # Case A: density_ice ~ 820 (specialized branch)
    mp_820 = GEMB.ModelParameters(
        density_ice=820.0,
        densification_method=:Ligtenberg,
        densification_coeffs_M01=:Gre_RACMO_GS_SW0,
    )
    (_, density_820) = GEMB.calculate_density(t_vec, copy(dz), copy(density), grain_radius, cfs, mp_820)

    # Case B: density_ice ~ 917 (standard branch)
    mp_917 = GEMB.ModelParameters(
        density_ice=917.0,
        densification_method=:Ligtenberg,
        densification_coeffs_M01=:Gre_RACMO_GS_SW0,
    )
    (_, density_917) = GEMB.calculate_density(t_vec, copy(dz), copy(density), grain_radius, cfs, mp_917)

    # Different density_ice branches should produce different results
    # Check only values well below 820
    check_idx = density .< 700
    @test density_820[check_idx] != density_917[check_idx]
end
