# Tests for calculate_albedo
# Translated from MATLAB test_calculate_albedo.m

@testset "Method None" begin
    n = 10
    temperature = fill(260.0, n)
    dz = fill(0.1, n)
    density = fill(400.0, n)
    water = zeros(n)
    grain_radius = fill(0.5, n)
    albedo_in = fill(0.8, n)
    albedo_diffuse_in = fill(0.8, n)
    evaporation_condensation = 0.0
    melt_surface = 0.0

    cfs = GEMB.ClimateForcingStep(
        3600.0, 260.0, 100000.0, 0.0, 5.0, 0.0, 0.0,
        100.0, 260.0, 5.0, 0.0, 2.0, 2.0,
        0.1, 0.1, 1.0, 60.0, 0.0, 0.5
    )

    mp = GEMB.ModelParameters(
        albedo_method="None",
        albedo_fixed=0.75,
        albedo_density_threshold=1023.0,
        density_ice=917.0,
        albedo_snow=0.85,
        albedo_ice=0.48,
        albedo_wet_snow_t0=15.0,
        albedo_dry_snow_t0=30.0,
        albedo_K=7.0
    )

    a_out, _ = GEMB.calculate_albedo(
        temperature, dz, density, water, grain_radius,
        albedo_in, albedo_diffuse_in, evaporation_condensation, melt_surface, cfs, mp)

    @test a_out[1] ≈ 0.75 atol=1e-10
end

@testset "Density threshold override" begin
    n = 10
    temperature = fill(260.0, n)
    dz = fill(0.1, n)
    density = fill(400.0, n)
    density[1] = 350.0  # Higher than threshold
    water = zeros(n)
    grain_radius = fill(0.5, n)
    albedo_in = fill(0.8, n)
    albedo_diffuse_in = fill(0.8, n)
    evaporation_condensation = 0.0
    melt_surface = 0.0

    cfs = GEMB.ClimateForcingStep(
        3600.0, 260.0, 100000.0, 0.0, 5.0, 0.0, 0.0,
        100.0, 260.0, 5.0, 0.0, 2.0, 2.0,
        0.1, 0.1, 1.0, 60.0, 0.0, 0.5
    )

    mp = GEMB.ModelParameters(
        albedo_method="GardnerSharp",
        albedo_fixed=0.4,
        albedo_density_threshold=300.0,  # Low threshold
        density_ice=917.0,
        albedo_snow=0.85,
        albedo_ice=0.48,
        albedo_wet_snow_t0=15.0,
        albedo_dry_snow_t0=30.0,
        albedo_K=7.0
    )

    a_out, _ = GEMB.calculate_albedo(
        temperature, dz, density, water, grain_radius,
        albedo_in, albedo_diffuse_in, evaporation_condensation, melt_surface, cfs, mp)

    @test a_out[1] ≈ 0.4 atol=1e-10
end

@testset "GardnerSharp valid range" begin
    n = 10
    temperature = fill(260.0, n)
    dz = fill(0.1, n)
    density = fill(400.0, n)
    water = zeros(n)
    grain_radius = fill(0.5, n)
    albedo_in = fill(0.8, n)
    albedo_diffuse_in = fill(0.8, n)
    evaporation_condensation = 0.0
    melt_surface = 0.0

    cfs = GEMB.ClimateForcingStep(
        3600.0, 260.0, 100000.0, 0.0, 5.0, 0.0, 0.0,
        100.0, 260.0, 5.0, 0.0, 2.0, 2.0,
        0.1, 0.1, 1.0, 60.0, 0.0, 0.5
    )

    mp = GEMB.ModelParameters(
        albedo_method="GardnerSharp",
        albedo_fixed=0.8,
        albedo_density_threshold=Inf,
        density_ice=917.0,
        albedo_snow=0.85,
        albedo_ice=0.48,
        albedo_wet_snow_t0=15.0,
        albedo_dry_snow_t0=30.0,
        albedo_K=7.0
    )

    a_out, a_diff_out = GEMB.calculate_albedo(
        temperature, dz, density, water, grain_radius,
        albedo_in, albedo_diffuse_in, evaporation_condensation, melt_surface, cfs, mp)

    @test 0.0 <= a_out[1] <= 1.0
    @test 0.0 <= a_diff_out[1] <= 1.0
end

@testset "GreuellKonzelmann" begin
    n = 10
    temperature = fill(260.0, n)
    dz = fill(0.1, n)
    density = fill(400.0, n)
    density[1] = 450.0
    water = zeros(n)
    grain_radius = fill(0.5, n)
    albedo_in = fill(0.8, n)
    albedo_diffuse_in = fill(0.8, n)
    evaporation_condensation = 0.0
    melt_surface = 0.0

    cfs = GEMB.ClimateForcingStep(
        3600.0, 260.0, 100000.0, 0.0, 5.0, 0.0, 0.0,
        100.0, 260.0, 5.0, 0.0, 2.0, 2.0,
        0.1, 0.1, 1.0, 60.0, 0.0, 0.5
    )

    mp = GEMB.ModelParameters(
        albedo_method="GreuellKonzelmann",
        albedo_fixed=0.8,
        albedo_density_threshold=Inf,
        density_ice=900.0,
        albedo_snow=0.85,
        albedo_ice=0.5,
        albedo_wet_snow_t0=15.0,
        albedo_dry_snow_t0=30.0,
        albedo_K=7.0
    )

    a_out, _ = GEMB.calculate_albedo(
        temperature, dz, density, water, grain_radius,
        albedo_in, albedo_diffuse_in, evaporation_condensation, melt_surface, cfs, mp)

    # Expected: a_ice + (density - density_ice) * (a_snow - a_ice) / (300 - density_ice) + 0.05*(cloud - 0.5)
    expected = 0.5 + (450 - 900) * (0.85 - 0.5) / (300 - 900)
    @test a_out[1] ≈ expected atol=1e-5
end

@testset "Bougamont2005 decay" begin
    n = 10
    temperature = fill(260.0, n)
    temperature[1] = 250.0  # Cold snow (-23 C)
    dz = fill(0.1, n)
    density = fill(400.0, n)
    water = zeros(n)
    grain_radius = fill(0.5, n)
    albedo_in = fill(0.8, n)
    albedo_diffuse_in = fill(0.8, n)
    evaporation_condensation = 0.0
    melt_surface = 0.0

    cfs = GEMB.ClimateForcingStep(
        86400.0, 260.0, 100000.0, 0.0, 5.0, 0.0, 0.0,  # dt = 1 day
        100.0, 260.0, 5.0, 0.0, 2.0, 2.0,
        0.1, 0.1, 1.0, 60.0, 0.0, 0.5
    )

    mp = GEMB.ModelParameters(
        albedo_method="Bougamont2005",
        albedo_fixed=0.8,
        albedo_density_threshold=Inf,
        density_ice=917.0,
        albedo_snow=0.85,
        albedo_ice=0.48,
        albedo_wet_snow_t0=15.0,
        albedo_dry_snow_t0=30.0,
        albedo_K=7.0
    )

    a_out, _ = GEMB.calculate_albedo(
        temperature, dz, density, water, grain_radius,
        albedo_in, albedo_diffuse_in, evaporation_condensation, melt_surface, cfs, mp)

    # t0 = 10 * K + dry_snow_t0 = 10*7 + 30 = 100 days
    t0 = 100.0
    da = (0.8 - 0.48) / t0 * 1.0  # dt = 1 day
    expected = 0.8 - da
    @test a_out[1] ≈ expected atol=1e-5
end
