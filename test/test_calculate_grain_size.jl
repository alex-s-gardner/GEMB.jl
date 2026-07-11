# Tests for calculate_grain_size - translated from MATLAB test_calculate_grain_size.m

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

@testset "Albedo method skip" begin
    n = 5
    temperature = 260.0 * ones(n)
    dz = 0.1 * ones(n)
    density = 300.0 * ones(n)
    water = zeros(n)
    grain_radius = 0.5 * ones(n)
    grain_dendricity = 0.5 * ones(n)
    grain_sphericity = 0.5 * ones(n)
    cfs = _make_grain_cfs()

    # Method "None" should skip grain evolution
    mp_none = GEMB.ModelParameters(albedo_method="None")
    (gs_out, gdn_out, gsp_out) = GEMB.calculate_grain_size(
        temperature, dz, density, water, grain_radius,
        grain_dendricity, grain_sphericity, cfs, mp_none)

    @test gs_out == grain_radius
    @test gdn_out == grain_dendricity
    @test gsp_out == grain_sphericity

    # Method "GreuellKonzelmann" should also skip
    mp_gk = GEMB.ModelParameters(albedo_method="GreuellKonzelmann")
    (gs_out2, _, _) = GEMB.calculate_grain_size(
        temperature, dz, density, water, grain_radius,
        grain_dendricity, grain_sphericity, cfs, mp_gk)

    @test gs_out2 == grain_radius
end

@testset "Dendritic dry low gradient" begin
    mp = GEMB.ModelParameters(albedo_method="GardnerSharp")
    cfs = _make_grain_cfs(dt=86400.0)

    temperature = [260.0, 260.0, 260.0]
    dz = [0.1, 0.1, 0.1]
    density = [200.0, 200.0, 200.0]
    water = [0.0, 0.0, 0.0]
    grain_radius = [0.2, 0.2, 0.2]
    grain_dendricity = [0.8, 0.8, 0.8]
    grain_sphericity = [0.2, 0.2, 0.2]

    (gs_out, gdn_out, gsp_out) = GEMB.calculate_grain_size(
        temperature, dz, density, water, grain_radius,
        grain_dendricity, grain_sphericity, cfs, mp)

    # Dendricity should decrease (decay)
    @test all(gdn_out .< grain_dendricity)
    # Sphericity should increase
    @test all(gsp_out .> grain_sphericity)
    # Grain radius should change
    @test gs_out != grain_radius
end

@testset "Dendritic dry high gradient" begin
    mp = GEMB.ModelParameters(albedo_method="GardnerSharp")
    cfs = _make_grain_cfs(dt=86400.0)

    # High temperature gradient (> 5 K/m)
    temperature = [260.0, 270.0, 280.0]
    dz = [0.1, 0.1, 0.1]
    density = [200.0, 200.0, 200.0]
    water = [0.0, 0.0, 0.0]
    grain_radius = [0.2, 0.2, 0.2]
    grain_dendricity = [0.8, 0.8, 0.8]
    grain_sphericity = [0.2, 0.2, 0.2]

    (_, gdn_out, gsp_out) = GEMB.calculate_grain_size(
        temperature, dz, density, water, grain_radius,
        grain_dendricity, grain_sphericity, cfs, mp)

    # Under high gradient: dendricity decreases, sphericity decreases
    @test all(gdn_out .< grain_dendricity)
    @test all(gsp_out .< grain_sphericity)
end

@testset "Dendritic wet snow" begin
    mp = GEMB.ModelParameters(albedo_method="GardnerSharp")
    cfs = _make_grain_cfs(dt=86400.0)

    temperature = [GEMB.CtoK, GEMB.CtoK, GEMB.CtoK]
    dz = [0.1, 0.1, 0.1]
    density = [250.0, 250.0, 250.0]
    water = [1.0, 1.0, 1.0]
    grain_radius = [0.2, 0.2, 0.2]
    grain_dendricity = [0.8, 0.8, 0.8]
    grain_sphericity = [0.2, 0.2, 0.2]

    (_, gdn_out, gsp_out) = GEMB.calculate_grain_size(
        temperature, dz, density, water, grain_radius,
        grain_dendricity, grain_sphericity, cfs, mp)

    # Wet snow causes rapid rounding
    @test all(gdn_out .< grain_dendricity)
    @test all(gsp_out .> grain_sphericity)
end

@testset "Nondendritic dry (Marbouty)" begin
    mp = GEMB.ModelParameters(albedo_method="GardnerSharp")
    cfs = _make_grain_cfs(dt=86400.0)

    # Moderate gradient (~20 K/m)
    temperature = [250.0, 252.0, 254.0]
    dz = [0.1, 0.1, 0.1]
    density = [300.0, 300.0, 300.0]  # must be < 400 for growth
    water = [0.0, 0.0, 0.0]
    grain_radius = [0.5, 0.5, 0.5]
    grain_dendricity = [0.0, 0.0, 0.0]
    grain_sphericity = [0.5, 0.5, 0.5]

    (gs_out, gdn_out, _) = GEMB.calculate_grain_size(
        temperature, dz, density, water, grain_radius,
        grain_dendricity, grain_sphericity, cfs, mp)

    # Expect growth in grain size
    @test all(gs_out .> grain_radius)
    # Dendricity stays at 0
    @test gdn_out == grain_dendricity
end

@testset "Marbouty density limit" begin
    mp = GEMB.ModelParameters(albedo_method="GardnerSharp")
    cfs = _make_grain_cfs(dt=86400.0)

    temperature = [250.0, 252.0, 254.0]
    dz = [0.1, 0.1, 0.1]
    density = [450.0, 450.0, 450.0]  # > 400 threshold
    water = [0.0, 0.0, 0.0]
    grain_radius = [0.5, 0.5, 0.5]
    grain_dendricity = [0.0, 0.0, 0.0]
    grain_sphericity = [0.5, 0.5, 0.5]

    (gs_out, _, _) = GEMB.calculate_grain_size(
        temperature, dz, density, water, grain_radius,
        grain_dendricity, grain_sphericity, cfs, mp)

    # No growth expected above density threshold
    @test gs_out ≈ grain_radius atol = 1e-10
end

@testset "Nondendritic wet snow (Brun)" begin
    mp = GEMB.ModelParameters(albedo_method="GardnerSharp")
    cfs = _make_grain_cfs(dt=86400.0)

    temperature = [GEMB.CtoK, GEMB.CtoK, GEMB.CtoK]
    dz = [0.1, 0.1, 0.1]
    density = [300.0, 300.0, 300.0]
    water = [1.5, 1.5, 1.5]
    grain_radius = [0.5, 0.5, 0.5]
    grain_dendricity = [0.0, 0.0, 0.0]
    grain_sphericity = [0.5, 0.5, 0.5]

    (gs_out, gdn_out, _) = GEMB.calculate_grain_size(
        temperature, dz, density, water, grain_radius,
        grain_dendricity, grain_sphericity, cfs, mp)

    # Expect growth via wet snow mechanism
    @test all(gs_out .> grain_radius)
    # Dendricity stays at 0
    @test gdn_out == grain_dendricity
end

@testset "Clamping limits" begin
    mp = GEMB.ModelParameters(albedo_method="GardnerSharp")
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
    mp = GEMB.ModelParameters(albedo_method="GardnerSharp")
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
