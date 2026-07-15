# Tests for calculate_accumulation - matches MATLAB test_calculate_accumulation.m

# Helper to create a ClimateForcingStep for accumulation tests
function _make_accum_cfs(; precipitation=0.0, temperature_air=270.0, wind_speed=5.0,
    precipitation_mean=200.0, temperature_air_mean=260.0, wind_speed_mean=5.0)
    return GEMB.ClimateForcingStep(
        86400.0,          # dt [s]
        temperature_air,  # temperature_air
        80000.0,          # pressure_air
        precipitation,    # precipitation
        wind_speed,       # wind_speed
        100.0,            # shortwave_downward
        250.0,            # longwave_downward
        300.0,            # vapor_pressure
        temperature_air_mean,   # temperature_air_mean
        wind_speed_mean,        # wind_speed_mean
        precipitation_mean,     # precipitation_mean
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

@testset "No precipitation (no changes)" begin
    n = 5
    t_vec = 260.0 * ones(n)
    dz = 0.1 * ones(n)
    density = 400.0 * ones(n)
    water = zeros(n)
    grain_radius = 0.5 * ones(n)
    grain_dendricity = 0.5 * ones(n)
    grain_sphericity = 0.5 * ones(n)
    albedo_in = 0.7 * ones(n)
    albedo_diffuse_in = 0.7 * ones(n)

    mp = GEMB.ModelParameters(
        new_snow_method=Symbol("150kgm2"),
        column_dzmin=0.05,
        density_ice=917.0,
        albedo_snow=0.85,
        albedo_method=:GardnerSharp,
        rain_temperature_threshold=273.15,
    )
    cfs = _make_accum_cfs(precipitation=0.0, temperature_air=270.0)

    (t_out, dz_out, d_out, _, _, _, _, _, _, ra_out) = GEMB.calculate_accumulation(
        t_vec, dz, density, water, grain_radius,
        grain_dendricity, grain_sphericity, albedo_in, albedo_diffuse_in,
        cfs, mp, false)

    @test dz_out == dz
    @test d_out == density
    @test t_out == t_vec
    @test ra_out == 0.0
end

@testset "Large snow event (new layer)" begin
    n = 5
    t_vec = 260.0 * ones(n)
    dz = 0.1 * ones(n)
    density = 400.0 * ones(n)
    water = zeros(n)
    grain_radius = 0.5 * ones(n)
    grain_dendricity = 0.5 * ones(n)
    grain_sphericity = 0.5 * ones(n)
    albedo_in = 0.7 * ones(n)
    albedo_diffuse_in = 0.7 * ones(n)

    density_snow = 150.0

    mp = GEMB.ModelParameters(
        new_snow_method=Symbol("150kgm2"),
        column_dzmin=0.05,
        density_ice=917.0,
        albedo_snow=0.85,
        albedo_method=:GardnerSharp,
        rain_temperature_threshold=273.15,
    )
    cfs = _make_accum_cfs(precipitation=50.0, temperature_air=260.0)

    (t_out, dz_out, d_out, _, _, gdn_out, gsp_out, a_out, _, _) = GEMB.calculate_accumulation(
        t_vec, dz, density, water, grain_radius,
        grain_dendricity, grain_sphericity, albedo_in, albedo_diffuse_in,
        cfs, mp, false)

    expected_dz = 50.0 / density_snow

    # Large snow should add a layer
    @test length(dz_out) == n + 1

    # Top layer properties
    @test d_out[1] == density_snow
    @test dz_out[1] ≈ expected_dz atol = 1e-10
    @test t_out[1] == 260.0
    @test a_out[1] == 0.85

    # Default microstructure for new snow
    @test gdn_out[1] == 1.0
    @test gsp_out[1] == 0.5
end

@testset "Small snow event (merge with top layer)" begin
    n = 5
    t_vec = 260.0 * ones(n)
    dz = 0.1 * ones(n)
    density = 400.0 * ones(n)
    water = zeros(n)
    grain_radius = 0.5 * ones(n)
    grain_dendricity = 0.5 * ones(n)
    grain_sphericity = 0.5 * ones(n)
    albedo_in = 0.7 * ones(n)
    albedo_diffuse_in = 0.7 * ones(n)

    density_snow = 150.0

    mp = GEMB.ModelParameters(
        new_snow_method=Symbol("150kgm2"),
        column_dzmin=0.05,
        density_ice=917.0,
        albedo_snow=0.85,
        albedo_method=:GardnerSharp,
        rain_temperature_threshold=273.15,
    )
    cfs = _make_accum_cfs(precipitation=2.0, temperature_air=260.0)

    old_mass = density[1] * dz[1]
    old_dz1 = dz[1]
    old_t1 = t_vec[1]
    old_a1 = albedo_in[1]

    (t_out, dz_out, d_out, _, _, _, _, a_out, _, _) = GEMB.calculate_accumulation(
        t_vec, dz, density, water, grain_radius,
        grain_dendricity, grain_sphericity, albedo_in, albedo_diffuse_in,
        cfs, mp, false)

    # Small snow should merge, no new layer
    @test length(dz_out) == n

    # Mass conservation and mixing
    new_mass = old_mass + 2.0
    expected_dz = old_dz1 + 2.0 / density_snow
    expected_d = new_mass / expected_dz

    @test dz_out[1] ≈ expected_dz atol = 1e-10
    @test d_out[1] ≈ expected_d atol = 1e-10

    # Temperature weighting
    expected_t = (260.0 * 2.0 + old_t1 * old_mass) / new_mass
    @test t_out[1] ≈ expected_t atol = 1e-10

    # Albedo weighting
    expected_a = (0.85 * 2.0 + old_a1 * old_mass) / new_mass
    @test a_out[1] ≈ expected_a atol = 1e-10
end

@testset "Rain event" begin
    n = 5
    t_vec = 260.0 * ones(n)
    dz = 0.1 * ones(n)
    density = 400.0 * ones(n)
    water = zeros(n)
    grain_radius = 0.5 * ones(n)
    grain_dendricity = 0.5 * ones(n)
    grain_sphericity = 0.5 * ones(n)
    albedo_in = 0.7 * ones(n)
    albedo_diffuse_in = 0.7 * ones(n)

    mp = GEMB.ModelParameters(
        new_snow_method=Symbol("150kgm2"),
        column_dzmin=0.05,
        density_ice=917.0,
        albedo_snow=0.85,
        albedo_method=:GardnerSharp,
        rain_temperature_threshold=273.15,
    )
    # temperature_air > 273.15 -> Rain
    cfs = _make_accum_cfs(precipitation=10.0, temperature_air=275.0, wind_speed=0.0)

    # Constants from GEMB
    lf = GEMB.LF
    ci = GEMB.C_ICE

    old_mass = density[1] * dz[1]
    old_dz1 = dz[1]
    old_t1 = t_vec[1]

    (t_out, dz_out, d_out, _, _, _, _, _, _, ra_out) = GEMB.calculate_accumulation(
        t_vec, dz, density, water, grain_radius,
        grain_dendricity, grain_sphericity, albedo_in, albedo_diffuse_in,
        cfs, mp, false)

    # Rain output flag
    @test ra_out == 10.0

    # Mass update: thickness stays same, density increases
    new_mass = old_mass + 10.0
    @test d_out[1] ≈ new_mass / old_dz1 atol = 1e-10
    @test dz_out[1] ≈ old_dz1 atol = 1e-10

    # Temperature update includes latent heat logic
    # T(1) = (precipitation * (temperature_air + LF/CI) + T(1) * mInit(1)) / mass
    term_rain = 10.0 * (275.0 + lf / ci)
    term_snow = old_t1 * old_mass
    expected_t = (term_rain + term_snow) / new_mass

    @test t_out[1] ≈ expected_t atol = 1e-8
end

@testset "Rain density cap" begin
    # Huge rain event causing density to exceed ice density -> clamp and expand
    n = 5
    t_vec = 260.0 * ones(n)
    dz = 0.1 * ones(n)
    density = 400.0 * ones(n)
    density[1] = 900.0  # Near ice density
    dz[1] = 0.1         # Mass = 90
    water = zeros(n)
    grain_radius = 0.5 * ones(n)
    grain_dendricity = 0.5 * ones(n)
    grain_sphericity = 0.5 * ones(n)
    albedo_in = 0.7 * ones(n)
    albedo_diffuse_in = 0.7 * ones(n)

    mp = GEMB.ModelParameters(
        new_snow_method=Symbol("150kgm2"),
        column_dzmin=0.05,
        density_ice=917.0,
        albedo_snow=0.85,
        albedo_method=:GardnerSharp,
        rain_temperature_threshold=273.15,
    )
    # Huge rain
    cfs = _make_accum_cfs(precipitation=500.0, temperature_air=275.0, wind_speed=0.0)

    (_, dz_out, d_out, _, _, _, _, _, _, _) = GEMB.calculate_accumulation(
        t_vec, dz, density, water, grain_radius,
        grain_dendricity, grain_sphericity, albedo_in, albedo_diffuse_in,
        cfs, mp, false)

    # Density should be capped at ice density
    @test d_out[1] ≈ 917.0 atol = 1e-10

    # Thickness should expand to conserve mass
    total_mass = 90.0 + 500.0
    expected_dz = total_mass / 917.0
    @test dz_out[1] ≈ expected_dz atol = 1e-10
end

@testset "New snow density methods" begin
    n = 5
    t_vec = 260.0 * ones(n)
    dz = 0.1 * ones(n)
    density = 400.0 * ones(n)
    water = zeros(n)
    grain_radius = 0.5 * ones(n)
    grain_dendricity = 0.5 * ones(n)
    grain_sphericity = 0.5 * ones(n)
    albedo_in = 0.7 * ones(n)
    albedo_diffuse_in = 0.7 * ones(n)

    temperature_air_mean = 260.0
    precipitation_mean = 200.0
    wind_speed_mean = 5.0

    # Common forcing (large precip to ensure new layer)
    cfs = _make_accum_cfs(precipitation=50.0, temperature_air=250.0, wind_speed=5.0,
        precipitation_mean=precipitation_mean, temperature_air_mean=temperature_air_mean,
        wind_speed_mean=wind_speed_mean)

    # Method 1: "350kgm2"
    mp1 = GEMB.ModelParameters(
        new_snow_method=Symbol("350kgm2"),
        column_dzmin=0.05,
        density_ice=917.0,
        albedo_snow=0.85,
        albedo_method=:GardnerSharp,
        rain_temperature_threshold=273.15,
    )
    (_, _, d1, _, _, _, _, _, _, _) = GEMB.calculate_accumulation(
        t_vec, dz, density, water, grain_radius,
        grain_dendricity, grain_sphericity, albedo_in, albedo_diffuse_in,
        cfs, mp1, false)
    @test d1[1] == 350.0

    # Method 2: "Fausto"
    mp2 = GEMB.ModelParameters(
        new_snow_method=:Fausto,
        column_dzmin=0.05,
        density_ice=917.0,
        albedo_snow=0.85,
        albedo_method=:GardnerSharp,
        rain_temperature_threshold=273.15,
    )
    (_, _, d2, _, _, _, _, _, _, _) = GEMB.calculate_accumulation(
        t_vec, dz, density, water, grain_radius,
        grain_dendricity, grain_sphericity, albedo_in, albedo_diffuse_in,
        cfs, mp2, false)
    @test d2[1] == 315.0

    # Method 3: "Kaspers"
    mp3 = GEMB.ModelParameters(
        new_snow_method=:Kaspers,
        column_dzmin=0.05,
        density_ice=917.0,
        albedo_snow=0.85,
        albedo_method=:GardnerSharp,
        rain_temperature_threshold=273.15,
    )
    expected_3 = (7.36e-2 + 1.06e-3 * temperature_air_mean + 6.69e-2 * precipitation_mean / 1000.0 + 4.77e-3 * wind_speed_mean) * 1000.0
    (_, _, d3, _, _, _, _, _, _, _) = GEMB.calculate_accumulation(
        t_vec, dz, density, water, grain_radius,
        grain_dendricity, grain_sphericity, albedo_in, albedo_diffuse_in,
        cfs, mp3, false)
    @test d3[1] ≈ expected_3 atol = 1e-5

    # Method 4: "KuipersMunneke"
    mp4 = GEMB.ModelParameters(
        new_snow_method=:KuipersMunneke,
        column_dzmin=0.05,
        density_ice=917.0,
        albedo_snow=0.85,
        albedo_method=:GardnerSharp,
        rain_temperature_threshold=273.15,
    )
    expected_4 = 481.0 + 4.834 * (temperature_air_mean - 273.15)
    (_, _, d4, _, _, _, _, _, _, _) = GEMB.calculate_accumulation(
        t_vec, dz, density, water, grain_radius,
        grain_dendricity, grain_sphericity, albedo_in, albedo_diffuse_in,
        cfs, mp4, false)
    @test d4[1] ≈ expected_4 atol = 1e-5
end
