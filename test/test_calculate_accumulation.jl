# Tests for calculate_accumulation

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
    age = zeros(n)

    mp = GEMB.ModelParameters(
        new_snow_method=:Constant150,
        column_dzmin=0.05,
        density_ice=917.0,
        albedo_snow=0.85,
        albedo_method=:GardnerSharp,
        rain_temperature_threshold=273.15,
    )
    cfs = _make_accum_cfs(precipitation=0.0, temperature_air=270.0)

    (t_out, dz_out, d_out, _, _, _, _, _, ra_out) = GEMB.calculate_accumulation(
        t_vec, dz, density, water, grain_radius,
        grain_dendricity, grain_sphericity, age,
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
    age = zeros(n)

    density_snow = 150.0

    mp = GEMB.ModelParameters(
        new_snow_method=:Constant150,
        column_dzmin=0.05,
        density_ice=917.0,
        albedo_snow=0.85,
        albedo_method=:GardnerSharp,
        rain_temperature_threshold=273.15,
    )
    cfs = _make_accum_cfs(precipitation=50.0, temperature_air=260.0)

    (t_out, dz_out, d_out, _, _, gdn_out, gsp_out, _, _) = GEMB.calculate_accumulation(
        t_vec, dz, density, water, grain_radius,
        grain_dendricity, grain_sphericity, age,
        cfs, mp, false)

    expected_dz = 50.0 / density_snow

    # Large snow should add a layer
    @test length(dz_out) == n + 1

    # Top layer properties
    @test d_out[1] == density_snow
    @test dz_out[1] ≈ expected_dz atol = 1e-10
    @test t_out[1] == 260.0

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
    age = zeros(n)

    density_snow = 150.0

    mp = GEMB.ModelParameters(
        new_snow_method=:Constant150,
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

    (t_out, dz_out, d_out, _, _, _, _, _, _) = GEMB.calculate_accumulation(
        t_vec, dz, density, water, grain_radius,
        grain_dendricity, grain_sphericity, age,
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
    age = zeros(n)

    mp = GEMB.ModelParameters(
        new_snow_method=:Constant150,
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
    ci = mp.heat_capacity_ice

    old_mass = density[1] * dz[1]
    old_dz1 = dz[1]
    old_t1 = t_vec[1]

    # `calculate_accumulation` mutates its column arguments in place, so the pristine state
    # is captured here for the second, `rain_heat_capacity=:ice` call below.
    pristine = map(copy, (t_vec, dz, density, water, grain_radius,
                          grain_dendricity, grain_sphericity, age))

    (t_out, dz_out, d_out, _, _, _, _, _, ra_out) = GEMB.calculate_accumulation(
        t_vec, dz, density, water, grain_radius,
        grain_dendricity, grain_sphericity, age,
        cfs, mp, false)

    # Rain output flag
    @test ra_out == 10.0

    # Mass update: thickness stays same, density increases
    new_mass = old_mass + 10.0
    @test d_out[1] ≈ new_mass / old_dz1 atol = 1e-10
    @test dz_out[1] ≈ old_dz1 atol = 1e-10

    # Temperature update includes latent heat logic. Rain arrives at 275 K, so it carries
    # sensible heat above the melting point as well as `LF`; under the default
    # `rain_heat_capacity = :water` that sensible part is at the *water* heat capacity, not
    # the ice one. Expressed as the equivalent temperature the mixing algebra sees:
    # `T_liq = CtoK + LF/c_ice + c_water*(T_air - CtoK)/c_ice`.
    cw = GEMB.HEAT_CAPACITY_WATER
    term_rain = 10.0 * (GEMB.CtoK + lf / ci + cw * (275.0 - GEMB.CtoK) / ci)
    term_snow = old_t1 * old_mass
    expected_t = (term_rain + term_snow) / new_mass

    @test t_out[1] ≈ expected_t atol = 1e-8

    # `rain_heat_capacity = :ice` reproduces the pre-Fix-5 convention exactly, so the old
    # reference is retained as the check that the gate is a faithful selector.
    mp_ice = GEMB.ModelParameters(
        new_snow_method=:Constant150, column_dzmin=0.05, density_ice=917.0,
        albedo_snow=0.85, albedo_method=:GardnerSharp,
        rain_temperature_threshold=273.15, rain_heat_capacity=:ice,
    )
    t_ice = GEMB.calculate_accumulation(pristine..., cfs, mp_ice, false)[1]
    @test t_ice[1] ≈ (10.0 * (275.0 + lf / ci) + term_snow) / new_mass atol = 1e-8
    # Water's larger heat capacity means more energy arrives with the same rain, so the
    # default warms the surface cell strictly more.
    @test t_out[1] > t_ice[1]
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
    age = zeros(n)

    mp = GEMB.ModelParameters(
        new_snow_method=:Constant150,
        column_dzmin=0.05,
        density_ice=917.0,
        albedo_snow=0.85,
        albedo_method=:GardnerSharp,
        rain_temperature_threshold=273.15,
    )
    # Huge rain
    cfs = _make_accum_cfs(precipitation=500.0, temperature_air=275.0, wind_speed=0.0)

    (_, dz_out, d_out, _, _, _, _, _, _) = GEMB.calculate_accumulation(
        t_vec, dz, density, water, grain_radius,
        grain_dendricity, grain_sphericity, age,
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
    age = zeros(n)

    temperature_air_mean = 260.0
    precipitation_mean = 200.0
    wind_speed_mean = 5.0

    # Common forcing (large precip to ensure new layer)
    cfs = _make_accum_cfs(precipitation=50.0, temperature_air=250.0, wind_speed=5.0,
        precipitation_mean=precipitation_mean, temperature_air_mean=temperature_air_mean,
        wind_speed_mean=wind_speed_mean)

    # Method 1: :Constant350
    mp1 = GEMB.ModelParameters(
        new_snow_method=:Constant350,
        column_dzmin=0.05,
        density_ice=917.0,
        albedo_snow=0.85,
        albedo_method=:GardnerSharp,
        rain_temperature_threshold=273.15,
    )
    (_, _, d1, _, _, _, _, _, _) = GEMB.calculate_accumulation(
        t_vec, dz, density, water, grain_radius,
        grain_dendricity, grain_sphericity, age,
        cfs, mp1, false)
    @test d1[1] == 350.0

    # :Constant315 — the bare constant, numerically equal to :Fausto but without the
    # Crocus wind-dependent fresh-grain properties that :Fausto also selects.
    mp1b = GEMB.ModelParameters(
        new_snow_method=:Constant315,
        column_dzmin=0.05,
        density_ice=917.0,
        albedo_snow=0.85,
        albedo_method=:GardnerSharp,
        rain_temperature_threshold=273.15,
    )
    (_, _, d1b, _, re1b, gdn1b, gsp1b, _, _) = GEMB.calculate_accumulation(
        t_vec, dz, density, water, grain_radius,
        grain_dendricity, grain_sphericity, age,
        cfs, mp1b, false)
    @test d1b[1] == 315.0
    # ...and the default grain properties, not the Crocus wind-dependent ones.
    @test re1b[1] == GEMB.RE_NEW_SNOW
    @test gdn1b[1] == GEMB.GDN_NEW_SNOW
    @test gsp1b[1] == GEMB.GSP_NEW_SNOW

    # Method 2: "Fausto"
    mp2 = GEMB.ModelParameters(
        new_snow_method=:Fausto,
        column_dzmin=0.05,
        density_ice=917.0,
        albedo_snow=0.85,
        albedo_method=:GardnerSharp,
        rain_temperature_threshold=273.15,
    )
    (_, _, d2, _, _, _, _, _, _) = GEMB.calculate_accumulation(
        t_vec, dz, density, water, grain_radius,
        grain_dendricity, grain_sphericity, age,
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
    (_, _, d3, _, _, _, _, _, _) = GEMB.calculate_accumulation(
        t_vec, dz, density, water, grain_radius,
        grain_dendricity, grain_sphericity, age,
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
    (_, _, d4, _, _, _, _, _, _) = GEMB.calculate_accumulation(
        t_vec, dz, density, water, grain_radius,
        grain_dendricity, grain_sphericity, age,
        cfs, mp4, false)
    @test d4[1] ≈ expected_4 atol = 1e-5
end

@testset "FaustoFit new snow density" begin
    _mp(method) = GEMB.ModelParameters(
        new_snow_method=method,
        column_dzmin=0.05,
        density_ice=917.0,
        albedo_snow=0.85,
        albedo_method=:GardnerSharp,
        rain_temperature_threshold=273.15,
    )
    mp_fit = _mp(:FaustoFit)
    mp_const = _mp(:Fausto)

    # The published Fausto et al. (2018) Greenland regression, written out independently of
    # the implementation. IMAU-FDM applies exactly this form (`initialise_model.f90`).
    rho0_fausto(T) = 362.1 + 2.78 * (T - 273.15)

    # `fresh_snow_density` reads the *fifth* argument (instantaneous air temperature) for
    # this method and ignores the three climatological means.
    for T in (233.15, 240.0, 256.2, 260.0, 273.15)
        @test GEMB.fresh_snow_density(mp_fit, 260.0, 200.0, 5.0, T) == rho0_fausto(T)
    end
    # With no fifth argument it falls back to the climatological mean, which is what
    # `steady_state_profile` relies on for a well-defined ρ₀.
    @test GEMB.fresh_snow_density(mp_fit, 248.0, 200.0, 5.0) == rho0_fausto(248.0)

    # `:Fausto`'s bare constant is this same fit evaluated at one fixed temperature. Recorded
    # as a test so the relationship between the two options cannot drift apart silently.
    T_equiv = 273.15 + (GEMB.DENSITY_NEW_SNOW_FAUSTO_CONSTANT - 362.1) / 2.78
    @test T_equiv ≈ 256.2 atol = 0.05
    @test GEMB.fresh_snow_density(mp_fit, 260.0, 200.0, 5.0, T_equiv) ≈
          GEMB.fresh_snow_density(mp_const, 260.0, 200.0, 5.0) rtol = 1e-14

    # The four older methods must be bit-identical with and without the new argument: none
    # of them reads it, so adding it cannot move an existing run.
    for method in (:Constant150, :Constant315, :Constant350, :Fausto, :Kaspers, :KuipersMunneke)
        mp = _mp(method)
        @test GEMB.fresh_snow_density(mp, 260.0, 200.0, 5.0) ===
              GEMB.fresh_snow_density(mp, 260.0, 200.0, 5.0, 240.0)
    end

    # Through `calculate_accumulation`: the new cell takes the instantaneous temperature's
    # density, not the climatological mean's (the two means differ here on purpose).
    function _fresh_cell(mp; temperature_air)
        n = 5
        cfs = _make_accum_cfs(precipitation=50.0, temperature_air=temperature_air,
            wind_speed=5.0, precipitation_mean=200.0, temperature_air_mean=260.0,
            wind_speed_mean=5.0)
        return GEMB.calculate_accumulation(260.0 * ones(n), 0.1 * ones(n), 400.0 * ones(n),
            zeros(n), 0.5 * ones(n), 0.5 * ones(n), 0.5 * ones(n), zeros(n), cfs, mp, false)
    end

    (_, dz_a, dens_a, _, re_a, gdn_a, gsp_a, _, _) = _fresh_cell(mp_fit; temperature_air=250.0)
    @test dens_a[1] ≈ rho0_fausto(250.0) rtol = 1e-14
    (_, _, dens_b, _, _, _, _, _, _) = _fresh_cell(mp_fit; temperature_air=270.0)
    @test dens_b[1] ≈ rho0_fausto(270.0) rtol = 1e-14
    @test dens_a[1] < dens_b[1]                        # colder snow is lighter
    @test dz_a[1] ≈ 50.0 / dens_a[1] rtol = 1e-14      # mass/density is consistent

    # `:FaustoFit` must take the Crocus wind-dependent grain branch that `:Fausto` takes —
    # the grain coupling is orthogonal to the density fit, and missing it would silently
    # change fresh-snow dendricity.
    (_, _, _, _, re_c, gdn_c, gsp_c, _, _) = _fresh_cell(mp_const; temperature_air=250.0)
    @test gdn_a[1] == gdn_c[1]
    @test gsp_a[1] == gsp_c[1]
    @test re_a[1] == re_c[1]
    # And those differ from the constants the other methods use, so the test above has teeth.
    (_, _, _, _, re_d, gdn_d, gsp_d, _, _) = _fresh_cell(_mp(:Constant350); temperature_air=250.0)
    @test gdn_d[1] == GEMB.GDN_NEW_SNOW
    @test gdn_a[1] != GEMB.GDN_NEW_SNOW

    # The regression is unbounded below (it reaches 0 at T ≈ 143 K), so the call site clamps.
    # Far outside any real forcing, but the clamp is what keeps `dz = mass/density` finite.
    (_, dz_e, dens_e, _, _, _, _, _, _) = _fresh_cell(mp_fit; temperature_air=100.0)
    @test rho0_fausto(100.0) < 0                        # the bare fit is unphysical here
    @test dens_e[1] == 1.0
    @test isfinite(dz_e[1]) && dz_e[1] > 0

    @test GEMB.validate_parameters(mp_fit) === nothing
end

@testset "Pahaut new snow density" begin
    _mp(method) = GEMB.ModelParameters(
        new_snow_method=method,
        column_dzmin=0.05,
        density_ice=917.0,
        albedo_snow=0.85,
        albedo_method=:GardnerSharp,
        rain_temperature_threshold=273.15,
    )
    mp_pahaut = _mp(:Pahaut)

    # Pahaut (1975) as Crocus implements it, via Lafaysse et al. (2026) eq. 35. Written out
    # independently of the implementation, in the paper's own form.
    rho0_pahaut(T, U) = max(50.0, 109.0 + 6.0 * (T - 273.15) + 26.0 * sqrt(U))

    # Reads both instantaneous arguments — the fifth (air temperature) and the sixth (wind
    # speed) — and ignores all three climatological means.
    for T in (243.15, 253.15, 263.15, 273.15), U in (0.0, 1.0, 5.0, 12.0)
        @test GEMB.fresh_snow_density(mp_pahaut, 260.0, 200.0, 3.0, T, U) == rho0_pahaut(T, U)
    end

    # Two published reference points, computed by hand from eq. 35.
    @test GEMB.fresh_snow_density(mp_pahaut, 260.0, 200.0, 3.0, 273.15, 0.0) ≈ 109.0 rtol = 1e-14
    # T = -10 C, U = 4 m/s: 109 - 60 + 52 = 101
    @test GEMB.fresh_snow_density(mp_pahaut, 260.0, 200.0, 3.0, 263.15, 4.0) ≈ 101.0 rtol = 1e-14

    # Monotone increasing in both arguments, which is the physical content of the fit:
    # warmer snowfall gives denser crystals, and wind fragments them.
    ρ_T = [GEMB.fresh_snow_density(mp_pahaut, 260.0, 200.0, 3.0, T, 8.0) for T in 258.0:2.0:273.15]
    @test all(diff(ρ_T) .> 0)
    ρ_U = [GEMB.fresh_snow_density(mp_pahaut, 260.0, 200.0, 3.0, 270.0, U) for U in 0.0:1.0:15.0]
    @test all(diff(ρ_U) .> 0)

    # The published 50 kg m-3 floor binds in cold, calm air and nowhere warmer.
    @test GEMB.fresh_snow_density(mp_pahaut, 260.0, 200.0, 3.0, 250.0, 0.0) == 50.0
    @test GEMB.fresh_snow_density(mp_pahaut, 260.0, 200.0, 3.0, 200.0, 0.0) == 50.0
    @test GEMB.fresh_snow_density(mp_pahaut, 260.0, 200.0, 3.0, 273.15, 0.0) > 50.0

    # With neither instantaneous argument it falls back to the climatological means, which is
    # what `steady_state_profile` relies on for a well-defined ρ₀.
    @test GEMB.fresh_snow_density(mp_pahaut, 265.0, 200.0, 4.0) == rho0_pahaut(265.0, 4.0)

    # Every method that does not read the sixth argument must be bit-identical with and
    # without it, so adding it cannot move an existing run.
    for method in (:Constant150, :Constant315, :Constant350, :Fausto, :FaustoFit,
                   :Kaspers, :KuipersMunneke)
        mp = _mp(method)
        @test GEMB.fresh_snow_density(mp, 260.0, 200.0, 5.0, 250.0) ===
              GEMB.fresh_snow_density(mp, 260.0, 200.0, 5.0, 250.0, 11.0)
    end

    # Much lighter than the polar fits over their common range: alpine snowfall is warmer,
    # wetter, and far less wind-packed than a katabatic-scoured ice-sheet surface. This is
    # the reason the option exists, so it is pinned rather than left implicit.
    ρ_pahaut = GEMB.fresh_snow_density(mp_pahaut, 260.0, 200.0, 5.0, 263.15, 5.0)
    ρ_fausto = GEMB.fresh_snow_density(_mp(:FaustoFit), 260.0, 200.0, 5.0, 263.15, 5.0)
    @test ρ_pahaut < ρ_fausto
    @test ρ_pahaut ≈ 107.14 atol = 0.01
    @test ρ_fausto ≈ 334.30 atol = 0.01

    # Through `calculate_accumulation`: the new cell takes the instantaneous forcing's
    # density, and `dz = mass/density` follows from it.
    function _fresh_cell(mp; temperature_air, wind_speed)
        n = 5
        cfs = _make_accum_cfs(precipitation=50.0, temperature_air=temperature_air,
            wind_speed=wind_speed, precipitation_mean=200.0, temperature_air_mean=260.0,
            wind_speed_mean=3.0)
        return GEMB.calculate_accumulation(260.0 * ones(n), 0.1 * ones(n), 400.0 * ones(n),
            zeros(n), 0.5 * ones(n), 0.5 * ones(n), 0.5 * ones(n), zeros(n), cfs, mp, false)
    end

    (_, dz_a, dens_a, _, re_a, gdn_a, gsp_a, _, _) =
        _fresh_cell(mp_pahaut; temperature_air=268.0, wind_speed=6.0)
    @test dens_a[1] ≈ rho0_pahaut(268.0, 6.0) rtol = 1e-14
    @test dz_a[1] ≈ 50.0 / dens_a[1] rtol = 1e-14

    # Wind, not just temperature, must reach the call site — the sixth argument is new, and
    # passing the climatological mean by mistake would be invisible without this.
    (_, _, dens_b, _, _, _, _, _, _) =
        _fresh_cell(mp_pahaut; temperature_air=268.0, wind_speed=12.0)
    @test dens_b[1] > dens_a[1]
    @test dens_b[1] ≈ rho0_pahaut(268.0, 12.0) rtol = 1e-14

    # `:Pahaut` takes the Crocus wind-dependent grain branch, being the fresh-snow density
    # Crocus itself pairs those properties with.
    (_, _, _, _, re_c, gdn_c, gsp_c, _, _) =
        _fresh_cell(_mp(:Fausto); temperature_air=268.0, wind_speed=6.0)
    @test gdn_a[1] == gdn_c[1]
    @test gsp_a[1] == gsp_c[1]
    @test re_a[1] == re_c[1]
    @test gdn_a[1] != GEMB.GDN_NEW_SNOW   # so the test above has teeth

    @test GEMB.validate_parameters(mp_pahaut) === nothing
end
