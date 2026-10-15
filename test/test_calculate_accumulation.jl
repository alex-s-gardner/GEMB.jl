# Tests for calculate_accumulation

# Helper to create a ClimateForcingStep for accumulation tests
function _make_accum_cfs(; precipitation=0.0, temperature_air=270.0, wind_speed=5.0,
    precipitation_mean=200.0, temperature_air_mean=260.0, wind_speed_mean=5.0, dt=86400.0)
    return GEMB.ClimateForcingStep(
        dt,               # dt [s]
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
    # The top cell has the fresh-snow grain state, so the snow is similar and merges.
    n = 5
    t_vec = 260.0 * ones(n)
    dz = 0.1 * ones(n)
    density = 400.0 * ones(n)
    water = zeros(n)
    grain_radius = GEMB.RE_NEW_SNOW * ones(n)
    grain_dendricity = GEMB.GDN_NEW_SNOW * ones(n)
    grain_sphericity = GEMB.GSP_NEW_SNOW * ones(n)
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
    # Rain joins water[1] at the melt-point enthalpy; its sensible heat above 0 °C warms the
    # matrix. Density and thickness do not change.
    n = 5
    t_vec = 260.0 * ones(n)
    dz = 0.1 * ones(n)
    density = 400.0 * ones(n)
    water = zeros(n)
    grain_radius = 0.5 * ones(n)
    grain_dendricity = 0.5 * ones(n)
    grain_sphericity = 0.5 * ones(n)
    age = zeros(n)
    pristine = map(copy, (t_vec, dz, density, water, grain_radius,
                          grain_dendricity, grain_sphericity, age))

    cfs = _make_accum_cfs(precipitation=10.0, temperature_air=275.0, wind_speed=0.0)
    # Matrix heat capacity, and its chord between two temperatures: 2102 under `:constant`;
    # Cuffey & Paterson (2010) eq. 9.1, c = 152.5 + 7.122·T, whose enthalpy change between
    # T1 and T2 is (T2 − T1) times c at the midpoint, under `:CuffeyPaterson`.
    c_matrix(hcm, T) = hcm === :constant ? 2102.0 : 152.5 + 7.122 * T
    c_chord(hcm, T1, T2) = c_matrix(hcm, (T1 + T2) / 2)
    for hcm in (:CuffeyPaterson, :constant), rain_heat_capacity in (:water, :ice)
        c_rain = rain_heat_capacity === :water ? GEMB.HEAT_CAPACITY_WATER :
                 c_matrix(hcm, GEMB.CtoK)
        mp = GEMB.ModelParameters(new_snow_method=:Constant150, column_dzmin=0.05,
            density_ice=917.0, rain_temperature_threshold=273.15,
            rain_heat_capacity=rain_heat_capacity, heat_capacity_method=hcm)
        (t_out, dz_out, d_out, w_out, _, _, _, _, ra_out) = GEMB.calculate_accumulation(
            map(copy, pristine)..., cfs, mp, true)
        @test ra_out == 10.0
        @test w_out[1] == 10.0
        @test d_out[1] == 400.0
        @test dz_out[1] == 0.1
        M1 = 400.0 * 0.1
        @test M1 * c_chord(hcm, 260.0, t_out[1]) * (t_out[1] - 260.0) ≈
              10.0 * c_rain * (275.0 - GEMB.CtoK) rtol = 1e-12
    end
end

@testset "Rain does not compact snow" begin
    # 5 kg of rain at 275.15 K on a 0.05 m, 350 kg m-3 cell at 0 °C, through accumulation
    # and melt. The cell keeps its density; its water gains the rain plus the ice its
    # sensible heat melts, 5·c_w·2/LF; nothing is refrozen; and the deep column below holds
    # all of the water, so nothing runs off.
    n = 20
    t_vec = fill(GEMB.CtoK, n)
    dz = fill(0.05, n)
    density = fill(350.0, n)
    water = zeros(n)
    gr, gd, gs, age = fill(0.5, n), fill(0.0, n), fill(0.5, n), zeros(n)
    mp = GEMB.ModelParameters(new_snow_method=:Constant150, density_ice=917.0,
        rain_temperature_threshold=273.15, melt_geometry=:thickness)
    cfs = _make_accum_cfs(precipitation=5.0, temperature_air=275.15, wind_speed=0.0)
    t1, dz1, ρ1, w1, gr1, gd1, gs1, a1, rain = GEMB.calculate_accumulation(t_vec, dz,
        density, water, gr, gd, gs, age, cfs, mp, true)
    @test ρ1[1] == 350.0
    @test dz1[1] == 0.05
    t2, dz2, ρ2, w2, _, _, _, _, melt, _, runoff, refreeze, _ = GEMB.calculate_melt(t1, dz1,
        ρ1, w1, gr1, gd1, gs1, a1, mp, true)
    melted = 5.0 * GEMB.HEAT_CAPACITY_WATER * 2.0 / GEMB.LF
    @test all(ρ2 .== 350.0)
    @test sum(w2) ≈ 5.0 + melted rtol = 1e-6
    @test melt ≈ melted rtol = 1e-9
    @test refreeze == 0.0
    @test runoff == 0.0
    @test sum(dz2 .* ρ2) ≈ n * 0.05 * 350.0 - melted rtol = 1e-12
end

@testset "Rain on cold snow refreezes and is counted" begin
    # 2 kg of rain at 274.15 K on 17.5 kg of snow at -10 °C. The rain's sensible heat warms
    # the cell, and the remaining cold content refreezes [c·M·10 − P·c_w·1]/LF of the rain,
    # with c the chord heat capacity between 263.15 and 273.15 K: under the default
    # `:CuffeyPaterson`, Cuffey & Paterson (2010) eq. 9.1 at the midpoint, 152.5 + 7.122·268.15.
    n = 10
    t_vec = fill(GEMB.CtoK - 10.0, n)
    dz = fill(0.05, n)
    density = fill(350.0, n)
    gr, gd, gs, age = fill(0.5, n), fill(0.0, n), fill(0.5, n), zeros(n)
    mp = GEMB.ModelParameters(new_snow_method=:Constant150, density_ice=917.0,
        rain_temperature_threshold=273.15)
    cfs = _make_accum_cfs(precipitation=2.0, temperature_air=274.15, wind_speed=0.0)
    out = GEMB.calculate_accumulation(t_vec, dz, density, zeros(n), gr, gd, gs, age, cfs, mp,
        true)
    M = 0.05 * 350.0
    w_before = sum(out[4])
    t2, dz2, ρ2, w2, _, _, _, _, melt, _, runoff, refreeze, _ =
        GEMB.calculate_melt(out[1:8]..., mp, true)
    c = 152.5 + 7.122 * (GEMB.CtoK - 5.0)
    expected = (c * M * 10.0 - 2.0 * GEMB.HEAT_CAPACITY_WATER * 1.0) / GEMB.LF
    @test refreeze ≈ expected rtol = 1e-9
    @test melt == 0.0
    @test t2[1] ≈ GEMB.CtoK atol = 1e-9
    # Liquid budget: melt + rain − refreeze − runoff = change in column water.
    @test melt + out[9] - refreeze - runoff ≈ sum(w2) - 0.0 atol = 1e-12
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


# Snow layering after Crocus (see `calculate_accumulation` and `snow_layer_distance`).
# `_brun_dendritic_radius` is the Brun et al. (1992) dendritic grain diameter, halved, as GEMB
# carries it; `_crocus_distance` is Lafaysse et al. (2026) eq. 7 with d_opt = 2r [m].
_brun_dendritic_radius(d, s) = 0.1 * (d / 0.99 + (1 - d / 0.99) * (3s / 0.99 + 4(1 - s / 0.99))) / 2
_crocus_distance(d1, s1, r1, w1, a1, d2, s2, r2, w2, a2) = abs(a1 - a2) < 300 ?
    50 * (abs(d1 - d2) + abs(s1 - s2) + 2000 * abs(2r1 / 1000 - 2r2 / 1000) + (w1 != w2)) : 200.0
# Mass-weighted r^0.07 (Gardner & Sharp 2010 eq. 7 conserved), in BigFloat.
_albedo_mean_radius(r1, M1, r2, M2) =
    Float64(((big(M1) * big(r1)^big"0.07" + big(M2) * big(r2)^big"0.07") / (big(M1) + big(M2)))^(1 / big"0.07"))

# A column of `n` cells: cell 1 as given, cells 2… 0.05 m of old snow (350 kg m-3, r 1 mm,
# sphericity 0.8, non-dendritic, 100 days old) over ice from cell 6.
function _layer_column(; dz1=0.05, ρ1=350.0, r1=1.0, d1=0.0, s1=0.8, age1=100.0, T1=260.0,
    w1=0.0, n=10, ρ_ice=917.0)
    ρ = [ρ1; fill(350.0, 4); fill(ρ_ice, n - 5)]
    return ([T1; fill(260.0, n - 1)], [dz1; fill(0.05, n - 1)], ρ, [w1; zeros(n - 1)],
        [r1; fill(1.0, 4); fill(GEMB.GRAIN_RADIUS_ICE, n - 5)], [d1; zeros(n - 1)],
        [s1; fill(0.8, n - 1)], [age1; fill(100.0, n - 1)])
end

@testset "Snow layer distance and distinctness (Crocus eq. 7, 11)" begin
    mp = GEMB.ModelParameters()
    for (d1, s1, r1, w1, a1, d2, s2, r2, w2, a2) in (
        (1.0, 0.5, 0.05, false, 0.0, 0.0, 0.8, 1.0, false, 100.0),
        (1.0, 0.5, 0.05, false, 0.0, 0.98, 0.55, 0.0611, true, 0.25),
        (0.0, 0.3, 0.2, true, 10.0, 0.0, 0.3, 0.2, true, 20.0),
        (1.0, 0.5, 0.05, false, 0.0, 0.0, 0.5, 0.06, false, 3.0),     # fine non-dendritic
        (1.0, 0.5, 0.05, false, 0.0, 0.0, 0.9, 0.4, false, 400.0))
        D = GEMB.snow_layer_distance(d1, s1, r1, w1, a1, 300.0, d2, s2, r2, w2, a2, 400.0)
        @test D ≈ _crocus_distance(d1, s1, r1, w1, a1, d2, s2, r2, w2, a2) rtol = 1e-14
        @test D == GEMB.snow_layer_distance(d2, s2, r2, w2, a2, 400.0, d1, s1, r1, w1, a1, 300.0)
    end
    @test GEMB.snow_layer_distance(1.0, 0.5, 0.05, false, 0.0, 300.0,
        1.0, 0.5, 0.05, false, 0.0, 300.0) == 0
    # Firn or ice on either side: never aggregated.
    for ρ in (GEMB.DENSITY_PORE_CLOSEOFF, 870.0, 917.0)
        @test GEMB.snow_layer_distance(1.0, 0.5, 0.05, false, 0.0, 300.0,
            1.0, 0.5, 0.05, false, 0.0, ρ) == Inf
        @test GEMB.snow_layers_distinct(1e-9, Inf, mp)
    end
    # z < 5·z*/D aggregates, z* = column_dztop.
    @test GEMB.snow_layers_distinct(0.01, 5 * mp.column_dztop / 0.01, mp)
    @test !GEMB.snow_layers_distinct(0.01, prevfloat(5 * mp.column_dztop / 0.01), mp)
    @test !GEMB.snow_layers_distinct(0.02, 0.0, mp)
end

@testset "Snowfall onto similar fresh snow merges by mass" begin
    # Top cell: 3-hour-old fresh snow, 0.02 m at 150 kg m-3 (3 kg m-2); 2 kg m-2 more at
    # 150 kg m-3 is 0.0133 m, under column_dzmin 0.05 m and not distinct from it.
    mp = GEMB.ModelParameters(new_snow_method=:Constant150, column_dzmin=0.05,
        heat_capacity_method=:constant)
    d1, s1 = 0.98, 0.52
    r1 = _brun_dendritic_radius(d1, s1)
    col = _layer_column(dz1=0.02, ρ1=150.0, r1=r1, d1=d1, s1=s1, age1=0.125, T1=262.0)
    before = deepcopy(col)
    P, Ta = 2.0, 255.0
    D = _crocus_distance(GEMB.GDN_NEW_SNOW, GEMB.GSP_NEW_SNOW, GEMB.RE_NEW_SNOW, false, 0.0,
        d1, s1, r1, false, 0.125)
    @test D * (P / 150.0) < 5 * mp.column_dztop
    out = GEMB.calculate_accumulation(col..., _make_accum_cfs(dt=10800.0, precipitation=P,
        temperature_air=Ta, wind_speed=0.0), mp, true)
    M = 0.02 * 150.0
    f = P / (M + P)
    @test length(out[2]) == 10
    @test out[2][1] ≈ 0.02 + P / 150.0 rtol = 1e-14
    @test out[3][1] ≈ (M + P) / (0.02 + P / 150.0) rtol = 1e-14
    @test out[1][1] ≈ (M * 262.0 + P * Ta) / (M + P) rtol = 1e-14       # c·T, :constant
    @test out[8][1] ≈ 0.125 * M / (M + P) rtol = 1e-14
    d = (M * d1 + P * GEMB.GDN_NEW_SNOW) / (M + P)
    s = (M * s1 + P * GEMB.GSP_NEW_SNOW) / (M + P)
    @test out[6][1] ≈ d rtol = 1e-14
    @test out[7][1] ≈ s rtol = 1e-14
    @test out[5][1] ≈ _brun_dendritic_radius(d, s) rtol = 1e-14
    for k in 1:8
        @test out[k][2:end] == before[k][2:end]                          # cells 2… untouched
    end
end

@testset "Snowfall onto dissimilar snow, firn or ice forms its own cell" begin
    # 2 kg m-2 at 350 kg m-3 is 5.7 mm: above the 3-hour minimum mass, below column_dzmin.
    mp = GEMB.ModelParameters(new_snow_method=:Constant350)
    P, Ta = 2.0, 258.0
    z = P / 350.0
    @test GEMB.snow_layer_min_mass(10800.0) <= P && z < mp.column_dzmin
    # Old coarse snow; fine-grained non-dendritic snow (dissimilar through dendricity alone);
    # firn; ice.
    for (ρ1, r1, s1) in ((350.0, 1.0, 0.8), (350.0, 0.06, 0.5), (870.0, 1.5, 0.9),
        (917.0, GEMB.GRAIN_RADIUS_ICE, 0.9))
        col = _layer_column(; ρ1, r1, s1)
        before = deepcopy(col)
        out = GEMB.calculate_accumulation(col..., _make_accum_cfs(dt=10800.0, precipitation=P,
            temperature_air=Ta, wind_speed=0.0), mp, true)
        @test length(out[2]) == 11
        @test (out[1][1], out[2][1], out[3][1], out[4][1]) == (Ta, z, 350.0, 0.0)
        @test (out[5][1], out[6][1], out[7][1], out[8][1]) ==
              (GEMB.RE_NEW_SNOW, GEMB.GDN_NEW_SNOW, GEMB.GSP_NEW_SNOW, 0.0)
        for k in 1:8
            @test out[k][2:end] == before[k]                                # old cell untouched
        end
        # Layer management keeps the thin cell (distinct, at least the minimum thickness) and
        # leaves the old cell as it was.
        lay = GEMB.manage_layer_thickness(out[1:8]..., mp, true; dt=10800.0, n_target=11)
        for k in 1:8
            @test lay[k][1] == out[k][1]
            @test lay[k][2] == before[k][1]
        end
    end
    # Below the minimum mass dissimilar snow merges by mass, into snow, firn or ice alike:
    # repeated light snowfall onto ice lowers its density until it is snow.
    P = 0.5
    @test P < GEMB.snow_layer_min_mass(10800.0)
    for (ρ1, r1) in ((350.0, 1.0), (870.0, 1.5), (917.0, 1.5))
        out = GEMB.calculate_accumulation(_layer_column(; ρ1, r1)...,
            _make_accum_cfs(dt=10800.0, precipitation=P, temperature_air=Ta, wind_speed=0.0), mp, true)
        M = 0.05 * ρ1
        @test length(out[2]) == 10
        @test out[2][1] ≈ 0.05 + P / 350.0 rtol = 1e-14
        @test out[3][1] ≈ (M + P) / (0.05 + P / 350.0) rtol = 1e-14
        @test out[6][1] == 0.0
        @test out[7][1] ≈ (M * 0.8 + P * GEMB.GSP_NEW_SNOW) / (M + P) rtol = 1e-14
        @test out[5][1] ≈ _albedo_mean_radius(r1, M, GEMB.RE_NEW_SNOW, P) rtol = 1e-13
    end
end

@testset "Repeated small snowfalls onto an old surface accumulate into a fresh layer" begin
    # Four 3-hourly snowfalls of 2 kg m-2 (5.7 mm each at 350 kg m-3) onto old non-dendritic
    # snow of 1 mm grains, with a grain-growth step and layer management between them. The
    # first forms a fresh layer, the others join it; the old cell is never touched.
    mp = GEMB.ModelParameters(new_snow_method=:Constant350)
    col = _layer_column()
    n = length(col[2])
    cfs_dry = GEMB.ClimateForcingStep(; dt=10800.0, temperature_air=258.0)
    for k in 1:4
        # The old cell as it stands before this snowfall (the grain steps age it).
        old_top = [c[k == 1 ? 1 : 2] for c in col]
        out = GEMB.calculate_accumulation(col..., _make_accum_cfs(dt=10800.0, precipitation=2.0,
            temperature_air=258.0, wind_speed=0.0), mp, true)
        out = GEMB.manage_layer_thickness(out[1:8]..., mp, true; dt=10800.0, n_target=n + 1)
        col = collect(out[1:8])
        @test length(col[2]) == n + 1
        @test [c[2] for c in col] == old_top
        r, d, s = GEMB.calculate_grain_size(col[1], col[2], col[3], col[4], col[5], col[6],
            col[7], cfs_dry, mp)
        col[5], col[6], col[7] = r, d, s
    end
    @test col[2][1] * col[3][1] ≈ 8.0 rtol = 1e-13
    @test col[6][1] > 0.9                        # still dendritic fresh snow
    @test col[5][1] < 0.06
    # The albedo is that of the fresh layer, not of the old surface.
    cfs = GEMB.ClimateForcingStep(; dt=10800.0, solar_zenith_angle=50.0)
    @test GEMB.calculate_albedo(col[2], col[3], col[5], cfs, mp)[1] ==
          GEMB._albedo_gardner(col[5][1], 0.0, 50.0, 0.0)
end

@testset "Trace and zero snowfall leave the surface albedo unchanged" begin
    # Old snow, a thin fresh layer over old snow, a thin snow layer over ice, firn and ice.
    mp = GEMB.ModelParameters(albedo_ice=0.35)
    cfs_alb = GEMB.ClimateForcingStep(; dt=3600.0, solar_zenith_angle=60.0,
        cloud_optical_thickness=2.0)
    columns = (_layer_column(), _layer_column(dz1=0.008, ρ1=300.0, r1=0.05, d1=1.0, s1=0.5),
        (c = _layer_column(dz1=0.006, ρ1=300.0, r1=0.3, s1=0.6); c[3][2:end] .= 917.0; c),
        _layer_column(ρ1=870.0, r1=1.5), _layer_column(ρ1=917.0, r1=GEMB.GRAIN_RADIUS_ICE))
    for col in columns, P in (5e-324, 0.0)
        before = deepcopy(col)
        α0 = GEMB.calculate_albedo(col[2], col[3], col[5], cfs_alb, mp)
        out = GEMB.calculate_accumulation(deepcopy(col)..., _make_accum_cfs(dt=10800.0, precipitation=P,
            temperature_air=250.0, wind_speed=0.0), mp, true)
        @test length(out[2]) == 10
        α1 = GEMB.calculate_albedo(out[2], out[3], out[5], cfs_alb, mp)
        @test abs(α1[1] - α0[1]) <= 1e-12
        @test abs(α1[2] - α0[2]) <= 1e-12
        # Every field: grain state bit for bit, the rest to rounding.
        @test (out[5], out[6], out[7]) == (before[5], before[6], before[7])
        for k in (1, 2, 3, 8)
            @test all(isapprox.(out[k], before[k]; rtol=1e-15))
        end
        P == 0 && @test all(out[k] == before[k] for k in 1:8)
    end
end

@testset "Snowfall conserves mass and energy on every path" begin
    # Column mass and enthalpy gain exactly P and P·h(T_air), with h written from Cuffey &
    # Paterson (2010) eq. 9.1 integrated (c = 152.5 + 7.122 T), or c·T under :constant. P
    # spans a merge below the minimum thickness, a thin new cell and a full new cell; the top
    # cell is old snow, fresh snow, or ice.
    for method in (:CuffeyPaterson, :constant), P in (0.01, 1.0, 2.0, 8.0, 20.0),
        top in (:old, :fresh, :ice)
        mp = GEMB.ModelParameters(new_snow_method=:Constant350, heat_capacity_method=method)
        h(T) = method === :constant ? mp.heat_capacity_ice * T : 152.5 * T + 7.122 * T^2 / 2
        col = top === :old ? _layer_column(dz1=0.04) :
              top === :fresh ? _layer_column(dz1=0.04, r1=0.05, d1=1.0, s1=0.5, age1=0.1) :
              _layer_column(dz1=0.04, ρ1=917.0, r1=GEMB.GRAIN_RADIUS_ICE)
        col[1][1:5] .= [262.0, 261.0, 259.0, 258.0, 257.0]
        mass(c) = sum(c[2] .* c[3])
        enth(c) = sum(c[2] .* c[3] .* h.(c[1]))
        M0, E0 = mass(col), enth(col)
        out = GEMB.calculate_accumulation(deepcopy(col)..., _make_accum_cfs(dt=10800.0, precipitation=P,
            temperature_air=248.0, wind_speed=0.0), mp, true)
        @test mass(out) ≈ M0 + P rtol = 1e-14
        @test enth(out) ≈ E0 + P * h(248.0) rtol = 1e-13
    end
end

@testset "Snow layering closes the budgets over a multi-event run" begin
    # Twenty hourly days over bare ice: snowfalls of a trace, sub-minimum, thin-layer and
    # full-cell size, separated by warm sunny spells that melt the snow back to ice. `verbose`
    # makes `gemb_core` check mass and energy every step and `gemb` the whole-run mass budget;
    # a violation throws. Both thermal solvers, both heat-capacity methods.
    nh = 20 * 24
    time = DateTime(2001, 7, 1) .+ Hour.(0:nh-1)
    hour = [Dates.hour(t) for t in time]
    day = [(Dates.value(t - time[1]) ÷ 86_400_000) for t in time]
    cold = [d % 5 < 2 for d in day]                        # 2 cold, snowy days in every 5
    P = zeros(nh)
    for (i, d) in pairs(day)
        cold[i] || continue
        h = hour[i]
        P[i] = h == 3 ? 5e-324 : h == 6 ? 0.5 : h in (9, 10) ? 2.0 : h == 15 ? 12.0 :
               h in (18, 19, 20) ? 0.3 : 0.0
    end
    T = [c ? 263.0 : 276.0 + 3 * sin(2π * (h - 9) / 24) for (c, h) in zip(cold, hour)]
    sza = [clamp(40.0 + 45.0 * abs(h - 12) / 12, 0.0, 95.0) for h in hour]
    sw = [s < 90 ? 600.0 * cosd(s) : 0.0 for s in sza]
    cf = initialize_forcing(time, T, fill(85000.0, nh), P, fill(4.0, nh), sw,
        fill(290.0, nh), fill(500.0, nh); temperature_air_mean=268.0, wind_speed_mean=4.0,
        precipitation_mean=sum(P) * 365 / 20, temperature_observation_height=2.0,
        wind_observation_height=10.0, solar_zenith_angle=sza, cloud_optical_thickness=1.0,
        shortwave_downward_diffuse=0.3 .* sw)
    for solver in (GEMB.ExplicitThermal(), GEMB.ImplicitThermal()),
        hc in (:CuffeyPaterson, :constant)
        mp = ModelParameters(output_frequency=:all, thermal_solver=solver,
            heat_capacity_method=hc, albedo_ice=0.40)
        profile = initialize_profile(mp, cf; constant_density=true, constant_temperature=true)
        out = gemb(profile, cf, mp; verbose=true)
        dz1 = parent(out[:dz])[1, :]
        ρ1 = parent(out[:density])[1, :]
        # The thin-layer path ran: snow cells over ice thinner than column_dzmin.
        @test any((dz1 .< mp.column_dzmin) .& (ρ1 .< GEMB.DENSITY_PORE_CLOSEOFF))
        @test all((dz1 .* ρ1 .>= GEMB.snow_layer_min_mass(3600.0) * (1 - 1e-12)) .||
                  (ρ1 .>= GEMB.DENSITY_PORE_CLOSEOFF) .|| (dz1 .>= mp.column_dzmin))
        @test sum(parent(out[:melt])) > 0
    end
end
