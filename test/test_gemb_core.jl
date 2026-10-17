# Tests for gemb_core

# Helper to create a standard ClimateForcingStep for gemb_core tests
function _make_core_cfs(; dt=3600.0, temperature_air=265.0, precipitation=0.0,
                          shortwave_downward=200.0, wind_speed=5.0)
    return GEMB.ClimateForcingStep(
        dt,                   # dt [s]
        temperature_air,      # temperature_air
        100000.0,             # pressure_air
        precipitation,        # precipitation
        wind_speed,           # wind_speed
        shortwave_downward,   # shortwave_downward
        300.0,                # longwave_downward
        400.0,                # vapor_pressure
        260.0,                # temperature_air_mean
        5.0,                  # wind_speed_mean
        200.0,                # precipitation_mean
        2.0,                  # temperature_observation_height
        2.0,                  # wind_observation_height
        0.0,                  # black_carbon_snow
        0.0,                  # black_carbon_ice
        0.0,                  # cloud_optical_thickness
        60.0,                 # solar_zenith_angle
        50.0,                 # shortwave_downward_diffuse
        0.0,                  # cloud_fraction
    )
end

# Helper to create ModelParameters for gemb_core tests
function _make_core_mp(; dt=3600.0, column_depth_max=0.9)
    return GEMB.ModelParameters(
        albedo_method=:GardnerSharp,
        albedo_ice=0.45,
        albedo_snow=0.85,
        albedo_fixed=0.7,
        albedo_density_threshold=1023.0,
        shortwave_subsurface_absorption=true,
        emissivity=0.98,
        emissivity_grain_radius_large=0.97,
        emissivity_method=:uniform,
        emissivity_grain_radius_threshold=10.0,
        thermal_conductivity_method=:Sturm,
        column_dzmin=0.05,
        column_dzmax=0.10,
        column_depth_max=column_depth_max,
        column_ztop=2.0,
        column_zy=1.1,
        new_snow_method=:Constant150,
        density_ice=917.0,
        water_irreducible_saturation=0.07,
        densification_method=:HerronLangway,
        densification_coeffs_M01=:Gre_RACMO_GS_SW0,
        surface_roughness_effective_ratio=0.1,
        rain_temperature_threshold=273.15,
        dt_divisors=GEMB.fast_divisors(round(Int, dt * 10000)) ./ 10000,
    )
end

# Helper to create standard initial state NamedTuple
function _make_core_state(; n=10)
    return (
        temperature = 260.0 * ones(n),
        dz = 0.08 * ones(n),
        density = 400.0 * ones(n),
        water = zeros(n),
        grain_radius = 0.5 * ones(n),
        grain_dendricity = 0.5 * ones(n),
        grain_sphericity = 0.5 * ones(n),
        age = zeros(n),
        evaporation_condensation = 0.0,
        melt_surface = 0.0,
    )
end

@testset "Pipeline execution (smoke test)" begin
    state = _make_core_state()
    cfs = _make_core_cfs()
    mp = _make_core_mp()

    new_state, flux = GEMB.gemb_core(state, cfs, mp, true)

    # The column length is fixed: `gemb_core` defaults `n_target` to the incoming
    # length, and the count controller restores it before returning.
    @test length(new_state.dz) == 10
    # Total depth is likewise pinned to the incoming depth (0.8 m) by the basal trim.
    @test sum(new_state.dz) ≈ 0.8 atol=1e-9

    # Output sizes should be consistent
    @test length(new_state.dz) == length(new_state.temperature)
    @test length(new_state.dz) == length(new_state.density)

    # Fluxes should be scalar
    @test isa(flux.shortwave_net, Float64)
    @test isa(flux.heat_flux_sensible, Float64)
    @test isa(flux.heat_flux_latent, Float64)
    @test isa(flux.longwave_upward, Float64)

    # Compaction variables should be scalar
    @test isa(flux.densification_from_compaction, Float64)
    @test isa(flux.densification_from_melt, Float64)

    # Non-local mass exchanges: basal (trim_bottom!) and lateral (horizontal strain).
    @test isa(flux.mass_added, Float64)
    @test isa(flux.mass_lateral, Float64)
    @test flux.mass_lateral == 0.0   # horizontal_strain_rate defaults to 0
end

@testset "Accumulation event" begin
    state = _make_core_state()
    cfs = _make_core_cfs(precipitation=10.0, temperature_air=260.0)
    mp = _make_core_mp(column_depth_max=100.0)

    new_state, flux = GEMB.gemb_core(state, cfs, mp, true)

    # Accumulation prepends a fresh-snow cell (10 -> 11), and the count controller
    # reclaims the slot by merging a deep pair, so the length returns to 10.
    @test length(new_state.dz) == 10

    # Top layer should be fresh snow density
    @test round(new_state.density[1]; digits=2) == 150.0
end

@testset "Melt event" begin
    dt = 3600.0 * 3  # 3 hours
    temperature = 260.0 * ones(10)
    temperature[1] = GEMB.CtoK
    state = (
        temperature = temperature,
        dz = 0.08 * ones(10),
        density = 400.0 * ones(10),
        water = zeros(10),
        grain_radius = 0.5 * ones(10),
        grain_dendricity = 0.5 * ones(10),
        grain_sphericity = 0.5 * ones(10),
        age = zeros(10),
        evaporation_condensation = 0.0,
        melt_surface = 0.0,
    )
    cfs = _make_core_cfs(dt=dt, temperature_air=280.0, shortwave_downward=800.0)
    mp = _make_core_mp(dt=dt)

    new_state, flux = GEMB.gemb_core(state, cfs, mp, true)

    @test flux.melt > 0.0
end

@testset "Densification compaction" begin
    state = (
        temperature = 260.0 * ones(10),
        dz = 0.08 * ones(10),
        density = 300.0 * ones(10),
        water = zeros(10),
        grain_radius = 0.5 * ones(10),
        grain_dendricity = 0.5 * ones(10),
        grain_sphericity = 0.5 * ones(10),
        age = zeros(10),
        evaporation_condensation = 0.0,
        melt_surface = 0.0,
    )
    cfs = _make_core_cfs()
    mp = _make_core_mp()

    new_state, flux = GEMB.gemb_core(state, cfs, mp, true)

    @test flux.densification_from_compaction > 0.0
    @test all(new_state.density .>= 300.0)
end

@testset "Vapor exchange uses the phase of its latent heat" begin
    mp = _make_core_mp()
    h_w = GEMB.specific_enthalpy_water(mp)
    h_i0 = GEMB.specific_enthalpy(mp, GEMB.CtoK)
    n = 4
    col() = (fill(GEMB.CtoK, n), fill(0.05, n), fill(400.0, n), [0.5; zeros(n - 1)])

    # Condensation at a melting surface: c kg enters water[1] at h_ice(0 °C) + LF.
    T, dz, ρ, w = col()
    c = 0.2
    E = GEMB.apply_vapor_exchange!(T, dz, ρ, w, c, c, mp)
    @test w[1] ≈ 0.7 rtol = 1e-14
    @test dz[1] == 0.05
    @test E ≈ c * (h_i0 + GEMB.LF) rtol = 1e-12
    @test h_w ≈ h_i0 + GEMB.LF rtol = 1e-12

    # Evaporation smaller than the pore water takes only water.
    T, dz, ρ, w = col()
    E = GEMB.apply_vapor_exchange!(T, dz, ρ, w, -0.3, -0.3, mp)
    @test w[1] ≈ 0.2 rtol = 1e-14
    @test dz[1] == 0.05 && T[1] == GEMB.CtoK
    @test E ≈ -0.3 * h_w rtol = 1e-12

    # Evaporation beyond the pore water takes the rest from ice, which costs the cell LF per kg
    # on top of the LV the flux charged: the cell cools by LF·s/(M₁·c).
    T, dz, ρ, w = col()
    s = 0.4
    M1 = 0.05 * 400.0
    E = GEMB.apply_vapor_exchange!(T, dz, ρ, w, -(0.5 + s), -(0.5 + s), mp)
    @test w[1] == 0.0
    @test dz[1] * ρ[1] ≈ M1 - s rtol = 1e-12
    @test (M1 - s) * GEMB.specific_enthalpy(mp, T[1]) ≈ M1 * h_i0 - s * (h_i0 + GEMB.LF) rtol = 1e-12
    @test E ≈ -0.5 * h_w - s * (h_i0 + GEMB.LF) rtol = 1e-12

    # Below the melting point the exchange is with ice at constant density.
    T, dz, ρ, w = col()
    T .= 260.0
    E = GEMB.apply_vapor_exchange!(T, dz, ρ, w, -0.1, 0.0, mp)
    @test dz[1] ≈ 0.05 - 0.1 / 400.0 rtol = 1e-12
    @test w[1] == 0.5
    @test E ≈ -0.1 * GEMB.specific_enthalpy(mp, 260.0) rtol = 1e-12

    # The flux charged LS per kg exchanged with ice and LV per kg with water.
    cfs = _make_core_cfs()
    for T_s in (260.0, GEMB.CtoK)
        _, lhf, L = GEMB.turbulent_heat_flux(T_s, 1.2, 1e-3, 1e-4, 1e-4, cfs)
        @test L == (T_s < GEMB.CtoK ? GEMB.LS : GEMB.LV)
    end
    temperature = fill(GEMB.CtoK, 10)
    _, _, _, lhf, _, ec, ec_liquid = GEMB.calculate_temperature(temperature, fill(0.08, 10),
        fill(400.0, 10), 0.5, fill(0.5, 10), [800.0; zeros(9)],
        _make_core_cfs(temperature_air=280.0, shortwave_downward=800.0), mp, true)
    @test ec_liquid == ec
    @test lhf * 3600.0 ≈ GEMB.LV * ec rtol = 1e-12
end

@testset "Condensation on a melting surface closes the column energy budget" begin
    # Humid warm air over a melting, wet surface condenses onto it. `verbose = true` makes
    # gemb_core check the column energy against the vapor enthalpy above.
    temperature = fill(GEMB.CtoK, 10)
    state = merge(_make_core_state(), (temperature=temperature, water=[0.3; zeros(9)]))
    cfs = GEMB.ClimateForcingStep(; dt=3600.0, temperature_air=280.0, pressure_air=100000.0,
        wind_speed=8.0, shortwave_downward=0.0, longwave_downward=320.0, vapor_pressure=900.0,
        temperature_air_mean=260.0, wind_speed_mean=5.0, precipitation_mean=200.0,
        temperature_observation_height=2.0, wind_observation_height=2.0)
    new_state, flux = GEMB.gemb_core(state, cfs, _make_core_mp(), true)
    @test new_state.evaporation_condensation > 0.0
    @test flux.heat_flux_latent > 0.0
end

@testset "Water the thermal solve freezes is counted" begin
    # The thin, water-heavy surface cell of `_wet_cell_handoff!` (W/M = 1.04) over a dry column at the
    # melting point, losing 100 W m-2: the air is at 273.15 K and saturated over water, so the
    # turbulent fluxes vanish, and there is no shortwave. In 1800 s the solve itself freezes
    # 100·1800/LF of the water at the melting point, in one cell.
    T0 = GEMB.CtoK
    state = merge(_make_core_state(), (temperature=fill(T0, 10), dz=[0.0061; fill(0.08, 9)],
        density=[93.4; fill(400.0, 9)], water=[0.597; zeros(9)]))
    cfs = GEMB.ClimateForcingStep(; dt=1800.0, temperature_air=T0, pressure_air=100000.0,
        wind_speed=5.0, shortwave_downward=0.0, longwave_downward=GEMB.SB * T0^4 - 100.0 / 0.98,
        vapor_pressure=GEMB.saturation_vapor_pressure_water(T0), temperature_air_mean=260.0,
        wind_speed_mean=5.0, precipitation_mean=200.0, temperature_observation_height=2.0,
        wind_observation_height=2.0)
    _, flux = GEMB.gemb_core(state, cfs, _make_core_mp(dt=1800.0), true)
    @test flux.refreeze_thermal ≈ 100.0 * 1800.0 / GEMB.LF rtol = 1e-9
    @test flux.refreeze_thermal_count === 1
    # It is part of refreeze, and calculate_melt finds nothing more to refreeze.
    @test flux.refreeze ≈ flux.refreeze_thermal rtol = 1e-12

    # An ordinary column never takes the branch.
    _, flux = GEMB.gemb_core(_make_core_state(), _make_core_cfs(), _make_core_mp(), true)
    @test flux.refreeze_thermal === 0.0
    @test flux.refreeze_thermal_count === 0
end
