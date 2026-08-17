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
