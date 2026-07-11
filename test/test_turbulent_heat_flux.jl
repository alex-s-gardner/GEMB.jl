# Tests for turbulent_heat_flux - matches MATLAB test_turbulent_heat_flux.m

@testset "Stable conditions" begin
    # Set up forcing (stable: T_air > T_surface)
    cfs = GEMB.ClimateForcingStep(
        10800.0,      # dt
        268.0,        # temperature_air (warmer than surface)
        80000.0,      # pressure_air
        0.0,          # precipitation
        5.0,          # wind_speed
        0.0,          # shortwave_downward
        300.0,        # longwave_downward
        300.0,        # vapor_pressure
        255.0,        # temperature_air_mean
        5.0,          # wind_speed_mean
        200.0,        # precipitation_mean
        2.0,          # temperature_observation_height
        10.0,         # wind_observation_height
        0.0, 0.0, 0.0, 0.0, 0.0, 0.1  # BC, COT, SZA, SWdiff, CF
    )

    T_surface = 265.0
    density_air = 1.225
    z0 = 0.00012
    zT = z0 * 0.10
    zQ = z0 * 0.10

    shf, lhf, lh = GEMB.turbulent_heat_flux(T_surface, density_air, z0, zT, zQ, cfs)

    # In stable conditions with T_air > T_surface, SHF should be positive (toward surface)
    @test shf > 0
    # Latent heat should be sublimation since T_surface < 273.15
    @test lh ≈ GEMB.LS atol = 1e-6
    # Values should be finite
    @test isfinite(shf)
    @test isfinite(lhf)
end

@testset "Unstable conditions" begin
    # Set up forcing (unstable: T_surface > T_air)
    cfs = GEMB.ClimateForcingStep(
        10800.0,      # dt
        260.0,        # temperature_air (colder than surface)
        80000.0,      # pressure_air
        0.0,          # precipitation
        5.0,          # wind_speed
        0.0,          # shortwave_downward
        300.0,        # longwave_downward
        300.0,        # vapor_pressure
        255.0,        # temperature_air_mean
        5.0,          # wind_speed_mean
        200.0,        # precipitation_mean
        2.0,          # temperature_observation_height
        10.0,         # wind_observation_height
        0.0, 0.0, 0.0, 0.0, 0.0, 0.1
    )

    T_surface = 272.0
    density_air = 1.225
    z0 = 0.00012
    zT = z0 * 0.10
    zQ = z0 * 0.10

    shf, lhf, lh = GEMB.turbulent_heat_flux(T_surface, density_air, z0, zT, zQ, cfs)

    # In unstable conditions with T_surface > T_air, SHF should be negative (away from surface)
    @test shf < 0
    @test isfinite(shf)
    @test isfinite(lhf)
end

@testset "Melting surface uses vaporization" begin
    cfs = GEMB.ClimateForcingStep(
        10800.0, 275.0, 80000.0, 0.0, 5.0, 0.0, 300.0, 600.0,
        255.0, 5.0, 200.0, 2.0, 10.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.1
    )

    T_surface = 273.15  # At melting point
    density_air = 1.225
    z0 = 0.0013
    zT = z0 * 0.10
    zQ = z0 * 0.10

    _, _, lh = GEMB.turbulent_heat_flux(T_surface, density_air, z0, zT, zQ, cfs)

    # At melting point, should use latent heat of vaporization
    @test lh ≈ GEMB.LV atol = 1e-6
end
