# Tests for calculate_temperature
# Translated from MATLAB test_calculate_temperature.m

@testset "Steady state" begin
    # If T_surf = temperature_air and LW balanced, T should remain roughly constant
    n = 10
    t_vec = fill(260.0, n)
    dz = fill(0.1, n)
    density = fill(400.0, n)
    water_surface = 0.0
    grain_radius = fill(0.5, n)
    shortwave_flux = zeros(n)

    sb = 5.67e-8
    cfs = GEMB.ClimateForcingStep(
        3600.0, 260.0, 100000.0, 0.0, 5.0, 0.0,
        sb * 260.0^4 * 0.97,  # longwave_downward balanced
        100.0, 260.0, 5.0, 0.0, 2.0, 2.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    )

    mp = GEMB.ModelParameters(
        density_ice=917.0,
        emissivity=0.97,
        emissivity_method=:uniform,
        emissivity_grain_radius_large=0.97,
        emissivity_grain_radius_threshold=10.0,
        surface_roughness_effective_ratio=0.1,
        thermal_conductivity_method=:Sturm,
        dt_divisors=Float64.(GEMB.fast_divisors(36000000)) ./ 10000
    )

    t_out, _, _, _, _, _ = GEMB.calculate_temperature(
        t_vec, dz, density, water_surface, grain_radius,
        shortwave_flux, cfs, mp, false)

    @test abs(t_out[1] - 260.0) < 1.0  # Allow some drift from turbulent fluxes
end

@testset "Solar heating" begin
    n = 10
    t_vec = fill(260.0, n)
    dz = fill(0.1, n)
    density = fill(400.0, n)
    water_surface = 0.0
    grain_radius = fill(0.5, n)
    shortwave_flux = zeros(n)
    shortwave_flux[1] = 200.0  # 200 W/m2 absorbed in top layer

    sb = 5.67e-8
    cfs = GEMB.ClimateForcingStep(
        3600.0, 260.0, 100000.0, 0.0, 5.0, 0.0,
        sb * 260.0^4 * 0.97,  # balance LW
        100.0, 260.0, 5.0, 0.0, 2.0, 2.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    )

    mp = GEMB.ModelParameters(
        density_ice=917.0,
        emissivity=0.97,
        emissivity_method=:uniform,
        emissivity_grain_radius_large=0.97,
        emissivity_grain_radius_threshold=10.0,
        surface_roughness_effective_ratio=0.1,
        thermal_conductivity_method=:Sturm,
        dt_divisors=Float64.(GEMB.fast_divisors(36000000)) ./ 10000
    )

    t_out, _, _, _, _, _ = GEMB.calculate_temperature(
        t_vec, dz, density, water_surface, grain_radius,
        shortwave_flux, cfs, mp, false)

    @test t_out[1] > 260.0  # Top layer should warm
end

@testset "Thermal diffusion" begin
    n = 10
    t_vec = fill(250.0, n)
    t_vec[1] = 273.0  # Hot surface

    dz = fill(0.1, n)
    density = fill(400.0, n)
    water_surface = 0.0
    grain_radius = fill(0.5, n)
    shortwave_flux = zeros(n)

    cfs = GEMB.ClimateForcingStep(
        10800.0, 273.0, 100000.0, 0.0, 0.1, 0.0,  # warm air, low wind to minimize turbulent cooling
        0.0,  # no longwave
        100.0, 260.0, 5.0, 0.0, 2.0, 2.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    )

    mp = GEMB.ModelParameters(
        density_ice=917.0,
        emissivity=0.0,  # disable radiative cooling
        emissivity_method=:uniform,
        emissivity_grain_radius_large=0.0,
        emissivity_grain_radius_threshold=10.0,
        surface_roughness_effective_ratio=0.1,
        thermal_conductivity_method=:Sturm,
        dt_divisors=Float64.(GEMB.fast_divisors(108000000)) ./ 10000
    )

    initial_gradient = t_vec[1] - t_vec[2]  # 23.0 (save before call mutates t_vec)

    t_out, _, _, _, _, _ = GEMB.calculate_temperature(
        t_vec, dz, density, water_surface, grain_radius,
        shortwave_flux, cfs, mp, false)

    # Diffusion should reduce the temperature gradient between layers 1 and 2
    final_gradient = t_out[1] - t_out[2]
    @test final_gradient < initial_gradient  # Gradient should decrease via diffusion
    @test t_out[2] > 250.0  # Layer 2 should warm from diffusion
end

@testset "Bottom boundary condition" begin
    n = 10
    t_vec = fill(260.0, n)
    t_vec[end] = 240.0  # Distinct bottom temp

    dz = fill(0.1, n)
    density = fill(400.0, n)
    water_surface = 0.0
    grain_radius = fill(0.5, n)
    shortwave_flux = zeros(n)

    cfs = GEMB.ClimateForcingStep(
        3600.0, 260.0, 100000.0, 0.0, 5.0, 0.0, 0.0,
        100.0, 260.0, 5.0, 0.0, 2.0, 2.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    )

    mp = GEMB.ModelParameters(
        density_ice=917.0,
        emissivity=0.97,
        emissivity_method=:uniform,
        emissivity_grain_radius_large=0.97,
        emissivity_grain_radius_threshold=10.0,
        surface_roughness_effective_ratio=0.1,
        thermal_conductivity_method=:Sturm,
        dt_divisors=Float64.(GEMB.fast_divisors(36000000)) ./ 10000
    )

    t_out, _, _, _, _, _ = GEMB.calculate_temperature(
        t_vec, dz, density, water_surface, grain_radius,
        shortwave_flux, cfs, mp, false)

    @test t_out[end] == 240.0  # Bottom T fixed (Dirichlet BC)
end

@testset "No NaN with large timestep" begin
    n = 10
    t_vec = fill(260.0, n)
    dz = fill(0.1, n)
    density = fill(400.0, n)
    water_surface = 0.0
    grain_radius = fill(0.5, n)
    shortwave_flux = zeros(n)

    cfs = GEMB.ClimateForcingStep(
        86400.0, 260.0, 100000.0, 0.0, 5.0, 0.0, 0.0,
        100.0, 260.0, 5.0, 0.0, 2.0, 2.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    )

    mp = GEMB.ModelParameters(
        density_ice=917.0,
        emissivity=0.97,
        emissivity_method=:uniform,
        emissivity_grain_radius_large=0.97,
        emissivity_grain_radius_threshold=10.0,
        surface_roughness_effective_ratio=0.1,
        thermal_conductivity_method=:Sturm,
        dt_divisors=Float64.(GEMB.fast_divisors(864000000)) ./ 10000
    )

    t_out, _, _, _, _, _ = GEMB.calculate_temperature(
        t_vec, dz, density, water_surface, grain_radius,
        shortwave_flux, cfs, mp, false)

    @test !any(isnan.(t_out))  # Should not explode
end
