# Tests for thermal_conductivity - matches MATLAB test_thermal_conductivity.m

@testset "Sturm method snow" begin
    mp = GEMB.ModelParameters(density_ice=917.0, thermal_conductivity_method="Sturm")
    density_snow = [300.0]
    temperature_in = [260.0]

    k_out = GEMB.thermal_conductivity(temperature_in, density_snow, mp)

    expected = 0.138 - 1.01e-3 * 300.0 + 3.233e-6 * 300.0^2
    @test k_out[1] ≈ expected atol = 1e-8
end

@testset "Calonne method snow" begin
    mp = GEMB.ModelParameters(density_ice=917.0, thermal_conductivity_method="Calonne")
    density_snow = [300.0]
    temperature_in = [260.0]

    k_out = GEMB.thermal_conductivity(temperature_in, density_snow, mp)

    expected = 0.024 - 1.23e-4 * 300.0 + 2.5e-6 * 300.0^2
    @test k_out[1] ≈ expected atol = 1e-8
end

@testset "Ice conductivity" begin
    mp = GEMB.ModelParameters(density_ice=917.0, thermal_conductivity_method="Sturm")
    density_ice = [917.0]

    # Cold ice
    k_cold = GEMB.thermal_conductivity([240.0], density_ice, mp)
    expected_cold = 9.828 * exp(-5.7e-3 * 240.0)
    @test k_cold[1] ≈ expected_cold atol = 1e-8

    # Warm ice
    k_warm = GEMB.thermal_conductivity([270.0], density_ice, mp)
    expected_warm = 9.828 * exp(-5.7e-3 * 270.0)
    @test k_warm[1] ≈ expected_warm atol = 1e-8

    # Temperature dependence
    @test k_cold[1] != k_warm[1]
end

@testset "Mixed profile" begin
    mp = GEMB.ModelParameters(density_ice=917.0, thermal_conductivity_method="Sturm")
    temperature_vec = [260.0, 250.0]
    density_vec = [400.0, 920.0]  # 400=Snow, 920=Ice

    k_vec = GEMB.thermal_conductivity(temperature_vec, density_vec, mp)

    # Expected Snow (Sturm)
    exp_snow = 0.138 - 1.01e-3 * 400.0 + 3.233e-6 * 400.0^2
    # Expected Ice
    exp_ice = 9.828 * exp(-5.7e-3 * 250.0)

    @test k_vec[1] ≈ exp_snow atol = 1e-8
    @test k_vec[2] ≈ exp_ice atol = 1e-8
end

@testset "Density threshold boundary" begin
    mp = GEMB.ModelParameters(density_ice=917.0, thermal_conductivity_method="Sturm")
    d_vals = [917.0 - 1e-10, 917.0]
    t_val = 260.0

    k_out = GEMB.thermal_conductivity([t_val, t_val], d_vals, mp)

    # Just below threshold should be snow
    exp_snow = 0.138 - 1.01e-3 * d_vals[1] + 3.233e-6 * d_vals[1]^2
    @test k_out[1] ≈ exp_snow atol = 1e-8

    # At threshold should be ice
    exp_ice = 9.828 * exp(-5.7e-3 * t_val)
    @test k_out[2] ≈ exp_ice atol = 1e-8
end
