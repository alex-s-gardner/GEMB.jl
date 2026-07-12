using Test
using GEMB
using Dates

@testset "gemb_spinup" begin
    @testset "Basic spinup execution" begin
        # Create parameters
        params = initialize_parameters(
            densification_method = "Arthern",
            output_frequency = "last"
        )

        # Create 1-year climatology
        n_days = 365
        start_time = DateTime(2020, 1, 1)
        time = start_time .+ Day.(0:n_days-1)

        # Simple seasonal forcing
        day_of_year = collect(1:n_days)
        temp_seasonal = 260.0 .+ 15.0 .* cos.(2π .* day_of_year ./ 365.0 .- π)

        forcing = initialize_forcing(
            time,
            temp_seasonal,  # temperature_air
            fill(85000.0, n_days),  # pressure_air
            fill(1.0, n_days),  # precipitation
            fill(5.0, n_days),  # wind_speed
            fill(100.0, n_days),  # shortwave_downward
            fill(200.0, n_days),  # longwave_downward
            fill(100.0, n_days),  # vapor_pressure
            temperature_observation_height = 2.0,
            wind_observation_height = 10.0
        )

        # Initialize profile
        profile = initialize_profile(params, forcing)

        # Run spinup with 3 cycles (fast test)
        n_cycles = 3
        output = gemb_spinup(profile, forcing, params, n_cycles, verbose=false)

        # Check output structure
        @test output isa DimStack
        @test haskey(output, :temperature)
        @test haskey(output, :density)
        @test haskey(output, :dz)

        # Check that we get output for final state
        @test length(dims(output, Ti)) > 0
        @test length(dims(output, Z)) > 0
    end

    @testset "Spinup convergence test" begin
        # Test that running more cycles leads to more stable profiles

        params = initialize_parameters(
            output_frequency = "last"
        )

        # Annual climatology
        n_days = 365
        start_time = DateTime(2020, 1, 1)
        time = start_time .+ Day.(0:n_days-1)

        forcing = initialize_forcing(
            time,
            fill(255.0, n_days),  # temperature_air - Cold, stable
            fill(85000.0, n_days),  # pressure_air
            fill(0.5, n_days),  # precipitation
            fill(3.0, n_days),  # wind_speed
            fill(50.0, n_days),  # shortwave_downward
            fill(180.0, n_days),  # longwave_downward
            fill(80.0, n_days),  # vapor_pressure
            temperature_observation_height = 2.0,
            wind_observation_height = 10.0
        )

        # Run with different numbers of cycles
        profile = initialize_profile(params, forcing)
        output_3 = gemb_spinup(profile, forcing, params, 3)
        output_5 = gemb_spinup(profile, forcing, params, 5)

        # Extract final profiles
        temp_3 = output_3[:temperature][Ti=1]
        temp_5 = output_5[:temperature][Ti=1]

        # Both should be valid
        @test all(isfinite.(temp_3[.!isnan.(temp_3)]))
        @test all(isfinite.(temp_5[.!isnan.(temp_5)]))

        # Longer spinup should produce deeper profiles
        n_layers_3 = sum(.!isnan.(temp_3))
        n_layers_5 = sum(.!isnan.(temp_5))
        @test n_layers_5 >= n_layers_3
    end

    @testset "Profile extraction after spinup" begin
        # Test that gemb_profile works after spinup

        params = initialize_parameters(
            output_frequency = "last"
        )

        n_days = 365
        start_time = DateTime(2020, 1, 1)
        time = start_time .+ Day.(0:n_days-1)

        forcing = initialize_forcing(
            time,
            fill(260.0, n_days),  # temperature_air
            fill(85000.0, n_days),  # pressure_air
            fill(1.0, n_days),  # precipitation
            fill(5.0, n_days),  # wind_speed
            fill(100.0, n_days),  # shortwave_downward
            fill(200.0, n_days),  # longwave_downward
            fill(100.0, n_days),  # vapor_pressure
            temperature_observation_height = 2.0,
            wind_observation_height = 10.0
        )

        profile = initialize_profile(params, forcing)
        output = gemb_spinup(profile, forcing, params, 3)

        # Extract the final profile
        profile = gemb_profile(output)

        # Check profile has required fields
        @test hasfield(typeof(profile), :temperature)
        @test hasfield(typeof(profile), :density)
        @test hasfield(typeof(profile), :dz)
        @test hasfield(typeof(profile), :grain_radius)

        # Profile arrays should have same length
        n_layers = length(profile.temperature)
        @test length(profile.density) == n_layers
        @test length(profile.dz) == n_layers
    end

    @testset "Spinup with zero accumulation" begin
        # Edge case: no precipitation

        params = initialize_parameters(
            output_frequency = "last"
        )

        n_days = 100  # Shorter for zero accumulation test
        start_time = DateTime(2020, 1, 1)
        time = start_time .+ Day.(0:n_days-1)

        forcing = initialize_forcing(
            time,
            fill(250.0, n_days),  # temperature_air
            fill(85000.0, n_days),  # pressure_air
            zeros(n_days),  # precipitation - Zero accumulation
            fill(3.0, n_days),  # wind_speed
            zeros(n_days),  # shortwave_downward
            fill(150.0, n_days),  # longwave_downward
            fill(50.0, n_days),  # vapor_pressure
            temperature_observation_height = 2.0,
            wind_observation_height = 10.0
        )

        # Should still work, just won't grow
        profile = initialize_profile(params, forcing)
        output = gemb_spinup(profile, forcing, params, 2)

        @test output isa DimStack
        @test haskey(output, :temperature)

        # Profile should still exist
        temps = output[:temperature][Ti=1]
        @test sum(.!isnan.(temps)) > 0
    end
end
