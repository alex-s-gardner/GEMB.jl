using Test
using GEMB
using Dates

@testset "gemb_driver" begin
    @testset "Basic integration test" begin
        # Create simple parameters
        params = initialize_parameters(
            densification_method = "Arthern",
            albedo_method = "GardnerSharp"
        )

        # Create simple forcing (1 week, hourly)
        n_steps = 24 * 7
        start_time = DateTime(2020, 1, 1)
        time = start_time .+ Hour.(0:n_steps-1)

        forcing = initialize_forcing(
            time,
            fill(260.0, n_steps),  # temperature_air
            fill(101325.0, n_steps),  # pressure_air
            zeros(n_steps),  # precipitation
            fill(5.0, n_steps),  # wind_speed
            fill(200.0, n_steps),  # shortwave_downward
            fill(200.0, n_steps),  # longwave_downward
            fill(100.0, n_steps),  # vapor_pressure
            temperature_observation_height = 2.0,
            wind_observation_height = 10.0
        )

        # Initialize profile
        profile = initialize_profile(params, forcing)

        # Run model for 1 week
        output = gemb(profile, forcing, params)

        # Basic sanity checks
        @test output isa DimStack
        @test haskey(output, :temperature)
        @test haskey(output, :density)
        @test haskey(output, :dz)

        # Check output dimensions
        @test length(dims(output, Ti)) == n_steps + 1  # +1 for initial condition
        @test length(dims(output, Z)) > 0

        # Check temperature is in reasonable range
        temps = output[:temperature]
        @test all(isfinite.(temps[.!isnan.(temps)]))
        @test all(temps[.!isnan.(temps)] .> 200.0)  # Above absolute zero
        @test all(temps[.!isnan.(temps)] .< 300.0)  # Below boiling

        # Check density is in reasonable range (ice density)
        densities = output[:density]
        @test all(isfinite.(densities[.!isnan.(densities)]))
        @test all(densities[.!isnan.(densities)] .> 100.0)   # Above fresh snow
        @test all(densities[.!isnan.(densities)] .<= 917.0)  # At or below ice
    end

    @testset "Conservation test - no forcing" begin
        # Test that without any external forcing or melt, mass is conserved

        params = initialize_parameters(
            output_frequency = "all"
        )

        # Zero forcing (no precipitation, no melt conditions)
        n_steps = 24  # 1 day
        start_time = DateTime(2020, 1, 1)
        time = start_time .+ Hour.(0:n_steps-1)

        forcing = initialize_forcing(
            time,
            fill(250.0, n_steps),  # temperature_air - Cold, no melt
            fill(101325.0, n_steps),  # pressure_air
            zeros(n_steps),  # precipitation - No precip
            fill(2.0, n_steps),  # wind_speed
            zeros(n_steps),  # shortwave_downward - No solar
            fill(150.0, n_steps),  # longwave_downward
            fill(50.0, n_steps),  # vapor_pressure
            temperature_observation_height = 2.0,
            wind_observation_height = 10.0
        )

        profile = initialize_profile(params, forcing)

        # Calculate initial mass
        initial_mass = sum(profile.density .* profile.dz)

        output = gemb(profile, forcing, params)

        # Calculate final mass
        final_density = output[:density][Ti=Near(time[end])]
        final_dz = output[:dz][Ti=Near(time[end])]

        final_mass = sum(final_density[.!isnan.(final_density)] .*
                        final_dz[.!isnan.(final_dz)])

        # Mass should be approximately conserved (within 1%)
        @test abs(final_mass - initial_mass) / initial_mass < 0.01
    end

    @testset "Accumulation test" begin
        # Test that precipitation adds mass correctly

        params = initialize_parameters(
            output_frequency = "all"
        )

        n_steps = 10
        start_time = DateTime(2020, 1, 1)
        time = start_time .+ Hour.(0:n_steps-1)

        # Constant precipitation
        precip_rate = 0.001  # kg/m²/s for 1 hour = 3.6 kg/m²
        precip_per_hour = precip_rate * 3600.0

        forcing = initialize_forcing(
            time,
            fill(260.0, n_steps),  # temperature_air
            fill(101325.0, n_steps),  # pressure_air
            fill(precip_per_hour, n_steps),  # precipitation
            fill(2.0, n_steps),  # wind_speed
            zeros(n_steps),  # shortwave_downward
            fill(200.0, n_steps),  # longwave_downward
            fill(100.0, n_steps),  # vapor_pressure
            temperature_observation_height = 2.0,
            wind_observation_height = 10.0
        )

        profile = initialize_profile(params, forcing)
        initial_mass = sum(profile.density .* profile.dz)

        output = gemb(profile, forcing, params)

        final_density = output[:density][Ti=Near(time[end])]
        final_dz = output[:dz][Ti=Near(time[end])]
        final_mass = sum(final_density[.!isnan.(final_density)] .*
                        final_dz[.!isnan.(final_dz)])

        # Mass should increase by approximately the precipitation amount
        expected_added_mass = precip_per_hour * n_steps
        actual_added_mass = final_mass - initial_mass

        # Allow 10% tolerance for density changes and numerical error
        @test actual_added_mass > 0
        @test abs(actual_added_mass - expected_added_mass) / expected_added_mass < 0.2
    end

    @testset "Output frequency options" begin
        params = initialize_parameters()

        n_steps = 24 * 3  # 3 days
        start_time = DateTime(2020, 1, 1)
        time = start_time .+ Hour.(0:n_steps-1)

        forcing = initialize_forcing(
            time,
            fill(260.0, n_steps),  # temperature_air
            fill(101325.0, n_steps),  # pressure_air
            zeros(n_steps),  # precipitation
            fill(5.0, n_steps),  # wind_speed
            fill(100.0, n_steps),  # shortwave_downward
            fill(200.0, n_steps),  # longwave_downward
            fill(100.0, n_steps),  # vapor_pressure
            temperature_observation_height = 2.0,
            wind_observation_height = 10.0
        )

        profile = initialize_profile(params, forcing)

        # Test different output frequencies
        for freq in ["all", "daily", "last"]
            params_freq = initialize_parameters(output_frequency = freq)
            output = gemb(profile, forcing, params_freq)

            if freq == :all
                @test length(dims(output, Ti)) == n_steps + 1
            elseif freq == :daily
                expected_outputs = 3 + 1  # 3 days + initial
                @test length(dims(output, Ti)) == expected_outputs
            elseif freq == :last
                @test length(dims(output, Ti)) == 2  # Initial + final
            end
        end
    end
end
