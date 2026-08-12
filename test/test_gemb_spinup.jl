using Test
using GEMB
using Dates

@testset "gemb_spinup" begin
    @testset "Basic spinup execution" begin
        # Create parameters
        params = initialize_parameters(
            densification_method = :Arthern,
            output_frequency = :last
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
        output = gemb_spinup(profile, forcing, params; max_iterations=3, verbose=false)

        # Check output structure
        @test output isa DimStack
        @test haskey(output, :temperature)
        @test haskey(output, :density)
        @test haskey(output, :dz)

        # Check that we get a profile with Z dimension
        @test length(output[:dz]) > 0
    end

    @testset "Spinup convergence test" begin
        # Test that running more cycles leads to more stable profiles

        params = initialize_parameters(
            output_frequency = :last
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
        output_3 = gemb_spinup(profile, forcing, params; max_iterations=3)
        output_5 = gemb_spinup(profile, forcing, params; max_iterations=5)

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
        # Test that gemb_spinup returns a valid profile

        params = initialize_parameters(
            output_frequency = :last
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
        spunup_profile = gemb_spinup(profile, forcing, params; max_iterations=3)

        # Check profile has required fields
        @test haskey(spunup_profile, :temperature)
        @test haskey(spunup_profile, :density)
        @test haskey(spunup_profile, :dz)
        @test haskey(spunup_profile, :grain_radius)

        # Profile arrays should have same length
        n_layers = length(spunup_profile[:dz])
        @test length(spunup_profile[:density]) == n_layers
        @test length(spunup_profile[:temperature]) == n_layers
    end

    @testset "Herron-Langway steady-state profile" begin
        zc = -collect(0.05:0.1:60.0)   # cell-center depths, negative below surface
        ρ0, ρi = 350.0, 910.0
        d = herron_langway_steady_state(zc, 253.0, 300.0, ρ0, ρi)

        # Surface ≈ ρ0, deep asymptote → ρ_ice, monotonic non-decreasing.
        @test isapprox(d[1], ρ0; atol=2.0)
        @test all(diff(d) .>= -1e-9)
        @test all(d .<= ρi + 1e-9)
        d_deep = herron_langway_steady_state(-collect(0.5:1.0:250.0), 253.0, 300.0, ρ0, ρi)
        @test d_deep[end] > 900.0

        # 550 kg/m³ crossover at a physically sensible depth (a few m to tens of m).
        i55 = findfirst(>=(550.0), d)
        crossover_depth = -zc[i55]
        @test 2.0 < crossover_depth < 40.0

        # Negligible accumulation → no firn forms → pure ice.
        d_ice = herron_langway_steady_state(zc, 253.0, 0.0, ρ0, ρi)
        @test all(d_ice .== ρi)
    end

    @testset "annual_pdd_melt" begin
        n = 365
        time = DateTime(2020, 1, 1) .+ Day.(0:n-1)
        # Constant +5 °C above the melt point.
        cf = initialize_forcing(
            time, fill(278.15, n), fill(85000.0, n), fill(1.0, n), fill(5.0, n),
            fill(100.0, n), fill(200.0, n), fill(100.0, n);
            temperature_air_mean=278.15, wind_speed_mean=5.0, precipitation_mean=365.25)
        melt = annual_pdd_melt(cf; ddf_snow=3.0)
        expected = 3.0 * (5.0 * n) * (365.25 / n)   # DDF * PDD, annualized
        @test isapprox(melt, expected; rtol=1e-9)

        # All sub-freezing → zero melt.
        cf_cold = initialize_forcing(
            time, fill(253.0, n), fill(85000.0, n), fill(1.0, n), fill(5.0, n),
            fill(100.0, n), fill(200.0, n), fill(100.0, n);
            temperature_air_mean=253.0, wind_speed_mean=5.0, precipitation_mean=365.25)
        @test annual_pdd_melt(cf_cold) == 0.0
    end

    @testset "Ablation regime → pure ice + clamped temperature" begin
        n = 365
        time = DateTime(2020, 1, 1) .+ Day.(0:n-1)
        # Warm, melt-dominated site.
        cf = initialize_forcing(
            time, fill(278.15, n), fill(85000.0, n), fill(1.0, n), fill(5.0, n),
            fill(100.0, n), fill(200.0, n), fill(100.0, n);
            temperature_air_mean=278.15, wind_speed_mean=5.0, precipitation_mean=365.25)
        mp = initialize_parameters()
        prof = initialize_profile(mp, cf)

        @test all(parent(prof[:density]) .== mp.density_ice)
        @test all(parent(prof[:temperature]) .<= 273.15 + 1e-9)
    end

    @testset "Steady-state init: firn profile + faster convergence" begin
        n = 365
        time = DateTime(2020, 1, 1) .+ Day.(0:n-1)
        doy = collect(1:n)
        # Cold, dry, accumulating site (dry-snow zone).
        temp_seasonal = 253.0 .+ 10.0 .* cos.(2π .* doy ./ 365.0 .- π)
        cf = initialize_forcing(
            time, temp_seasonal, fill(85000.0, n), fill(0.8, n), fill(4.0, n),
            fill(80.0, n), fill(190.0, n), fill(80.0, n);
            temperature_air_mean=253.0, wind_speed_mean=4.0,
            precipitation_mean=0.8 * 365.25)
        mp = initialize_parameters(densification_method=:HerronLangway,
                                   output_frequency=:last)

        prof_ss  = initialize_profile(mp, cf)                     # steady-state (default)
        prof_ice = initialize_profile(mp, cf; steady_state=false) # legacy all-ice

        # Steady-state start is a graded firn column, not pure ice.
        d_ss = parent(prof_ss[:density])
        @test !all(d_ss .== mp.density_ice)
        @test d_ss[1] < 400.0            # low-density surface firn
        @test all(parent(prof_ice[:density]) .== mp.density_ice)

        # Both converge to the same depth-averaged density; steady-state faster.
        function spin_cycles(prof0, thr, depth; maxit=120)
            prof = prof0
            prev = nothing
            cycles = maxit
            for c in 1:maxit
                prof = gemb_profile(gemb(prof, cf, mp; verbose=false))
                zc = -GEMB.dz2z(prof[:dz])
                rho = prof[:density]
                idx = zc .<= depth
                a = sum(rho[idx]) / count(idx)
                if prev !== nothing && abs(a - prev) < thr
                    cycles = c
                    break
                end
                prev = a
            end
            zc = -GEMB.dz2z(prof[:dz]); rho = prof[:density]
            idx = zc .<= depth
            return sum(rho[idx]) / count(idx), cycles
        end

        ρ_ss,  c_ss  = spin_cycles(prof_ss,  0.2, 40.0)
        ρ_ice, c_ice = spin_cycles(prof_ice, 0.2, 40.0)

        # Same equilibrium (loose tolerance: early-exit thresholds differ slightly).
        @test isapprox(ρ_ss, ρ_ice; atol=20.0)
        # Steady-state initialization converges in far fewer cycles.
        @test c_ss < c_ice
    end

    @testset "Spinup with zero accumulation" begin
        # Edge case: no precipitation

        params = initialize_parameters(
            output_frequency = :last
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
        output = gemb_spinup(profile, forcing, params; max_iterations=2)

        @test output isa DimStack
        @test haskey(output, :temperature)

        # Profile should still exist
        temps = output[:temperature][Ti=1]
        @test sum(.!isnan.(temps)) > 0
    end
end
