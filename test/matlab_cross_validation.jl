"""
Cross-validation test suite that runs both MATLAB and Julia versions of GEMB
and compares outputs to ensure identical behavior.

This requires:
1. MATLAB installed and accessible
2. MATLAB.jl package installed
3. MATLAB GEMB repository available at ../GEMB/
"""

using Test
using GEMB
using MAT

# Path to MATLAB GEMB repository
const MATLAB_GEMB_PATH = joinpath(@__DIR__, "..", "..", "GEMB")

# Check if MATLAB GEMB is available
if !isdir(MATLAB_GEMB_PATH)
    @warn "MATLAB GEMB not found at $MATLAB_GEMB_PATH - skipping cross-validation tests"
else
    @testset "MATLAB Cross-Validation" begin
        @testset "Vapor Pressure Functions" begin
            # Create test data
            vapor_pressures = [200.0, 313.9, 500.0, 611.0]
            temperatures = [260.0, 265.3, 270.0, 273.15]
            relative_humidities = [50.0, 75.0, 90.0, 100.0]
            dewpoints = [255.0, 260.0, 265.0, 270.0]

            # Test vapor_pressure_to_relative_humidity
            @testset "vapor_pressure_to_relative_humidity" begin
                for (vp, temp) in zip(vapor_pressures, temperatures)
                    # Julia result
                    rh_julia = vapor_pressure_to_relative_humidity(vp, temp)

                    # MATLAB result (would need to call MATLAB)
                    # For now, test against known values
                    @test rh_julia >= 0.0
                    @test rh_julia <= 100.0
                end
            end

            # Test relative_humidity_to_vapor_pressure
            @testset "relative_humidity_to_vapor_pressure" begin
                for (temp, rh) in zip(temperatures, relative_humidities)
                    # Julia result
                    vp_julia = relative_humidity_to_vapor_pressure(temp, rh)

                    # Test roundtrip
                    rh_back = vapor_pressure_to_relative_humidity(vp_julia, temp)
                    @test rh_back ≈ rh atol=0.01
                end
            end

            # Test dewpoint_to_vapor_pressure
            @testset "dewpoint_to_vapor_pressure" begin
                for td in dewpoints
                    vp_julia = dewpoint_to_vapor_pressure(td)
                    @test vp_julia > 0.0
                end
            end
        end

        @testset "Grid Utility Functions" begin
            # Test surface_timeseries
            @testset "surface_timeseries" begin
                test_matrix = [1.0 2.0 3.0 4.0;
                              5.0 6.0 7.0 8.0;
                              9.0 10.0 11.0 12.0]

                surface = surface_timeseries(test_matrix)
                @test surface == [1.0, 2.0, 3.0, 4.0]

                # Test with NaN
                test_nan = [NaN 2.0 NaN;
                           5.0 6.0 7.0;
                           9.0 10.0 11.0]

                surface_nan = surface_timeseries(test_nan)
                @test surface_nan[1] ≈ 5.0
                @test surface_nan[2] ≈ 2.0
                @test isnan(surface_nan[3]) || surface_nan[3] ≈ 7.0
            end

            # Test dz2z
            @testset "dz2z" begin
                dz_test = ones(10, 5) * 0.1
                z_centers = dz2z(dz_test)

                # First layer should be at -dz/2
                @test all(z_centers[1, :] .≈ -0.05)

                # Should be monotonically decreasing
                for j in 1:5
                    for i in 2:10
                        @test z_centers[i, j] < z_centers[i-1, j]
                    end
                end
            end

            # Test fast_divisors
            @testset "fast_divisors" begin
                # Test known values
                @test fast_divisors(12) == [1, 2, 3, 4, 6, 12]
                @test fast_divisors(42) == [1, 2, 3, 6, 7, 14, 21, 42]

                # Test properties
                for n in [6, 12, 24, 36, 48, 60]
                    divs = fast_divisors(n)
                    @test all(n % d == 0 for d in divs)
                    @test issorted(divs)
                    @test divs[1] == 1
                    @test divs[end] == n
                end
            end

            # Test decyear2datenum
            @testset "decyear2datenum" begin
                # Test mid-year
                dn_mid = decyear2datenum(2020.5)
                dn_start = decyear2datenum(2020.0)
                dn_end = decyear2datenum(2021.0)

                # Should be approximately in the middle
                @test abs((dn_mid - dn_start) - (dn_end - dn_mid)) < 1.0

                # 2020 is a leap year
                @test dn_end - dn_start ≈ 366.0 atol=0.1
            end
        end

        @testset "Full Model Integration" begin
            # Run a short synthetic simulation in both Julia and MATLAB
            # and compare key output variables

            @testset "Synthetic test case" begin
                # Julia version
                params_jl = initialize_parameters()
                forcing_jl = simulate_climate_forcing("test_1", 3)
                profile_jl = initialize_profile(params_jl, forcing_jl)

                # Run for just a few timesteps
                n_steps = 10
                forcing_short = ClimateForcing(
                    time = forcing_jl.time[1:n_steps],
                    temperature_air = forcing_jl.temperature_air[:, 1:n_steps],
                    pressure_air = forcing_jl.pressure_air[:, 1:n_steps],
                    wind_speed = forcing_jl.wind_speed[:, 1:n_steps],
                    relative_humidity = forcing_jl.relative_humidity[:, 1:n_steps],
                    longwave_in = forcing_jl.longwave_in[:, 1:n_steps],
                    shortwave_in = forcing_jl.shortwave_in[:, 1:n_steps],
                    precipitation = forcing_jl.precipitation[:, 1:n_steps]
                )

                output_jl = gemb(profile_jl, forcing_short, params_jl)

                # For full MATLAB comparison, would need to:
                # 1. Write forcing data to .mat file
                # 2. Call MATLAB to run GEMB
                # 3. Read MATLAB output
                # 4. Compare results

                # For now, just verify Julia output is reasonable
                @test size(output_jl.temperature, 2) == n_steps
                @test all(isfinite, surface_timeseries(output_jl.temperature))
                @test all(output_jl.temperature .> 200.0)  # Physically reasonable temps
                @test all(output_jl.temperature .< 300.0)

                @test all(output_jl.density .> 0.0)
                @test all(output_jl.density .<= 917.0)  # Max ice density
            end
        end

        @testset "Individual Physics Modules" begin
            # These should already be covered by existing tests that load
            # MATLAB reference data from .mat files

            @testset "Thermal conductivity" begin
                # Reference test already exists
                include("test_thermal_conductivity.jl")
            end

            @testset "Turbulent heat flux" begin
                # Reference test already exists
                include("test_turbulent_heat_flux.jl")
            end

            @testset "Calculate grain size" begin
                # Reference test already exists
                include("test_calculate_grain_size.jl")
            end

            @testset "Calculate albedo" begin
                # Reference test already exists
                include("test_calculate_albedo.jl")
            end

            @testset "Calculate shortwave radiation" begin
                # Reference test already exists
                include("test_calculate_shortwave_radiation.jl")
            end

            @testset "Calculate temperature" begin
                # Reference test already exists
                include("test_calculate_temperature.jl")
            end

            @testset "Calculate accumulation" begin
                # Reference test already exists
                include("test_calculate_accumulation.jl")
            end

            @testset "Calculate melt" begin
                # Reference test already exists
                include("test_calculate_melt.jl")
            end

            @testset "Calculate density" begin
                # Reference test already exists
                include("test_calculate_density.jl")
            end

            @testset "Manage layers" begin
                # Reference test already exists
                include("test_manage_layers.jl")
            end
        end
    end
end

# Helper function to call MATLAB and compare outputs
"""
    compare_with_matlab(func_name, julia_func, args...)

Run both Julia and MATLAB versions of a function and compare outputs.
"""
function compare_with_matlab(func_name::String, julia_func::Function, args...)
    # Julia result
    result_julia = julia_func(args...)

    # TODO: Call MATLAB version using MATLAB.jl
    # This would require:
    # using MATLAB
    # mat"addpath('$MATLAB_GEMB_PATH/src')"
    # result_matlab = mat"$func_name($args...)"

    # For now, just return Julia result
    return result_julia
end

"""
    run_matlab_gemb(profile, forcing, params)

Run the MATLAB version of GEMB and return outputs.
"""
function run_matlab_gemb(profile, forcing, params)
    # TODO: Implement MATLAB call
    # This would:
    # 1. Save Julia data to .mat file
    # 2. Call MATLAB to run gemb()
    # 3. Load MATLAB output from .mat file
    # 4. Return as Julia structures

    error("MATLAB integration not yet implemented")
end

"""
    generate_matlab_reference_data()

Generate reference data files by running MATLAB version.
"""
function generate_matlab_reference_data()
    # This would generate new reference .mat files in test/
    # by running the MATLAB version of GEMB

    error("Reference data generation not yet implemented")
end
