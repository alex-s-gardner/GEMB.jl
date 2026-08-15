# Tests for thermal_conductivity - matches MATLAB test_thermal_conductivity.m

# Try to load MATLAB reference data if available
const MATLAB_REF_THERMAL = let
    ref_file = joinpath(@__DIR__, "reference_data", "thermal_conductivity.mat")
    if isfile(ref_file)
        try
            using MAT
            matread(ref_file)
        catch
            nothing
        end
    else
        nothing
    end
end

@testset "Sturm method snow" begin
    mp = GEMB.ModelParameters(density_ice=917.0, thermal_conductivity_method=:Sturm)
    density_snow = [300.0]
    temperature_in = [260.0]

    k_out = GEMB.thermal_conductivity(temperature_in, density_snow, mp)

    expected = 0.138 - 1.01e-3 * 300.0 + 3.233e-6 * 300.0^2
    @test k_out[1] ≈ expected atol = 1e-8
end

@testset "Calonne method snow" begin
    mp = GEMB.ModelParameters(density_ice=917.0, thermal_conductivity_method=:Calonne)
    density_snow = [300.0]
    temperature_in = [260.0]

    k_out = GEMB.thermal_conductivity(temperature_in, density_snow, mp)

    expected = 0.024 - 1.23e-4 * 300.0 + 2.5e-6 * 300.0^2
    @test k_out[1] ≈ expected atol = 1e-8
end

@testset "Calonne2019 method" begin
    mp = GEMB.ModelParameters(density_ice=917.0, thermal_conductivity_method=:Calonne2019)

    # Hand-computed from Calonne et al. (2019) eq. 5, written out independently of the
    # implementation (the same expression the Community Firn Model uses).
    function k_reference(T, ρ)
        θ = 1 / (1 + exp(-2 * 0.02 * (ρ - 450.0)))
        k_snow = 0.024 - 1.23e-4 * ρ + 2.5e-6 * ρ^2
        k_firn = 2.107 + 0.003618 * (ρ - 917.0)
        K_ice = 9.828 * exp(-5.7e-3 * T)
        return (1 - θ) * (K_ice / 2.107) * k_snow + θ * (K_ice / 2.107) * k_firn
    end

    for T in (240.0, 273.15), ρ in (200.0, 450.0, 700.0, 917.0)
        @test GEMB.thermal_conductivity(T, ρ, mp) ≈ k_reference(T, ρ) rtol = 1e-14
    end

    # Continuous into ice by construction: no ice branch, and K(917) is K_ice(T). The
    # residual is the logistic weight's distance from 1 at ρ = 917 (~3e-10 relative).
    for T in (240.0, 273.15)
        @test GEMB.thermal_conductivity(T, 917.0, mp) ≈ 9.828 * exp(-5.7e-3 * T) rtol = 1e-8
    end

    # Temperature dependence, which the 2011 fit lacks
    @test GEMB.thermal_conductivity(240.0, 500.0, mp) > GEMB.thermal_conductivity(270.0, 500.0, mp)

    # Monotone in density over the full column range
    k_profile = [GEMB.thermal_conductivity(260.0, ρ, mp) for ρ in 50.0:5.0:916.0]
    @test all(diff(k_profile) .> 0)
    @test all(k_profile .> 0)

    # Vector method agrees with the scalar method
    @test GEMB.thermal_conductivity([260.0, 250.0], [300.0, 800.0], mp) ==
          [GEMB.thermal_conductivity(260.0, 300.0, mp), GEMB.thermal_conductivity(250.0, 800.0, mp)]
end

@testset "Marchenko2019 method" begin
    mp = GEMB.ModelParameters(density_ice=917.0, thermal_conductivity_method=:Marchenko2019)
    k_marchenko(ρ) = 0.301e-2 * ρ - 0.724
    k_calonne(ρ) = 0.024 - 1.23e-4 * ρ + 2.5e-6 * ρ^2

    # In the calibration range (350-900) the linear fit is the larger of the two and is used
    for ρ in (400.0, 600.0, 900.0)
        @test GEMB.thermal_conductivity(260.0, ρ, mp) ≈ k_marchenko(ρ) rtol = 1e-14
    end

    # Below the ~321 crossing the Calonne floor takes over, so K never goes negative even
    # though the bare fit crosses zero at ρ ≈ 241
    for ρ in (50.0, 150.0, 240.0, 300.0)
        @test GEMB.thermal_conductivity(260.0, ρ, mp) ≈ k_calonne(ρ) rtol = 1e-14
    end
    @test k_marchenko(240.0) < 0                      # the fit alone is unphysical here
    @test GEMB.thermal_conductivity(260.0, 240.0, mp) > 0

    # Positive and monotone across the whole firn range, and continuous at the handover
    k_profile = [GEMB.thermal_conductivity(260.0, ρ, mp) for ρ in 50.0:1.0:916.0]
    @test all(k_profile .> 0)
    @test all(diff(k_profile) .> 0)
    @test maximum(abs.(diff(k_profile))) < 1e-2       # no step at the floor crossover

    # Higher than Sturm and Calonne through the firn, the direction RetMIP's cold bias implies
    mp_sturm = GEMB.ModelParameters(density_ice=917.0, thermal_conductivity_method=:Sturm)
    mp_calonne = GEMB.ModelParameters(density_ice=917.0, thermal_conductivity_method=:Calonne)
    for ρ in (500.0, 700.0, 800.0)
        @test GEMB.thermal_conductivity(260.0, ρ, mp) > GEMB.thermal_conductivity(260.0, ρ, mp_sturm)
        @test GEMB.thermal_conductivity(260.0, ρ, mp) > GEMB.thermal_conductivity(260.0, ρ, mp_calonne)
    end

    # Temperature-independent in firn, but the ice branch still applies at density_ice
    @test GEMB.thermal_conductivity(240.0, 600.0, mp) == GEMB.thermal_conductivity(270.0, 600.0, mp)
    @test GEMB.thermal_conductivity(240.0, 917.0, mp) ≈ 9.828 * exp(-5.7e-3 * 240.0) atol = 1e-8
end

@testset "Unknown conductivity method errors" begin
    # `ModelParameters` is not validated on construction, so the guard in the dispatch helper
    # is what catches a typo reaching the physics.
    mp = GEMB.ModelParameters(thermal_conductivity_method=:NotAMethod)
    @test_throws ErrorException GEMB.thermal_conductivity(260.0, 400.0, mp)
    @test_throws AssertionError GEMB.validate_parameters(mp)
    for m in (:Sturm, :Calonne, :Calonne2019, :Marchenko2019)
        @test GEMB.validate_parameters(GEMB.ModelParameters(thermal_conductivity_method=m)) === nothing
    end
end

@testset "Ice conductivity" begin
    mp = GEMB.ModelParameters(density_ice=917.0, thermal_conductivity_method=:Sturm)
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
    mp = GEMB.ModelParameters(density_ice=917.0, thermal_conductivity_method=:Sturm)
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
    mp = GEMB.ModelParameters(density_ice=917.0, thermal_conductivity_method=:Sturm)
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

@testset "MATLAB reference validation" begin
    if !isnothing(MATLAB_REF_THERMAL)
        # Test case from generate_reference_data.m
        temperature = MATLAB_REF_THERMAL["temperature"][:]
        density = MATLAB_REF_THERMAL["density"][:]

        # Sturm method
        mp_sturm = GEMB.ModelParameters(density_ice=910.0, thermal_conductivity_method=:Sturm)
        k_sturm_julia = GEMB.thermal_conductivity(temperature, density, mp_sturm)
        k_sturm_matlab = MATLAB_REF_THERMAL["K_sturm"][:]

        @test k_sturm_julia ≈ k_sturm_matlab rtol=1e-12 atol=1e-14

        # Calonne method
        mp_calonne = GEMB.ModelParameters(density_ice=910.0, thermal_conductivity_method=:Calonne)
        k_calonne_julia = GEMB.thermal_conductivity(temperature, density, mp_calonne)
        k_calonne_matlab = MATLAB_REF_THERMAL["K_calonne"][:]

        @test k_calonne_julia ≈ k_calonne_matlab rtol=1e-12 atol=1e-14

        @info "✓ MATLAB reference validation passed for thermal_conductivity"
    else
        @info "⊘ MATLAB reference data not available - skipping cross-validation"
        @info "  Run test/generate_reference_data.m in MATLAB to enable validation"
    end
end
