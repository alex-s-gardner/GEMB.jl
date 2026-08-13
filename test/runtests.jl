using GEMB
using GEMB_ClimateForcing
using Test
using Dates
using GEMB: DimensionalData, DimArray, DimStack, Ti, Z, At, Near, dims

# Include test utilities for MATLAB validation
include("test_utils.jl")

@testset "GEMB.jl" begin
    @testset "Thermal Conductivity" begin
        include("test_thermal_conductivity.jl")
    end
    @testset "Turbulent Heat Flux" begin
        include("test_turbulent_heat_flux.jl")
    end
    @testset "Initialize Profile" begin
        include("test_initialize_profile.jl")
    end
    @testset "Shortwave Radiation" begin
        include("test_calculate_shortwave_radiation.jl")
    end
    @testset "Temperature" begin
        include("test_calculate_temperature.jl")
    end
    @testset "Albedo" begin
        include("test_calculate_albedo.jl")
    end
    @testset "Density" begin
        include("test_calculate_density.jl")
    end
    @testset "Accumulation" begin
        include("test_calculate_accumulation.jl")
    end
    @testset "Melt" begin
        include("test_calculate_melt.jl")
    end
    @testset "Grain Size" begin
        include("test_calculate_grain_size.jl")
    end
    @testset "Manage Layers" begin
        include("test_manage_layer_thickness.jl")
    end
    @testset "GEMB Core" begin
        include("test_gemb_core.jl")
    end
    @testset "GEMB Driver" begin
        include("test_gemb_driver.jl")
    end
    @testset "GEMB Spinup" begin
        include("test_gemb_spinup.jl")
    end
    @testset "Ablation Regime" begin
        include("test_ablation_regime.jl")
    end
    @testset "GEMB_ClimateForcing Integration" begin
        # Verify simulate → convert → run model end-to-end
        ds = GEMB_ClimateForcing.simulate_climate_forcing("test_1", 3)
        @test ds isa DimStack
        cf = GEMB.initialize_forcing(ds)
        @test cf isa GEMB.ClimateForcing
        mp = GEMB.ModelParameters()
        profile, mp = GEMB.initialize_profile(mp, cf)
        @test length(profile.dz) > 0
    end
    @testset "Grid Utilities" begin
        include("test_grid_utilities.jl")
    end

    @testset "Synthetic Regression" begin
        # Full model validation with relaxed tolerances for cross-platform reproducibility
        include("test_synthetic_regression.jl")
    end
end
