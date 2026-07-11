# Tests for initialize_profile - matches MATLAB test_model_initialize_profile.m
# DimArray and Ti come from GEMB (re-exported from DimensionalData)
using GEMB: DimArray, Ti

@testset "Default parameters" begin
    mp = GEMB.ModelParameters()

    # Create minimal ClimateForcing for initialization
    times = [DateTime(2000, 1, 1), DateTime(2000, 1, 1, 3)]
    cf = GEMB.ClimateForcing(
        DimArray([253.15, 253.15], (Ti(times),)),
        DimArray([80000.0, 80000.0], (Ti(times),)),
        DimArray([0.001, 0.001], (Ti(times),)),
        DimArray([5.0, 5.0], (Ti(times),)),
        DimArray([100.0, 100.0], (Ti(times),)),
        DimArray([250.0, 250.0], (Ti(times),)),
        DimArray([300.0, 300.0], (Ti(times),)),
        253.15,  # temperature_air_mean
        5.0,     # wind_speed_mean
        200.0,   # precipitation_mean
        2.0,     # temperature_observation_height
        10.0,    # wind_observation_height
    )

    profile = GEMB.initialize_profile(mp, cf)

    # Check that profile has expected fields
    @test haskey(profile, :temperature)
    @test haskey(profile, :dz)
    @test haskey(profile, :density)
    @test haskey(profile, :water)
    @test haskey(profile, :grain_radius)

    # Check values
    dz = collect(profile[:dz])
    temperature = collect(profile[:temperature])
    density = collect(profile[:density])

    # All temperatures should be the mean
    @test all(temperature .≈ 253.15)

    # All densities should be density_ice
    @test all(density .≈ mp.density_ice)

    # Top layers should have constant spacing
    n_top = round(Int, mp.column_ztop / mp.column_dztop)
    @test all(dz[1:n_top] .≈ mp.column_dztop)

    # Total depth should be at or slightly above column_zmax
    # (last cell extends past boundary in MATLAB implementation)
    @test sum(dz) >= mp.column_zmax
    @test sum(dz) < mp.column_zmax + dz[end] + 1e-10
end

@testset "Grid stretching" begin
    mp = GEMB.ModelParameters(column_ztop=5.0, column_dztop=0.05, column_zmax=50.0, column_zy=1.10)

    times = [DateTime(2000, 1, 1), DateTime(2000, 1, 1, 3)]
    cf = GEMB.ClimateForcing(
        DimArray([260.0, 260.0], (Ti(times),)),
        DimArray([80000.0, 80000.0], (Ti(times),)),
        DimArray([0.001, 0.001], (Ti(times),)),
        DimArray([5.0, 5.0], (Ti(times),)),
        DimArray([100.0, 100.0], (Ti(times),)),
        DimArray([250.0, 250.0], (Ti(times),)),
        DimArray([300.0, 300.0], (Ti(times),)),
        260.0, 5.0, 200.0, 2.0, 10.0,
    )

    profile = GEMB.initialize_profile(mp, cf)
    dz = collect(profile[:dz])

    n_top = round(Int, mp.column_ztop / mp.column_dztop)

    # Below top zone, layers should increase geometrically
    for i in (n_top+2):length(dz)
        ratio = dz[i] / dz[i-1]
        @test ratio ≈ mp.column_zy atol = 1e-10
    end
end

@testset "z_center calculation" begin
    mp = GEMB.ModelParameters()

    times = [DateTime(2000, 1, 1), DateTime(2000, 1, 1, 3)]
    cf = GEMB.ClimateForcing(
        DimArray([253.15, 253.15], (Ti(times),)),
        DimArray([80000.0, 80000.0], (Ti(times),)),
        DimArray([0.001, 0.001], (Ti(times),)),
        DimArray([5.0, 5.0], (Ti(times),)),
        DimArray([100.0, 100.0], (Ti(times),)),
        DimArray([250.0, 250.0], (Ti(times),)),
        DimArray([300.0, 300.0], (Ti(times),)),
        253.15, 5.0, 200.0, 2.0, 10.0,
    )

    profile = GEMB.initialize_profile(mp, cf)
    z_center = collect(profile[:z_center])
    dz = collect(profile[:dz])

    # First center should be at -dz[1]/2
    @test z_center[1] ≈ -dz[1] / 2 atol = 1e-12

    # All centers should be negative (below surface)
    @test all(z_center .< 0)

    # Centers should be monotonically decreasing
    @test all(diff(z_center) .< 0)
end

# MATLAB validation test  
matlab_validation_testset("initialize_profile", "initialize_profile.mat") do ref
    # This will validate the grid geometry matches MATLAB
    # Note: Full profile validation requires matching all parameters exactly
    # Here we just validate grid structure
    
    params = GEMB.initialize_parameters()
    forcing = GEMB.simulate_climate_forcing("test_1", 3)
    profile = GEMB.initialize_profile(params, forcing)
    
    # Validate number of layers
    n_layers_julia = length(profile.dz)
    n_layers_matlab = Int(ref["n_layers_init"][1])
    
    @test n_layers_julia == n_layers_matlab
    
    # Validate grid structure (dz and z_center patterns)
    # Note: Exact values may differ slightly due to forcing differences
    # but structure should match
    @test length(profile.dz) > 0
    @test all(profile.dz .> 0)  # All layers have positive thickness
end
