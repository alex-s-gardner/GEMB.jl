# Full synthetic data regression test
# Validates against MATLAB GEMB reference output from GEMB_example_synthetic.m

using GEMB: Statistics

@testset "Synthetic regression (MATLAB reference)" begin
    # Generate 3-hourly synthetic climate forcing
    cf = simulate_climate_forcing("test_1", 3)

    # Initialize model parameters
    mp = ModelParameters(output_frequency=:daily)

    # Initialize profile
    profile = initialize_profile(mp, cf)

    # Create climatological forcing and spin up
    cf_climatology = forcing_climatology(cf)
    mp_spinup = ModelParameters(output_frequency=:last)
    profile_spunup = gemb_spinup(profile, cf_climatology, mp_spinup, 75)

    # Run GEMB with spun-up profile
    output = gemb(profile_spunup, cf, mp)

    # MATLAB reference values (from GEMB_example_synthetic.m)
    mean_albedo = Statistics.mean(parent(output[:albedo_surface]))
    total_melt = sum(parent(output[:melt]))
    total_runoff = sum(parent(output[:runoff]))

    @test mean_albedo ≈ 0.821303 atol=1e-4
    @test total_melt ≈ 11504.085424 atol=0.1
    @test total_runoff ≈ 5217.635140 atol=0.1
end
