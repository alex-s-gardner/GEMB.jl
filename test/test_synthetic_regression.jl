# Full synthetic data regression test
# Validates against MATLAB GEMB reference output from GEMB_example_synthetic.m
#
# NOTE: This test validates full GEMB workflow with 75-cycle spinup.
# Tolerances are relaxed compared to unit physics tests (rtol=1e-12) because:
# - 75 spinup cycles amplify floating-point differences
# - Non-linear physics creates chaotic sensitivity
# - Platform/compiler differences affect numerical libraries
# Results vary across Julia versions but remain physically consistent.

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

    # MATLAB reference values (generated with MATLAB GEMB)
    # Note: Full 75-cycle spinup creates chaotic sensitivity to platform/version differences
    # Tolerances reflect realistic expectations for accumulated numerical error

    if VERSION >= v"1.11"
        # Julia 1.11+: Moderate tolerances (modern versions show good agreement with platform variations)
        # Note: atol reflects that scalar-loop rewrites change FP evaluation order vs the original
        # broadcast chain, causing tiny per-step differences that compound over 75 spinup cycles
        # (7.9M iterations). arm64 and x86_64 diverge by up to ~13 kg/m² for runoff/melt and
        # ~5e-4 for albedo relative to the MATLAB reference values below.
        @test mean_albedo ≈ 0.821303 atol=5e-4       # 0.06% relative
        @test total_melt ≈ 11504.085424 atol=15.0    # 0.13% relative
        @test total_runoff ≈ 5217.635140 atol=15.0   # 0.3% relative
    else
        # Julia 1.10: Relaxed tolerances (significant platform/version differences observed)
        @test mean_albedo ≈ 0.821303 atol=1e-3       # 0.1% relative
        @test total_melt ≈ 11504.085424 atol=100.0   # 0.9% relative (CI variations: up to ±47 kg/m²)
        @test total_runoff ≈ 5217.635140 atol=500.0  # 10% relative (CI variations: up to ±413 kg/m²)
    end
end
