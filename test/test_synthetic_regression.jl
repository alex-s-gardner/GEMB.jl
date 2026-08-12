# Full synthetic data regression test (Julia self-consistency pin)
#
# IMPORTANT: This is NOT a MATLAB cross-validation test. The reference values
# below are GEMB.jl's own deterministic output, used to catch unintended changes
# to the full-model result. Per-module MATLAB fidelity is validated separately
# by the *.mat reference tests (thermal conductivity, melt, density, etc. at
# ~1e-12). The FULL synthetic run cannot be compared bit-for-bit to MATLAB for
# two structural reasons:
#
#   1. RNG streams differ. simulate_climate_forcing seeds a Mersenne Twister
#      with seed 42 in both languages, but MATLAB `randn` and Julia
#      `randn(::MersenneTwister)` draw different sequences (different Gaussian
#      sampling + MT seeding/tempering). The synthetic forcing is therefore a
#      different realization in each language (mean climatology matches to
#      ~0.03%, but noise-driven fields differ: precip ~5%, wind ~4%), which
#      propagates to different melt/runoff totals.
#   2. The Julia spinup was optimized in a way not present in the MATLAB version
#      (see commit 5afce9f, "profile and spinup optimization ... not implemented
#      in the Matlab version"), so the full workflow is intentionally not the
#      same algorithm.
#
# The reference values are the deterministic Julia output (bit-for-bit identical
# across runs on a given platform). Small cross-platform/version spread remains
# because the 75-cycle spinup (~7.9M iterations) amplifies floating-point
# evaluation-order differences (arm64 vs x86_64), so tolerances are modest
# rather than 1e-12.

using GEMB: Statistics
using GEMB_ClimateForcing

@testset "Synthetic regression (Julia self-consistency pin)" begin
    # Generate 3-hourly synthetic climate forcing (returns DimStack)
    ds = simulate_climate_forcing("test_1", 3)
    cf = GEMB.initialize_forcing(ds)

    # Initialize model parameters
    mp = ModelParameters(output_frequency=:daily)

    # Initialize profile. steady_state=false uses the pure-ice initialization
    # that the reference values below were generated with.
    profile = initialize_profile(mp, cf; steady_state=false)

    # Create climatological forcing and spin up
    cf_climatology = forcing_climatology(cf)
    mp_spinup = ModelParameters(output_frequency=:last)
    profile_spunup = gemb_spinup(profile, cf_climatology, mp_spinup; max_iterations=75)

    # Run GEMB with spun-up profile
    output = gemb(profile_spunup, cf, mp)

    # Julia reference values (GEMB.jl's own deterministic output — NOT MATLAB).
    mean_albedo = Statistics.mean(parent(output[:albedo_surface]))
    total_melt = sum(parent(output[:melt]))
    total_runoff = sum(parent(output[:runoff]))

    # Reference = deterministic Julia output on arm64 (macOS). Tolerances bracket
    # the cross-platform FP evaluation-order spread amplified over the 75-cycle
    # spinup (~7.9M iterations): x86_64 CI vs arm64 differs by ~28-32 kg/m² for
    # melt/runoff and ~6e-4 for albedo. atols are set ~2x that spread.
    @test mean_albedo ≈ 0.822099 atol=2e-3      # ~0.24% relative
    @test total_melt ≈ 11384.264341 atol=60.0   # ~0.53% relative
    @test total_runoff ≈ 5024.274200 atol=60.0  # ~1.2% relative
end
