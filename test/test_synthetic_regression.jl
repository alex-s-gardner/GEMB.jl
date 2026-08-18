# Full synthetic data regression test (self-consistency pin)
#
# The reference values below are GEMB.jl's own deterministic output, used to catch
# *unintended* changes to the full-model result. They are bit-for-bit reproducible
# across runs on a given platform; small cross-platform spread remains because the
# 75-cycle spinup (~7.9M iterations) amplifies floating-point evaluation-order
# differences (arm64 vs x86_64), so the tolerances are modest rather than 1e-12.
#
# This is a change detector, not a physics validation: an intended change to the
# defaults or to a physics scheme is expected to move these numbers, and re-pinning
# them is part of making such a change. Every re-pin is recorded below with what
# moved it and by how much, so the history of the fingerprint is auditable.

using GEMB: Statistics
using GEMB_ClimateForcing

@testset "Synthetic regression (self-consistency pin)" begin
    # Generate 3-hourly synthetic climate forcing (returns DimStack)
    ds = simulate_climate_forcing("test_1", 3)
    cf = initialize_forcing(ds)

    # Initialize model parameters
    mp = ModelParameters(output_frequency=:daily)

    # Initialize profile. Both escape-hatch flags give the pure-ice initialization,
    # and the deep column, that the reference values below were generated with.
    profile = initialize_profile(mp, cf; constant_density=true, constant_temperature=true)

    # Create climatological forcing and spin up
    cf_climatology = GEMB.forcing_climatology(cf)
    mp_spinup = ModelParameters(output_frequency=:last)
    profile_spunup = gemb_spinup(profile, cf_climatology, mp_spinup; max_iterations=75)

    # Run GEMB with spun-up profile
    output = gemb(profile_spunup, cf, mp)

    # The three pinned fields.
    mean_albedo = Statistics.mean(parent(output[:albedo_broadband]))
    total_melt = sum(parent(output[:melt]))
    total_runoff = sum(parent(output[:runoff]))

    # Reference = deterministic output on arm64 (macOS). Tolerances bracket
    # the cross-platform FP evaluation-order spread amplified over the 75-cycle
    # spinup (~7.9M iterations): measured x86_64 CI vs arm64 is ~51 kg/m² for melt.
    # atols are set ~2x that spread.
    #
    # The spread widened from ~30 to ~51 kg/m² with the fixed-length grid: merge and
    # split are threshold comparisons against the per-cell bands, so a 1-ulp band
    # difference can flip a single merge and slightly change the deep discretization.
    # This is evaluation-order sensitivity in the regridding, not a conservation
    # problem — the per-timestep mass/energy checks run under `verbose=true` and the
    # whole-run budget closes to ~1e-11 kg m-2 on both platforms.
    # Re-centered when internal energy became the enthalpy integral ∫c_p dT rather than
    # M·T·c_p. For the default :constant heat capacity the two are algebraically identical,
    # so the shift is arithmetic reordering only (~1e-13 per operation) — but the 75-cycle
    # spinup is not converged and amplifies it to ~0.3% here. The previous values
    # (0.821798 / 11398.551015 / 5017.057989) still fell inside these atols; they are
    # re-centered so the tolerance budget stays available for the cross-platform spread it
    # was sized for rather than being half-consumed by this change.
    # Re-centered again for two physics fixes, both of which change output by design:
    # densification now compacts against mean *snowfall* rather than mean precipitation
    # (`densification_accumulation = :accumulation`), and rain's sensible heat above the
    # melting point is carried at the water heat capacity rather than the ice one
    # (`rain_heat_capacity = :water`). Melt rises 48.7 kg m-2 (+0.43%) — the expected sign
    # for the second, which adds real energy at a site where 3.2% of precipitation falls as
    # rain. The previous values (0.821914 / 11368.809855 / 4997.438019) still fall inside
    # these atols; re-centered so the tolerance budget stays available for the
    # cross-platform spread it was sized for.
    #
    # Verified inert: with `densification_accumulation=:precipitation` and
    # `rain_heat_capacity=:ice` the run reproduces the pre-fix result to the last bit on all
    # three fields — 0.8219139818790556 / 11368.809854556823 / 4997.438018517874 — so the
    # shift below is the two physics changes and nothing else. (The albedo reference had
    # drifted to 0.822103 in transcription; melt and runoff were exact.)
    #
    # Re-pinned for v2.0.0, which changed four defaults on the recommendation of the two
    # firn-model intercomparisons and the shipped configurations of the CFM and IMAU-FDM:
    # `density_ice` 910 -> 917, `thermal_conductivity_method` :Sturm -> :Calonne2019,
    # `water_irreducible_method` :constant -> :ColeouLesaffre, and `grain_growth_method`
    # :Marbouty -> :Arthern. This is the first re-pin that moves the numbers by more than
    # the tolerance, which is what the major version records. Melt falls 453.6 kg m-2
    # (-4.0%) and runoff 255.1 (-5.1%); albedo moves 2e-5. The conductivity change is
    # responsible for most of it — the higher firn conductivity conducts winter cold
    # deeper, so the column enters the melt season colder and refreezes more of what melts
    # (measured in isolation: melt -570, minimum column temperature +4.25 K). The atols
    # are unchanged, since the cross-platform spread they were sized for has not moved.
    #
    # Re-centered for the thermal sub-step stability limit, which now comes from the scheme
    # actually being solved — `dt ≤ ρᵢcᵢdzᵢ/(Gᵢ + Gᵢ₋₁)` over the harmonic-mean face
    # conductances, excluding the Dirichlet bottom cell — rather than the textbook
    # uniform-grid form `0.5·ρᵢcᵢdzᵢ²/Kᵢ`. On the graded grid the old form was not
    # conservative (per-cell ratio 0.66–1.82); the correct limit is looser at the cell that
    # binds, so mean sub-steps per timestep fall 40.6 → 29.4. Melt falls 14.6 kg m-2 (-0.13%)
    # and runoff 8.7 (-0.18%); albedo moves 9e-5. The previous values (0.821926 /
    # 10963.944470 / 4760.650556) still fall inside these atols; re-centered so the tolerance
    # budget stays available for the cross-platform spread it was sized for. See
    # `GEMB._max_safe_dt`.
    #
    # Re-pinned for the corrected unstable-branch integrated stability functions. Paulson's
    # (1970) `Ψ` takes the *inverse* profile function as its argument, so the exponents are
    # positive, not negative; the previous form also carried Högström's (1988) `0.95`
    # `κ_H/κ_M` ratio inside the integral, which put `Ψ_h(0) = -0.0999` instead of 0 and left
    # the fluxes discontinuous at neutral stability. Both branches now reproduce
    # `Ψ(ζ) = ∫₀^ζ (1 - φ)/z dz` under numerical integration of their own published `φ`. The
    # corrected `Ψ` is larger, so `coef` shrinks and unstable-side turbulent exchange
    # strengthens: melt falls 576.4 kg m-2 (-5.3%) and runoff 283.8 (-6.0%), the expected sign
    # (more turbulent cooling and more sublimation at a melting surface). Albedo moves 9e-4.
    # This exceeds the atols, which is why it is a re-pin rather than a re-centering; the atols
    # themselves are unchanged, the cross-platform spread they were sized for having not moved.
    # See `GEMB._turbulent_heat_flux` and `GEMB.ZETA_UNSTABLE_MIN`.
    @test mean_albedo ≈ 0.822755 atol=3e-3      # ~0.36% relative
    @test total_melt ≈ 10372.968525 atol=120.0  # ~1.2% relative
    @test total_runoff ≈ 4468.185615 atol=120.0 # ~2.7% relative
end
