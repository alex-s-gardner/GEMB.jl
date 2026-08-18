# Full synthetic data regression test (self-consistency pin)
#
# The reference values below are GEMB.jl's own deterministic output, used to catch
# *unintended* changes to the full-model result. They are bit-for-bit reproducible
# across runs on a given platform; small cross-platform spread remains because the
# 75-cycle spinup (~7.9M iterations) amplifies floating-point evaluation-order
# differences (arm64 vs x86_64), so the tolerances are modest rather than 1e-12.
#
# Note this pin depends on the `test_1` forcing from GEMB_ClimateForcing.jl, so a change to
# that site's parameters moves these numbers without anything in GEMB.jl changing. The most
# recent re-pin (see the history above the assertions) was exactly that.
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

    # The three pinned fields. `refreeze` replaced `runoff` when `test_1` moved to 2200 m:
    # the site now retains all of its meltwater, so `runoff` is identically zero and pinning
    # it would assert nothing. `refreeze` covers the same melt/percolation path with signal
    # in it (402.6 kg m-2), which is the property the third pin was there for.
    mean_albedo = Statistics.mean(parent(output[:albedo_broadband]))
    total_melt = sum(parent(output[:melt]))
    total_refreeze = sum(parent(output[:refreeze]))

    # Guard the substitution's premise: if a change ever puts runoff back at this site, the
    # third pin is watching the wrong field and this fails rather than passing quietly.
    @test sum(parent(output[:runoff])) == 0.0

    # Reference = deterministic output on arm64 (macOS). Tolerances bracket
    # the cross-platform FP evaluation-order spread amplified over the 75-cycle
    # spinup (~7.9M iterations): measured x86_64 CI vs arm64 is ~51 kg/m² for melt,
    # against the ~10000 kg/m² totals of the time — i.e. ~0.49% — and the windows are
    # set ~2.4x that. The mass windows are expressed relatively for exactly this
    # reason, since the spread scales with the totals; see the note at the assertions.
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
    # Re-pinned when the `test_1` synthetic site moved from 700 m to 2200 m elevation in
    # GEMB_ClimateForcing.jl. This is a change to the *forcing*, not to GEMB's physics: at
    # 700 m the site melted hard (332 kg m-2 yr-1, 3.9% of steps at or above freezing), and
    # melt is convex in temperature, so `forcing_climatology` — which averages the 32 years
    # into one — cancelled the warm excursions that carried it and produced a climatology
    # with *zero* melt. The spinup therefore equilibrated a melt-free column (FAC ~4.06 m)
    # and the transient run spent its whole length relaxing to the real attractor
    # (FAC ~0.74 m), so the docs' diagnostic figure showed a 7 m elevation drift that was an
    # artifact of the spinup/transient mismatch rather than a property of the site. At
    # 2200 m melt survives averaging (12.5 kg m-2 yr-1) and a spun-up column holds its
    # elevation: FAC drifts +0.17 m over 32 years against -5.57 m from a cold start.
    #
    # Every number below moves by far more than any tolerance, because it is a different
    # climate: melt 10372.97 -> 397.22 kg m-2 (-96%), albedo 0.822755 -> 0.835208 (the colder,
    # drier surface keeps finer grains), and runoff 4468.19 -> 0.0, which is why the third
    # pin is now `refreeze`. Note the fingerprint is now a *colder* site than before, so it
    # exercises the melt/refreeze path more lightly; `test_ablation_regime.jl` remains the
    # melt-dominated end-to-end case.
    #
    # The mass tolerances became relative-with-a-floor as a consequence of this re-pin: an
    # absolute window sized against ~10000 kg m-2 totals does not transfer to ~400. See the
    # note at the assertions.
    # The mass tolerances are `max(rtol·reference, floor)` rather than a bare `atol`, because
    # the quantity they have to bracket — cross-platform FP evaluation-order spread amplified
    # over the 75-cycle spinup — scales with the totals, while a fixed `atol` does not. The
    # inherited `atol = 120.0` was sized against ~10000 kg m-2 totals, where it was 1.16%
    # relative and 2.35x the measured arm64/x86_64 spread of ~51 kg m-2 (0.49%). Carried
    # unchanged onto the ~400 kg m-2 totals this site produces, the same window is ~30% —
    # loose enough to miss a real regression. The relative form keeps the ratio the pin was
    # actually sized for: `MASS_RTOL = 0.012` is 2.4x the spread expressed as a fraction, so
    # it stays correctly sized if a future re-pin moves these totals again in either
    # direction.
    #
    # The floor is what a pure `rtol` would get wrong at the other end. `refreeze` is one
    # substitution away from having been `runoff`, which this site drove to exactly zero;
    # should a future change take either total near zero, `rtol·reference` collapses with it
    # and the pin becomes unfalsifiable-by-arithmetic — asserting agreement to a window
    # narrower than the FP noise it exists to tolerate. `MASS_ATOL_FLOOR = 2.0` kg m-2 holds a
    # usable window there; it binds only below ~167 kg m-2, so it is inert at the current
    # references (4.77 and 4.83 kg m-2 respectively).
    MASS_RTOL = 0.012
    MASS_ATOL_FLOOR = 2.0
    mass_tol(reference) = max(MASS_RTOL * abs(reference), MASS_ATOL_FLOOR)

    @test mean_albedo ≈ 0.835208 atol=3e-3      # ~0.36% relative
    @test total_melt ≈ 397.216343 atol=mass_tol(397.216343)
    @test total_refreeze ≈ 402.605404 atol=mass_tol(402.605404)

    # Pin the tolerance policy itself, so a later edit cannot quietly reintroduce the failure
    # mode either half exists to prevent: a window that does not scale, or one that collapses
    # to zero. These are assertions about the *test*, which is why they are cheap and here.
    @test mass_tol(397.216343) ≈ 4.7666 atol=1e-3     # rtol binds at the current references
    @test mass_tol(0.0) == MASS_ATOL_FLOOR            # floor binds if a total ever goes to 0
    @test mass_tol(1e5) > MASS_ATOL_FLOOR             # and scales up with a larger total

    # Still worth measuring rather than inferring: 0.012 is the *old* site's spread-to-window
    # ratio carried across, not a spread measured at 2200 m. The references above are exact on
    # arm64, so an x86_64 CI run gives the real figure for this site directly; if it comes in
    # materially below 0.49%, this can tighten.
end
