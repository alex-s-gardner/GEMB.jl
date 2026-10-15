# Full-model regression tests (self-consistency pins)
#
# The reference values below are GEMB.jl's own deterministic output, used to catch *unintended*
# changes to the full-model result. They are change detectors, not physics validations: an
# intended change to the defaults or to a physics scheme moves them, and re-pinning them, with
# the old and new values and the reason in the commit message, is part of making that change.
#
# Two sites:
#
# - `test_1` from GEMB_ClimateForcing.jl (3-hourly, 1994-2025, 75-cycle spin-up on its
#   climatology): a cold, dry accumulation site, T̄ = 249.4 K, where melt (≈ 3.6 kg m-2 yr-1)
#   is the small residual of the surface energy balance and all of it refreezes. Its runoff is
#   exactly zero and its rain is 0.05% of precipitation, so it sees the surface energy balance
#   and the radiation inputs but not the liquid-water paths. A change to `test_1` in
#   GEMB_ClimateForcing.jl moves these values without anything in GEMB.jl changing.
# - The maritime site of `synthetic_maritime_forcing.jl` (hourly, 4 years, no spin-up):
#   seasonal snow over ice with 1300 kg m-2 yr-1 of rain, 1600 of melt, 700 of refreeze and
#   2200 of runoff. It sees the rain, pore-water, refreeze, runoff and vapor-phase paths. The
#   forcing is closed-form and lives in this repository.
#
# The tolerances come from measured variability, on arm64 macOS:
#
# - Both sites are bit-for-bit reproducible run to run, and Julia 1.11.9 and 1.12.7 give
#   bit-identical values at both.
# - Raising every nonzero value of one forcing series by one ulp moves the pinned mass totals by
#   at most 4e-12 kg m-2 at `test_1` and 1.5e-10 at the maritime site, and the mean albedo by at
#   most 7e-16: neither site amplifies rounding differences.
# - The response stays linear up to a relative perturbation of the whole temperature series
#   of 1e-7 at `test_1`, and 1e-11 at the maritime site, where discrete events (cell merges,
#   phase switches) begin to flip at 1e-9 and add at most 7e-4 kg m-2 there.
#
# `mass_tol` is therefore max(1e-5 × total, 0.01 kg m-2) and the albedo tolerance 1e-5: more than
# ten times the largest of these, while a change of 0.01 kg m-2 to a `test_1` total, or of 1e-5
# relative to a maritime one, still fails.

using GEMB: Statistics
using GEMB_ClimateForcing

include("synthetic_maritime_forcing.jl")

const MASS_RTOL = 1e-5
const MASS_ATOL_FLOOR = 0.01    # kg m-2
const ALBEDO_ATOL = 1e-5
mass_tol(reference) = max(MASS_RTOL * abs(reference), MASS_ATOL_FLOOR)

@testset "Synthetic regression: test_1 (cold accumulation site)" begin
    ds = simulate_climate_forcing("test_1", 3)
    cf = initialize_forcing(ds)
    mp = ModelParameters(output_frequency=:daily)

    # Both escape-hatch flags give the pure-ice initialization and the deep column the
    # references were generated with.
    profile = initialize_profile(mp, cf; constant_density=true, constant_temperature=true)
    cf_climatology = GEMB.forcing_climatology(cf)
    mp_spinup = ModelParameters(output_frequency=:last)
    profile_spunup = gemb_spinup(profile, cf_climatology, mp_spinup;
        simulation_years_maximum=75 * GEMB._years_per_cycle(cf_climatology))
    output = gemb(profile_spunup, cf, mp)

    mean_albedo = Statistics.mean(parent(output[:albedo_broadband]))
    total_melt = sum(parent(output[:melt]))
    total_refreeze = sum(parent(output[:refreeze]))

    # The site retains all of its meltwater; `refreeze` stands in for `runoff`. If runoff ever
    # appears here this fails rather than leaving the pins watching the wrong field.
    @test sum(parent(output[:runoff])) == 0.0

    @test mean_albedo ≈ 0.870836 atol = ALBEDO_ATOL
    @test total_melt ≈ 114.669238 atol = mass_tol(114.669238)
    @test total_refreeze ≈ 119.610019 atol = mass_tol(119.610019)
end

@testset "Synthetic regression: maritime site (rain, refreeze, runoff)" begin
    cf = synthetic_maritime_forcing()
    mp = ModelParameters(output_frequency=:daily)
    profile = initialize_profile(mp, cf)
    output = gemb(profile, cf, mp)
    total(k) = sum(parent(output[k]))

    # The premise of this site: rain and runoff are each a large part of the budget.
    @test total(:rain) > 0.4 * total(:precipitation)
    @test total(:runoff) > 0.5 * total(:precipitation)
    @test total(:refreeze) > 0.2 * total(:precipitation)

    @test Statistics.mean(parent(output[:albedo_broadband])) ≈ 0.798615 atol = ALBEDO_ATOL
    @test total(:melt) ≈ 6443.665400 atol = mass_tol(6443.665400)
    @test total(:rain) ≈ 5195.191007 atol = mass_tol(5195.191007)
    @test total(:refreeze) ≈ 2868.168168 atol = mass_tol(2868.168168)
    @test total(:runoff) ≈ 8849.323336 atol = mass_tol(8849.323336)
    @test total(:evaporation_condensation) ≈ 139.513451 atol = mass_tol(139.513451)
end

@testset "Synthetic regression: tolerance policy" begin
    # The window scales with the total and keeps a floor if a total goes to zero.
    @test mass_tol(114.669238) == MASS_ATOL_FLOOR
    @test mass_tol(6443.665400) ≈ 0.0644366540 rtol = 1e-12
    @test mass_tol(0.0) == MASS_ATOL_FLOOR
end
