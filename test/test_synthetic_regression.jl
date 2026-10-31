# Full-model regression tests (self-consistency pins)
#
# The reference values below are GEMB.jl's own deterministic output, used to catch *unintended*
# changes to the full-model result. They are change detectors, not physics validations: an
# intended change to the defaults or to a physics scheme moves them, and re-pinning them, with
# the old and new values and the reason in the commit message, is part of making that change.
#
# Three sites:
#
# - `test_1` from GEMB_ClimateForcing.jl (3-hourly, 1994-2025, 75-cycle spin-up on its
#   climatology): a cold, dry accumulation site, T̄ = 249.4 K, where melt (≈ 10.9 kg m-2 yr-1)
#   is the small residual of the surface energy balance and all of it refreezes. Its runoff is
#   exactly zero and its rain is 0.05% of precipitation, so it sees the surface energy balance
#   and the radiation inputs but not the liquid-water paths. A change to `test_1` in
#   GEMB_ClimateForcing.jl moves these values without anything in GEMB.jl changing.
# - The maritime site of `synthetic_maritime_forcing.jl` (hourly, 4 years, no spin-up):
#   seasonal snow over ice with 1300 kg m-2 yr-1 of rain, 1386 of melt, 673 of refreeze and
#   2033 of runoff, under a varying solar zenith angle and cloud optical thickness. It sees the
#   rain, pore-water, refreeze, runoff and vapor-phase paths and the albedo's zenith and cloud
#   terms. The forcing is closed-form and lives in this repository.
# - The thin-snow site of `synthetic_thin_snow_forcing.jl` (hourly, 90 spring days at 72° N, no
#   spin-up, `albedo_ice = 0.30`): repeated light snowfall on bare dark ice that coarsens and
#   melts off within days, with the sun at least 50° from the zenith and a cloud optical
#   thickness of at least 3. Thin coarse snow over ice puts the albedo on the `A < 1` branch of
#   Gardner & Sharp (2010) eq. 13, which the other two sites barely reach. 1058 kg m-2 of melt
#   and 1111 of runoff from 253 of precipitation.
#
# The tolerances come from measured variability, on arm64 macOS with Julia 1.12.7:
#
# - All three sites are bit-for-bit reproducible run to run.
# - Raising every value of the air-temperature series by one ulp moves the pinned mass totals
#   by at most 5e-12 kg m-2 at `test_1`, 1.8e-10 at the maritime site and 4e-12 at the
#   thin-snow site, and the mean albedo by at most 7e-16: no site amplifies rounding differences.
# - A relative perturbation of the whole temperature series of 1e-9 moves the totals by at most
#   6.7e-6 kg m-2 at `test_1`, 3.2e-4 at the maritime site and 5.8e-6 at the thin-snow site,
#   linearly at the last two (at the thin-snow site 1e-7 still moves them by only 6.4e-4). At 1e-7, discrete events (new layers, merges, phase switches) flip at
#   `test_1` and move its melt by 0.55 kg m-2.
#
# `mass_tol` is therefore max(1e-5 × total, 0.01 kg m-2) and the albedo tolerance 1e-5: more than
# ten times the largest of these, while a change of 0.01 kg m-2 to a `test_1` total, or of 1e-5
# relative to a maritime one, still fails.

using GEMB: Statistics
using GEMB_ClimateForcing

include("synthetic_maritime_forcing.jl")
include("synthetic_thin_snow_forcing.jl")

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

    @test mean_albedo ≈ 0.851263 atol = ALBEDO_ATOL
    @test total_melt ≈ 348.555912 atol = mass_tol(348.555912)
    @test total_refreeze ≈ 347.259360 atol = mass_tol(347.259360)
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

    @test Statistics.mean(parent(output[:albedo_broadband])) ≈ 0.856312 atol = ALBEDO_ATOL
    @test total(:melt) ≈ 5538.731369 atol = mass_tol(5538.731369)
    @test total(:rain) ≈ 5195.191007 atol = mass_tol(5195.191007)
    @test total(:refreeze) ≈ 2699.522068 atol = mass_tol(2699.522068)
    @test total(:runoff) ≈ 8122.407532 atol = mass_tol(8122.407532)
    @test total(:evaporation_condensation) ≈ 158.119231 atol = mass_tol(158.119231)
end

@testset "Synthetic regression: thin snow over dark ice at low sun" begin
    cf = synthetic_thin_snow_forcing()
    mp = ModelParameters(output_frequency=:daily, albedo_ice=0.30)
    # The mean air temperature is above the melting point, so the ice starts at 0 °C.
    profile = @test_logs (:warn, r"above the melt point") initialize_profile(mp, cf;
        constant_density=true, constant_temperature=true)
    output = gemb(profile, cf, mp)
    total(k) = sum(parent(output[k]))

    # The premise of this site: the snow melts off and the ice beneath melts too.
    @test total(:melt) > total(:precipitation)

    @test Statistics.mean(parent(output[:albedo_broadband])) ≈ 0.736489 atol = ALBEDO_ATOL
    @test total(:melt) ≈ 1058.360946 atol = mass_tol(1058.360946)
    @test total(:refreeze) ≈ 100.362076 atol = mass_tol(100.362076)
    @test total(:runoff) ≈ 1111.273253 atol = mass_tol(1111.273253)
end

@testset "Synthetic regression: tolerance policy" begin
    # The window scales with the total and keeps a floor if a total goes to zero.
    @test mass_tol(127.969237) == MASS_ATOL_FLOOR
    @test mass_tol(5376.685582) ≈ 0.05376685582 rtol = 1e-12
    @test mass_tol(0.0) == MASS_ATOL_FLOOR
end
