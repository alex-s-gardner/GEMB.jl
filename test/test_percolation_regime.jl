# Percolation-regime integration test.
#
# `initialize_climate_summary` documents the sign of `balance` as "the only regime
# discriminator", but only the negative branch had end-to-end coverage
# (`test_ablation_regime.jl`). This covers the case that is neither branch's easy end: a site
# that is **net accumulating** and yet melts, so meltwater has somewhere to go and the
# percolation/refreeze path runs inside a growing firn column rather than on bare ice.
#
# That combination is what the unit tests cannot reach. `test_calculate_melt.jl` and
# `test_hydrology.jl` cover the bucket scheme itself, including a multi-timestep aquifer, but
# they drive it directly. Here the water comes from a surface energy balance, percolates
# through firn whose density the run itself produced, and refreezes against cold content the
# thermal solver is tracking — with both grid controllers active throughout.
#
# `verbose=true` throughout, so the per-timestep mass/energy checks, the grid invariants, and
# the whole-run mass budget are enforced rather than merely exercised.

@testset "Percolation regime" begin
    # Warm-but-accumulating site: 2 kg m-2 d-1 of precipitation against a mean of 263 K with a
    # 12 K seasonal amplitude, so summer crosses the melt point while the annual balance stays
    # firmly positive. Incoming longwave of 291 W m-2 is what makes the surface actually reach
    # melting — see the note in `test_ablation_regime.jl` on why it has to be realistic.
    n = 365 * 3
    time = DateTime(2000, 1, 1) .+ Day.(0:n-1)
    doy = collect(1:n)
    temp = 263.0 .+ 12.0 .* cos.(2π .* doy ./ 365.0 .- π)

    cf = initialize_forcing(
        time, temp, fill(75000.0, n), fill(2.0, n), fill(4.0, n),
        fill(140.0, n), fill(291.0, n), fill(350.0, n);
        temperature_air_mean=263.0, wind_speed_mean=4.0,
        precipitation_mean=2.0 * 365.25)

    mp = ModelParameters(output_frequency=:daily)
    profile = initialize_profile(mp, cf)

    # --- The regime, as the summary sees it ----------------------------------------
    # Net accumulating: this is the *positive* branch of the discriminator, unlike
    # `test_ablation_regime.jl`. Firn is buried, so the column is graded rather than pure ice.
    cs = GEMB.initialize_climate_summary(cf, mp)
    @test cs.balance > 0.0
    @test cs.accumulation > 0.0
    @test any(<(mp.density_ice), parent(profile[:density]))

    # Some of the precipitation falls as rain, which is the other liquid source the
    # percolation path has to handle alongside melt.
    @test cs.rainfall > 0.0

    N = length(profile[:dz])
    Z_fixed = sum(profile[:dz])

    out = gemb(profile, cf, mp; verbose=true)

    # --- Melt happens, and most of it is retained ---------------------------------
    # The point of the regime: liquid water is produced *and* the column can absorb it.
    # Refreeze exceeds runoff here (measured: 486 vs 236 kg m-2), which is the property that
    # distinguishes percolation into cold firn from melt running off bare ice.
    melt = sum(parent(out[:melt]))
    refreeze = sum(parent(out[:refreeze]))
    runoff = sum(parent(out[:runoff]))
    rain = sum(parent(out[:rain]))
    @test melt > 0.0
    @test rain > 0.0
    @test refreeze > 0.0
    @test refreeze > runoff
    # Refreeze cannot exceed the liquid available to it. It legitimately exceeds *melt* alone
    # (486 > 369) because rain is the other source, so the bound is on the sum.
    @test refreeze <= melt + rain + 1e-6

    # --- The water is in the column, at depth --------------------------------------
    # `percolation_depth` is a diagnostic that takes no part in the budget, so it is pinned
    # here as a claim about where water went rather than assumed from the mass balance.
    water = parent(out[:water])
    perc_depth = parent(out[:percolation_depth])
    @test any(>(0.0), water)
    @test maximum(perc_depth) > 1.0            # measured 9.49 m; well past the surface cell
    @test all(d -> 0.0 <= d <= Z_fixed, perc_depth)

    # Water sits in pores, so it can never exceed the pore volume of the cell holding it.
    dz_out = parent(out[:dz])
    density = parent(out[:density])
    pore_capacity = @. dz_out * (mp.density_ice - min(density, mp.density_ice))
    @test all(water .<= pore_capacity .+ 1e-6)

    # --- Firn air content falls, and melt is part of why ---------------------------
    # Under accumulation alone a column loses air only by compaction. Here the melt term is
    # separately reported and non-zero, so the FAC loss is not compaction misread as melt.
    fac = parent(out[:firn_air_content])
    @test fac[end] < fac[1]                                        # 3.69 -> 2.78 m
    @test sum(parent(out[:densification_from_melt])) > 0.0
    @test sum(parent(out[:densification_from_compaction])) > 0.0

    # --- Fixed cell count and fixed total depth ------------------------------------
    @test size(dz_out, 1) == N
    @test !any(isnan, dz_out)
    depths = vec(sum(dz_out; dims=1))
    @test all(d -> abs(d - Z_fixed) < 1e-9, depths)

    # --- The elevation identity, with pore-water storage ---------------------------
    # The ablation form of the identity,
    #
    #     -cumsum(ice_flux) = SMB / density_ice + Δ(firn_air_content)
    #
    # is a *special case*: it assumes the column stores no liquid water, which is true of the
    # bare-ice site in `test_ablation_regime.jl` and false here. Water retained in pores adds
    # mass to `SMB` without adding thickness — `firn_air_content` counts only the ice
    # fraction, so a pore filling with water leaves both `Σdz` and FAC unchanged. The mass
    # term must therefore be the SMB the *matrix* received, i.e. less the change in column
    # water storage:
    #
    #     -cumsum(ice_flux) = (SMB − ΔW) / density_ice + Δ(firn_air_content)
    #
    # Without the ΔW term the residual here is -0.0395 m against 37.6 kg m-2 of stored water;
    # with it, 0.0016 m, the same half-interval offset the ablation test tolerates.
    water_stored = sum(@view water[:, end]) - sum(@view water[:, 1])
    @test water_stored > 0.0
    smb_ie = (sum(parent(out[:precipitation])) +
              sum(parent(out[:evaporation_condensation])) -
              sum(parent(out[:runoff])) - water_stored) / mp.density_ice
    cum_flux = cumsum(parent(out[:ice_flux]))
    @test all(isfinite, parent(out[:ice_flux]))
    # `firn_air_content` is an interval mean and `ice_flux` an interval sum, so the two are
    # offset by half an output interval; the tolerance covers that, not a physical slop.
    @test -cum_flux[end] ≈ smb_ie + (fac[end] - fac[1]) atol = 1e-2

    # --- The column stays physical -------------------------------------------------
    T = parent(out[:temperature])
    @test all(isfinite, T)
    @test all(T .<= GEMB.CtoK + 1e-6)          # refreezing cannot push firn past melting
    @test all(density .> 0.0)
    @test all(density .<= mp.density_ice + 1e-6)
    @test all(dz_out .> 0.0)

    # A percolating column reaches the melt point somewhere — otherwise refreeze would be
    # limited by cold content that never ran out, and the bucket scheme's saturated branch
    # would go untested.
    @test maximum(T) ≈ GEMB.CtoK atol = 1e-6

    # --- Spinup holds the same invariants -----------------------------------------
    # Cycling the forcing exposes a controller or storage interaction that drifts too slowly
    # per cycle to show up in the single pass above.
    spun = gemb_spinup(profile, cf, mp; max_iterations=4, verbose=true)
    @test length(spun[:dz]) == N
    @test sum(parent(spun[:dz])) ≈ Z_fixed atol = 1e-9
    @test all(parent(spun[:water]) .>= 0.0)
end
