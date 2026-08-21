# Accumulation-regime integration test — the mirror image of `test_ablation_regime.jl`.
#
# Cold and dry: no melt at all, so the only processes moving the column are burial and
# compaction. That inverts every sign the ablation test pins:
#
#   * mass arrives at the surface, so the *count* controller must **merge** deep cells to
#     hold the fixed count (it splits under ablation);
#   * the column grows, so the *mass* controller must **shorten** the bottom cell — mass
#     leaving at the base — giving a **negative** cumulative basal flux (positive under
#     ablation), at every single step rather than only in the net.
#
# This is the branch `test_synthetic_regression.jl` used to cover incidentally, and covers
# less well since `test_1` moved to a colder site; it was never asserted *as* a regime. The
# value of the dry limit is that the melt, percolation, and refreeze terms are all exactly
# zero, so any mass or energy the column gains or loses has nowhere to hide.
#
# `verbose=true` throughout, so the per-timestep mass/energy checks, the grid invariants, and
# the whole-run mass budget are enforced rather than merely exercised.

@testset "Accumulation regime" begin
    # Cold interior site: mean 248 K with a 15 K seasonal amplitude, 0.6 kg m-2 d-1 of
    # snowfall, and a low longwave and shortwave load, so the surface never approaches
    # melting even at the summer peak (measured maximum column temperature: 252.99 K).
    n = 365 * 3
    time = DateTime(2000, 1, 1) .+ Day.(0:n-1)
    doy = collect(1:n)
    temp = 248.0 .+ 15.0 .* cos.(2π .* doy ./ 365.0 .- π)

    cf = initialize_forcing(
        time, temp, fill(75000.0, n), fill(0.6, n), fill(4.0, n),
        fill(110.0, n), fill(185.0, n), fill(60.0, n);
        temperature_air_mean=248.0, wind_speed_mean=4.0,
        precipitation_mean=0.6 * 365.25)

    mp = ModelParameters(output_frequency=:daily)
    profile = initialize_profile(mp, cf)

    # --- The regime, as the summary sees it ----------------------------------------
    # The positive branch of the discriminator in its clean form: all accumulation, and the
    # summary's own energy balance finds no melt to refreeze.
    cs = GEMB.initialize_climate_summary(cf, mp)
    @test cs.balance > 0.0
    @test cs.accumulation > 0.0
    @test cs.melt == 0.0
    @test cs.refreeze == 0.0
    @test cs.rainfall == 0.0                   # every step is well below the rain threshold
    @test cs.balance ≈ cs.accumulation

    # Everything precipitates as snow and is buried, so the column is graded firn, not the
    # pure ice the ablation site initializes to.
    @test any(<(mp.density_ice), parent(profile[:density]))
    @test all(parent(profile[:water]) .== 0.0)

    N = length(profile[:dz])
    Z_fixed = sum(profile[:dz])

    out = gemb(profile, cf, mp; verbose=true)

    # --- The dry limit, exactly ----------------------------------------------------
    # Not "small" — zero. Every liquid-water term is off, which is what makes the mass budget
    # below a statement about burial and compaction alone.
    @test sum(parent(out[:melt])) == 0.0
    @test sum(parent(out[:runoff])) == 0.0
    @test sum(parent(out[:refreeze])) == 0.0
    @test sum(parent(out[:rain])) == 0.0
    @test all(==(0.0), parent(out[:water]))
    @test all(==(0.0), parent(out[:percolation_depth]))
    @test sum(parent(out[:densification_from_melt])) == 0.0
    # Compaction is therefore the *only* densification mechanism acting.
    @test sum(parent(out[:densification_from_compaction])) > 0.0

    # --- Fixed cell count and fixed total depth ------------------------------------
    dz_out = parent(out[:dz])
    @test size(dz_out, 1) == N
    @test !any(isnan, dz_out)
    depths = vec(sum(dz_out; dims=1))
    @test all(d -> abs(d - Z_fixed) < 1e-9, depths)

    # --- Basal flux is negative at every step, not merely in net --------------------
    # The inverse of the ablation test's assertion, and stronger than its cumulative form:
    # this site accumulates on every timestep, so the mass controller must remove mass at the
    # base on every timestep. A single positive entry would mean the column momentarily
    # shrank, which under monotone burial with no melt it cannot.
    flux = parent(out[:ice_flux])
    @test all(isfinite, flux)
    @test all(<(0.0), flux)
    cum_flux = cumsum(flux)
    @test cum_flux[end] < 0.0
    @test issorted(cum_flux; rev=true)         # monotonically decreasing, by the above

    # Net surface accumulation, confirming the column really is gaining mass at the top.
    net_surface = sum(parent(out[:precipitation])) +
                  sum(parent(out[:evaporation_condensation])) -
                  sum(parent(out[:runoff]))
    @test net_surface > 0.0

    # Mass leaves at the base to make room for it — the intended consequence of pinning the
    # depth, not a leak. The whole-run budget closes under `verbose=true` above.
    #
    # The two do not cancel, and the *sign* of the residual is set by which way the column's
    # firn air content moved, not by a fixed inequality. Rearranging the elevation identity
    # below (`-cum_flux = SMB/ρᵢ + ΔFAC`) and multiplying through by `ρᵢ`:
    #
    #     |basal_mass| − net_surface = ΔFAC · density_ice
    #
    # so a column whose FAC shrinks loses *less* at the base than it gains at the top, and one
    # whose FAC grows loses more. Both are physical. Measured here: ΔFAC = −0.0502 m, i.e.
    # −46.0 kg m-2, against |basal_mass| − net_surface = −45.7 — which is the real invariant.
    # Asserting `abs(basal_mass) > net_surface` instead was pinning the sign of ΔFAC for this
    # particular forcing, which is not a property of the accumulation regime.
    basal_mass = cum_flux[end] * mp.density_ice
    fac = parent(out[:firn_air_content])
    @test basal_mass < 0.0
    @test abs(basal_mass) - net_surface ≈ (fac[end] - fac[1]) * mp.density_ice atol = 1.0

    # --- The elevation identity ----------------------------------------------------
    # Same identity as the ablation test, and it takes its simple form here because the column
    # stores no liquid water (see `test_percolation_regime.jl` for the ΔW term a wet column
    # needs):
    #
    #     -cumsum(ice_flux) = SMB / density_ice + Δ(firn_air_content)
    #
    # A dry accumulating column is the cleanest case for it: the surface *rises* here
    # (-cum_flux > 0), and the mass term carries most of it (+0.760 m of the +0.710 m total,
    # with the FAC term contributing -0.050 m as the column compacts).
    smb_ie = net_surface / mp.density_ice
    @test -cum_flux[end] ≈ smb_ie + (fac[end] - fac[1]) atol = 1e-2
    @test -cum_flux[end] > 0.0                 # the surface rises

    # --- The column stays physical -------------------------------------------------
    T = parent(out[:temperature])
    density = parent(out[:density])
    @test all(isfinite, T)
    @test all(T .< GEMB.CtoK)                  # strictly below melting — nothing ever melts
    @test all(density .> 0.0)
    @test all(density .<= mp.density_ice + 1e-6)
    @test all(dz_out .> 0.0)

    # Burial with no melt gives a densifying column: every cell is at least as dense as the one
    # above it, at the final step.
    #
    # The tolerance is not slop. Fresh snow arrives at a density set by the air temperature at
    # the moment of snowfall, which swings 30 K over the seasonal cycle here, so cells buried in
    # different seasons enter the column at different densities and the compaction that follows
    # does not always erase the ordering. Measured: 6 inversions in 264 cells, worst
    # -0.043 kg m-3 — a relative 1.1e-4, three orders of magnitude below the ~50 kg m-3 the
    # fresh-snow density itself varies by across a year. A strict `issorted` asserts that
    # seasonal density banding is fully overprinted by depth, which is not a property of the
    # regime.
    final_density = density[:, end]
    @test all(diff(final_density) .>= -0.1)
    # And densifying overall by far more than the banding amplitude, which is the physical
    # claim the strict sort was standing in for.
    @test final_density[end] - final_density[1] > 100.0

    # --- Spinup holds the same invariants -----------------------------------------
    spun = gemb_spinup(profile, cf, mp; max_iterations=4, verbose=true)
    @test length(spun[:dz]) == N
    @test sum(parent(spun[:dz])) ≈ Z_fixed atol = 1e-9
    @test all(parent(spun[:water]) .== 0.0)
end
