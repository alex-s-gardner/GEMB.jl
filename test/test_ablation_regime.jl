# Ablation-regime integration test.
#
# The rest of the suite exercises accumulation-dominated columns (`test_synthetic_regression`
# uses the accumulating `test_1` forcing). Ablation inverts the sign of both grid controllers
# and is a distinct code path worth covering end to end:
#
#   * mass leaves at the surface, so cells melt out and the *count* controller must **split**
#     deep cells to restore the fixed count (under accumulation it merges instead);
#   * the column shrinks, so the *mass* controller must **grow** the bottom cell — basal
#     accretion — giving a **positive** cumulative basal mass flux (negative under
#     accumulation).
#
# `verbose=true` throughout, so every per-timestep mass/energy check, the grid invariants,
# and the whole-run mass budget are enforced rather than merely exercised.

@testset "Ablation regime" begin
    # Warm, melt-dominated site with a seasonal cycle: enough melt that the net annual
    # mass balance is negative, with modest snowfall so accumulation cannot keep up.
    # `initialize_profile` therefore never buries anything and returns the exposed ice
    # column on a grid sized to the thermal wave rather than to nonexistent firn.
    n = 365 * 3
    time = DateTime(2000, 1, 1) .+ Day.(0:n-1)
    doy = collect(1:n)
    temp = 276.15 .+ 8.0 .* cos.(2π .* doy ./ 365.0 .- π)

    # Incoming longwave has to be realistic for a melting surface: at 260 W m-2 the
    # surface energy balance never reaches the melt point (εσT⁴ ≈ 306 W m-2 there),
    # so this site had *zero* melt despite 276 K air, and only registered as
    # ablating through a since-corrected rainfall term in the balance.
    cf = initialize_forcing(
        time, temp, fill(85000.0, n), fill(0.3, n), fill(4.0, n),
        fill(220.0, n), fill(310.0, n), fill(400.0, n);
        temperature_air_mean=276.15, wind_speed_mean=4.0,
        precipitation_mean=0.3 * 365.25)

    mp = ModelParameters(output_frequency=:daily)
    profile = initialize_profile(mp, cf)

    # The regime emerges as the degenerate limit of the marcher, not from a threshold:
    # nothing is buried, so the column is ice throughout.
    @test all(parent(profile[:density]) .== mp.density_ice)
    @test all(parent(profile[:water]) .== 0.0)

    # No firn to resolve, so the derived depth collapses to the annual-thermal-wave
    # floor — far shallower than the configured default, but not degenerate. The
    # derived depth is read off the grid, which is where it lives.
    @test sum(profile[:dz]) < 0.2 * mp.column_depth_max
    @test sum(profile[:dz]) >= mp.column_ztop

    # The escape hatch still yields the old exact pure-ice column on the *configured*
    # grid, which is what the MATLAB-fidelity sites rely on.
    # This site's mean annual air temperature is above freezing, so the fidelity path
    # clamps it: no initialized cell may start above the melt point, since the column
    # has no way to carry that enthalpy.
    prof_const = @test_logs (:warn, r"above the melt point") match_mode = :any (
        initialize_profile(mp, cf; constant_density=true, constant_temperature=true))
    @test parent(prof_const[:dz]) == GEMB.initialize_grid(mp)
    @test all(parent(prof_const[:density]) .== mp.density_ice)
    @test cf.temperature_air_mean > GEMB.CtoK
    @test all(parent(prof_const[:temperature]) .== GEMB.CtoK)

    N = length(profile[:dz])
    Z_fixed = sum(profile[:dz])

    # Net annual balance must genuinely be negative, or this is not an ablation test.
    cs = GEMB.initialize_climate_summary(cf, mp)
    @test cs.balance < 0.0
    @test cs.melt > cs.accumulation

    out = gemb(profile, cf, mp; verbose=true)

    # --- Fixed cell count, exactly ------------------------------------------------
    dz_out = parent(out[:dz])
    @test size(dz_out, 1) == N
    @test all(parent(out[:valid_profile_length]) .== N)
    @test !any(isnan, dz_out)

    # --- Fixed total depth, every output step -------------------------------------
    depths = vec(sum(dz_out; dims=1))
    @test all(d -> abs(d - Z_fixed) < 1e-9, depths)

    # --- Basal flux is positive (accretion) and plateaus ---------------------------
    # `thickness_cumulative` is the running basal mass flux divided by density_ice, so its
    # sign is the sign of cumulative `mass_added`: positive means mass entering at the base,
    # which is what a net-ablating column requires.
    thick = parent(out[:thickness_cumulative])
    @test all(isfinite, thick)
    @test thick[end] > 0.0

    # Net surface ablation, confirming the column really is losing mass at the top.
    net_surface = sum(parent(out[:precipitation])) - sum(parent(out[:runoff]))
    @test net_surface < 0.0

    # Basal inflow exceeds the net surface loss here (measured: +10106 vs -9502 kg m-2), the
    # intended consequence of pinning the depth rather than a leak — the budget closes under
    # `verbose=true` above. Pinned as a checked property, not just a docstring claim.
    basal_mass = thick[end] * mp.density_ice
    @test basal_mass > 0.0
    @test basal_mass > abs(net_surface)

    # --- The column stays physical -------------------------------------------------
    T = parent(out[:temperature])
    ρ = parent(out[:density])
    @test all(isfinite, T)
    @test all(T .<= GEMB.CtoK + 1e-6)          # nothing above the melt point
    @test all(ρ .> 0.0)
    @test all(ρ .<= mp.density_ice + 1e-6)     # an ablating column cannot exceed ice
    @test all(dz_out .> 0.0)

    # --- Spinup holds the same invariants -----------------------------------------
    # Spinup cycles the forcing repeatedly, so a controller interaction that drifts only
    # slowly per cycle shows up here and not in the single pass above.
    spun = gemb_spinup(profile, cf, mp; max_iterations=4, verbose=true)
    @test length(spun[:dz]) == N
    @test sum(parent(spun[:dz])) ≈ Z_fixed atol = 1e-9
end
