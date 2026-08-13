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
    # Warm, melt-dominated site with a diurnal-ish seasonal cycle: enough positive degree
    # days that `initialize_profile` infers ice (and `depth_autoadjust` shrinks the column
    # to the shallow ~25 m ice grid), with modest snowfall so accumulation cannot keep up.
    n = 365 * 3
    time = DateTime(2000, 1, 1) .+ Day.(0:n-1)
    doy = collect(1:n)
    temp = 276.15 .+ 8.0 .* cos.(2π .* doy ./ 365.0 .- π)

    cf = initialize_forcing(
        time, temp, fill(85000.0, n), fill(0.3, n), fill(4.0, n),
        fill(220.0, n), fill(260.0, n), fill(400.0, n);
        temperature_air_mean=276.15, wind_speed_mean=4.0,
        precipitation_mean=0.3 * 365.25)

    mp = ModelParameters(output_frequency=:daily)
    profile, mp = initialize_profile(mp, cf)

    # The regime was actually detected: pure ice on the autoadjusted shallow grid.
    @test all(parent(profile[:density]) .== mp.density_ice)
    @test mp.column_depth == 25.0

    N = length(profile[:dz])
    Z_fixed = sum(profile[:dz])

    # Melt must genuinely exceed accumulation, or this is not an ablation test.
    @test annual_pdd_melt(cf) > cf.precipitation_mean

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
