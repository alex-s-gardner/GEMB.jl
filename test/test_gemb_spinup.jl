using Test
using GEMB
using Dates

@testset "gemb_spinup" begin
    @testset "Basic spinup execution" begin
        # Create parameters
        params = initialize_parameters(
            densification_method = :Arthern,
            output_frequency = :last
        )

        # Create 1-year climatology
        n_days = 365
        start_time = DateTime(2020, 1, 1)
        time = start_time .+ Day.(0:n_days-1)

        # Simple seasonal forcing
        day_of_year = collect(1:n_days)
        temp_seasonal = 260.0 .+ 15.0 .* cos.(2π .* day_of_year ./ 365.0 .- π)

        forcing = initialize_forcing(
            time,
            temp_seasonal,  # temperature_air
            fill(85000.0, n_days),  # pressure_air
            fill(1.0, n_days),  # precipitation
            fill(5.0, n_days),  # wind_speed
            fill(100.0, n_days),  # shortwave_downward
            fill(200.0, n_days),  # longwave_downward
            fill(100.0, n_days),  # vapor_pressure
            temperature_observation_height = 2.0,
            wind_observation_height = 10.0
        )

        # Initialize profile
        profile = initialize_profile(params, forcing)

        # Run spinup with 3 cycles (fast test)
        output = gemb_spinup(profile, forcing, params; max_iterations=3, verbose=false)

        # Check output structure
        @test output isa DimStack
        @test haskey(output, :temperature)
        @test haskey(output, :density)
        @test haskey(output, :dz)

        # Check that we get a profile with Z dimension
        @test length(output[:dz]) > 0
    end

    @testset "Spinup convergence test" begin
        # Test that running more cycles leads to more stable profiles

        params = initialize_parameters(
            output_frequency = :last
        )

        # Annual climatology
        n_days = 365
        start_time = DateTime(2020, 1, 1)
        time = start_time .+ Day.(0:n_days-1)

        forcing = initialize_forcing(
            time,
            fill(255.0, n_days),  # temperature_air - Cold, stable
            fill(85000.0, n_days),  # pressure_air
            fill(0.5, n_days),  # precipitation
            fill(3.0, n_days),  # wind_speed
            fill(50.0, n_days),  # shortwave_downward
            fill(180.0, n_days),  # longwave_downward
            fill(80.0, n_days),  # vapor_pressure
            temperature_observation_height = 2.0,
            wind_observation_height = 10.0
        )

        # Run with different numbers of cycles
        profile = initialize_profile(params, forcing)
        output_3 = gemb_spinup(profile, forcing, params; max_iterations=3)
        output_5 = gemb_spinup(profile, forcing, params; max_iterations=5)

        # Extract final profiles
        temp_3 = output_3[:temperature][Ti=1]
        temp_5 = output_5[:temperature][Ti=1]

        # Both should be valid
        @test all(isfinite.(temp_3[.!isnan.(temp_3)]))
        @test all(isfinite.(temp_5[.!isnan.(temp_5)]))

        # Longer spinup should produce deeper profiles
        n_layers_3 = sum(.!isnan.(temp_3))
        n_layers_5 = sum(.!isnan.(temp_5))
        @test n_layers_5 >= n_layers_3
    end

    @testset "Profile extraction after spinup" begin
        # Test that gemb_spinup returns a valid profile

        params = initialize_parameters(
            output_frequency = :last
        )

        n_days = 365
        start_time = DateTime(2020, 1, 1)
        time = start_time .+ Day.(0:n_days-1)

        forcing = initialize_forcing(
            time,
            fill(260.0, n_days),  # temperature_air
            fill(85000.0, n_days),  # pressure_air
            fill(1.0, n_days),  # precipitation
            fill(5.0, n_days),  # wind_speed
            fill(100.0, n_days),  # shortwave_downward
            fill(200.0, n_days),  # longwave_downward
            fill(100.0, n_days),  # vapor_pressure
            temperature_observation_height = 2.0,
            wind_observation_height = 10.0
        )

        profile = initialize_profile(params, forcing)
        spunup_profile = gemb_spinup(profile, forcing, params; max_iterations=3)

        # Check profile has required fields
        @test haskey(spunup_profile, :temperature)
        @test haskey(spunup_profile, :density)
        @test haskey(spunup_profile, :dz)
        @test haskey(spunup_profile, :grain_radius)

        # Profile arrays should have same length
        n_layers = length(spunup_profile[:dz])
        @test length(spunup_profile[:density]) == n_layers
        @test length(spunup_profile[:temperature]) == n_layers
    end

    @testset "steady_state_density: Herron-Langway" begin
        mp_hl = initialize_parameters(densification_method=:HerronLangway)
        zc = -collect(0.05:0.1:60.0)   # cell-center depths, negative below surface
        ρ0, ρi = 350.0, 910.0
        d = steady_state_density(zc, 253.0, 300.0, ρ0, ρi, mp_hl)

        # Surface ≈ ρ0, deep asymptote → ρ_ice, monotonic non-decreasing.
        @test isapprox(d[1], ρ0; atol=2.0)
        @test all(diff(d) .>= -1e-9)
        @test all(d .<= ρi + 1e-9)
        d_deep = steady_state_density(-collect(0.5:1.0:250.0), 253.0, 300.0, ρ0, ρi, mp_hl)
        @test d_deep[end] > 900.0

        # 550 kg/m³ crossover at a physically sensible depth (a few m to tens of m).
        i55 = findfirst(>=(550.0), d)
        crossover_depth = -zc[i55]
        @test 2.0 < crossover_depth < 40.0

        # Negligible accumulation → no firn forms → pure ice.
        d_ice = steady_state_density(zc, 253.0, 0.0, ρ0, ρi, mp_hl)
        @test all(d_ice .== ρi)

        # The densification method is genuinely plumbed through: Arthern (2010) is a
        # different law and must give a different profile at identical forcing.
        mp_ar = initialize_parameters(densification_method=:Arthern)
        d_ar = steady_state_density(zc, 253.0, 300.0, ρ0, ρi, mp_ar)
        @test !isapprox(d_ar, d; rtol=1e-3)
        @test all(diff(d_ar) .>= -1e-9)
        @test all(d_ar .<= ρi + 1e-9)
    end

    @testset "Ablation regime → ice column + clamped temperature" begin
        n = 365
        time = DateTime(2020, 1, 1) .+ Day.(0:n-1)
        # Warm, melt-dominated site: net annual balance is negative, so nothing is
        # ever buried and the column is the ice it exposes.
        cf = initialize_forcing(
            time, fill(278.15, n), fill(85000.0, n), fill(1.0, n), fill(5.0, n),
            fill(100.0, n), fill(200.0, n), fill(100.0, n);
            temperature_air_mean=278.15, wind_speed_mean=5.0, precipitation_mean=365.25)
        mp = initialize_parameters()

        cs = GEMB.initialize_climate_summary(cf, mp)
        @test cs.balance <= 0.0

        prof = initialize_profile(mp, cf)
        @test all(parent(prof[:density]) .== mp.density_ice)
        @test all(parent(prof[:temperature]) .<= 273.15 + 1e-9)

        # No firn to resolve, so the derived column collapses to the thermal floor:
        # deep enough for the annual wave, far shallower than the 250 m default. The
        # derived depth is read off the grid, which is where it lives.
        @test sum(prof[:dz]) < 0.2 * mp.column_depth_max
        @test sum(prof[:dz]) >= mp.column_ztop
    end

    @testset "Steady-state init: firn profile + faster convergence" begin
        n = 365
        time = DateTime(2020, 1, 1) .+ Day.(0:n-1)
        doy = collect(1:n)
        # Cold, dry, accumulating site (dry-snow zone).
        temp_seasonal = 253.0 .+ 10.0 .* cos.(2π .* doy ./ 365.0 .- π)
        cf = initialize_forcing(
            time, temp_seasonal, fill(85000.0, n), fill(0.8, n), fill(4.0, n),
            fill(80.0, n), fill(190.0, n), fill(80.0, n);
            temperature_air_mean=253.0, wind_speed_mean=4.0,
            precipitation_mean=0.8 * 365.25)
        mp = initialize_parameters(densification_method=:HerronLangway,
                                   output_frequency=:last)

        prof_ss = initialize_profile(mp, cf)       # steady-state (default)

        # Both escape-hatch flags give the legacy all-ice column, built on the
        # *configured* column_depth_max. Pin that to the depth the steady-state column
        # derived so the two starts share a grid — the convergence claim below is
        # about the initial state, not the geometry. (Asking for exactly `sum(dz)`
        # reproduces the grid cell for cell: `initialize_grid` stops at the first
        # depth at or past its target, and that depth is already `sum(dz)`.)
        mp_ice = initialize_parameters(densification_method=:HerronLangway,
                                       output_frequency=:last,
                                       column_depth_max=sum(prof_ss[:dz]))
        prof_ice = initialize_profile(mp_ice, cf; constant_density=true, constant_temperature=true)
        @test length(prof_ice[:dz]) == length(prof_ss[:dz])

        # Steady-state start is a graded firn column, not pure ice.
        d_ss = parent(prof_ss[:density])
        @test !all(d_ss .== mp.density_ice)
        @test d_ss[1] < 400.0            # low-density surface firn
        @test all(parent(prof_ice[:density]) .== mp.density_ice)

        # Both converge to the same depth-averaged density; steady-state faster.
        #
        # Convergence is measured as distance to the common equilibrium, not as the first
        # cycle where the change between consecutive cycles falls below a threshold. The
        # consecutive-delta proxy is not a convergence measure here: with multi-year
        # forcing the column settles into a small seasonal limit cycle, so once both
        # starts are in the asymptotic tail their per-cycle deltas are indistinguishable
        # regardless of how much sooner one of them arrived.
        function spin_trajectory(prof0, depth; maxit=20)
            prof = prof0
            rho_avg = Float64[]
            for _ in 1:maxit
                prof = gemb_profile(gemb(prof, cf, mp; verbose=false))
                zc = -GEMB.dz2z(prof[:dz])
                rho = prof[:density]
                idx = zc .<= depth
                push!(rho_avg, sum(rho[idx]) / count(idx))
            end
            return rho_avg
        end

        # 40 cycles, not 20: the all-ice start has to densify ~410 kg m-3 of gap away, and at
        # cycle 20 it is still en route (gap 23 kg m-3, closing to 8.5 by cycle 30 and 1.7 by
        # cycle 60). Sampling before the slower trajectory has arrived tests the *rate* of
        # convergence, not the claim being made here, which is that both starts reach a
        # common equilibrium.
        traj_ss  = spin_trajectory(prof_ss,  40.0; maxit=40)
        traj_ice = spin_trajectory(prof_ice, 40.0; maxit=40)
        ρ_eq = traj_ss[end]

        # Same equilibrium from either start.
        @test isapprox(traj_ss[end], traj_ice[end]; atol=20.0)

        # Steady-state initialization starts far closer to equilibrium and reaches it in
        # fewer cycles. (Observed: initial gap ~14 vs ~338 kg m-3; 1 cycle vs 4 to get
        # within 20 kg m-3.)
        @test abs(traj_ss[1] - ρ_eq) < abs(traj_ice[1] - ρ_eq)
        cycles_within(traj, tol) = something(findfirst(a -> abs(a - ρ_eq) < tol, traj),
                                             length(traj) + 1)
        @test cycles_within(traj_ss, 20.0) < cycles_within(traj_ice, 20.0)
    end

    @testset "Steady-state init converges at three contrasting sites" begin
        # The testset above establishes the payoff at one dry-cold site. The claim is
        # that it holds across regimes, so repeat it at a dry-snow site, a
        # percolation site with substantial melt, and an ablation site — and in
        # particular that all three reach the *same* fixed point as an ice-block
        # start. A better guess must be a shortcut to the same answer, not a
        # different answer.
        function make_site(; T_mean, T_amp, precip, sw, lw, vapor)
            n = 365
            time = DateTime(2020, 1, 1) .+ Day.(0:n-1)
            temp = T_mean .+ T_amp .* cos.(2π .* collect(1:n) ./ 365.0 .- π)
            return initialize_forcing(
                time, temp, fill(85000.0, n), fill(precip, n), fill(4.0, n),
                fill(sw, n), fill(lw, n), fill(vapor, n);
                temperature_air_mean=T_mean, wind_speed_mean=4.0,
                precipitation_mean=precip * 365.25)
        end

        sites = (
            # name, forcing, expected sign of the net annual balance
            ("dry-cold",    make_site(T_mean=250.0, T_amp=10.0, precip=0.8,
                                      sw=80.0,  lw=190.0, vapor=60.0),  :positive),
            ("percolation", make_site(T_mean=262.0, T_amp=8.0,  precip=1.5,
                                      sw=140.0, lw=270.0, vapor=200.0), :positive),
            ("ablation",    make_site(T_mean=276.0, T_amp=8.0,  precip=0.3,
                                      sw=220.0, lw=310.0, vapor=400.0), :negative),
        )

        mp_base = initialize_parameters(output_frequency=:last)

        for (name, cf, sign) in sites
            cs = GEMB.initialize_climate_summary(cf, mp_base)
            if sign === :positive
                @test cs.balance > 0.0
            else
                @test cs.balance < 0.0
            end

            prof_ss = initialize_profile(mp_base, cf)
            # Ice-block start on the *same* grid, so this compares initial states. The
            # escape hatch builds on the configured column_depth_max, so pin it to the
            # depth the steady-state column derived.
            mp_ice = initialize_parameters(output_frequency=:last,
                                           column_depth_max=sum(prof_ss[:dz]))
            prof_ice = initialize_profile(mp_ice, cf;
                constant_density=true, constant_temperature=true)
            @test length(prof_ice[:dz]) == length(prof_ss[:dz])

            function traj(prof0; maxit=25, depth=20.0)
                prof = prof0
                out = Float64[]
                for _ in 1:maxit
                    prof = gemb_profile(gemb(prof, cf, mp_base; verbose=false))
                    zc = -GEMB.dz2z(prof[:dz])
                    idx = zc .<= depth
                    push!(out, sum(prof[:density][idx]) / count(idx))
                end
                return out
            end

            t_ss = traj(prof_ss)
            t_ice = traj(prof_ice)
            ρ_eq = t_ss[end]

            # Same fixed point from either start.
            @test isapprox(t_ss[end], t_ice[end]; atol=20.0)

            # The steady-state guess starts no further from it, and converges no
            # slower. Measured initial gaps (ss vs ice): dry-cold 5 vs 391,
            # percolation 61 vs 253, ablation 0.6 vs 0.3 — the ablation column *is*
            # the ice block (both are solid ice at `density_ice`), so the two coincide
            # there to within a fraction of a kg m-3 rather than one winning.
            #
            # The tolerance is what makes that third case meaningful. The two ablation columns
            # differ only in their initial temperature: the escape-hatch path fills exactly
            # `CtoK` while the steady-state march returns the mean surface temperature the
            # energy balance gives (273.106 K here), and that 0.04 K seeds a sub-kg m-3
            # difference in how the first cycles compact. Requiring the march to win that
            # comparison outright would be asserting a coincidence of the clamp, not the claim
            # this testset makes.
            @test abs(t_ss[1] - ρ_eq) <= abs(t_ice[1] - ρ_eq) + 1.0
            within(t, tol) = something(findfirst(a -> abs(a - ρ_eq) < tol, t),
                                       length(t) + 1)
            @test within(t_ss, 20.0) <= within(t_ice, 20.0)
        end
    end

    @testset "Climatology provenance" begin
        # Multi-year forcing so forcing_climatology has complete years to average.
        n = 365 * 3
        time = DateTime(1950, 1, 1) .+ Day.(0:n-1)
        cf = initialize_forcing(
            time, fill(255.0, n), fill(85000.0, n), fill(0.5, n), fill(3.0, n),
            fill(50.0, n), fill(180.0, n), fill(80.0, n);
            temperature_air_mean=255.0, wind_speed_mean=3.0, precipitation_mean=182.6)

        window = (DateTime(1950, 1, 1), DateTime(1952, 12, 31))
        clim = forcing_climatology(cf, window)
        m = DimensionalData.metadata(clim)

        # Requested window is recorded verbatim.
        @test m[:climatology_window_start] == window[1]
        @test m[:climatology_window_stop] == window[2]
        # Three complete non-leap years averaged.
        @test m[:climatology_n_years] == 3
        @test m[:climatology_steps_per_year] == 365
        # Accessible via the cf.<field> interface too.
        @test clim.climatology_n_years == 3

        # No-window form falls back to the extent of the averaged years.
        clim2 = forcing_climatology(cf)
        m2 = DimensionalData.metadata(clim2)
        @test m2[:climatology_n_years] == 3
        @test m2[:climatology_window_start] isa DateTime
    end

    @testset "Spinup provenance" begin
        n = 365
        time = DateTime(2020, 1, 1) .+ Day.(0:n-1)
        forcing = initialize_forcing(
            time, fill(255.0, n), fill(85000.0, n), fill(0.5, n), fill(3.0, n),
            fill(50.0, n), fill(180.0, n), fill(80.0, n);
            temperature_air_mean=255.0, wind_speed_mean=3.0, precipitation_mean=182.6)
        params = initialize_parameters(output_frequency=:last)
        profile = initialize_profile(params, forcing)

        # No convergence check → runs the full max_iterations, converged=false.
        prof_max = gemb_spinup(profile, forcing, params; max_iterations=3)
        pm = DimensionalData.metadata(prof_max)
        @test pm[:spinup_cycles] == 3
        @test pm[:spinup_converged] == false
        @test pm[:spinup_max_iterations] == 3
        @test isnan(pm[:spinup_final_delta_density])
        @test isnan(pm[:spinup_final_drift_density])

        # Loose tolerance → converges early, before max_iterations.
        prof_conv = gemb_spinup(profile, forcing, params;
                                max_iterations=50, convergence_delta_density=1e3)
        pc = DimensionalData.metadata(prof_conv)
        @test pc[:spinup_converged] == true
        @test pc[:spinup_cycles] < 50
        @test pc[:spinup_convergence_delta_density] == 1e3
        @test isfinite(pc[:spinup_final_delta_density])
        # Drift was not requested, so it is recorded as unset rather than as passed.
        @test pc[:spinup_convergence_drift_density] === nothing
        @test isnan(pc[:spinup_final_drift_density])
        @test pc[:spinup_drift_window] == GEMB.SPINUP_DRIFT_WINDOW

        # --- Mean SMB rate over the final cycle -----------------------------------
        # SMB is a surface-flux quantity: precipitation in, sublimation out, runoff away.
        # Recomputed here from those three outputs rather than pinned to a number, so the
        # test fixes the definition and the sign convention.
        r = pm[:spinup_smb_rate]
        @test isfinite(r)
        # Accumulating site (precipitation_mean > 0, 255 K, no melt) → positive.
        @test r > 0.0

        # Measured over the *final* cycle, so rerunning that cycle from the cycle-2 profile
        # must reproduce it exactly.
        prof_2 = gemb_spinup(profile, forcing, params; max_iterations=2)
        last_cycle = gemb(prof_2, forcing, params)
        # One year of daily forcing: 365 steps of 86400 s. The last step integrates a full
        # step of its own, so the length is n*dt, not the first-to-last time span.
        years = n * forcing.time_step / (365.25 * 86400)
        tot(v) = sum(parent(last_cycle[v]))
        smb_mass = tot(:precipitation) + tot(:evaporation_condensation) - tot(:runoff)
        @test r ≈ smb_mass / params.density_ice / years rtol = 1e-12

        # NOT the basal ice flux. `-cumsum(ice_flux)` is surface *elevation* change, which is
        # SMB plus a compaction term, so it equals SMB only where compaction has equilibrated
        # and can even carry the opposite sign. SMB is a climate quantity and `ice_flux` a
        # model construct; pinned so a future refactor cannot quietly conflate them again.
        elevation_rate = -sum(parent(last_cycle[:ice_flux])) / years
        @test !isapprox(r, elevation_rate; rtol=1e-3)

        # A rate must be reported even from a single-cycle spinup.
        @test isfinite(DimensionalData.metadata(
            gemb_spinup(profile, forcing, params; max_iterations=1))[:spinup_smb_rate])

        # Ice density is recorded on the run output so a consumer can convert the mass-flux
        # outputs to metres of ice with the value the run actually used.
        @test DimensionalData.metadata(last_cycle)["density_ice"] == params.density_ice
    end

    @testset "Column-mean density is mass-weighted and grid-independent" begin
        mk(dz, rho) = DimStack((dz=DimArray(dz, (Z(1:length(dz)),)),
                                density=DimArray(rho, (Z(1:length(dz)),))))

        # Mass-weighted, not cell-count-weighted: one thick dense cell must dominate
        # three thin light ones. Arithmetic mean would be 550.
        s = mk([0.1, 0.1, 0.1, 9.7], [400.0, 400.0, 400.0, 900.0])
        @test GEMB._column_mean_density(s) ≈ (3 * 0.1 * 400.0 + 9.7 * 900.0) / 10.0
        @test GEMB._column_mean_density(s) > 850.0

        # Splitting a cell in two without moving mass leaves the mean untouched — this
        # is the grid-independence the removed spline used to buy over a partial depth.
        coarse = mk([1.0, 1.0], [500.0, 700.0])
        fine   = mk([0.5, 0.5, 0.25, 0.75], [500.0, 500.0, 700.0, 700.0])
        @test GEMB._column_mean_density(coarse) ≈ GEMB._column_mean_density(fine)
        @test GEMB._column_mean_density(coarse) ≈ 600.0
    end

    @testset "Drift: least-squares slope" begin
        # Quantity-agnostic: the same regression serves the density and FAC criteria, so this
        # exercises it on bare numbers rather than on either quantity's units.
        drift = GEMB._least_squares_slope

        # Exact recovery of a known slope, and sign follows the trend direction.
        rising = [100.0 + 3.0 * k for k in 0:9]
        @test drift(rising, 10) ≈ 3.0
        @test drift(reverse(rising), 10) ≈ -3.0

        # A settled column oscillating about a fixed mean has ~zero slope even though
        # its consecutive deltas are large — the case the step test cannot detect.
        osc = [500.0 + (isodd(k) ? 5.0 : -5.0) for k in 1:10]
        @test abs(drift(osc, 10)) < 1.5
        @test maximum(abs.(diff(osc))) ≈ 10.0

        # Window truncates to the trailing cycles: the early transient is excluded.
        history = vcat([100.0, 200.0, 300.0], [500.0 + 0.1 * k for k in 0:9])
        @test drift(history, 10) ≈ 0.1 atol = 1e-9
        # Fewer points than the window is allowed (fits what is there); fewer than 2 is not.
        @test drift(rising, 100) ≈ 3.0
        @test isnan(drift([1.0], 10))
    end

    @testset "Drift convergence criterion" begin
        n = 365
        time = DateTime(2020, 1, 1) .+ Day.(0:n-1)
        forcing = initialize_forcing(
            time, fill(255.0, n), fill(85000.0, n), fill(0.5, n), fill(3.0, n),
            fill(50.0, n), fill(180.0, n), fill(80.0, n);
            temperature_air_mean=255.0, wind_speed_mean=3.0, precipitation_mean=182.6)
        params = initialize_parameters(output_frequency=:last)
        profile = initialize_profile(params, forcing)

        # Drift alone, loose tolerance. It cannot fire before the window is full, so the
        # earliest possible exit is at cycle == drift_window.
        prof_d = gemb_spinup(profile, forcing, params;
                             max_iterations=12, convergence_drift_density=1e3,
                             drift_window=3)
        pd = DimensionalData.metadata(prof_d)
        @test pd[:spinup_converged] == true
        @test pd[:spinup_cycles] == 3
        @test pd[:spinup_drift_window] == 3
        @test pd[:spinup_convergence_drift_density] == 1e3
        @test isfinite(pd[:spinup_final_drift_density])
        @test pd[:spinup_final_drift_density] >= 0.0   # magnitude, not signed slope

        # An unreachable drift tolerance cannot converge, however loose the delta is.
        prof_and = gemb_spinup(profile, forcing, params;
                               max_iterations=4, convergence_delta_density=1e3,
                               convergence_drift_density=0.0, drift_window=2)
        pa = DimensionalData.metadata(prof_and)
        @test pa[:spinup_converged] == false
        @test pa[:spinup_cycles] == 4
        # ...and the delta criterion did pass on its own, so AND is what blocked it.
        @test pa[:spinup_final_delta_density] < 1e3

        # Symmetrically: a loose drift cannot rescue an unreachable delta.
        prof_and2 = gemb_spinup(profile, forcing, params;
                                max_iterations=4, convergence_delta_density=0.0,
                                convergence_drift_density=1e3, drift_window=2)
        @test DimensionalData.metadata(prof_and2)[:spinup_converged] == false

        # Both loose → converges as soon as both are computable (cycle 2 here).
        prof_both = gemb_spinup(profile, forcing, params;
                                max_iterations=12, convergence_delta_density=1e3,
                                convergence_drift_density=1e3, drift_window=2)
        @test DimensionalData.metadata(prof_both)[:spinup_cycles] == 2

        # A window too short to fit a slope is rejected up front, not silently ignored.
        @test_throws ErrorException gemb_spinup(profile, forcing, params;
            max_iterations=1, convergence_drift_density=1.0, drift_window=1)
    end

    @testset "FAC convergence criterion" begin
        # Why FAC is a criterion in its own right rather than a density tolerance in disguise:
        # `FAC = Σ dz (ρᵢ − ρ)/ρᵢ`, so at fixed column depth `Z`, `∂FAC/∂ρ̄ = −Z/ρᵢ`. The FAC
        # residual a given density tolerance admits therefore scales with column depth, and the
        # deep cold columns that altimetry work cares about are exactly where it is loosest
        # (measured across a 21-site fleet: 0.02 mm of FAC per 1e-3 kg/m³ at 14 m against
        # 0.23 mm at 212 m — a 12x spread). The first testset below pins that scaling directly,
        # since it is the whole justification for the criterion existing.
        n = 365
        time = DateTime(2020, 1, 1) .+ Day.(0:n-1)
        forcing = initialize_forcing(
            time, fill(255.0, n), fill(85000.0, n), fill(0.5, n), fill(3.0, n),
            fill(50.0, n), fill(180.0, n), fill(80.0, n);
            temperature_air_mean=255.0, wind_speed_mean=3.0, precipitation_mean=182.6)
        params = initialize_parameters(output_frequency=:last)
        profile = initialize_profile(params, forcing)

        @testset "FAC sensitivity to mean density scales with depth" begin
            # Two columns of the same uniform density, differing only in depth. Perturbing the
            # density by the same amount moves the deeper column's FAC proportionally more.
            zdim(m) = Z(1:m)
            uniform(m, dz, rho) = DimStack((
                dz=DimArray(fill(dz, m), (zdim(m),)),
                density=DimArray(fill(rho, m), (zdim(m),))))

            ρ, Δρ = 600.0, 1.0
            shallow, deep = uniform(10, 1.0, ρ), uniform(100, 1.0, ρ)      # 10 m vs 100 m
            dfac(col_a, col_b) = abs(GEMB._column_fac(col_a, params) -
                                     GEMB._column_fac(col_b, params))

            d_shallow = dfac(shallow, uniform(10, 1.0, ρ + Δρ))
            d_deep    = dfac(deep,    uniform(100, 1.0, ρ + Δρ))

            # Ten times the depth, ten times the FAC movement for the same Δρ̄.
            @test d_deep ≈ 10 * d_shallow rtol = 1e-12
            # And the analytic value, Z·Δρ/ρᵢ.
            @test d_deep ≈ 100.0 * Δρ / params.density_ice rtol = 1e-12
        end

        # FAC alone converges, and reports its own measures while leaving the density ones
        # untouched at NaN — an unrequested quantity must be distinguishable from a measured
        # zero.
        prof_f = gemb_spinup(profile, forcing, params;
                             max_iterations=6, convergence_delta_fac=1e3)
        pf = DimensionalData.metadata(prof_f)
        @test pf[:spinup_converged] == true
        @test pf[:spinup_cycles] == 2                  # earliest a step test can fire
        @test isfinite(pf[:spinup_final_delta_fac])
        @test pf[:spinup_final_delta_fac] >= 0.0        # magnitude
        @test pf[:spinup_convergence_delta_fac] == 1e3
        @test isnan(pf[:spinup_final_delta_density])
        @test pf[:spinup_convergence_delta_density] === nothing

        # FAC drift, like density drift, cannot fire before the window is full.
        prof_fd = gemb_spinup(profile, forcing, params;
                              max_iterations=12, convergence_drift_fac=1e3,
                              drift_window=3)
        pfd = DimensionalData.metadata(prof_fd)
        @test pfd[:spinup_converged] == true
        @test pfd[:spinup_cycles] == 3
        @test isfinite(pfd[:spinup_final_drift_fac])
        @test pfd[:spinup_final_drift_fac] >= 0.0
        @test isnan(pfd[:spinup_final_drift_density])

        # The criteria are conjunctive *across quantities*, not only within one: an
        # unreachable FAC tolerance blocks convergence however loose the density one is, and
        # vice versa. This is the property that makes adding a criterion always a tightening.
        pb1 = DimensionalData.metadata(gemb_spinup(profile, forcing, params;
            max_iterations=4, convergence_delta_density=1e3, convergence_delta_fac=0.0))
        @test pb1[:spinup_converged] == false
        @test pb1[:spinup_final_delta_density] < 1e3    # density alone would have passed

        pb2 = DimensionalData.metadata(gemb_spinup(profile, forcing, params;
            max_iterations=4, convergence_delta_density=0.0, convergence_delta_fac=1e3))
        @test pb2[:spinup_converged] == false
        @test pb2[:spinup_final_delta_fac] < 1e3        # FAC alone would have passed

        # Both loose → converges as soon as both are computable.
        pb3 = DimensionalData.metadata(gemb_spinup(profile, forcing, params;
            max_iterations=12, convergence_delta_density=1e3, convergence_delta_fac=1e3))
        @test pb3[:spinup_cycles] == 2
        @test pb3[:spinup_converged] == true

        # `drift_window` is validated for the FAC drift criterion too, not only the density one.
        @test_throws ErrorException gemb_spinup(profile, forcing, params;
            max_iterations=1, convergence_drift_fac=1.0, drift_window=1)

        # The default call requests nothing, so every measure is NaN/nothing and the spinup
        # runs its full iteration count — unchanged by this feature existing.
        pn = DimensionalData.metadata(gemb_spinup(profile, forcing, params; max_iterations=2))
        @test pn[:spinup_converged] == false
        @test pn[:spinup_cycles] == 2
        @test isnan(pn[:spinup_final_delta_fac])
        @test isnan(pn[:spinup_final_drift_fac])
        @test pn[:spinup_convergence_delta_fac] === nothing
        @test pn[:spinup_convergence_drift_fac] === nothing

        # `_column_fac` is the same integral the `gemb` output reports, so a criterion and a
        # diagnostic cannot disagree about what FAC means.
        @test GEMB._column_fac(profile, params) ≈
              firn_air_content(parent(profile[:dz]), parent(profile[:density]),
                               params.density_ice)
    end

    @testset "Spinup with zero accumulation" begin
        # Edge case: no precipitation

        params = initialize_parameters(
            output_frequency = :last
        )

        n_days = 100  # Shorter for zero accumulation test
        start_time = DateTime(2020, 1, 1)
        time = start_time .+ Day.(0:n_days-1)

        forcing = initialize_forcing(
            time,
            fill(250.0, n_days),  # temperature_air
            fill(85000.0, n_days),  # pressure_air
            zeros(n_days),  # precipitation - Zero accumulation
            fill(3.0, n_days),  # wind_speed
            zeros(n_days),  # shortwave_downward
            fill(150.0, n_days),  # longwave_downward
            fill(50.0, n_days),  # vapor_pressure
            temperature_observation_height = 2.0,
            wind_observation_height = 10.0
        )

        # Should still work, just won't grow
        profile = initialize_profile(params, forcing)
        output = gemb_spinup(profile, forcing, params; max_iterations=2)

        @test output isa DimStack
        @test haskey(output, :temperature)

        # Profile should still exist
        temps = output[:temperature][Ti=1]
        @test sum(.!isnan.(temps)) > 0
    end
end
