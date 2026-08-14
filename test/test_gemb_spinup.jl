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
        profile, params = initialize_profile(params, forcing)

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
        profile, params = initialize_profile(params, forcing)
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

        profile, params = initialize_profile(params, forcing)
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

        prof, mp_out = initialize_profile(mp, cf)
        @test all(parent(prof[:density]) .== mp.density_ice)
        @test all(parent(prof[:temperature]) .<= 273.15 + 1e-9)

        # No firn to resolve, so the derived column collapses to the thermal floor:
        # deep enough for the annual wave, far shallower than the 250 m default.
        @test mp_out.column_depth < 0.2 * mp.column_depth
        @test mp_out.column_depth >= mp.column_ztop
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

        prof_ss, mp = initialize_profile(mp, cf)   # steady-state (default); rebind the
                                                   # derived column depth
        # Both escape-hatch flags give the legacy all-ice column. Built from the
        # *rebound* mp so the two starts share a grid — the convergence claim below is
        # about the initial state, not the geometry.
        prof_ice, _ = initialize_profile(mp, cf; constant_density=true, constant_temperature=true)
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

            prof_ss, mp = initialize_profile(mp_base, cf)
            # Ice-block start on the *same* grid, so this compares initial states.
            prof_ice, _ = initialize_profile(mp, cf;
                constant_density=true, constant_temperature=true)

            function traj(prof0; maxit=25, depth=20.0)
                prof = prof0
                out = Float64[]
                for _ in 1:maxit
                    prof = gemb_profile(gemb(prof, cf, mp; verbose=false))
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
            # percolation 61 vs 253, ablation 0.6 vs 0.6 — the ablation column *is*
            # the ice block, so the two coincide there rather than one winning.
            @test abs(t_ss[1] - ρ_eq) <= abs(t_ice[1] - ρ_eq) + 1e-9
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
        profile, params = initialize_profile(params, forcing)

        # No convergence check → runs the full max_iterations, converged=false.
        prof_max = gemb_spinup(profile, forcing, params; max_iterations=3)
        pm = DimensionalData.metadata(prof_max)
        @test pm[:spinup_cycles] == 3
        @test pm[:spinup_converged] == false
        @test pm[:spinup_max_iterations] == 3
        @test isnan(pm[:spinup_final_delta_density])

        # Loose tolerance → converges early, before max_iterations.
        prof_conv = gemb_spinup(profile, forcing, params;
                                max_iterations=50, convergence_delta_density=1e3)
        pc = DimensionalData.metadata(prof_conv)
        @test pc[:spinup_converged] == true
        @test pc[:spinup_cycles] < 50
        @test pc[:spinup_convergence_delta_density] == 1e3
        @test isfinite(pc[:spinup_final_delta_density])
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
        profile, params = initialize_profile(params, forcing)
        output = gemb_spinup(profile, forcing, params; max_iterations=2)

        @test output isa DimStack
        @test haskey(output, :temperature)

        # Profile should still exist
        temps = output[:temperature][Ti=1]
        @test sum(.!isnan.(temps)) > 0
    end
end
