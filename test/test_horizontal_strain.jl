using Test
using GEMB

# Horizontal (ice-dynamic) strain thinning. There is no MATLAB reference for this term —
# MATLAB GEMB has no ice-dynamic term at all — so the tests pin it against its own closed
# form, `dz *= exp(-D dt)`, and against the mass/energy the column reports losing.

_cols(; dz=[1.0, 2.0], water=[5.0, 0.0], density=[400.0, 600.0],
    temperature=[260.0, 265.0]) = (
    temperature=copy(temperature), dz=copy(dz), density=copy(density),
    water=copy(water), grain_radius=[0.1, 0.2], grain_dendricity=[0.0, 0.0],
    grain_sphericity=[0.5, 0.5])

@testset "apply_horizontal_strain!" begin
    yr = GEMB.SECONDS_PER_YEAR

    @testset "analytic thinning, exactly exp(-D dt)" begin
        mp = initialize_parameters(horizontal_strain_rate=0.01)
        cols = _cols(dz=[1.0, 2.0])
        GEMB.apply_horizontal_strain!(cols, yr, mp)
        @test cols.dz[1] ≈ exp(-0.01) atol = 1e-14
        @test cols.dz[2] ≈ 2 * exp(-0.01) atol = 1e-14

        # The exponential composes, so sub-stepping gives the same answer. Forward Euler
        # would not: ten steps of D dt/10 would land ~5e-5 short.
        cols10 = _cols(dz=[1.0, 2.0])
        for _ in 1:10
            GEMB.apply_horizontal_strain!(cols10, yr / 10, mp)
        end
        @test cols10.dz ≈ cols.dz atol = 1e-13
    end

    @testset "density and temperature untouched; water scales with dz" begin
        mp = initialize_parameters(horizontal_strain_rate=0.05)
        cols = _cols()
        d0 = copy(cols.density)
        t0 = copy(cols.temperature)
        dz0 = copy(cols.dz)
        w0 = copy(cols.water)
        GEMB.apply_horizontal_strain!(cols, yr, mp)
        @test cols.density == d0          # bit-for-bit: strain is at constant density
        @test cols.temperature == t0
        @test cols.water[1] / w0[1] ≈ cols.dz[1] / dz0[1] atol = 1e-15
        @test cols.grain_radius == [0.1, 0.2]   # intensive, unchanged
    end

    @testset "sign convention" begin
        cols_div = _cols()
        mass_div, _ = GEMB.apply_horizontal_strain!(cols_div,
            yr, initialize_parameters(horizontal_strain_rate=0.01))
        @test cols_div.dz[1] < 1.0        # divergence thins
        @test mass_div < 0.0              # mass leaves laterally

        cols_con = _cols()
        mass_con, _ = GEMB.apply_horizontal_strain!(cols_con,
            yr, initialize_parameters(horizontal_strain_rate=-0.01))
        @test cols_con.dz[1] > 1.0        # convergence thickens
        @test mass_con > 0.0

        # Zero is an exact no-op: this is what keeps the default bit-identical.
        cols_zero = _cols()
        m, e = GEMB.apply_horizontal_strain!(cols_zero, yr,
            initialize_parameters(horizontal_strain_rate=0.0))
        @test (m, e) === (0.0, 0.0)
        @test cols_zero.dz == [1.0, 2.0]
        @test cols_zero.water == [5.0, 0.0]
    end

    @testset "reported mass and energy match the column change" begin
        mp = initialize_parameters(horizontal_strain_rate=0.02)
        cols = _cols()
        mass_before = GEMB.column_mass(cols)
        _, energy_before = GEMB.column_mass_energy(cols, mp)
        mass_lateral, E_lateral = GEMB.apply_horizontal_strain!(cols, yr, mp)
        mass_after = GEMB.column_mass(cols)
        _, energy_after = GEMB.column_mass_energy(cols, mp)
        @test mass_lateral ≈ mass_after - mass_before atol = 1e-12
        @test E_lateral ≈ energy_after - energy_before atol =
            GEMB.energy_tolerance(energy_before)
    end

    @testset "validator bounds" begin
        @test_throws AssertionError initialize_parameters(horizontal_strain_rate=2.0)
        @test_throws AssertionError initialize_parameters(horizontal_strain_rate=-2.0)
        @test_throws AssertionError initialize_parameters(horizontal_strain_rate=NaN)
    end
end

@testset "horizontal strain through gemb_core and gemb" begin
    ds = GEMB_ClimateForcing.simulate_climate_forcing("test_1", 3)
    cf = GEMB.initialize_forcing(ds)

    mp0 = initialize_parameters(output_frequency=:monthly)
    profile = initialize_profile(mp0, cf)

    # Off by default: `strain_thinning` is present and identically zero, and the run is
    # unaffected. `verbose=true` exercises the per-step and whole-run mass budgets, both of
    # which now carry the lateral term.
    out0 = gemb(profile, cf, mp0; verbose=true)
    @test haskey(out0, :strain_thinning)
    @test all(iszero, out0[:strain_thinning])

    # On: the budget must still close (this is the test that catches an incomplete wiring —
    # omitting `mass_lateral` from either check makes `gemb` error), the grid invariants
    # must hold, and the column must genuinely thin.
    #
    # `2e-2 yr-1` rather than a more typical rate because the firn-air comparison below has
    # to clear this run's chaotic noise floor: a 1e-14 relative perturbation to a single
    # cell's initial temperature moves summed `firn_air_content` by ~0.7 m, so bit-level
    # rounding differences between platforms are not small here. At 2e-2 the strain signal
    # is ~12 m, a ~17x margin; at 5e-3 it is ~1.9 m and the ordering flips between x86 and
    # arm64. This is a test-robustness bound, not a physical one.
    mp = initialize_parameters(output_frequency=:monthly, horizontal_strain_rate=2e-2)
    out = gemb(profile, cf, mp; verbose=true)
    @test all(>(0.0), out[:strain_thinning])
    @test sum(out[:strain_thinning]) > sum(out0[:strain_thinning])

    # `trim_bottom!` refills the column, so the fixed-depth/fixed-count invariants hold
    # exactly as they do without strain.
    z_target = sum(profile[:dz])
    for i in axes(out[:dz], 2)
        @test sum(view(out[:dz], :, i)) ≈ z_target atol = 1e-9
    end

    # Thinning exports firn air, so the column carries less of it than the D = 0 run.
    @test sum(out[:firn_air_content]) < sum(out0[:firn_air_content])

    # Convergence has the opposite sign in both the diagnostic and the air content.
    mp_con = initialize_parameters(output_frequency=:monthly, horizontal_strain_rate=-2e-2)
    out_con = gemb(profile, cf, mp_con; verbose=true)
    @test all(<(0.0), out_con[:strain_thinning])
    @test sum(out_con[:firn_air_content]) > sum(out0[:firn_air_content])
end
