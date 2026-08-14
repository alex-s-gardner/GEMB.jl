# Tests for heat_capacity.jl — specific heat capacity and the enthalpy it integrates to.
# No MATLAB counterpart: MATLAB has only the constant C_ICE = 2102.

_hc_const() = GEMB.ModelParameters(heat_capacity_method=:constant)
_hc_cp() = GEMB.ModelParameters(heat_capacity_method=:CuffeyPaterson)

@testset "heat_capacity values" begin
    mp = _hc_const()
    # :constant is temperature-independent by construction.
    for T in (180.0, 240.0, 273.15)
        @test GEMB.heat_capacity(mp, T) === GEMB.HEAT_CAPACITY_ICE_DEFAULT
    end
    mp_alt = GEMB.ModelParameters(heat_capacity_method=:constant, heat_capacity_ice=2000.0)
    @test GEMB.heat_capacity(mp_alt, 250.0) === 2000.0

    # Cuffey & Paterson (2010) eq. 9.1: c_p = 152.5 + 7.122 T.
    cp = _hc_cp()
    @test GEMB.heat_capacity(cp, GEMB.CtoK) ≈ 2097.87 atol = 0.01
    @test GEMB.heat_capacity(cp, 240.0) ≈ 1861.78 atol = 0.01
    @test GEMB.heat_capacity(cp, 210.0) ≈ 1648.12 atol = 0.01

    # `heat_capacity_ice` is ignored under :CuffeyPaterson.
    cp_alt = GEMB.ModelParameters(heat_capacity_method=:CuffeyPaterson, heat_capacity_ice=1600.0)
    @test GEMB.heat_capacity(cp_alt, 250.0) == GEMB.heat_capacity(cp, 250.0)

    # The constant value over-permits refreezing in cold firn by these margins:
    # +12.9% at 240 K, +27.5% at 210 K.
    @test GEMB.HEAT_CAPACITY_ICE_DEFAULT / GEMB.heat_capacity(cp, 240.0) ≈ 1.129 atol = 1e-3
    @test GEMB.HEAT_CAPACITY_ICE_DEFAULT / GEMB.heat_capacity(cp, 210.0) ≈ 1.275 atol = 1e-3
end

@testset "specific_enthalpy" begin
    mp = _hc_const()
    # Exact equality, not approximate: this is what keeps the default path
    # bit-identical to the MATLAB `M*T*C_ICE` accounting.
    for T in (180.0, 213.7, 250.0, GEMB.CtoK)
        @test GEMB.specific_enthalpy(mp, T) == mp.heat_capacity_ice * T
    end

    cp = _hc_cp()
    a, b = GEMB.HEAT_CAPACITY_CUFFEY_A, GEMB.HEAT_CAPACITY_CUFFEY_B
    for T in (180.0, 250.0, GEMB.CtoK)
        @test GEMB.specific_enthalpy(cp, T) ≈ a * T + (b / 2) * T^2 rtol = 1e-15
    end

    # Against numerical quadrature of ∫₀ᵀ c_p dT′ (trapezoid, fine grid).
    for T_end in (200.0, GEMB.CtoK)
        grid = range(0.0, T_end, length=200_001)
        c = [GEMB.heat_capacity(cp, T) for T in grid]
        h_quad = sum((c[1:end-1] .+ c[2:end]) ./ 2) * step(grid)
        @test GEMB.specific_enthalpy(cp, T_end) ≈ h_quad rtol = 1e-9
    end

    # Strictly increasing (c_p > 0 over the physical range), both methods.
    for m in (mp, cp)
        h = [GEMB.specific_enthalpy(m, T) for T in 150.0:0.5:GEMB.CtoK]
        @test all(diff(h) .> 0.0)
    end
end

@testset "specific_enthalpy is not T*c_p(T)" begin
    # Regression guard for the bug this module exists to prevent. For c_p = a + bT,
    # `T*c_p(T)` = aT + bT² overstates the enthalpy aT + (b/2)T² by (b/2)T² — which at
    # the melting point is 2.66e5 J/kg, or 0.79 x the latent heat of fusion. That is a
    # dominant error, not a correction, so it must never be "simplified" back.
    cp = _hc_cp()
    T = GEMB.CtoK
    naive = T * GEMB.heat_capacity(cp, T)
    correct = GEMB.specific_enthalpy(cp, T)
    @test naive - correct > 2.6e5
    @test (naive - correct) / GEMB.LF > 0.79
    @test (GEMB.HEAT_CAPACITY_CUFFEY_B / 2) * T^2 ≈ naive - correct rtol = 1e-12

    # Under :constant the two coincide — which is exactly why the error hid.
    mp = _hc_const()
    @test T * GEMB.heat_capacity(mp, T) == GEMB.specific_enthalpy(mp, T)
end

@testset "temperature_from_specific_enthalpy round-trips" begin
    for m in (_hc_const(), _hc_cp())
        for T in 180.0:0.5:GEMB.CtoK
            h = GEMB.specific_enthalpy(m, T)
            @test GEMB.temperature_from_specific_enthalpy(m, h) ≈ T atol = 1e-11
        end
        # Monotone, and well-conditioned near h -> 0 (the cancellation-free root).
        T_small = GEMB.temperature_from_specific_enthalpy(m, 1e-6)
        @test 0.0 < T_small < 1e-6
        @test GEMB.temperature_from_specific_enthalpy(m, 0.0) == 0.0
    end
end

@testset "specific_enthalpy_water" begin
    mp = _hc_const()
    # The MATLAB budgets' `LF + CtoK * C_ICE`.
    @test GEMB.specific_enthalpy_water(mp) == GEMB.LF + GEMB.CtoK * mp.heat_capacity_ice
    cp = _hc_cp()
    @test GEMB.specific_enthalpy_water(cp) == GEMB.LF + GEMB.specific_enthalpy(cp, GEMB.CtoK)
    # Liquid water is warmer than ice at the same temperature, by exactly LF.
    @test GEMB.specific_enthalpy_water(cp) - GEMB.specific_enthalpy(cp, GEMB.CtoK) == GEMB.LF
end

@testset "mix_temperature conserves enthalpy" begin
    for m in (_hc_const(), _hc_cp())
        T1, M1, T2, M2 = 250.0, 3.0, 270.0, 7.0
        T = GEMB.mix_temperature(m, T1, M1, T2, M2)
        # Bracketed by the inputs...
        @test T1 < T < T2
        # ...and enthalpy-conserving to round-off.
        E_in = M1 * GEMB.specific_enthalpy(m, T1) + M2 * GEMB.specific_enthalpy(m, T2)
        E_out = (M1 + M2) * GEMB.specific_enthalpy(m, T)
        @test E_out ≈ E_in rtol = 1e-14

        # Degenerate cases.
        @test GEMB.mix_temperature(m, T1, M1, T2, 0.0) ≈ T1 atol = 1e-11
        @test GEMB.mix_temperature(m, T1, 0.0, T2, 0.0) ≈ T1 atol = 1e-11
        @test GEMB.mix_temperature(m, T1, M1, T1, M2) ≈ T1 atol = 1e-11

        # The enthalpy-argument form agrees with the temperature form.
        @test GEMB.mix_temperature(m, GEMB.specific_enthalpy(m, T1), M1,
            GEMB.specific_enthalpy(m, T2), M2, Val(:enthalpy)) ≈ T atol = 1e-12
    end
end

@testset "mass-weighted mean temperature loses energy" begin
    # Negative control for `mix_temperature`: under :constant the mass-weighted mean is
    # exact, but under :CuffeyPaterson it loses (b/2)·M·Var_M(T) joules because h is
    # convex. This is the quantity that breaks cell merging if mixing is done in
    # temperature space.
    T1, M1, T2, M2 = 250.0, 5.0, 255.0, 5.0
    M = M1 + M2
    T_mean = (M1 * T1 + M2 * T2) / M

    mp = _hc_const()
    @test GEMB.mix_temperature(mp, T1, M1, T2, M2) ≈ T_mean atol = 1e-11

    cp = _hc_cp()
    T_mix = GEMB.mix_temperature(cp, T1, M1, T2, M2)
    @test T_mix > T_mean          # the correct mixture is warmer than the mean
    E_true = M1 * GEMB.specific_enthalpy(cp, T1) + M2 * GEMB.specific_enthalpy(cp, T2)
    E_mean = M * GEMB.specific_enthalpy(cp, T_mean)
    var_M = (M1 * (T1 - T_mean)^2 + M2 * (T2 - T_mean)^2) / M
    @test E_true - E_mean ≈ (GEMB.HEAT_CAPACITY_CUFFEY_B / 2) * M * var_M rtol = 1e-10
    # 223 J for this 10 kg pair, and it scales with mass: a 0.5 m cell of 400 kg m-3 firn
    # is 200 kg, so a realistic deep merge loses ~4.5 kJ. Against E_TOLERANCE = 1e-3 J,
    # mixing in temperature space would fail the conservation check immediately.
    @test E_true - E_mean ≈ 222.56 atol = 0.05
end

@testset "add_energy_temperature" begin
    for m in (_hc_const(), _hc_cp())
        T, M, Q = 250.0, 4.0, 1.0e5
        T_new = GEMB.add_energy_temperature(m, T, M, Q)
        @test T_new > T
        @test M * (GEMB.specific_enthalpy(m, T_new) - GEMB.specific_enthalpy(m, T)) ≈ Q rtol = 1e-13
        # Cooling, zero energy, and zero mass.
        @test GEMB.add_energy_temperature(m, T, M, -Q) < T
        @test GEMB.add_energy_temperature(m, T, M, 0.0) ≈ T atol = 1e-12
        @test GEMB.add_energy_temperature(m, T, 0.0, Q) == T
    end
    # Under :constant this is exactly T + Q/(M*c).
    mp = _hc_const()
    @test GEMB.add_energy_temperature(mp, 250.0, 4.0, 1.0e5) ≈
        250.0 + 1.0e5 / (4.0 * mp.heat_capacity_ice) rtol = 1e-15
end

@testset "heat capacity validator" begin
    @test_throws AssertionError initialize_parameters(heat_capacity_method=:Bogus)
    @test_throws AssertionError initialize_parameters(heat_capacity_ice=1000.0)
    @test_throws AssertionError initialize_parameters(heat_capacity_ice=3000.0)
    @test initialize_parameters(heat_capacity_method=:CuffeyPaterson).heat_capacity_method ===
        :CuffeyPaterson
    @test initialize_parameters().heat_capacity_method === :constant
    @test initialize_parameters().heat_capacity_ice == GEMB.HEAT_CAPACITY_ICE_DEFAULT
end

@testset "melt-equation energy helpers" begin
    for m in (_hc_const(), _hc_cp())
        # Refreezing exactly the cold content lands the cell at the melting point: this is
        # the consistency `cold_content_mass` and `refreeze_temperature` must share.
        T, M = 255.0, 300.0
        cc = GEMB.cold_content_mass(m, T, M)
        @test cc > 0.0
        @test GEMB.refreeze_temperature(m, T, M + cc, cc) ≈ GEMB.CtoK atol = 1e-9

        # No cold content at or above the melting point.
        @test GEMB.cold_content_mass(m, GEMB.CtoK, M) == 0.0
        @test GEMB.cold_content_mass(m, GEMB.CtoK + 1.0, M) == 0.0

        # The (density, dz) form agrees with the mass form up to grouping.
        @test GEMB.cold_content_mass(m, T, 400.0, 0.75) ≈
              GEMB.cold_content_mass(m, T, 400.0 * 0.75) rtol = 1e-14

        # Excess enthalpy round-trips through its own inverse.
        for dT in (0.0, 1e-6, 0.5, 5.0, 200.0)
            e = GEMB.excess_specific_enthalpy(m, dT)
            @test GEMB.excess_temperature_from_specific_enthalpy(m, e) ≈ dT atol = 1e-10
        end

        # Excess enthalpy is the enthalpy difference across the melting point.
        for dT in (0.5, 5.0, 50.0)
            @test GEMB.excess_specific_enthalpy(m, dT) ≈
                  GEMB.specific_enthalpy(m, GEMB.CtoK + dT) -
                  GEMB.specific_enthalpy(m, GEMB.CtoK) rtol = 1e-12
        end

        # At the full-melt threshold the excess enthalpy is exactly the latent heat, so a
        # cell that hot melts out entirely and has nothing left to pass down.
        T_fm = GEMB.full_melt_excess_temperature(m)
        @test GEMB.excess_specific_enthalpy(m, T_fm) ≈ GEMB.LF rtol = 1e-12
        @test GEMB.surplus_energy(m, T_fm, 0.0, 100.0) == 0.0
        # Twice the threshold leaves one latent heat per kg of surplus.
        dT2 = GEMB.excess_temperature_from_specific_enthalpy(m, 2 * GEMB.LF)
        @test GEMB.surplus_energy(m, dT2, dT2 - T_fm, 100.0) ≈ 100.0 * GEMB.LF rtol = 1e-12

        # Melt mass is the excess enthalpy divided by the latent heat.
        @test GEMB.melt_mass_from_excess(m, 2.0, 400.0, 0.5) ≈
              GEMB.excess_specific_enthalpy(m, 2.0) * 400.0 * 0.5 / GEMB.LF rtol = 1e-14
        @test GEMB.melt_mass_from_excess(m, 0.0, 400.0, 0.5) == 0.0
    end

    # Under :constant every helper reduces to the reference's exact expression, grouping
    # included — this is the bit-parity gate for the default path.
    mp = _hc_const()
    c = mp.heat_capacity_ice
    @test GEMB.full_melt_excess_temperature(mp) === GEMB.LF / c
    @test GEMB.cold_content_mass(mp, 260.0, 200.0) ===
          max(0.0, -((260.0 - GEMB.CtoK) * 200.0 * c) / GEMB.LF)
    @test GEMB.cold_content_mass(mp, 260.0, 400.0, 0.1) ===
          max(0.0, -((260.0 - GEMB.CtoK) * 400.0 * 0.1 * c) / GEMB.LF)
    @test GEMB.melt_mass_from_excess(mp, 1.5, 400.0, 0.1) === 1.5 * 400.0 * 0.1 * c / GEMB.LF
    @test GEMB.surplus_energy(mp, 200.0, 40.0, 30.0) === 40.0 * c * 30.0
    @test GEMB.refreeze_temperature(mp, 260.0, 210.0, 10.0) ===
          260.0 + (10.0 * (GEMB.LF + (GEMB.CtoK - 260.0) * c) / (210.0 * c))
end

@testset "column_enthalpy" begin
    for m in (_hc_const(), _hc_cp())
        M = [100.0, 250.0, 400.0]
        T = [240.0, 260.0, 271.0]
        water = [0.0, 2.0, 5.0]
        expect = sum(M[i] * GEMB.specific_enthalpy(m, T[i]) for i in eachindex(M))
        @test GEMB.column_enthalpy(m, M, T) ≈ expect rtol = 1e-14
        @test GEMB.column_enthalpy(m, M, T, water) ≈
              expect + sum(water) * GEMB.specific_enthalpy_water(m) rtol = 1e-14
        # The water-free form is the three-argument form with no water.
        @test GEMB.column_enthalpy(m, M, T, zeros(3)) ≈ GEMB.column_enthalpy(m, M, T) rtol = 1e-14
    end
end

@testset "calculate_melt stays type-stable and non-allocating" begin
    # Regression guard. The verbose energy budgets were originally written as generators
    # (`sum(M[i] * specific_enthalpy(...) for i in ...)`). Inside a function as large as
    # `calculate_melt` that closure pushed inference past its limit, so *every* local
    # inferred as `Any` and the whole hot path boxed — an 18x allocation regression that was
    # invisible to correctness tests because the numbers were unchanged. `column_enthalpy`
    # exists to keep these budgets as plain loops.
    mp = _hc_const()
    n = 200
    args = (fill(272.0, n), fill(0.1, n), collect(range(350.0, 800.0, n)), fill(0.5, n),
        fill(0.5, n), fill(0.5, n), fill(0.5, n))

    # No local may infer as Any.
    ir = code_typed(GEMB.calculate_melt,
        (typeof.(args)..., Float64, GEMB.ModelParameters, Bool); optimize=false)[1][1]
    slot_types = ir.slottypes
    @test !any(t -> t === Any, slot_types)

    # And the allocation count must stay near the array-count floor rather than scaling with
    # the cell count, which is what boxed arithmetic does.
    GEMB.calculate_melt(copy.(args)..., 0.0, mp, false)
    inputs = copy.(args)
    allocated = @allocated GEMB.calculate_melt(inputs..., 0.0, mp, false)
    @test allocated < 60_000
end

@testset "thermal solver conserves energy under both heat-capacity methods" begin
    # The solver carries enthalpy as its prognostic and diffuses by applying one flux per
    # face with opposite signs, so the column total is conserved independently of c_p. This
    # is the property the old normalized-coefficient form could not have under a
    # temperature-dependent c_p: it froze c_p per sub-step, injecting (b/2)ΔT² J kg-1 each
    # time. Run with verbose=true so the solver's own budget check is the assertion.
    n = 40
    dz = fill(0.05, n); dz[10:end] .= 0.2
    density = collect(range(320.0, 800.0, length=n))
    grain_radius = fill(0.2, n)
    shortwave_flux = [200.0 * exp(-i * 0.3) for i in 1:n]
    cfs = GEMB.ClimateForcingStep(10800.0, 265.0, 90000.0, 0.0, 5.0, 300.0, 250.0,
        200.0, 265.0, 5.0, 0.0, 2.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

    for method in (:constant, :CuffeyPaterson)
        mp = GEMB.ModelParameters(density_ice=917.0, emissivity=0.97,
            emissivity_method=:uniform, emissivity_grain_radius_large=0.97,
            emissivity_grain_radius_threshold=10.0, surface_roughness_effective_ratio=0.1,
            thermal_conductivity_method=:Sturm, heat_capacity_method=method,
            dt_divisors=Float64.(GEMB.fast_divisors(36000000)) ./ 10000)
        T0 = [250.0 + 20.0 * exp(-i / 8) for i in 1:n]

        # verbose=true makes the solver assert its own closure every sub-step, including the
        # exact `!=` check that the Dirichlet bottom cell never moved.
        T, _, _, _, _, _ = GEMB.calculate_temperature(copy(T0), copy(dz), copy(density),
            0.0, grain_radius, copy(shortwave_flux), cfs, mp, true)

        @test !any(isnan, T)
        @test T[end] == T0[end]                  # Dirichlet bottom, bit-exact
        @test T[1] > T0[1]                       # surface warms under this forcing
        # The solver does not cap at the melting point — excess enthalpy above CtoK is what
        # `calculate_melt` later converts to melt — but it must stay physically bounded.
        @test all(200.0 .< T .< 290.0)
    end
end

@testset "solver agrees between methods where c_p agrees" begin
    # A near-isothermal column at the melting point, where c_p(CtoK) = 2097.87 differs from
    # the constant 2102 by 0.2%: the two methods must land within that of each other. This
    # pins the enthalpy inverse against the constant path — a sign error or a missing factor
    # of 2 in `T = 2h/(a + √(a² + 2bh))` would show here as a gross disagreement, which
    # per-method conservation checks cannot see because each closes against its own map.
    n = 20
    dz = fill(0.1, n)
    density = fill(500.0, n)
    grain_radius = fill(0.3, n)
    shortwave_flux = zeros(n)
    cfs = GEMB.ClimateForcingStep(3600.0, 272.0, 100000.0, 0.0, 3.0, 0.0,
        GEMB.SB * 272.0^4 * 0.97, 400.0, 272.0, 3.0, 0.0, 2.0, 2.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    _mp(method) = GEMB.ModelParameters(density_ice=917.0, emissivity=0.97,
        emissivity_method=:uniform, emissivity_grain_radius_large=0.97,
        emissivity_grain_radius_threshold=10.0, surface_roughness_effective_ratio=0.1,
        thermal_conductivity_method=:Sturm, heat_capacity_method=method,
        dt_divisors=Float64.(GEMB.fast_divisors(36000000)) ./ 10000)

    T0 = fill(272.0, n)
    Ta, = GEMB.calculate_temperature(copy(T0), dz, density, 0.0, grain_radius,
        shortwave_flux, cfs, _mp(:constant), false)
    Tb, = GEMB.calculate_temperature(copy(T0), dz, density, 0.0, grain_radius,
        shortwave_flux, cfs, _mp(:CuffeyPaterson), false)
    @test maximum(abs.(Ta .- Tb)) < 0.05
end

@testset "full model run closes its budgets under :CuffeyPaterson" begin
    # The end-to-end gate for the enthalpy refactor. `verbose=true` turns on all five
    # per-timestep conservation checks (thermal, melt, accumulation, layer management,
    # gemb_core) plus the whole-run budget in the driver, so any site still computing energy
    # as `M·T·c_p` instead of `M·∫c_p dT` throws rather than silently drifting. Two months of
    # forcing and no spinup: this is a correctness gate, not a physics pin.
    ds = simulate_climate_forcing("test_1", 3)
    cf_full = GEMB.initialize_forcing(ds)
    t = dims(cf_full, Ti)
    cf = cf_full[Ti=1:min(480, length(t))]      # ~2 months at 3-hourly

    for method in (:constant, :CuffeyPaterson)
        mp = ModelParameters(output_frequency=:daily, heat_capacity_method=method)
        profile = initialize_profile(mp, cf; constant_density=true,
            constant_temperature=true)
        output = gemb(profile, cf, mp; verbose=true)

        T = parent(output[:temperature])
        @test !any(isnan, T)
        @test all(200.0 .< T .<= GEMB.CtoK)
        @test sum(parent(output[:melt])) >= 0.0
        @test all(0.0 .<= parent(output[:albedo_broadband]) .<= 1.0)
    end
end
