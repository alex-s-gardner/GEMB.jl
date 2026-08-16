# Tests for calculate_density

# Helper to create a ClimateForcingStep for density tests
function _make_density_cfs(; dt=86400.0 * 30, precipitation_mean=200.0, temperature_air_mean=250.0)
    return GEMB.ClimateForcingStep(
        dt,             # dt [s]
        260.0,          # temperature_air
        80000.0,        # pressure_air
        0.0,            # precipitation
        5.0,            # wind_speed
        100.0,          # shortwave_downward
        250.0,          # longwave_downward
        300.0,          # vapor_pressure
        temperature_air_mean,  # temperature_air_mean
        5.0,            # wind_speed_mean
        precipitation_mean,    # precipitation_mean
        2.0,            # temperature_observation_height
        10.0,           # wind_observation_height
        0.0,            # black_carbon_snow
        0.0,            # black_carbon_ice
        0.0,            # cloud_optical_thickness
        0.0,            # solar_zenith_angle
        0.0,            # shortwave_downward_diffuse
        0.1,            # cloud_fraction
    )
end

@testset "Mass conservation (HerronLangway)" begin
    n = 10
    t_vec = 260.0 * ones(n)
    dz = 0.5 * ones(n)
    density = collect(range(300.0, 700.0, length=n))
    grain_radius = 0.5 * ones(n)

    mp = GEMB.ModelParameters(
        density_ice=917.0,
        densification_method=:HerronLangway,
        densification_coeffs_M01=:Ant_RACMO_GS_SW0,
    )
    cfs = _make_density_cfs()

    mass_initial = density .* dz
    density_before = copy(density)
    dz_before = copy(dz)

    (dz_out, density_out) = GEMB.calculate_density(t_vec, copy(dz), copy(density), grain_radius, zeros(length(dz)), cfs, mp)

    mass_final = density_out .* dz_out

    # Mass must be conserved during densification
    @test mass_final ≈ mass_initial atol = 1e-10

    # Density should increase over time
    @test all(density_out .> density_before)

    # Thickness should decrease over time
    @test all(dz_out .< dz_before)
end

@testset "Density clamping at ice density" begin
    n = 10
    t_vec = 260.0 * ones(n)
    dz = 0.5 * ones(n)
    grain_radius = 0.5 * ones(n)

    # Set density very close to ice density
    density_near_ice = fill(917.0 - 0.1, n)

    # Long timestep to force overshoot
    cfs = _make_density_cfs(dt=86400.0 * 365 * 100)

    mp = GEMB.ModelParameters(
        density_ice=917.0,
        densification_method=:HerronLangway,
        densification_coeffs_M01=:Ant_RACMO_GS_SW0,
    )

    (_, density_out) = GEMB.calculate_density(t_vec, copy(dz), copy(density_near_ice), grain_radius, zeros(length(dz)), cfs, mp)

    # Density must be clamped at density_ice
    @test all(density_out .<= 917.0 + 1e-10)
end

@testset "HerronLangway densification" begin
    n = 10
    t_vec = 260.0 * ones(n)
    dz = 0.5 * ones(n)
    density = collect(range(300.0, 700.0, length=n))
    grain_radius = 0.5 * ones(n)

    mp = GEMB.ModelParameters(
        density_ice=917.0,
        densification_method=:HerronLangway,
        densification_coeffs_M01=:Ant_RACMO_GS_SW0,
    )
    cfs = _make_density_cfs()

    density_before = copy(density)
    (_, density_out) = GEMB.calculate_density(t_vec, copy(dz), copy(density), grain_radius, zeros(length(dz)), cfs, mp)

    @test all(density_out .> density_before)
end

@testset "Arthern densification" begin
    n = 10
    t_vec = 260.0 * ones(n)
    dz = 0.5 * ones(n)
    density = collect(range(300.0, 700.0, length=n))
    grain_radius = 0.5 * ones(n)

    mp = GEMB.ModelParameters(
        density_ice=917.0,
        densification_method=:Arthern,
        densification_coeffs_M01=:Ant_RACMO_GS_SW0,
    )
    cfs = _make_density_cfs()

    density_before = copy(density)
    (_, density_out) = GEMB.calculate_density(t_vec, copy(dz), copy(density), grain_radius, zeros(length(dz)), cfs, mp)

    @test all(density_out .> density_before)
end

@testset "ArthernB grain size sensitivity" begin
    n = 10
    t_vec = 260.0 * ones(n)
    dz = 0.5 * ones(n)
    density = collect(range(300.0, 700.0, length=n))
    grain_radius = 0.5 * ones(n)

    mp = GEMB.ModelParameters(
        density_ice=917.0,
        densification_method=:ArthernB,
        densification_coeffs_M01=:Ant_RACMO_GS_SW0,
    )
    cfs = _make_density_cfs()

    # Standard grain size
    (_, density_std) = GEMB.calculate_density(t_vec, copy(dz), copy(density), grain_radius, zeros(length(dz)), cfs, mp)

    # Larger grain size -> slower densification
    grain_radius_large = grain_radius * 2
    (_, density_large) = GEMB.calculate_density(t_vec, copy(dz), copy(density), grain_radius_large, zeros(length(dz)), cfs, mp)

    # ArthernB: rate is proportional to 1/r^2
    # Larger grains -> slower densification -> lower final density
    # Exclude top layer (overburden = 0) and values near ice density
    valid_mask = (density_std .< 917.0 - 1.0)
    check_indices = findall(valid_mask)
    filter!(i -> i != 1, check_indices)

    if !isempty(check_indices)
        @test all(density_std[check_indices] .> density_large[check_indices])
    end
end

@testset "Ligtenberg with different coefficients" begin
    n = 10
    t_vec = 260.0 * ones(n)
    dz = 0.5 * ones(n)
    density = collect(range(300.0, 700.0, length=n))
    grain_radius = 0.5 * ones(n)
    cfs = _make_density_cfs()

    density_before = copy(density)

    # Test Case 1: Standard RACMO
    mp1 = GEMB.ModelParameters(
        density_ice=917.0,
        densification_method=:Ligtenberg,
        densification_coeffs_M01=:Ant_RACMO_GS_SW0,
    )
    (_, d1) = GEMB.calculate_density(t_vec, copy(dz), copy(density), grain_radius, zeros(length(dz)), cfs, mp1)

    # Test Case 2: ERA5 variant (different M coefficients)
    mp2 = GEMB.ModelParameters(
        density_ice=917.0,
        densification_method=:Ligtenberg,
        densification_coeffs_M01=:Ant_ERA5_BF_SW1,
    )
    (_, d2) = GEMB.calculate_density(t_vec, copy(dz), copy(density), grain_radius, zeros(length(dz)), cfs, mp2)

    # Ensure densification occurred
    @test all(d1 .> density_before)
    @test all(d2 .> density_before)

    # Different coefficients should produce different results
    @test d1 != d2
end

@testset "KuipersMunneke coefficients" begin
    n = 10
    t_vec = 260.0 * ones(n)
    dz = 0.5 * ones(n)
    density = collect(range(300.0, 700.0, length=n))
    grain_radius = 0.5 * ones(n)
    cfs = _make_density_cfs()

    density_before = copy(density)
    mp = GEMB.ModelParameters(
        density_ice=917.0,
        densification_method=:Ligtenberg,
        densification_coeffs_M01=:Gre_KuipersMunneke,
    )

    (_, density_out) = GEMB.calculate_density(t_vec, copy(dz), copy(density), grain_radius, zeros(length(dz)), cfs, mp)

    @test all(density_out .> density_before)
end

@testset "Zero time step (no change)" begin
    n = 10
    t_vec = 260.0 * ones(n)
    dz = 0.5 * ones(n)
    density = collect(range(300.0, 700.0, length=n))
    grain_radius = 0.5 * ones(n)

    cfs = _make_density_cfs(dt=0.0)

    mp = GEMB.ModelParameters(
        density_ice=917.0,
        densification_method=:HerronLangway,
        densification_coeffs_M01=:Ant_RACMO_GS_SW0,
    )

    density_before = copy(density)
    dz_before = copy(dz)
    (dz_out, density_out) = GEMB.calculate_density(t_vec, copy(dz), copy(density), grain_radius, zeros(length(dz)), cfs, mp)

    @test density_out ≈ density_before
    @test dz_out ≈ dz_before
end

@testset "Ligtenberg bare ice logic (density_ice = 820 vs 917)" begin
    n = 10
    t_vec = 260.0 * ones(n)
    dz = 0.5 * ones(n)
    density = collect(range(300.0, 700.0, length=n))
    grain_radius = 0.5 * ones(n)
    cfs = _make_density_cfs()

    # Case A: density_ice ~ 820 (specialized branch)
    mp_820 = GEMB.ModelParameters(
        density_ice=820.0,
        densification_method=:Ligtenberg,
        densification_coeffs_M01=:Gre_RACMO_GS_SW0,
    )
    (_, density_820) = GEMB.calculate_density(t_vec, copy(dz), copy(density), grain_radius, zeros(length(dz)), cfs, mp_820)

    # Case B: density_ice ~ 917 (standard branch)
    mp_917 = GEMB.ModelParameters(
        density_ice=917.0,
        densification_method=:Ligtenberg,
        densification_coeffs_M01=:Gre_RACMO_GS_SW0,
    )
    (_, density_917) = GEMB.calculate_density(t_vec, copy(dz), copy(density), grain_radius, zeros(length(dz)), cfs, mp_917)

    # Different density_ice branches should produce different results
    # Check only values well below 820
    check_idx = density .< 700
    @test density_820[check_idx] != density_917[check_idx]
end

@testset "ArthernB (Arthern et al. 2010 eq. B1)" begin
    # `:ArthernB` is the one stress-driven scheme: `c = kc·exp(-Ec/RT)·σ/r²` with
    # `σ = Σ(ρⱼdzⱼ + waterⱼ)·g`. The MATLAB reference omits gravity and the
    # per-second→per-year conversion (a combined factor of ~3.1e8, leaving the scheme
    # inert) and applies the overlying cell's density to the whole overlying depth
    # instead of integrating — see upstream GEMB issue #200.
    mp = GEMB.ModelParameters(density_ice=910.0, densification_method=:ArthernB)
    cfs = _make_density_cfs(dt=86400.0)     # 1 day

    dz = [0.05, 0.05, 0.05]
    density = [300.0, 350.0, 400.0]
    t_vec = [253.0, 253.0, 253.0]
    grain_radius = [0.5, 0.7, 0.9]          # mm
    water = zeros(3)

    (_, density_out) = GEMB.calculate_density(t_vec, copy(dz), copy(density),
        grain_radius, water, cfs, mp)
    Δρ = density_out .- density

    # The overburden integral is exclusive, so the top cell carries no load and
    # therefore cannot densify. This is the signature of a genuinely stress-driven law.
    @test Δρ[1] == 0.0
    @test all(Δρ[2:end] .> 0.0)

    # Hand-computed increment for cell 3, from σ = (ρ₁dz₁ + ρ₂dz₂)·g.
    load = density[1] * dz[1] + density[2] * dz[2]       # 32.5 kg m-2
    σ = load * GEMB.GRAVITY                              # 318.825 Pa
    r = grain_radius[3] / 1000                           # m
    H = exp(-60000.0 / (t_vec[3] * 8.314)) * σ / r^2
    c = 9.2e-9 * (365.0 * 86400.0) * H                   # ρ = 400 <= 550 -> kc1 [yr-1]
    # rtol reflects float reassociation between this expression and the loop's ordering,
    # not any physical difference.
    @test Δρ[3] ≈ c * (mp.density_ice - density[3]) / 365.0 rtol = 1e-9

    # Pore water adds to the load, so it must increase densification below it.
    (_, density_wet) = GEMB.calculate_density(t_vec, copy(dz), copy(density),
        grain_radius, [5.0, 5.0, 0.0], cfs, mp)
    @test (density_wet .- density)[3] > Δρ[3]

    # Mass conservation, as for the other schemes.
    (dz_out, d_out) = GEMB.calculate_density(t_vec, copy(dz), copy(density),
        grain_radius, water, cfs, mp)
    @test dz_out .* d_out ≈ dz .* density rtol = 1e-14
end

@testset "ArthernB is the same order as Arthern in deep firn" begin
    # Regression guard on the unit factors. Before the fix `:ArthernB` was ~3.1e8 too
    # small and produced no measurable densification; any sane bound catches that.
    n = 120
    dz = fill(0.5, n)
    z = cumsum(dz) .- 0.25
    density = [min(910.0, 350.0 + 560.0 * (1 - exp(-zz / 15))) for zz in z]
    t_vec = fill(253.0, n)
    grain_radius = [0.5 + 1.5 * (1 - exp(-zz / 25)) for zz in z]
    water = zeros(n)
    cfs = _make_density_cfs(dt=86400.0)

    rates = map((:Arthern, :ArthernB)) do meth
        mp = GEMB.ModelParameters(density_ice=910.0, densification_method=meth)
        (_, d) = GEMB.calculate_density(t_vec, copy(dz), copy(density),
            copy(grain_radius), copy(water), cfs, mp)
        (d .- density) .* 365.0     # per year
    end
    arthern, arthernb = rates

    # Below ~10 m the two agree to within an order of magnitude (measured ratio
    # 0.13-0.36 over 10-60 m). The pre-fix code gave ~1e-8 of Arthern.
    deep = 20:n
    @test all(arthernb[deep] .> 0.0)
    @test all(arthernb[deep] ./ arthern[deep] .> 0.05)
    @test all(arthernb[deep] ./ arthern[deep] .< 1.0)
end

@testset "ArthernB matches the Community Firn Model" begin
    # Independent cross-check against CFM's `Arthern2010T`, which uses the same
    # `kc·(ρᵢ-ρ)·exp(-Ec/RT)·σ/r²` form with `σ = cumsum((mass + LWC·ρ_w)·g)` and
    # applies `drho_dt` in kg m-3 s-1. Agreement is far tighter than the ~0.2%
    # g/year-length tolerance `:Arthern` meets, because neither factor differs here.
    mp = GEMB.ModelParameters(density_ice=910.0, densification_method=:ArthernB)
    cfs = _make_density_cfs(dt=86400.0)

    dz = [0.05, 0.05, 0.05]
    density = [300.0, 350.0, 400.0]
    t_vec = fill(253.0, 3)
    grain_radius = [0.5, 0.7, 0.9]

    (_, density_out) = GEMB.calculate_density(t_vec, copy(dz), copy(density),
        grain_radius, zeros(3), cfs, mp)

    # CFM, evaluated at the same (ρ, T, σ, r²) for cell 3:
    #   drho_dt = kc1*(RHO_I-rho)*exp(-Ec/(R*Tz))*sigma/r2   [kg m-3 s-1]
    σ = (density[1] * dz[1] + density[2] * dz[2]) * GEMB.GRAVITY
    cfm_drho_dt = 9.2e-9 * (mp.density_ice - density[3]) *
        exp(-60.0e3 / (8.314 * t_vec[3])) * σ / (0.9e-3)^2
    @test (density_out[3] - density[3]) ≈ cfm_drho_dt * 86400.0 rtol = 1e-8
end

@testset "GSFC2020 (Medley et al. 2022 eq. 18)" begin
    # Recalibrated Arthern: same `dρ/dt = c(ρᵢ − ρ)` relaxation, but the accumulation is
    # raised to a fitted exponent (α0 = 0.91, α1 = 0.644) and each stage carries its own
    # activation energy (59500 / 56870 rather than Arthern's single 60000).
    mp = GEMB.ModelParameters(density_ice=910.0, densification_method=:GSFC2020)
    pm, tam = 300.0, 250.0
    cfs = _make_density_cfs(dt=86400.0, precipitation_mean=pm, temperature_air_mean=tam)

    dz = fill(0.1, 4)
    density = [300.0, 500.0, 700.0, 850.0]
    t_vec = fill(250.0, 4)

    (dz_out, density_out) = GEMB.calculate_density(t_vec, copy(dz), copy(density),
        fill(0.5, 4), zeros(4), cfs, mp)

    # Densifies everywhere, and mass is conserved cell by cell.
    @test all(density_out .> density)
    @test all(density_out .< mp.density_ice)
    for i in 1:4
        @test density_out[i] * dz_out[i] ≈ density[i] * dz[i] rtol = 1e-14
    end

    # Hand-evaluated eq. 18 for cell 1 (stage 1, ρ <= 550) and cell 3 (stage 2).
    R = 8.314
    Eg = 42400.0 / (R * tam)
    c_stage1 = 0.07 * pm^0.91 * GEMB.GRAVITY * exp(-59500.0 / (R * 250.0) + Eg)
    c_stage2 = 0.03 * pm^0.644 * GEMB.GRAVITY * exp(-56870.0 / (R * 250.0) + Eg)
    # `_densify_cell!` applies c as yr-1 via c/365*dt with dt in days.
    # rtol 1e-11, not tighter: the code hoists `pm^α · g` into a single term while this
    # test groups it as `pm^α * g`, which differs in the last bit or two.
    @test density_out[1] - density[1] ≈
        c_stage1 * (mp.density_ice - density[1]) / 365 rtol = 1e-11
    @test density_out[3] - density[3] ≈
        c_stage2 * (mp.density_ice - density[3]) / 365 rtol = 1e-11

    # The 550 threshold is where the stage switches, and it is inclusive of 550 (stage 1),
    # matching the `<= 550 + D_TOLERANCE` convention the other schemes use.
    d_lo = [550.0]; d_hi = [550.0 + 1e-6]
    (_, out_lo) = GEMB.calculate_density([250.0], [0.1], copy(d_lo), [0.5], [0.0], cfs, mp)
    (_, out_hi) = GEMB.calculate_density([250.0], [0.1], copy(d_hi), [0.5], [0.0], cfs, mp)
    @test out_lo[1] - d_lo[1] ≈ c_stage1 * (mp.density_ice - 550.0) / 365 rtol = 1e-11
    @test (out_hi[1] - d_hi[1]) / (out_lo[1] - d_lo[1]) < 0.5   # stage 2 is much slower here
end

@testset "GSFC2020 matches the Community Firn Model" begin
    # Cross-check against CFM's `GSFC2020` (CFM_main/physics.py), which applies
    # `drho_dt = c*(RHO_I - rho)` in kg m-3 s-1 with the same eq. 18 coefficients.
    # CFM uses g = 9.8 against GEMB's 9.81, so agreement is expected at the 0.102%
    # level — the same known g offset `:Arthern` and `:HerronLangway` are validated to.
    mp = GEMB.ModelParameters(density_ice=917.0, densification_method=:GSFC2020)
    pm, tam = 300.0, 250.0
    cfs = _make_density_cfs(dt=86400.0, precipitation_mean=pm, temperature_air_mean=tam)

    # CFM reference, computed from physics.py with g = 9.8 (see test comment above):
    #   rho=400, T=250, Tmean=250, A=300  ->  c = 3.292501944424200e-02 yr-1
    #   rho=700, same                     ->  c = 1.096890463217807e-02 yr-1
    for (ρ0, c_cfm) in ((400.0, 3.292501944424200e-02), (700.0, 1.096890463217807e-02))
        (_, out) = GEMB.calculate_density([250.0], [0.1], [ρ0], [0.5], [0.0], cfs, mp)
        c_gemb = (out[1] - ρ0) * 365 / (mp.density_ice - ρ0)
        # Once g is matched the two agree to round-off (the reference is a decimal literal
        # and the groupings differ, so ~1e-11, not bit-exact).
        @test c_gemb ≈ c_cfm * (GEMB.GRAVITY / 9.8) rtol = 1e-10
        @test c_gemb ≈ c_cfm rtol = 1.1e-3   # 0.102% apart on g alone
    end
end

@testset "GSFC2020 differs from Arthern in the expected direction" begin
    # The refit is not cosmetic: at low accumulation the fractional exponent makes GSFC2020
    # densify faster than Arthern (whose accumulation enters linearly), and at high
    # accumulation slower. The crossover is where pm^α·g == pm·g, i.e. pm = 1 kg m-2 yr-1
    # for stage 1 — so across all realistic accumulations stage 1 is *slower* than Arthern.
    # This test pins the direction rather than a magnitude, since only the sign is robust.
    _c(method, pm) = begin
        mp = GEMB.ModelParameters(density_ice=910.0, densification_method=method)
        cfs = _make_density_cfs(dt=86400.0, precipitation_mean=pm, temperature_air_mean=250.0)
        (_, out) = GEMB.calculate_density([250.0], [0.1], [400.0], [0.5], [0.0], cfs, mp)
        (out[1] - 400.0) * 365 / (910.0 - 400.0)
    end
    for pm in (50.0, 200.0, 800.0)
        @test _c(:GSFC2020, pm) < _c(:Arthern, pm)      # α = 0.91 < 1 and Ec is larger
    end
    # Accumulation dependence is sublinear, unlike Arthern's exactly linear one.
    @test _c(:GSFC2020, 400.0) / _c(:GSFC2020, 100.0) ≈ 4.0^0.91 rtol = 1e-10
    @test _c(:Arthern, 400.0) / _c(:Arthern, 100.0) ≈ 4.0 rtol = 1e-10
end

@testset "GSFC2020 zero and negative mean accumulation" begin
    # `pm^0.91` is NaN for pm < 0 and the steady-state marcher can be handed a
    # non-positive effective accumulation on an ablation column, so the power is guarded.
    mp = GEMB.ModelParameters(density_ice=910.0, densification_method=:GSFC2020)
    for pm in (0.0, -50.0)
        cfs = _make_density_cfs(dt=86400.0, precipitation_mean=pm, temperature_air_mean=250.0)
        (_, out) = GEMB.calculate_density([250.0], [0.1], [400.0], [0.5], [0.0], cfs, mp)
        @test !isnan(out[1])
        @test out[1] == 400.0     # no accumulation, no densification
    end
end

@testset "GSFC2020 validator and steady-state march" begin
    @test initialize_parameters(densification_method=:GSFC2020).densification_method ===
        :GSFC2020
    @test_throws AssertionError initialize_parameters(densification_method=:GSFC2021)

    # Unlike :ArthernB, GSFC2020 is accumulation-driven, so the steady-state marcher uses
    # the real law rather than falling back to :Arthern. Assert the dispatch actually
    # reaches `_gsfc2020_c` rather than silently taking the :Arthern branch.
    mp = GEMB.ModelParameters(density_ice=910.0, densification_method=:GSFC2020)
    p = GEMB.DensificationCoeffs(mp, 300.0, 250.0)
    @test !isnan(p.k0)
    @test p.k0 ≈ 300.0^0.91 * GEMB.GRAVITY rtol = 1e-14
    @test GEMB._densification_rate(p, 400.0, 250.0) ≈
        GEMB._gsfc2020_c(400.0, 250.0, 250.0, p.k0, p.k1) rtol = 1e-15
    # ...and that it is not the Arthern rate.
    mp_a = GEMB.ModelParameters(density_ice=910.0, densification_method=:Arthern)
    p_a = GEMB.DensificationCoeffs(mp_a, 300.0, 250.0)
    @test GEMB._densification_rate(p, 400.0, 250.0) !=
        GEMB._densification_rate(p_a, 400.0, 250.0)
end

@testset "Simonsen2013 (Simonsen et al. 2013)" begin
    # Arthern's form and activation energies (60000 / 42400) retuned for Greenland:
    # a constant F0 = 0.8 below 550 kg m-3, and F1·γ above with γ = 61.7/sqrt(A)·exp(-3800/RT̄).
    mp = GEMB.ModelParameters(density_ice=910.0, densification_method=:Simonsen2013)
    pm, tam = 300.0, 250.0
    cfs = _make_density_cfs(dt=86400.0, precipitation_mean=pm, temperature_air_mean=tam)

    dz = fill(0.1, 4)
    density = [300.0, 500.0, 700.0, 850.0]
    t_vec = fill(250.0, 4)

    (dz_out, density_out) = GEMB.calculate_density(t_vec, copy(dz), copy(density),
        fill(0.5, 4), zeros(4), cfs, mp)

    @test all(density_out .> density)
    @test all(density_out .< mp.density_ice)
    for i in 1:4
        @test density_out[i] * dz_out[i] ≈ density[i] * dz[i] rtol = 1e-14
    end

    # Hand-evaluated for cell 1 (stage 1) and cell 3 (stage 2). Note γ multiplies stage 2
    # only — stage 1 carries F0 alone.
    R = 8.314
    E = exp(-60000.0 / (R * 250.0) + 42400.0 / (R * tam))
    γ = 61.7 / sqrt(pm) * exp(-3800.0 / (R * tam))
    c_stage1 = 0.8 * 0.07 * pm * GEMB.GRAVITY * E
    c_stage2 = 1.25 * γ * 0.03 * pm * GEMB.GRAVITY * E
    # rtol 1e-11, not tighter: the code hoists `F1·γ` and `pm·g` into single terms while this
    # test groups the factors differently, which differs in the last bit or two.
    @test density_out[1] - density[1] ≈
        c_stage1 * (mp.density_ice - density[1]) / 365 rtol = 1e-11
    @test density_out[3] - density[3] ≈
        c_stage2 * (mp.density_ice - density[3]) / 365 rtol = 1e-11

    # 550 is inclusive of stage 1, matching the other schemes' `<= 550 + D_TOLERANCE`.
    d_lo = [550.0]; d_hi = [550.0 + 1e-6]
    (_, out_lo) = GEMB.calculate_density([250.0], [0.1], copy(d_lo), [0.5], [0.0], cfs, mp)
    (_, out_hi) = GEMB.calculate_density([250.0], [0.1], copy(d_hi), [0.5], [0.0], cfs, mp)
    @test out_lo[1] - d_lo[1] ≈ c_stage1 * (mp.density_ice - 550.0) / 365 rtol = 1e-11
    @test (out_hi[1] - d_hi[1]) / (out_lo[1] - d_lo[1]) < 0.5   # stage 2 is slower here
end

@testset "Simonsen2013 matches the FirnMICE eq. A36-A37 form" begin
    # Independent structural check against Lundin et al. (2017) FirnMICE, which states the
    # scheme as scalar tunings of the Arthern coefficients:
    #   c0(SIM) = f0 · c0(ART)
    #   c1(SIM) = f1 · (61.7/ḃ^0.5) · exp(-3800/(R·Ta)) · c1(ART)
    # This pins the two properties the appendix makes unambiguous and that CFM's vectorized
    # masks could plausibly have transposed: stage 1 is a bare scalar multiple of Arthern with
    # no γ, and γ is a function of the *mean annual* temperature, not the layer temperature.
    R = 8.314
    _c(method, ρ0, pm, T_layer, tam) = begin
        mp = GEMB.ModelParameters(density_ice=910.0, densification_method=method)
        cfs = _make_density_cfs(dt=86400.0, precipitation_mean=pm, temperature_air_mean=tam)
        (_, out) = GEMB.calculate_density([T_layer], [0.1], [ρ0], [0.5], [0.0], cfs, mp)
        (out[1] - ρ0) * 365 / (910.0 - ρ0)
    end

    # Stage 1 carries no γ: the ratio to Arthern is f0 alone, invariant to both the
    # accumulation and the mean temperature that γ depends on.
    for pm in (100.0, 500.0), tam in (240.0, 260.0)
        @test _c(:Simonsen2013, 400.0, pm, 250.0, tam) /
            _c(:Arthern, 400.0, pm, 250.0, tam) ≈ 0.8 rtol = 1e-10
    end

    # Stage 2 γ uses the mean annual temperature, not the layer temperature. Hold the layer
    # temperature fixed and vary `tam`: the ratio to Arthern must track exp(-3800/(R·tam))
    # exactly (Arthern's own Eg/(R·tam) cancels in the ratio).
    ratio(tam) = _c(:Simonsen2013, 700.0, 300.0, 250.0, tam) /
                 _c(:Arthern, 700.0, 300.0, 250.0, tam)
    @test ratio(260.0) / ratio(240.0) ≈
        exp(-3800.0 / (R * 260.0)) / exp(-3800.0 / (R * 240.0)) rtol = 1e-10
    # ...and conversely the ratio is independent of the layer temperature.
    @test _c(:Simonsen2013, 700.0, 300.0, 230.0, 250.0) /
        _c(:Arthern, 700.0, 300.0, 230.0, 250.0) ≈
        _c(:Simonsen2013, 700.0, 300.0, 265.0, 250.0) /
        _c(:Arthern, 700.0, 300.0, 265.0, 250.0) rtol = 1e-10
end

@testset "Simonsen2013 matches the Community Firn Model" begin
    # Cross-check against CFM's `Simonsen2013` (CFM_main/physics.py, `bdot_type = 'mean'`),
    # which applies `drho_dt = c*(RHO_I - rho)` in kg m-3 s-1. CFM uses g = 9.8 against
    # GEMB's 9.81, so agreement is expected at the 0.102% level.
    mp = GEMB.ModelParameters(density_ice=917.0, densification_method=:Simonsen2013)
    pm, tam = 300.0, 250.0
    cfs = _make_density_cfs(dt=86400.0, precipitation_mean=pm, temperature_air_mean=tam)

    # CFM reference, computed from physics.py with g = 9.8:
    #   rho=400, T=250, Tmean=250, A=300  ->  c = 3.460062047878578e-02 yr-1  (stage 1, F0)
    #   rho=700, same                     ->  c = 1.326344865244790e-02 yr-1  (stage 2, F1·γ)
    for (ρ0, c_cfm) in ((400.0, 3.460062047878578e-02), (700.0, 1.326344865244790e-02))
        (_, out) = GEMB.calculate_density([250.0], [0.1], [ρ0], [0.5], [0.0], cfs, mp)
        c_gemb = (out[1] - ρ0) * 365 / (mp.density_ice - ρ0)
        @test c_gemb ≈ c_cfm * (GEMB.GRAVITY / 9.8) rtol = 1e-10
        @test c_gemb ≈ c_cfm rtol = 1.1e-3   # 0.102% apart on g alone
    end
end

@testset "Simonsen2013 is Arthern scaled by the stage factors" begin
    # Simonsen shares Arthern's activation energies and linear accumulation dependence, so
    # the ratio to `:Arthern` is exactly F0 in stage 1 and exactly F1·γ in stage 2 — an
    # identity that fails if either factor is applied to the wrong stage.
    _c(method, ρ0, pm) = begin
        mp = GEMB.ModelParameters(density_ice=910.0, densification_method=method)
        cfs = _make_density_cfs(dt=86400.0, precipitation_mean=pm, temperature_air_mean=250.0)
        (_, out) = GEMB.calculate_density([250.0], [0.1], [ρ0], [0.5], [0.0], cfs, mp)
        (out[1] - ρ0) * 365 / (910.0 - ρ0)
    end
    γ(pm) = 61.7 / sqrt(pm) * exp(-3800.0 / (8.314 * 250.0))
    for pm in (50.0, 200.0, 800.0)
        @test _c(:Simonsen2013, 400.0, pm) / _c(:Arthern, 400.0, pm) ≈ 0.8 rtol = 1e-10
        @test _c(:Simonsen2013, 700.0, pm) / _c(:Arthern, 700.0, pm) ≈
            1.25 * γ(pm) rtol = 1e-10
    end
    # Accumulation enters linearly in both stages (γ's A^-1/2 makes stage 2 sqrt-like).
    @test _c(:Simonsen2013, 400.0, 400.0) / _c(:Simonsen2013, 400.0, 100.0) ≈ 4.0 rtol = 1e-10
    @test _c(:Simonsen2013, 700.0, 400.0) / _c(:Simonsen2013, 700.0, 100.0) ≈
        sqrt(4.0) rtol = 1e-10
end

@testset "Simonsen2013 zero and negative mean accumulation" begin
    # γ ∝ 1/sqrt(pm) is Inf at pm = 0 and NaN below it, and the rate also carries a factor
    # of pm, so the stage-2 product would be 0·Inf = NaN without the guard.
    #
    # Only γ is guarded, not `precip_force`: for pm < 0 the rate stays exactly 0.8 (stage 1)
    # or 1.25·γ (stage 2) times `:Arthern`, which is itself negative there. Clamping the
    # accumulation would break that identity and make Simonsen disagree with the scheme it
    # rescales. A negative mean accumulation is an ablation column, not a valid climatology
    # for either scheme.
    mp = GEMB.ModelParameters(density_ice=910.0, densification_method=:Simonsen2013)
    mp_a = GEMB.ModelParameters(density_ice=910.0, densification_method=:Arthern)
    for ρ0 in (400.0, 700.0)
        # pm = 0: no accumulation, no densification, and γ's Inf does not leak through.
        cfs0 = _make_density_cfs(dt=86400.0, precipitation_mean=0.0, temperature_air_mean=250.0)
        (_, out0) = GEMB.calculate_density([250.0], [0.1], [ρ0], [0.5], [0.0], cfs0, mp)
        @test !isnan(out0[1])
        @test out0[1] == ρ0

        # pm < 0: finite, and still the Arthern rate times the stage factor.
        cfs_n = _make_density_cfs(dt=86400.0, precipitation_mean=-50.0, temperature_air_mean=250.0)
        (_, out_n) = GEMB.calculate_density([250.0], [0.1], [ρ0], [0.5], [0.0], cfs_n, mp)
        (_, out_a) = GEMB.calculate_density([250.0], [0.1], [ρ0], [0.5], [0.0], cfs_n, mp_a)
        @test isfinite(out_n[1])
        expected = ρ0 <= 550.0 ? 0.8 : 0.0     # γ = 0 for non-positive pm
        @test (out_n[1] - ρ0) ≈ expected * (out_a[1] - ρ0) rtol = 1e-10
    end
end

@testset "Simonsen2013 validator and steady-state march" begin
    @test initialize_parameters(densification_method=:Simonsen2013).densification_method ===
        :Simonsen2013
    @test_throws AssertionError initialize_parameters(densification_method=:Simonsen2014)

    # Accumulation-driven, so the steady-state marcher uses the real law rather than falling
    # back to :Arthern the way :ArthernB does.
    mp = GEMB.ModelParameters(density_ice=910.0, densification_method=:Simonsen2013)
    p = GEMB.DensificationCoeffs(mp, 300.0, 250.0)
    @test p.k0 == 0.8
    @test p.k1 ≈
        1.25 * 61.7 / sqrt(300.0) * exp(-3800.0 / (8.314 * 250.0)) rtol = 1e-14
    @test GEMB._densification_rate(p, 700.0, 250.0) ≈
        GEMB._arthern_scaled_c(700.0, 250.0, 250.0, p.precip_force, p.k0, p.k1) rtol = 1e-15
    mp_a = GEMB.ModelParameters(density_ice=910.0, densification_method=:Arthern)
    p_a = GEMB.DensificationCoeffs(mp_a, 300.0, 250.0)
    @test GEMB._densification_rate(p, 700.0, 250.0) !=
        GEMB._densification_rate(p_a, 700.0, 250.0)
end

@testset "Crocus viscosity (Vionnet et al. 2012 eqs. 7-9)" begin
    # `η = f1·f2·η₀(ρ/cη)·exp(aη(T_fus−T) + bη·ρ)`, evaluated against the paper's
    # coefficients directly rather than against the loop.
    ρ, T, D = 300.0, 263.15, 0.05
    r_fine = 0.05                              # mm; f2 = 1 at gs <= 0.2 mm
    η_dry = GEMB._crocus_viscosity(ρ, T, 0.0, D, r_fine)
    @test η_dry ≈ 7.62237e6 * (300.0 / 250.0) * exp(0.1 * 10.0 + 0.023 * 300.0) rtol = 1e-14

    # Eq. 8: f1 = 1/(1 + 60·W/(ρ_w·D)). 2.5 kg m-2 in a 0.05 m cell is a fully
    # saturated pore space by volume, so f1 bottoms out near 1/61.
    W = 2.5
    f1 = 1.0 / (1.0 + 60.0 * W / (1000.0 * D))
    @test f1 ≈ 1 / 4 rtol = 1e-14
    @test GEMB._crocus_viscosity(ρ, T, W, D, r_fine) ≈ f1 * η_dry rtol = 1e-14

    # Monotone in water: wetter is always weaker, never stronger.
    ηs = [GEMB._crocus_viscosity(ρ, T, w, D, r_fine) for w in 0.0:0.25:5.0]
    @test issorted(ηs, rev=true)

    # Eq. 9: the exponential is 1 at gs = 0.2 mm (a grain *diameter*, so radius 0.1 mm),
    # rises with grain size, and is capped at 4.
    @test GEMB._crocus_viscosity(ρ, T, 0.0, D, 0.1) ≈ η_dry rtol = 1e-14
    @test GEMB._crocus_viscosity(ρ, T, 0.0, D, 0.12) ≈ exp(0.4) * η_dry rtol = 1e-14
    @test GEMB._crocus_viscosity(ρ, T, 0.0, D, 2.5) ≈ 4.0 * η_dry rtol = 1e-14

    # Below gs = 0.2 mm eq. 9 falls under 1 and would *soften* the snow. The floor keeps it
    # a stiffening correction across eq. 9's valid domain and inert below it — see
    # `_crocus_viscosity` for why the paper never has to bound this.
    @test GEMB._crocus_viscosity(ρ, T, 0.0, D, GEMB.RE_NEW_SNOW) ≈ η_dry rtol = 1e-14
    @test GEMB._crocus_viscosity(ρ, T, 0.0, D, 0.01) ≈ η_dry rtol = 1e-14
    # Pin what the floor is suppressing: at the fresh-snow radius the raw eq. 9 is exp(-1),
    # a 2.7× softening. If the floor is ever removed this is the factor that appears.
    f2_raw(gs) = exp(min(GEMB.CROCUS_G1, gs - GEMB.CROCUS_G2) / GEMB.CROCUS_G3)
    @test f2_raw(2 * GEMB.RE_NEW_SNOW / 1000) ≈ exp(-1.0) rtol = 1e-14
    @test f2_raw(2 * GEMB.RE_NEW_SNOW / 1000) ≈ 0.3679 atol = 1e-4
    # And that the paper's own non-dendritic range never enters that regime: eq. 9 is
    # 2.7 at gs = 0.3 mm and at its cap of 4 by 0.339 mm.
    @test f2_raw(3.0e-4) ≈ exp(1.0) rtol = 1e-14
    @test f2_raw(3.39e-4) > GEMB.CROCUS_F2_MAX

    # Monotone non-decreasing in grain size across the whole GEMB range.
    η_gs = [GEMB._crocus_viscosity(ρ, T, 0.0, D, r) for r in 0.01:0.01:GEMB.GRAIN_RADIUS_ICE]
    @test issorted(η_gs)

    # Colder and denser snow is stiffer, in both cases by the exponential.
    @test GEMB._crocus_viscosity(ρ, 233.15, 0.0, D, r_fine) >
          GEMB._crocus_viscosity(ρ, 273.15, 0.0, D, r_fine)
    @test GEMB._crocus_viscosity(600.0, T, 0.0, D, r_fine) > η_dry
end

@testset "CrocusPure densification" begin
    mp = GEMB.ModelParameters(density_ice=910.0, densification_method=:CrocusPure)
    cfs = _make_density_cfs(dt=86400.0)     # 1 day

    dz = [0.05, 0.05, 0.05]
    density = [300.0, 350.0, 400.0]
    t_vec = fill(263.15, 3)
    grain_radius = [0.05, 0.07, 0.09]       # mm, all below eq. 9's f2 = 1 threshold
    water = zeros(3)

    (dz_out, density_out) = GEMB.calculate_density(t_vec, copy(dz), copy(density),
        grain_radius, water, cfs, mp)
    Δρ = density_out .- density

    # Unlike `:ArthernB`, the surface cell densifies: eq. 6's half-own-weight rule gives
    # it a nonzero stress.
    @test all(Δρ .> 0.0)

    # Hand-computed increment for cell 2 from the paper's equations:
    # σ = (ρ₁dz₁ + 0.5ρ₂dz₂)·g, dz' = dz·exp(-σ/η·dt).
    σ = (density[1] * dz[1] + 0.5 * density[2] * dz[2]) * GEMB.GRAVITY
    η = GEMB._crocus_viscosity(density[2], t_vec[2], 0.0, dz[2], grain_radius[2])
    dz_expected = dz[2] * exp(-σ / η * 86400.0)
    @test dz_out[2] ≈ dz_expected rtol = 1e-12
    @test density_out[2] ≈ density[2] * dz[2] / dz_expected rtol = 1e-12

    # Mass conservation.
    @test dz_out .* density_out ≈ dz .* density rtol = 1e-14
end

@testset "CrocusPure liquid water accelerates densification" begin
    # The reason this scheme exists: no other GEMB scheme lets pore water weaken the
    # matrix. `:ArthernB` reads `water` only as overburden, so it densifies the cells
    # *below* the wet one; Crocus densifies the wet cell itself.
    dz = fill(0.05, 3)
    density = [400.0, 400.0, 400.0]
    t_vec = fill(273.15, 3)                 # wet snow is at the melting point
    grain_radius = fill(0.05, 3)
    cfs = _make_density_cfs(dt=86400.0)

    rates = map((zeros(3), [0.0, 2.0, 0.0])) do water
        mp = GEMB.ModelParameters(density_ice=910.0, densification_method=:CrocusPure)
        (_, d) = GEMB.calculate_density(t_vec, copy(dz), copy(density),
            copy(grain_radius), water, cfs, mp)
        d .- density
    end
    dry, wet = rates

    # Cell 2 holds the water: f1 = 1/(1+60·2/(1000·0.05)) = 1/3.4 softens it by 3.4×, and
    # its own half-weight rises with the water from 30g to 31g. In the small-strain limit
    # the strain rate `σ/η` therefore rises by 3.4·31/30.
    @test wet[2] > dry[2]
    @test wet[2] / dry[2] ≈ 3.4 * 31 / 30 rtol = 1e-3

    # Cell 3 sees more overburden from the water above, so it too speeds up — but only
    # slightly, since the load rises by 2 against 20 kg m-2.
    @test wet[3] > dry[3]
    @test wet[3] / dry[3] < 1.2

    # Cell 1 is above the water and unaffected.
    @test wet[1] == dry[1]
end

@testset "CrocusPure is physically plausible in deep firn" begin
    # Guard on the unit chain (η in kg s-1 m-1, σ in Pa, dt in seconds): a missing
    # seconds-per-year factor either way is a 3.15e7 error and shows up immediately.
    n = 120
    dz = fill(0.5, n)
    z = cumsum(dz) .- 0.25
    density = [min(910.0, 350.0 + 560.0 * (1 - exp(-zz / 15))) for zz in z]
    t_vec = fill(253.0, n)
    grain_radius = [0.5 + 1.5 * (1 - exp(-zz / 25)) for zz in z]
    cfs = _make_density_cfs(dt=86400.0)

    rates = map((:Arthern, :CrocusPure)) do meth
        mp = GEMB.ModelParameters(density_ice=910.0, densification_method=meth)
        (_, d) = GEMB.calculate_density(t_vec, copy(dz), copy(density),
            copy(grain_radius), zeros(n), cfs, mp)
        (d .- density) .* 365.0     # per year
    end
    arthern, crocus = rates

    # Both must be positive and of a comparable order through the firn column.
    @test all(crocus .> 0.0)
    @test all(crocus[20:n] ./ arthern[20:n] .> 0.01)
    @test all(crocus[20:n] ./ arthern[20:n] .< 100.0)

    # Pin the published law's domain limit rather than leave it implicit: `exp(bη·ρ)` is
    # ~1e9 at ρ = 900, so Crocus compacts deep firn far more slowly than `:Arthern`. This
    # is the law as fitted (to a 1-2 m alpine snowpack), not a unit error, and it is why
    # `:Crocus` hands the deep column to `:GSFC2020` instead of using this branch.
    shallow = 1:6                               # top ~3 m, ρ <= 450
    @test all(0.3 .< crocus[shallow] ./ arthern[shallow] .< 1.0)
    deep = 40:n                                 # below 20 m, ρ >= 760
    @test all(crocus[deep] ./ arthern[deep] .< 0.1)
end

@testset "Crocus validator and steady-state fallback" begin
    for meth in (:Crocus, :CrocusPure)
        @test initialize_parameters(densification_method=meth).densification_method === meth
    end
    @test_throws AssertionError initialize_parameters(densification_method=:Crocus2012)

    # The steady-state marcher carries neither stress, grain radius, nor water, so both
    # Crocus variants fall back there. `:CrocusPure` falls back to `:Arthern`; `:Crocus`
    # falls back to `:GSFC2020`, which is already its own above-threshold branch.
    mp_a = GEMB.ModelParameters(density_ice=910.0, densification_method=:Arthern)
    mp_g = GEMB.ModelParameters(density_ice=910.0, densification_method=:GSFC2020)
    mp_c = GEMB.ModelParameters(density_ice=910.0, densification_method=:Crocus)
    mp_p = GEMB.ModelParameters(density_ice=910.0, densification_method=:CrocusPure)
    ps = map(m -> GEMB.DensificationCoeffs(m, 300.0, 250.0), (mp_a, mp_g, mp_c, mp_p))
    p_a, p_g, p_c, p_p = ps

    @test GEMB._densification_rate(p_p, 400.0, 250.0) ==
        GEMB._densification_rate(p_a, 400.0, 250.0)
    @test GEMB._densification_rate(p_c, 400.0, 250.0) ==
        GEMB._densification_rate(p_g, 400.0, 250.0)
    # ...and the two fallbacks are genuinely different laws, so the dispatch above is not
    # vacuously true.
    @test GEMB._densification_rate(p_c, 400.0, 250.0) !=
        GEMB._densification_rate(p_p, 400.0, 250.0)

    # `:Crocus` therefore needs GSFC2020's hoisted accumulation powers; `:CrocusPure` does
    # not use `k0`/`k1` at all.
    @test p_c.k0 ≈ 300.0^0.91 * GEMB.GRAVITY rtol = 1e-14
    @test p_c.k1 ≈ 300.0^0.644 * GEMB.GRAVITY rtol = 1e-14
    @test p_p.k0 == 1.0 && p_p.k1 == 1.0
end

@testset "CrocusPure never exceeds the density of ice" begin
    # `D·exp(-σ/η·dt)` is unconditionally positive, so the clamp is the only thing
    # standing between a long step and ρ > ρᵢ.
    mp = GEMB.ModelParameters(density_ice=910.0, densification_method=:CrocusPure)
    cfs = _make_density_cfs(dt=86400.0 * 365 * 100)     # absurdly long step
    n = 5
    dz = fill(1.0, n)
    density = fill(890.0, n)
    (dz_out, d_out) = GEMB.calculate_density(fill(273.15, n), copy(dz), copy(density),
        fill(0.05, n), fill(1.0, n), cfs, mp)
    @test all(d_out .<= mp.density_ice)
    @test all(dz_out .> 0.0)
    @test dz_out .* d_out ≈ dz .* density rtol = 1e-14
end

@testset "Crocus hybrid handover at CROCUS_HYBRID_DENSITY" begin
    # `:Crocus` is Crocus viscosity below the threshold and `:GSFC2020` at or above it.
    # Assert each cell equals the scheme it should have been handed to — exactly, since the
    # handover is a branch, not a blend.
    ρ_c = GEMB.CROCUS_HYBRID_DENSITY
    @test ρ_c == 450.0

    dz = fill(0.05, 4)
    density = [ρ_c - 50, ρ_c - 1e-9, ρ_c, ρ_c + 50]
    t_vec = fill(263.15, 4)
    grain_radius = fill(0.05, 4)
    water = zeros(4)
    cfs = _make_density_cfs(dt=86400.0)

    outs = map((:Crocus, :CrocusPure, :GSFC2020)) do meth
        mp = GEMB.ModelParameters(density_ice=910.0, densification_method=meth)
        GEMB.calculate_density(t_vec, copy(dz), copy(density),
            copy(grain_radius), copy(water), cfs, mp)
    end
    (dz_h, d_h), (dz_p, d_p), (dz_g, d_g) = outs

    # Below the threshold: the hybrid is the pure Crocus law.
    @test d_h[1:2] == d_p[1:2]
    @test dz_h[1:2] == dz_p[1:2]
    # At and above: the hybrid is GSFC2020. `:GSFC2020` ignores overburden, so its result
    # for these cells does not depend on what the cells above did.
    @test d_h[3:4] == d_g[3:4]
    @test dz_h[3:4] == dz_g[3:4]
    # The two branches genuinely differ here, so the equalities above are not vacuous.
    @test d_p[3:4] != d_g[3:4]

    # Mass conservation across the handover.
    @test dz_h .* d_h ≈ dz .* density rtol = 1e-14
end

@testset "Crocus hybrid keeps the wet-snow effect and fixes the deep column" begin
    n = 120
    dz = fill(0.5, n)
    z = cumsum(dz) .- 0.25
    density = [min(910.0, 350.0 + 560.0 * (1 - exp(-zz / 15))) for zz in z]
    t_vec = fill(253.0, n)
    grain_radius = [0.5 + 1.5 * (1 - exp(-zz / 25)) for zz in z]
    cfs = _make_density_cfs(dt=86400.0)

    rates = map((:Arthern, :Crocus, :CrocusPure)) do meth
        mp = GEMB.ModelParameters(density_ice=910.0, densification_method=meth)
        (_, d) = GEMB.calculate_density(t_vec, copy(dz), copy(density),
            copy(grain_radius), zeros(n), cfs, mp)
        (d .- density) .* 365.0
    end
    arthern, hybrid, pure = rates

    # The point of the hybrid: the deep column no longer stalls. `:CrocusPure` falls to
    # <0.1 of `:Arthern` below 20 m; the hybrid stays within a factor of a few.
    deep = 40:n
    @test all(hybrid[deep] .> pure[deep])
    @test all(0.2 .< hybrid[deep] ./ arthern[deep] .< 5.0)
    @test all(hybrid .> 0.0)

    # And the surface (ρ < 450, here the top ~2.5 m) is untouched by the handover, so the
    # wet-snow physics is still in force there.
    shallow = findall(<(GEMB.CROCUS_HYBRID_DENSITY), density)
    @test !isempty(shallow)
    @test hybrid[shallow] == pure[shallow]

    # Liquid water still weakens those surface cells under the hybrid.
    mp = GEMB.ModelParameters(density_ice=910.0, densification_method=:Crocus)
    water = zeros(n); water[2] = 5.0
    (_, d_wet) = GEMB.calculate_density(t_vec, copy(dz), copy(density),
        copy(grain_radius), water, cfs, mp)
    @test (d_wet[2] - density[2]) * 365.0 > hybrid[2]
end

@testset "Barnola f(ρ) (Barnola et al. 1991)" begin
    ρ_i = 910.0

    # Open-pore branch: the fitted polynomial, in g cm-3. Pinned against a direct
    # evaluation so a transcribed coefficient cannot drift silently.
    poly(ρ) = exp10(-37.455 * (ρ / 1000)^3 + 99.743 * (ρ / 1000)^2 -
                    95.027 * (ρ / 1000) + 30.673)
    for ρ in (560.0, 600.0, 700.0, 800.0)
        @test GEMB._barnola_f(ρ, ρ_i) ≈ poly(ρ) rtol = 1e-14
    end

    # Closed-pore branch: the analytic isolated-spherical-pore form.
    closed(ρ, ρi) = (3 / 16) * (1 - ρ / ρi) / (1 - (1 - ρ / ρi)^(1 / 3))^3
    for ρ in (810.0, 830.0, 870.0, 900.0)
        @test GEMB._barnola_f(ρ, ρ_i) ≈ closed(ρ, ρ_i) rtol = 1e-12
    end

    # The two branches cross at ρ_i ≈ 919.96, so the polynomial was fitted against an ice
    # density of ~920: at 920 the handover is continuous to 0.1%, and the gap at any other
    # `density_ice` measures the mismatch with the fit's own. Pin both ends so the
    # discontinuity at GEMB's default cannot grow unnoticed.
    f_open = GEMB._barnola_f(GEMB.BARNOLA_CLOSEOFF, ρ_i)
    @test abs(f_open - closed(GEMB.BARNOLA_CLOSEOFF, 920.0)) / f_open < 0.001
    @test abs(f_open - closed(GEMB.BARNOLA_CLOSEOFF, ρ_i)) / f_open < 0.15
    # The step is downward at the default ice density.
    @test closed(GEMB.BARNOLA_CLOSEOFF, ρ_i) < f_open
    # The paper states the polynomial was constructed to make both functions "and their first
    # derivatives equal for ρ = 0.8 g cm-3". Both conditions independently identify the fit's
    # ice density as ~920, which is what makes the 14% step at GEMB's 910 a ρ_i mismatch rather
    # than a transcription error. Solve each for the ρ_i that satisfies it.
    dρ = 1e-4
    d_poly = (poly(GEMB.BARNOLA_CLOSEOFF + dρ) - poly(GEMB.BARNOLA_CLOSEOFF - dρ)) / (2dρ)
    d_closed(ρi) = (closed(GEMB.BARNOLA_CLOSEOFF + dρ, ρi) -
                    closed(GEMB.BARNOLA_CLOSEOFF - dρ, ρi)) / (2dρ)
    function bisect(g, lo, hi)
        for _ in 1:200
            mid = (lo + hi) / 2
            g(mid) * g(lo) <= 0 ? (hi = mid) : (lo = mid)
        end
        return (lo + hi) / 2
    end
    ρi_value = bisect(ρi -> f_open - closed(GEMB.BARNOLA_CLOSEOFF, ρi), 900.0, 950.0)
    ρi_slope = bisect(ρi -> d_poly - d_closed(ρi), 900.0, 950.0)
    @test ρi_value ≈ 919.96 atol = 0.01
    @test ρi_slope ≈ 920.06 atol = 0.01
    # The paper's C1 claim holds only if the two agree; they do, to 0.1 kg m-3.
    @test abs(ρi_value - ρi_slope) < 0.15

    # 800, deliberately not GEMB's 830 pore-close-off constant. See BARNOLA_CLOSEOFF.
    @test GEMB.BARNOLA_CLOSEOFF == 800.0
    @test GEMB.BARNOLA_CLOSEOFF != GEMB.DENSITY_PORE_CLOSEOFF

    # Beyond the fit the polynomial diverges from the analytic form by ~5x, which is why
    # the handover exists at all rather than running the polynomial to ρ_i.
    @test poly(900.0) / closed(900.0, ρ_i) > 4.0

    # f vanishes with porosity, so the law self-limits as ρ → ρ_i without needing the clamp.
    # Swept over `density_ice` rather than pinned at the default: the property holds only
    # because `BARNOLA_MIN_DENSITY_ICE` keeps the handover strictly below ρ_i, so the
    # closed-pore branch is the one reached near ρ_i. Testing a single value hid that.
    for ρi in (GEMB.BARNOLA_MIN_DENSITY_ICE, 910.0, 917.0, 920.0, 950.0)
        @test GEMB._barnola_f(ρi, ρi) == 0.0
        @test GEMB._barnola_f(ρi - 1e-3, ρi) > 0.0
        @test GEMB._barnola_f(ρi - 1e-3, ρi) < 1e-5
        # Monotone decreasing in the closed-pore regime.
        fs = [GEMB._barnola_f(ρ, ρi) for ρ in 810.0:10.0:(ρi-10.0)]
        @test all(diff(fs) .< 0)
    end
    # Below the gate the property fails: at ρ_i = 800 the handover coincides with ρ_i, the
    # polynomial covers the whole column, and `f` never reaches 0. This is what the
    # `:Barnola1991` validator gate exists to exclude.
    @test GEMB._barnola_f(800.0, 800.0) > 0.2
    @test_throws AssertionError initialize_parameters(densification_method=:Barnola1991,
        density_ice=800.0)
    @test_throws AssertionError initialize_parameters(densification_method=:Barnola1991,
        density_ice=GEMB.BARNOLA_MIN_DENSITY_ICE - 1.0)
    # The gate is specific to `:Barnola1991`; the rest of the model still accepts 800.
    @test initialize_parameters(densification_method=:Arthern,
        density_ice=800.0).density_ice == 800.0
    @test initialize_parameters(densification_method=:Barnola1991,
        density_ice=GEMB.BARNOLA_MIN_DENSITY_ICE).density_ice ==
          GEMB.BARNOLA_MIN_DENSITY_ICE
    # The closed-pore branch follows the configured ice density; the fitted one cannot.
    @test GEMB._barnola_f(850.0, 917.0) != GEMB._barnola_f(850.0, 910.0)
    @test GEMB._barnola_f(700.0, 917.0) == GEMB._barnola_f(700.0, 910.0)

    # n = 3 throughout, as in the paper: A0 was fitted against n = 3 over 0.55-0.8 g cm-3, so
    # exponent and prefactor are one calibration. The paper's n = 1 remark applies to bulk ice
    # below close-off, not to low-stress firn. Pinned with the reasoning in BARNOLA_N.
    @test GEMB.BARNOLA_N == 3.0
    # Why a continuity-matched n=1 branch is not merely a missing option: A1 = A0*(1e5)^2, so
    # the rate would jump by (1e5/σ)² below the threshold — ~41x at the 15 kPa a shallow column
    # actually sees, against coefficients calibrated to reproduce observed profiles with n = 3.
    @test (1.0e5 / 1.55e4)^2 > 40.0
end

@testset "Barnola1991 handover at 550 is depth-dependent" begin
    # The 550 boundary joins a depth-independent accumulation-driven rate to a σ³ one, so the
    # step across it is set by the overburden, not by the densities. Pinned because the sign of
    # the step reverses with depth and a future change to the stress integral would move the
    # crossover silently. See the comment on the `:Barnola1991` branch.
    mp = GEMB.ModelParameters(density_ice=910.0, densification_method=:Barnola1991)
    cfs = GEMB.ClimateForcingStep(; dt=10800.0, precipitation_mean=300.0,
        temperature_air_mean=250.0)

    function rate_at(depth, ρ)
        dz = [depth, 0.1]
        dens = [530.0, ρ]
        _, d = GEMB.calculate_density(fill(250.0, 2), dz, dens, fill(1.0, 2), zeros(2), cfs, mp)
        return (d[2] - ρ) / cfs.dt
    end
    ratio(depth) = rate_at(depth, 550.1) / rate_at(depth, 549.9)

    # Sintering is far slower than Herron-Langway in the top metre, and faster deep.
    @test ratio(1.0) < 1e-3
    @test ratio(5.0) < 0.1
    @test ratio(20.0) > 1.0
    # Monotone in depth, so there is exactly one crossover.
    rs = [ratio(z) for z in (0.5, 1.0, 2.0, 5.0, 10.0, 20.0)]
    @test all(diff(rs) .> 0)
end

@testset "Barnola1991 stage 1 is exactly Herron-Langway" begin
    # Below DENSITY_STAGE_TRANSITION the two schemes share `_hl_c0`, so they must agree
    # bit-for-bit — not merely closely. This identity was verified against the Community
    # Firn Model's Barnola1991 zone 1 before being relied on here.
    n = 6
    dz = fill(0.1, n)
    density = [300.0, 350.0, 400.0, 450.0, 500.0, 540.0]
    temperature = [245.0, 248.0, 250.0, 252.0, 255.0, 258.0]
    grain_radius = fill(0.5, n)
    cfs = GEMB.ClimateForcingStep(; dt=10800.0, precipitation_mean=300.0,
        temperature_air_mean=250.0)

    dz_b, d_b = GEMB.calculate_density(copy(temperature), copy(dz), copy(density),
        copy(grain_radius), zeros(n), cfs,
        GEMB.ModelParameters(density_ice=910.0, densification_method=:Barnola1991))
    dz_h, d_h = GEMB.calculate_density(copy(temperature), copy(dz), copy(density),
        copy(grain_radius), zeros(n), cfs,
        GEMB.ModelParameters(density_ice=910.0, densification_method=:HerronLangway))

    @test d_b == d_h
    @test dz_b == dz_h
    @test all(d_b .> density)      # not vacuous
end

@testset "Barnola1991 densification above the transition" begin
    n = 8
    dz = fill(0.5, n)
    density = [560.0, 620.0, 680.0, 740.0, 790.0, 830.0, 870.0, 900.0]
    temperature = fill(250.0, n)
    grain_radius = fill(1.0, n)
    mp = GEMB.ModelParameters(density_ice=910.0, densification_method=:Barnola1991)
    cfs = GEMB.ClimateForcingStep(; dt=10800.0, precipitation_mean=300.0,
        temperature_air_mean=250.0)

    dz_out, d_out = GEMB.calculate_density(copy(temperature), copy(dz), copy(density),
        copy(grain_radius), zeros(n), cfs, mp)

    # Densifies everywhere, conserves mass, never exceeds ice.
    @test all(d_out .>= density)
    @test all(d_out .<= mp.density_ice)
    @test all(abs.(d_out .* dz_out .- density .* dz) .< 1e-9)

    # Hand-computed rate for cell 2, to pin the whole assembled law (stress integral,
    # Arrhenius, f, σ³) rather than just f.
    load = density[1] * dz[1]                       # cells above cell 2
    σ = (load + 0.5 * density[2] * dz[2]) * GEMB.GRAVITY
    expected = density[2] + density[2] * GEMB.BARNOLA_A0 *
                            exp(-GEMB.BARNOLA_Q / (temperature[2] * GEMB.R_GAS)) *
                            GEMB._barnola_f(density[2], mp.density_ice) *
                            σ^GEMB.BARNOLA_N * cfs.dt
    @test d_out[2] ≈ expected rtol = 1e-13

    # The surface cell compacts under half its own weight, not zero — the same midpoint
    # convention `:Crocus` uses. (Cell 1 here is above the transition, so it takes the
    # pressure-sintering branch.)
    @test d_out[1] > density[1]
end

@testset "Barnola1991 does not stall at the firn-ice transition" begin
    # The reason this scheme exists. `c·(ρᵢ−ρ)` schemes must go to zero as ρ → ρᵢ because
    # the driving term does; Barnola's rate is set by porosity-dependent creep instead, so
    # it stays finite deep in the column where the others have died off.
    n = 300
    dz = fill(0.2, n)
    z = cumsum(dz) .- 0.1
    density = [min(910.0, 350.0 + 560.0 * (1 - exp(-zz / 25))) for zz in z]
    temperature = fill(250.0, n)
    grain_radius = fill(0.5, n)
    cfs = GEMB.ClimateForcingStep(; dt=10800.0, precipitation_mean=300.0,
        temperature_air_mean=250.0)

    rate(method) = begin
        (_, d) = GEMB.calculate_density(copy(temperature), copy(dz), copy(density),
            copy(grain_radius), zeros(n), cfs,
            GEMB.ModelParameters(density_ice=910.0, densification_method=method))
        (d .- density) ./ (cfs.dt / (365.0 * 86400.0))     # kg m-3 yr-1
    end
    barnola = rate(:Barnola1991)
    arthern = rate(:Arthern)

    @test all(barnola .> 0.0)

    # Deepest cells (ρ > 830): Arthern is monotonically dying toward zero, Barnola is not.
    deep = findall(>(GEMB.DENSITY_PORE_CLOSEOFF), density)
    @test !isempty(deep)
    # Arthern's rate at the very bottom is below its rate at close-off; Barnola's is above.
    @test arthern[deep[end]] < arthern[deep[1]]
    @test barnola[deep[end]] > barnola[deep[1]]
    # Barnola stays below Arthern in absolute terms on this profile — the point is not that
    # it is faster, but that it recovers with depth where Arthern monotonically decays. So
    # the ratio closes with depth: 0.40 at close-off against 0.65 at the base.
    ratio = barnola[deep] ./ arthern[deep]
    @test ratio[end] > 1.5 * ratio[1]
    @test ratio[end] < 1.0
    @test argmin(ratio) == 1
end

@testset "Barnola1991 validator and steady-state fallback" begin
    @test initialize_parameters(densification_method=:Barnola1991).densification_method ===
          :Barnola1991
    @test_throws AssertionError initialize_parameters(densification_method=:Barnola)

    # The steady-state marcher has no stress, so `:Barnola1991` falls back to
    # `:HerronLangway` — whose stage 1 is its own below-transition branch, making the
    # fallback exact there and an approximation above.
    p_b = GEMB.DensificationCoeffs(
        GEMB.ModelParameters(densification_method=:Barnola1991), 300.0, 250.0)
    p_h = GEMB.DensificationCoeffs(
        GEMB.ModelParameters(densification_method=:HerronLangway), 300.0, 250.0)
    for ρ in (300.0, 450.0, 540.0, 600.0, 800.0)
        @test GEMB._densification_rate(p_b, ρ, 250.0) ==
              GEMB._densification_rate(p_h, ρ, 250.0)
    end
    # k0/k1 are unused by this scheme.
    @test p_b.k0 == 1.0
    @test p_b.k1 == 1.0
end

@testset "Barnola1991 reads overburden water" begin
    # Stress-driven, so pore water in the overlying cells adds load. Third scheme with
    # this property, after `:ArthernB` and `:Crocus`.
    n = 3
    dz = fill(0.5, n)
    density = fill(700.0, n)
    temperature = fill(250.0, n)
    grain_radius = fill(1.0, n)
    mp = GEMB.ModelParameters(density_ice=910.0, densification_method=:Barnola1991)
    cfs = GEMB.ClimateForcingStep(; dt=10800.0, precipitation_mean=300.0,
        temperature_air_mean=250.0)

    (_, d_dry) = GEMB.calculate_density(copy(temperature), copy(dz), copy(density),
        copy(grain_radius), zeros(n), cfs, mp)
    water = [20.0, 0.0, 0.0]
    (_, d_wet) = GEMB.calculate_density(copy(temperature), copy(dz), copy(density),
        copy(grain_radius), water, cfs, mp)

    # Cell 1's own load is unchanged (water is in cell 1, and only half its own weight
    # counts — but that half does include its water), cells 2-3 carry more.
    @test d_wet[2] > d_dry[2]
    @test d_wet[3] > d_dry[3]
    # σ³ dependence makes this strongly superlinear, not a rounding effect.
    @test (d_wet[3] - density[3]) / (d_dry[3] - density[3]) > 1.05
end

@testset "densification_accumulation and mean_temperature_method" begin
    # Both gates select which of two `ClimateForcingStep` scalars `calculate_density` reads.
    # `:precipitation`/`:arithmetic` must reproduce the pre-refinement behaviour exactly, so
    # the checks below are `==`, not `≈`.
    n = 4
    dz = fill(0.5, n)
    density = [350.0, 450.0, 600.0, 750.0]
    temperature = fill(250.0, n)
    grain_radius = fill(1.0, n)

    pm, am = 400.0, 300.0     # 25% of the precipitation falls as rain
    tam, teff = 250.0, 246.0  # T_eff is always the colder of the two (Arrhenius convexity)

    cfs = GEMB.ClimateForcingStep(; dt=10800.0,
        precipitation_mean=pm, temperature_air_mean=tam,
        accumulation_mean=am, temperature_air_effective=teff)

    _run(mp) = GEMB.calculate_density(copy(temperature), copy(dz), copy(density),
        copy(grain_radius), zeros(n), cfs, mp)[2]

    # Defaults: the refined scalars, which is why Fix 1 and Fix 2 are correctness fixes
    # rather than opt-in preferences.
    mp_default = GEMB.ModelParameters()
    @test mp_default.densification_accumulation === :accumulation
    @test mp_default.mean_temperature_method === :arithmetic

    # Every accumulation-driven scheme must respond to the accumulation gate, and the
    # Arrhenius-scaled subset to the temperature gate. `:Barnola1991` and `:Crocus` read
    # neither scalar, so they are the control: they must be invariant to both.
    # `:ArthernB` and `:CrocusPure` are the controls: both are purely stress-driven and read
    # neither scalar, so they must be invariant to both gates. `:Barnola1991` is *not* a
    # control despite being stress-driven — it falls back to Herron-Langway rates below the
    # stage transition, so it inherits that scheme's accumulation dependence.
    for method in (:HerronLangway, :Arthern, :Barnola1991, :Crocus, :GSFC2020,
                   :Simonsen2013, :Ligtenberg)
        kw = (; densification_method=method)
        d_p = _run(GEMB.ModelParameters(; kw..., densification_accumulation=:precipitation))
        d_a = _run(GEMB.ModelParameters(; kw..., densification_accumulation=:accumulation))
        # Less burial ⇒ less densification, at every cell the scheme touches.
        @test all(d_a .<= d_p)
        @test any(d_a .< d_p)

        # Reading `pm` from the accumulation slot must be *identical* to having been handed
        # that number as `precipitation_mean` — the gate is a selector, nothing more.
        cfs_as_precip = GEMB.ClimateForcingStep(; dt=10800.0,
            precipitation_mean=am, temperature_air_mean=tam)
        d_ref = GEMB.calculate_density(copy(temperature), copy(dz), copy(density),
            copy(grain_radius), zeros(n), cfs_as_precip,
            GEMB.ModelParameters(; kw..., densification_accumulation=:precipitation))[2]
        @test d_a == d_ref
    end

    for method in (:ArthernB, :CrocusPure)
        kw = (; densification_method=method)
        @test _run(GEMB.ModelParameters(; kw..., densification_accumulation=:precipitation)) ==
              _run(GEMB.ModelParameters(; kw..., densification_accumulation=:accumulation))
        @test _run(GEMB.ModelParameters(; kw..., mean_temperature_method=:arithmetic)) ==
              _run(GEMB.ModelParameters(; kw..., mean_temperature_method=:arrhenius))
    end

    # The temperature gate: every Arrhenius-scaled scheme densifies *faster* at the colder
    # effective temperature. For `:Arthern`/`:GSFC2020` this is direct — `exp(+Eg/(R·tam))`
    # rises as `tam` falls. This is the bias Fix 2 corrects: the arithmetic mean
    # under-densifies at any site with a real seasonal cycle.
    for method in (:Arthern, :Crocus, :GSFC2020, :Simonsen2013, :Ligtenberg)
        kw = (; densification_method=method)
        d_arith = _run(GEMB.ModelParameters(; kw..., mean_temperature_method=:arithmetic))
        d_arrh = _run(GEMB.ModelParameters(; kw..., mean_temperature_method=:arrhenius))
        @test all(d_arrh .>= d_arith)
        @test any(d_arrh .> d_arith)

        cfs_as_tam = GEMB.ClimateForcingStep(; dt=10800.0,
            precipitation_mean=am, temperature_air_mean=teff)
        d_ref = GEMB.calculate_density(copy(temperature), copy(dz), copy(density),
            copy(grain_radius), zeros(n), cfs_as_tam,
            GEMB.ModelParameters(; kw..., densification_accumulation=:precipitation))[2]
        @test _run(GEMB.ModelParameters(; kw..., mean_temperature_method=:arrhenius)) == d_ref
    end

    @test_throws AssertionError initialize_parameters(densification_accumulation=:bogus)
    @test_throws AssertionError initialize_parameters(mean_temperature_method=:bogus)
end

@testset "accumulation_mean / temperature_air_effective derivation" begin
    # Both scalars are derived from the record by `initialize_forcing` when not supplied.
    n = 400
    t = collect(DateTime(2000,1,1):Hour(3):DateTime(2000,1,1) + Hour(3*(n-1)))
    # Half the record above freezing, half below, with uniform precipitation, so the
    # snow fraction is exactly 1/2 and the expected value is closed-form.
    T = [i <= n ÷ 2 ? 260.0 : 280.0 for i in 1:n]
    precip = fill(1.0, n)
    cf = GEMB.initialize_forcing(t, T, fill(90000.0, n), precip, fill(3.0, n),
        fill(100.0, n), fill(250.0, n), fill(200.0, n);
        temperature_air_mean=270.0, wind_speed_mean=3.0, precipitation_mean=1000.0,
        temperature_observation_height=2.0, wind_observation_height=10.0)

    # A *fraction* of the supplied `precipitation_mean`, not an independent sum, so a
    # caller-supplied climatological mean keeps its own scale.
    @test cf.accumulation_mean == 500.0
    @test cf.accumulation_mean <= cf.precipitation_mean

    Eg, R = GEMB.GRAIN_GROWTH_EG, GEMB.R_GAS
    @test cf.temperature_air_effective ≈
          Eg / (R * log(GEMB.Statistics.mean(exp.(Eg ./ (R .* T))))) rtol=1e-14
    # Convexity: the Arrhenius mean is bounded by the coldest temperature and the
    # arithmetic mean of the record, and is strictly colder than the latter.
    @test minimum(T) < cf.temperature_air_effective < GEMB.Statistics.mean(T)

    # A rain-free record leaves `accumulation_mean` equal to `precipitation_mean`, so the
    # gate is inert at a dry-snow site rather than silently biased.
    cf_cold = GEMB.initialize_forcing(t, fill(250.0, n), fill(90000.0, n), precip,
        fill(3.0, n), fill(100.0, n), fill(250.0, n), fill(200.0, n);
        temperature_air_mean=250.0, wind_speed_mean=3.0, precipitation_mean=1000.0,
        temperature_observation_height=2.0, wind_observation_height=10.0)
    @test cf_cold.accumulation_mean == cf_cold.precipitation_mean
    # Isothermal record ⇒ the two temperature means coincide exactly.
    @test cf_cold.temperature_air_effective ≈ 250.0 rtol=1e-12

    # Explicitly supplied values are taken verbatim, never recomputed.
    cf_given = GEMB.initialize_forcing(t, T, fill(90000.0, n), precip, fill(3.0, n),
        fill(100.0, n), fill(250.0, n), fill(200.0, n);
        temperature_air_mean=270.0, wind_speed_mean=3.0, precipitation_mean=1000.0,
        temperature_observation_height=2.0, wind_observation_height=10.0,
        accumulation_mean=123.0, temperature_air_effective=234.0)
    @test cf_given.accumulation_mean == 123.0
    @test cf_given.temperature_air_effective == 234.0

    # `forcing_climatology` carries both over from the full record rather than recomputing
    # them from the averaged year, which would smooth away the variability they summarize.
    n2 = 8 * 366 * 2
    t2 = collect(DateTime(2000,1,1):Hour(3):DateTime(2000,1,1) + Hour(3*(n2-1)))
    T2 = 265.0 .+ 15.0 .* sin.(2π .* (0:n2-1) ./ (8*365.25))
    cf2 = GEMB.initialize_forcing(t2, T2, fill(90000.0, n2), fill(0.05, n2),
        fill(3.0, n2), fill(100.0, n2), fill(250.0, n2), fill(200.0, n2);
        temperature_air_mean=265.0, wind_speed_mean=3.0, precipitation_mean=150.0,
        temperature_observation_height=2.0, wind_observation_height=10.0)
    clim = forcing_climatology(cf2)
    @test clim.accumulation_mean == cf2.accumulation_mean
    @test clim.temperature_air_effective == cf2.temperature_air_effective
end

@testset "ClimateForcingStep positional ABI" begin
    # The two new scalars were appended at the *end* of the struct, with a 19-argument
    # positional constructor defaulting them, so the many positional call sites in this
    # suite keep binding their arguments to the fields they were written for.
    cfs = GEMB.ClimateForcingStep(10800.0, 260.0, 90000.0, 0.1, 3.0, 100.0, 250.0, 200.0,
        255.0, 4.0, 300.0, 2.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1)
    @test cfs.temperature_air_mean == 255.0
    @test cfs.precipitation_mean == 300.0
    @test cfs.cloud_fraction == 0.1
    # Defaulted to the scalars they refine, so a 19-argument step makes both gates inert.
    @test cfs.accumulation_mean == cfs.precipitation_mean
    @test cfs.temperature_air_effective == cfs.temperature_air_mean
end
