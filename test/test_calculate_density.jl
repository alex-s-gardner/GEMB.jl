# Tests for calculate_density - matches MATLAB test_calculate_density.m

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
