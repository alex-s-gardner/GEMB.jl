# Tests for calculate_temperature
# Translated from MATLAB test_calculate_temperature.m

# Tolerances and shared setup values for this file. Physical constants come from
# GEMB (`GEMB.SB`) rather than being restated, so a change there cannot leave the
# balanced-longwave forcing below silently out of balance.
const T_STEADY_DRIFT_MAX = 1.0   # allowed surface drift from turbulent fluxes [K]
const T_BOTTOM_ATOL = 0.0        # Dirichlet BC: bottom cell must be bit-unchanged [K]

# Shared column geometry / state for every case below.
const N_CELLS = 10
const DZ_UNIFORM = 0.1           # [m]
const DENSITY_UNIFORM = 400.0    # [kg m-3]
const GRAIN_RADIUS_UNIFORM = 0.5 # [mm]

"""
    _temperature_mp(dt_total_ms; emissivity=0.97)

`ModelParameters` for the temperature tests. `dt_total_ms` is the forcing timestep in
units of 1e-4 s, matching how `dt_divisors` is built from `fast_divisors` (integer
divisors of the timestep, rescaled by 1e-4 to get seconds).
"""
function _temperature_mp(dt_total_ms::Int; emissivity::Float64=0.97,
    thermal_solver::GEMB.AbstractThermalSolver=GEMB.ExplicitThermal(),
    heat_capacity_method::Symbol=:constant)
    GEMB.ModelParameters(
        density_ice=917.0,
        emissivity=emissivity,
        emissivity_method=:uniform,
        emissivity_grain_radius_large=emissivity,
        emissivity_grain_radius_threshold=10.0,
        surface_roughness_effective_ratio=0.1,
        thermal_conductivity_method=:Sturm,
        heat_capacity_method=heat_capacity_method,
        thermal_solver=thermal_solver,
        dt_divisors=Float64.(GEMB.fast_divisors(dt_total_ms)) ./ 10000
    )
end

"""
    _balanced_longwave(T)

Downward longwave that exactly balances the upward emission of a gray surface at `T`, so a
steady-state column neither warms nor cools radiatively: the surface absorbs `ε·LW↓` and emits
`εσT⁴`, so the balance `LW↓ = σT⁴` holds for every emissivity.
"""
_balanced_longwave(T::Float64) = GEMB.SB * T^4

@testset "Steady state" begin
    # If T_surf = temperature_air and LW balanced, T should remain roughly constant
    n = N_CELLS
    t_vec = fill(260.0, n)
    dz = fill(DZ_UNIFORM, n)
    density = fill(DENSITY_UNIFORM, n)
    water_surface = 0.0
    grain_radius = fill(GRAIN_RADIUS_UNIFORM, n)
    shortwave_flux = zeros(n)

    cfs = GEMB.ClimateForcingStep(
        3600.0, 260.0, 100000.0, 0.0, 5.0, 0.0,
        _balanced_longwave(260.0),  # longwave_downward balanced
        100.0, 260.0, 5.0, 0.0, 2.0, 2.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    )

    mp = _temperature_mp(36000000)

    t_out, _, _, _, _, _ = GEMB.calculate_temperature(
        t_vec, dz, density, water_surface, grain_radius,
        shortwave_flux, cfs, mp, false)

    @test abs(t_out[1] - 260.0) < T_STEADY_DRIFT_MAX  # Allow some drift from turbulent fluxes
end

@testset "Solar heating" begin
    n = N_CELLS
    t_vec = fill(260.0, n)
    dz = fill(DZ_UNIFORM, n)
    density = fill(DENSITY_UNIFORM, n)
    water_surface = 0.0
    grain_radius = fill(GRAIN_RADIUS_UNIFORM, n)
    shortwave_flux = zeros(n)
    shortwave_flux[1] = 200.0  # 200 W/m2 absorbed in top layer

    cfs = GEMB.ClimateForcingStep(
        3600.0, 260.0, 100000.0, 0.0, 5.0, 0.0,
        _balanced_longwave(260.0),  # balance LW
        100.0, 260.0, 5.0, 0.0, 2.0, 2.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    )

    mp = _temperature_mp(36000000)

    t_out, _, _, _, _, _ = GEMB.calculate_temperature(
        t_vec, dz, density, water_surface, grain_radius,
        shortwave_flux, cfs, mp, false)

    @test t_out[1] > 260.0  # Top layer should warm
end

@testset "Thermal diffusion" begin
    n = N_CELLS
    t_vec = fill(250.0, n)
    t_vec[1] = 273.0  # Hot surface

    dz = fill(DZ_UNIFORM, n)
    density = fill(DENSITY_UNIFORM, n)
    water_surface = 0.0
    grain_radius = fill(GRAIN_RADIUS_UNIFORM, n)
    shortwave_flux = zeros(n)

    cfs = GEMB.ClimateForcingStep(
        10800.0, 273.0, 100000.0, 0.0, 0.1, 0.0,  # warm air, low wind to minimize turbulent cooling
        0.0,  # no longwave
        100.0, 260.0, 5.0, 0.0, 2.0, 2.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    )

    mp = _temperature_mp(108000000; emissivity=0.0)  # disable radiative cooling

    initial_gradient = t_vec[1] - t_vec[2]  # 23.0 (save before call mutates t_vec)

    t_out, _, _, _, _, _ = GEMB.calculate_temperature(
        t_vec, dz, density, water_surface, grain_radius,
        shortwave_flux, cfs, mp, false)

    # Diffusion should reduce the temperature gradient between layers 1 and 2
    final_gradient = t_out[1] - t_out[2]
    @test final_gradient < initial_gradient  # Gradient should decrease via diffusion
    @test t_out[2] > 250.0  # Layer 2 should warm from diffusion
end

@testset "Bottom boundary condition" begin
    n = N_CELLS
    t_vec = fill(260.0, n)
    t_vec[end] = 240.0  # Distinct bottom temp

    dz = fill(DZ_UNIFORM, n)
    density = fill(DENSITY_UNIFORM, n)
    water_surface = 0.0
    grain_radius = fill(GRAIN_RADIUS_UNIFORM, n)
    shortwave_flux = zeros(n)

    cfs = GEMB.ClimateForcingStep(
        3600.0, 260.0, 100000.0, 0.0, 5.0, 0.0, 0.0,
        100.0, 260.0, 5.0, 0.0, 2.0, 2.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    )

    mp = _temperature_mp(36000000)

    t_out, _, _, _, _, _ = GEMB.calculate_temperature(
        t_vec, dz, density, water_surface, grain_radius,
        shortwave_flux, cfs, mp, false)

    @test t_out[end] ≈ 240.0 atol = T_BOTTOM_ATOL  # Bottom T fixed (Dirichlet BC)
end

@testset "No NaN with large timestep" begin
    n = N_CELLS
    t_vec = fill(260.0, n)
    dz = fill(DZ_UNIFORM, n)
    density = fill(DENSITY_UNIFORM, n)
    water_surface = 0.0
    grain_radius = fill(GRAIN_RADIUS_UNIFORM, n)
    shortwave_flux = zeros(n)

    cfs = GEMB.ClimateForcingStep(
        86400.0, 260.0, 100000.0, 0.0, 5.0, 0.0, 0.0,
        100.0, 260.0, 5.0, 0.0, 2.0, 2.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    )

    mp = _temperature_mp(864000000)

    t_out, _, _, _, _, _ = GEMB.calculate_temperature(
        t_vec, dz, density, water_surface, grain_radius,
        shortwave_flux, cfs, mp, false)

    @test !any(isnan.(t_out))  # Should not explode
end

@testset "Stability limit (_max_safe_dt)" begin
    # On a uniform grid with uniform K the face-based limit reduces analytically to the
    # textbook form it replaced: G_i = 1/(dz/2K + dz/2K) = K/dz, so for an interior cell
    # rho*c*dz/(G_i + G_{i-1}) = rho*c*dz/(2K/dz) = 0.5*rho*c*dz^2/K.
    n = 10
    dz = fill(DZ_UNIFORM, n)
    density = fill(DENSITY_UNIFORM, n)
    temperature = fill(260.0, n)
    mp = _temperature_mp(864000000)
    K = GEMB.thermal_conductivity(temperature, density, mp)
    @test all(≈(K[1]), K)   # uniform state => uniform K, precondition for the algebra below

    c = GEMB.heat_capacity(mp, temperature[1])
    analytic = 0.5 * DENSITY_UNIFORM * c * DZ_UNIFORM^2 / K[1]

    # Binds on an interior cell: cell 1 has only one face so its limit is 2x, and cell n is
    # excluded as the Dirichlet reservoir.
    @test GEMB._max_safe_dt(temperature, dz, density, K, mp) ≈ analytic rtol = 1e-12

    # The surface cell has a single face, so its limit is exactly twice the interior value —
    # verified by shrinking every other cell out of contention.
    dz_surf = copy(dz); dz_surf[2:end] .= 10.0
    G1 = 1.0 / (dz_surf[2] / (2 * K[2]) + dz_surf[1] / (2 * K[1]))
    @test GEMB._max_safe_dt(temperature, dz_surf, density, K, mp) ≈
          DENSITY_UNIFORM * c * dz_surf[1] / G1 rtol = 1e-12

    # The Dirichlet bottom cell's thermal mass imposes no constraint. Driving it to near zero
    # would dominate any per-cell minimum that included cell n; the limit must not move.
    # `density` is perturbed rather than `dz` because `dz[n]` legitimately enters face n-1's
    # resistance, whereas the bottom cell's mass enters nothing but its own excluded limit.
    # (`K` is supplied explicitly here, so it does not track the perturbed density.)
    density_light_bottom = copy(density); density_light_bottom[n] = 1e-6
    @test GEMB._max_safe_dt(temperature, dz, density_light_bottom, K, mp) ≈ analytic rtol = 1e-12

    # On a graded grid the face-based limit departs from the textbook form in BOTH directions,
    # which is why the latter was not conservative. A thin cell between thick neighbours has
    # G_i + G_{i-1} < 2K_i/dz_i, so the correct limit is looser...
    textbook(t, d, r, Kv, i) = 0.5 * r[i] * GEMB.heat_capacity(mp, t[i]) * d[i]^2 / Kv[i]
    dz_graded = [0.30, 0.30, 0.02, 0.30, 0.30, 0.30, 0.30, 0.30, 0.30, 0.30]
    @test GEMB._max_safe_dt(temperature, dz_graded, density, K, mp) >
          minimum(textbook(temperature, dz_graded, density, K, i) for i in 1:n)

    # ...while a thick cell between thin ones approaches 4K_i/dz_i, making the correct limit
    # up to 2x STRICTER than the textbook form — the unsound direction, since the safety
    # factor cannot absorb it.
    dz_thick = [0.02, 0.02, 0.40, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02]
    i = 3
    Gm = 1.0 / (dz_thick[i] / (2 * K[i]) + dz_thick[i-1] / (2 * K[i-1]))
    Gp = 1.0 / (dz_thick[i+1] / (2 * K[i+1]) + dz_thick[i] / (2 * K[i]))
    face_i = DENSITY_UNIFORM * c * dz_thick[i] / (Gm + Gp)
    @test face_i < textbook(temperature, dz_thick, density, K, i)
    # Beyond the margin: even the loosest admissible safety factor would not have covered it.
    @test face_i / textbook(temperature, dz_thick, density, K, i) <
          GEMB.THERMAL_EXPLICIT_SAFETY_FACTOR
end

@testset "Explicit safety factor is exposed and covers the surface term" begin
    # What the factor is FOR. `_max_safe_dt` bounds diffusion only; cell 1 additionally carries
    # the surface energy balance, whose slope `Λ <= 0` enters its own-temperature coefficient
    # exactly as a face conductance does. So the coefficient the scheme actually realizes is
    # `1 - dt*(G_1 + |Λ|)/(rho*c*dz)`, and the diffusive limit alone does not keep it
    # non-negative — the factor is the margin that does.
    @test 0.0 < GEMB.THERMAL_EXPLICIT_SAFETY_FACTOR <= 1.0
    @test GEMB.ModelParameters().thermal_explicit_safety_factor ===
          GEMB.THERMAL_EXPLICIT_SAFETY_FACTOR

    # Bounded in (0, 1]. Above 1 would scale up a limit that is already an upper bound, i.e.
    # run the explicit scheme knowingly unstable; 0 would collapse `dt` to the smallest divisor.
    @test_throws AssertionError initialize_parameters(thermal_explicit_safety_factor=1.5)
    @test_throws AssertionError initialize_parameters(thermal_explicit_safety_factor=0.0)
    @test_throws AssertionError initialize_parameters(thermal_explicit_safety_factor=-0.5)
    @test initialize_parameters(thermal_explicit_safety_factor=1.0) isa GEMB.ModelParameters

    # The factor reaches the sub-step: a small value must force more sub-steps than the default
    # on the same column. Checked through `calculate_temperature`'s own arithmetic rather than
    # by inspection, so a rename that failed to thread it through would fail here.
    n = 40
    dz = [0.02 * 1.08^(i - 1) for i in 1:n]
    density = [350.0 + 500.0 * (i - 1) / (n - 1) for i in 1:n]
    temperature = [250.0 + 20.0 * (i - 1) / (n - 1) for i in 1:n]
    mp_default = _temperature_mp(108000000)
    K = GEMB.thermal_conductivity(temperature, density, mp_default)
    limit = GEMB._max_safe_dt(temperature, dz, density, K, mp_default)
    dt_default = GEMB._find_dt_divisor(
        limit * mp_default.thermal_explicit_safety_factor, mp_default.dt_divisors)
    mp_tight = GEMB.ModelParameters(
        thermal_conductivity_method=:Sturm, thermal_explicit_safety_factor=0.05,
        dt_divisors=mp_default.dt_divisors)
    dt_tight = GEMB._find_dt_divisor(
        limit * mp_tight.thermal_explicit_safety_factor, mp_tight.dt_divisors)
    @test dt_tight < dt_default

    # The factor changes the sub-step count, so it necessarily changes the answer by the
    # scheme's own truncation error — under strong forcing that is ~2 K in the near-surface
    # cells over one 3-hour step, which is the explicit scheme's accuracy, not a defect. The
    # claim worth pinning is therefore made where the exact answer is known: an isothermal
    # column with air at its temperature, balanced longwave, and no shortwave is an exact
    # steady state, which every admissible factor must reproduce to the same drift bound.
    grain_radius = fill(GRAIN_RADIUS_UNIFORM, n)
    dz_u = fill(DZ_UNIFORM, N_CELLS)
    density_u = fill(DENSITY_UNIFORM, N_CELLS)
    t_steady = fill(260.0, N_CELLS)
    cfs_steady = GEMB.ClimateForcingStep(
        3600.0, 260.0, 100000.0, 0.0, 5.0, 0.0, _balanced_longwave(260.0),
        100.0, 260.0, 5.0, 0.0, 2.0, 2.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    )
    for sf in (0.05, 0.5, GEMB.THERMAL_EXPLICIT_SAFETY_FACTOR, 1.0)
        mp_sf = GEMB.ModelParameters(
            thermal_conductivity_method=:Sturm, thermal_explicit_safety_factor=sf,
            dt_divisors=Float64.(GEMB.fast_divisors(36000000)) ./ 10000)
        t_sf, _, _, _, _, _ = GEMB.calculate_temperature(
            copy(t_steady), dz_u, density_u, 0.0, fill(GRAIN_RADIUS_UNIFORM, N_CELLS),
            zeros(N_CELLS), cfs_steady, mp_sf, true)
        @test all(abs.(t_sf .- 260.0) .< T_STEADY_DRIFT_MAX)
        @test t_sf[end] === 260.0
    end

    # Under strong forcing the two differ only by truncation: neither rings below the initial
    # minimum, and `verbose=true` closes each budget internally. No upper bound is asserted —
    # the melting-point clamp is on the *flux*, so a cell may legitimately exceed `CtoK` and be
    # handed to `calculate_melt`.
    cfs = GEMB.ClimateForcingStep(
        10800.0, 265.0, 100000.0, 0.0, 5.0, 200.0, 300.0,
        200.0, 260.0, 5.0, 0.0, 2.0, 2.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    )
    shortwave_flux = [200.0 * exp(-3.0 * (i - 1)) for i in 1:n]
    t_default, _, _, _, _, _ = GEMB.calculate_temperature(
        copy(temperature), dz, density, 0.0, grain_radius, shortwave_flux, cfs, mp_default, true)
    t_tight, _, _, _, _, _ = GEMB.calculate_temperature(
        copy(temperature), dz, density, 0.0, grain_radius, shortwave_flux, cfs, mp_tight, true)
    for t_sf in (t_default, t_tight)
        @test !any(isnan, t_sf)
        @test all(isfinite, t_sf)
        @test minimum(t_sf) >= minimum(temperature) - T_STEADY_DRIFT_MAX
        @test t_sf[end] === temperature[end]
    end

    # `ImplicitThermal` must ignore the field entirely: it never consults `_max_safe_dt`, so
    # varying the factor cannot move a bit of its answer.
    mp_imp_a = GEMB.ModelParameters(thermal_solver=GEMB.ImplicitThermal(),
        thermal_conductivity_method=:Sturm, dt_divisors=Float64[])
    mp_imp_b = GEMB.ModelParameters(thermal_solver=GEMB.ImplicitThermal(),
        thermal_conductivity_method=:Sturm, dt_divisors=Float64[],
        thermal_explicit_safety_factor=0.05)
    t_imp_a, _, _, _, _, _ = GEMB.calculate_temperature(
        copy(temperature), dz, density, 0.0, grain_radius, shortwave_flux, cfs, mp_imp_a, true)
    t_imp_b, _, _, _, _, _ = GEMB.calculate_temperature(
        copy(temperature), dz, density, 0.0, grain_radius, shortwave_flux, cfs, mp_imp_b, true)
    @test t_imp_a == t_imp_b
end

@testset "Stability limit keeps the solve positive and stable" begin
    # The condition `_max_safe_dt(:face)` enforces is that every cell's own-temperature
    # coefficient stays non-negative. Verify it directly on a graded, thermally varied
    # column, and confirm the resulting solve neither explodes nor overshoots.
    n = 40
    dz = [0.02 * 1.08^(i - 1) for i in 1:n]          # graded: 0.02 m to ~0.37 m
    density = [350.0 + 500.0 * (i - 1) / (n - 1) for i in 1:n]
    temperature = [250.0 + 20.0 * (i - 1) / (n - 1) for i in 1:n]
    grain_radius = fill(GRAIN_RADIUS_UNIFORM, n)
    shortwave_flux = zeros(n)

    mp = _temperature_mp(864000000)

    K = GEMB.thermal_conductivity(temperature, density, mp)
    dt = GEMB._find_dt_divisor(
        GEMB._max_safe_dt(temperature, dz, density, K, mp) *
        mp.thermal_explicit_safety_factor,
        mp.dt_divisors)

    G_prev = 0.0
    for i in 1:n-1
        G = 1.0 / (dz[i+1] / (2 * K[i+1]) + dz[i] / (2 * K[i]))
        coef = 1.0 - dt * (G + G_prev) / (density[i] * GEMB.heat_capacity(mp, temperature[i]) * dz[i])
        @test coef >= 0.0
        G_prev = G
    end

    # Air at the surface temperature with balanced longwave and no shortwave, so the surface
    # is forced neutrally and the only active process is diffusion down the column.
    T_air = temperature[1]
    cfs = GEMB.ClimateForcingStep(
        86400.0, T_air, 100000.0, 0.0, 5.0, 0.0,
        _balanced_longwave(T_air),
        100.0, T_air, 5.0, 0.0, 2.0, 2.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    )
    t_in = copy(temperature)
    t_out, _, _, _, _, _ = GEMB.calculate_temperature(
        copy(temperature), dz, density, 0.0, grain_radius, shortwave_flux, cfs, mp, true)

    @test !any(isnan, t_out)
    @test all(t_out .<= GEMB.CtoK)
    # Diffusion cannot create a new extremum. This is the observable consequence of the
    # positivity condition checked above: a violated limit shows up as oscillatory
    # over/undershoot beyond the initial range, well before it becomes an outright blow-up.
    @test minimum(t_out) >= minimum(t_in) - T_STEADY_DRIFT_MAX
    @test maximum(t_out) <= maximum(t_in) + T_STEADY_DRIFT_MAX
    @test t_out[end] == t_in[end]   # Dirichlet bottom cell untouched
end

# ---------------------------------------------------------------------------------------
# ImplicitThermal
#
# The two testsets above are explicit-specific: they exercise `_max_safe_dt` and the Von
# Neumann positivity condition, neither of which exists on the implicit path. Everything
# below is the implicit solver's own contract.
# ---------------------------------------------------------------------------------------

"""
    _implicit_column(n)

Graded, thermally varied column shared by the implicit testsets: `dz` from 0.02 m to ~0.37 m
at n = 40, density increasing with depth, temperature increasing with depth.
"""
function _implicit_column(n::Int=40)
    dz = [0.02 * 1.08^(i - 1) for i in 1:n]
    density = [350.0 + 500.0 * (i - 1) / max(1, n - 1) for i in 1:n]
    temperature = [250.0 + 20.0 * (i - 1) / max(1, n - 1) for i in 1:n]
    grain_radius = fill(GRAIN_RADIUS_UNIFORM, n)
    return dz, density, temperature, grain_radius
end

@testset "Implicit solver: Dirichlet cell is bit-unchanged" begin
    # The invariant `AbstractThermalSolver` binds every implementation to, and the one the
    # `!=` check inside the solver and `gemb_core`'s drift check both rely on. Bit-exact,
    # not approximate: cell m is never an unknown of the tridiagonal system.
    for n in (2, 3, 40, 265), hcm in (:constant, :CuffeyPaterson)
        dz, density, temperature, grain_radius = _implicit_column(n)
        mp = _temperature_mp(108000000; thermal_solver=GEMB.ImplicitThermal(),
            heat_capacity_method=hcm)
        cfs = GEMB.ClimateForcingStep(
            10800.0, 265.0, 100000.0, 0.0, 5.0, 200.0, 300.0,
            200.0, 260.0, 5.0, 0.0, 2.0, 2.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0
        )
        t_in = copy(temperature)
        shortwave_flux = [200.0 * exp(-3.0 * (i - 1)) for i in 1:n]
        t_out, _, _, _, _, _ = GEMB.calculate_temperature(
            copy(temperature), dz, density, 0.0, grain_radius, shortwave_flux, cfs, mp, true)
        @test t_out[end] === t_in[end]
    end
end

@testset "Radiative equilibrium holds for one day (both solvers)" begin
    # Isothermal column at T_eq with LW↓ = σT_eq⁴, no shortwave, air at T_eq and saturated
    # over ice, so the sensible and latent fluxes vanish at T_s = T_eq. A gray surface absorbs
    # ε·LW↓ = εσT_eq⁴ and emits the same, so the column is an exact steady state.
    T_eq = 260.0
    n = N_CELLS
    for solver in (GEMB.ExplicitThermal(), GEMB.ImplicitThermal()), ε in (0.97, 0.8)
        mp = _temperature_mp(864000000; emissivity=ε, thermal_solver=solver)
        cfs = GEMB.ClimateForcingStep(
            86400.0, T_eq, 100000.0, 0.0, 5.0, 0.0, GEMB.SB * T_eq^4,
            GEMB.saturation_vapor_pressure_ice(T_eq), T_eq, 5.0, 0.0, 2.0, 2.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        t_out, lw_up, shf, lhf, _, _ = GEMB.calculate_temperature(
            fill(T_eq, n), fill(DZ_UNIFORM, n), fill(DENSITY_UNIFORM, n), 0.0,
            fill(GRAIN_RADIUS_UNIFORM, n), zeros(n), cfs, mp, true)
        @test maximum(abs.(t_out .- T_eq)) < 1e-3
        # Upward longwave is emitted plus reflected, which at equilibrium equals LW↓.
        @test lw_up ≈ GEMB.SB * T_eq^4 rtol = 1e-9
    end
end

@testset "Energy budget against fluxes recomputed outside the solver" begin
    # One sub-step, so the surface fluxes are those of a single surface state that the test can
    # recompute itself: on the implicit path (backward Euler) the balance at the new T₁; on the
    # explicit path the balance at the initial T₁ linearized in the change of T₁,
    # Q(T₁⁰) + Q'(T₁⁰)·ΔT₁. The expected enthalpy change is the textbook surface balance
    # (1 − α)SW + ε·LW↓ − εσT_s⁴ + H + LE plus the basal conduction, with H and LE from
    # `turbulent_heat_flux`, Q' by central difference, and the basal flux from the
    # harmonic-mean face conductance.
    n = N_CELLS
    dz = fill(DZ_UNIFORM, n)
    density = fill(DENSITY_UNIFORM, n)
    t_in = [255.0 + (i - 1) for i in 1:n]
    shortwave_flux = [150.0 * exp(-2.0 * (i - 1)) for i in 1:n]
    ε = 0.9
    for (solver, dt) in ((GEMB.ExplicitThermal(), 10.0), (GEMB.ImplicitThermal(), 900.0))
        mp = GEMB.ModelParameters(emissivity=ε, emissivity_method=:uniform,
            surface_roughness_effective_ratio=0.1, thermal_conductivity_method=:Sturm,
            thermal_solver=solver, dt_divisors=[dt])
        cfs = GEMB.ClimateForcingStep(
            dt, 265.0, 90000.0, 0.0, 6.0, 300.0, 280.0,
            150.0, 265.0, 6.0, 0.0, 2.0, 10.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        t_out, = GEMB.calculate_temperature(copy(t_in), dz, density, 0.0,
            fill(GRAIN_RADIUS_UNIFORM, n), shortwave_flux, cfs, mp, true)

        T_s = solver isa GEMB.ExplicitThermal ? t_in[1] : t_out[1]
        T_prof = solver isa GEMB.ExplicitThermal ? t_in : t_out
        z0 = GEMB.surface_roughness(mp, density[1], 0.0)
        ρa = GEMB.air_density(cfs.pressure_air, cfs.temperature_air)
        Q(T) = ε * cfs.longwave_downward - ε * GEMB.SB * T^4 +
               sum(GEMB.turbulent_heat_flux(T, ρa, z0, 0.1 * z0, 0.1 * z0, cfs)[1:2])
        Q_applied = Q(T_s)
        if solver isa GEMB.ExplicitThermal
            Q_applied += (Q(T_s + 1e-3) - Q(T_s - 1e-3)) / 2e-3 * (t_out[1] - t_in[1])
        end
        K = GEMB.thermal_conductivity(t_in, density, mp)
        G = (T_prof[n] - T_prof[n-1]) / (dz[n] / (2 * K[n]) + dz[n-1] / (2 * K[n-1]))
        expected = (sum(shortwave_flux) + Q_applied + G) * dt
        M = density .* dz
        ΔH = sum(M[i] * (GEMB.specific_enthalpy(mp, t_out[i]) -
                         GEMB.specific_enthalpy(mp, t_in[i])) for i in 1:n)
        @test ΔH ≈ expected atol = GEMB.energy_tolerance(GEMB.column_enthalpy(mp, M, t_in))
    end
end

@testset "Emissivity melt switch and reflected longwave per sub-step" begin
    # Two explicit sub-steps of 1 s on a surface held at the melting point. The first sub-step
    # uses `mp.emissivity`; the switch then selects `emissivity_grain_radius_large` for the
    # second. Each sub-step's upward longwave is εσ·273.15⁴ + (1 − ε)·LW↓ with its own ε.
    n = N_CELLS
    ε_small, ε_large, lw_down = 0.90, 0.98, 310.0
    mp = GEMB.ModelParameters(emissivity=ε_small, emissivity_grain_radius_large=ε_large,
        emissivity_method=:grain_radius_w_threshold, emissivity_grain_radius_threshold=10.0,
        thermal_conductivity_method=:Sturm, dt_divisors=[1.0])
    cfs = GEMB.ClimateForcingStep(
        2.0, 275.0, 100000.0, 0.0, 5.0, 800.0, lw_down,
        700.0, 275.0, 5.0, 0.0, 2.0, 2.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    shortwave_flux = zeros(n); shortwave_flux[1] = 800.0
    _, lw_up, = GEMB.calculate_temperature(fill(GEMB.CtoK, n), fill(DZ_UNIFORM, n),
        fill(DENSITY_UNIFORM, n), 0.0, fill(GRAIN_RADIUS_UNIFORM, n), shortwave_flux,
        cfs, mp, true)
    σT4 = GEMB.SB * GEMB.CtoK^4
    expected = 0.5 * ((ε_small * σT4 + (1 - ε_small) * lw_down) +
                      (ε_large * σT4 + (1 - ε_large) * lw_down))
    @test lw_up ≈ expected rtol = 1e-12
end

@testset "Implicit solver: energy budget closes" begin
    # `verbose=true` asserts conservation internally, but "didn't throw" is a weak
    # statement. Close the budget here explicitly against the same tolerance the solver
    # uses, so a widened tolerance inside cannot hide a leak.
    for hcm in (:constant, :CuffeyPaterson)
        n = 40
        dz, density, temperature, grain_radius = _implicit_column(n)
        mp = _temperature_mp(108000000; thermal_solver=GEMB.ImplicitThermal(),
            heat_capacity_method=hcm)
        shortwave_flux = [200.0 * exp(-3.0 * (i - 1)) for i in 1:n]
        cfs = GEMB.ClimateForcingStep(
            10800.0, 265.0, 100000.0, 0.0, 5.0, 200.0, 300.0,
            200.0, 260.0, 5.0, 0.0, 2.0, 2.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0
        )

        M = density .* dz
        E_initial = GEMB.column_enthalpy(mp, M, temperature)
        t_out, lw_up, shf, lhf, basal, _ = GEMB.calculate_temperature(
            copy(temperature), dz, density, 0.0, grain_radius, shortwave_flux, cfs, mp, true)
        E_final = GEMB.column_enthalpy(mp, M, t_out)

        # Supplied over the whole forcing step. `lw_up` is returned as a positive magnitude
        # of an outgoing flux, hence the minus sign.
        E_supplied = (sum(shortwave_flux) + cfs.longwave_downward - lw_up + shf + lhf +
                      basal) * cfs.dt
        @test E_final - E_initial ≈ E_supplied atol = GEMB.energy_tolerance(E_initial)
    end
end

@testset "Implicit solver: no new extremum (maximum principle)" begin
    # The system matrix is an irreducibly diagonally dominant M-matrix, so A^-1 >= 0 and
    # the solve cannot manufacture an over/undershoot. Forced neutrally (air at the surface
    # temperature, balanced longwave, no shortwave) so diffusion is the only active process
    # and the bound is tight — 1e-9, against the explicit path's 1.0 K allowance.
    n = 40
    dz, density, temperature, grain_radius = _implicit_column(n)
    mp = _temperature_mp(864000000; thermal_solver=GEMB.ImplicitThermal())
    T_air = temperature[1]
    cfs = GEMB.ClimateForcingStep(
        86400.0, T_air, 100000.0, 0.0, 5.0, 0.0,
        _balanced_longwave(T_air),
        100.0, T_air, 5.0, 0.0, 2.0, 2.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    )
    t_in = copy(temperature)
    t_out, _, _, _, _, _ = GEMB.calculate_temperature(
        copy(temperature), dz, density, 0.0, grain_radius, zeros(n), cfs, mp, true)

    @test !any(isnan, t_out)
    @test minimum(t_out) >= minimum(t_in) - 1e-9
    @test maximum(t_out) <= maximum(t_in) + 1e-9
    @test all(t_out .<= GEMB.CtoK)
    @test t_out[end] === t_in[end]
end

@testset "Implicit solver: m == 2, the one-unknown system" begin
    # `calculate_temperature` admits m == 2, which leaves exactly one unknown — the `n == 1`
    # branch, where cell 1 is simultaneously the surface row and the row that carries the
    # Dirichlet reservoir term. Neither the elimination loop nor the back-substitution runs.
    dz = [0.05, 0.20]
    density = [400.0, 700.0]
    temperature = [255.0, 265.0]
    grain_radius = fill(GRAIN_RADIUS_UNIFORM, 2)
    mp = _temperature_mp(36000000; thermal_solver=GEMB.ImplicitThermal())
    cfs = GEMB.ClimateForcingStep(
        3600.0, 265.0, 100000.0, 0.0, 5.0, 100.0, 250.0,
        200.0, 260.0, 5.0, 0.0, 2.0, 2.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    )
    t_out, _, _, _, basal, _ = GEMB.calculate_temperature(
        copy(temperature), dz, density, 0.0, grain_radius, [100.0, 0.0], cfs, mp, true)
    @test !any(isnan, t_out)
    @test t_out[end] === 265.0
    @test t_out[1] > 255.0        # warmed by the flux and the warmer reservoir below
    @test basal > 0.0             # heat flows up from the warmer bottom cell
end

@testset "Implicit solver: unaffected by a thin lens" begin
    # The structural claim. A 1e-4 m cell collapses the explicit stability limit, while the
    # implicit path neither consults the limit nor changes cost. Also pins that
    # `dt_divisors = Float64[]` is accepted here, which would `BoundsError` under
    # `ExplicitThermal`.
    n = 40
    dz, density, temperature, grain_radius = _implicit_column(n)
    dz[n÷2] = 1e-4
    cfs = GEMB.ClimateForcingStep(
        10800.0, 265.0, 100000.0, 0.0, 5.0, 200.0, 300.0,
        200.0, 260.0, 5.0, 0.0, 2.0, 2.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    )
    shortwave_flux = [200.0 * exp(-3.0 * (i - 1)) for i in 1:n]

    # No `dt_divisors` at all, which the explicit path cannot run without.
    mp = GEMB.ModelParameters(thermal_solver=GEMB.ImplicitThermal(),
        thermal_conductivity_method=:Sturm, dt_divisors=Float64[])
    t_out, _, _, _, _, _ = GEMB.calculate_temperature(
        copy(temperature), dz, density, 0.0, grain_radius, shortwave_flux, cfs, mp, true)
    @test !any(isnan, t_out)
    @test t_out[end] === temperature[end]

    # The explicit limit really has collapsed on this column, so the contrast is real:
    # measured 3.79 s with the lens against 832 s without it, a factor of 219 in sub-step
    # count for one thin cell. The implicit run above took `n_sub` from
    # `THERMAL_IMPLICIT_DT_TARGET` alone and never saw either number.
    K = GEMB.thermal_conductivity(temperature, density, mp)
    dz_no_lens = [0.02 * 1.08^(i - 1) for i in 1:n]
    dt_lens = GEMB._max_safe_dt(temperature, dz, density, K, mp)
    dt_no_lens = GEMB._max_safe_dt(temperature, dz_no_lens, density, K, mp)
    @test dt_lens ≈ 3.791 atol = 0.01
    @test dt_no_lens / dt_lens > 200.0
end

@testset "Implicit solver: melting-point clamp and melt switch" begin
    # The surface balance is non-smooth at CtoK — the `min(CtoK, T1)` clamp flattens every
    # flux and every derivative — and `emissivity_method = :grain_radius_w_threshold` adds a
    # second discontinuity there. Held outside the Newton loop; this pins that the solve
    # still lands, conserves, and never pushes a cell above the melting point.
    n = 40
    dz, density, temperature, grain_radius = _implicit_column(n)
    fill!(temperature, GEMB.CtoK - 1e-6)      # sitting on the clamp
    mp = GEMB.ModelParameters(
        thermal_solver=GEMB.ImplicitThermal(),
        emissivity_method=:grain_radius_w_threshold,
        emissivity_grain_radius_threshold=10.0,
        thermal_conductivity_method=:Sturm,
        dt_divisors=Float64[],
    )
    cfs = GEMB.ClimateForcingStep(
        10800.0, 275.0, 100000.0, 0.0, 5.0, 600.0, 320.0,
        400.0, 274.0, 5.0, 0.0, 2.0, 2.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    )
    t_in = copy(temperature)
    t_out, _, _, _, _, _ = GEMB.calculate_temperature(
        copy(temperature), dz, density, 0.5, grain_radius,
        [600.0 * exp(-3.0 * (i - 1)) for i in 1:n], cfs, mp, true)
    @test !any(isnan, t_out)
    @test t_out[end] === t_in[end]
    # Strong melt-season forcing drives the surface up, but the clamp is on the *flux*, not
    # on the temperature: the cell may exceed CtoK, which is what `calculate_melt` consumes.
    @test t_out[1] > t_in[1]

    # The slope is exactly zero once clamped — the flux no longer sees T1 at all.
    sfc = GEMB._ThermalSurface(1.3, 0.97, false, 0.16 * 5.0, GEMB.GRAVITY * 2.0 / 5.0^2,
        GEMB._stability_table(2.0, 2.0, GEMB.Z0_SNOW_DRY, 0.1 * GEMB.Z0_SNOW_DRY))
    @test GEMB._surface_energy_balance_and_slope(GEMB.CtoK + 5.0, 0.97, sfc, cfs)[6] === 0.0
    # Below the clamp it is strictly negative — the property that strengthens the diagonal.
    slope = GEMB._surface_energy_balance_and_slope(260.0, 0.97, sfc, cfs)[6]
    @test slope < 0.0

    # The analytic slope is the derivative of the balance it accompanies, transfer-coefficient
    # changes included.
    _Q(T) = let (lw2, s, l, _, _) = GEMB._surface_energy_balance(T, 0.97, sfc, cfs)
        lw2 + s + l
    end
    fd = (_Q(260.0 + 1e-3) - _Q(260.0 - 1e-3)) / 2e-3
    @test slope ≈ fd rtol = 1e-5
end

@testset "Implicit solver: steady state is a fixed point" begin
    # An isothermal column with air at its temperature, balanced longwave, and no shortwave
    # is an exact steady state of the continuous problem, so backward Euler — which is exact
    # on steady states regardless of step size — must reproduce it to within the turbulent
    # flux residual alone. Tighter than the explicit path's equivalent for the same reason.
    n = N_CELLS
    dz = fill(DZ_UNIFORM, n)
    density = fill(DENSITY_UNIFORM, n)
    temperature = fill(260.0, n)
    grain_radius = fill(GRAIN_RADIUS_UNIFORM, n)
    mp = _temperature_mp(36000000; thermal_solver=GEMB.ImplicitThermal())
    cfs = GEMB.ClimateForcingStep(
        3600.0, 260.0, 100000.0, 0.0, 5.0, 0.0, _balanced_longwave(260.0),
        100.0, 260.0, 5.0, 0.0, 2.0, 2.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    )
    t_out, _, _, _, _, _ = GEMB.calculate_temperature(
        copy(temperature), dz, density, 0.0, grain_radius, zeros(n), cfs, mp, false)
    @test abs(t_out[1] - 260.0) < T_STEADY_DRIFT_MAX
    @test all(abs.(t_out .- 260.0) .< T_STEADY_DRIFT_MAX)
end

@testset "Implicit solver: chord capacity conserves under :CuffeyPaterson" begin
    # The reason `chord_heat_capacity` exists. With a temperature-dependent c_p the lagged
    # form `c_p(T_old)` injects (b/2)*dT^2 per cell per step; the chord form does not. Drive
    # a large surface swing (the case where dT is biggest) and require the column enthalpy
    # to close far tighter than that spurious term would allow.
    n = 40
    dz, density, temperature, grain_radius = _implicit_column(n)
    mp = _temperature_mp(108000000; thermal_solver=GEMB.ImplicitThermal(),
        heat_capacity_method=:CuffeyPaterson)
    cfs = GEMB.ClimateForcingStep(
        10800.0, 275.0, 100000.0, 0.0, 10.0, 800.0, 330.0,
        400.0, 274.0, 10.0, 0.0, 2.0, 2.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    )
    shortwave_flux = [800.0 * exp(-3.0 * (i - 1)) for i in 1:n]
    M = density .* dz
    E_initial = GEMB.column_enthalpy(mp, M, temperature)
    t_out, lw_up, shf, lhf, basal, _ = GEMB.calculate_temperature(
        copy(temperature), dz, density, 0.0, grain_radius, shortwave_flux, cfs, mp, true)
    E_final = GEMB.column_enthalpy(mp, M, t_out)
    E_supplied = (sum(shortwave_flux) + cfs.longwave_downward - lw_up + shf + lhf +
                  basal) * cfs.dt
    @test E_final - E_initial ≈ E_supplied atol = GEMB.energy_tolerance(E_initial)

    # The surface really did swing far enough for the lagged form to have mattered: the
    # spurious term (b/2)*dT^2 per kg would exceed the tolerance many times over.
    dT_surface = abs(t_out[1] - temperature[1])
    @test dT_surface > 5.0
    @test M[1] * (GEMB.HEAT_CAPACITY_CUFFEY_B / 2) * dT_surface^2 >
          100 * GEMB.energy_tolerance(E_initial)
end

@testset "Implicit solver: Newton damping breaks limit cycles" begin
    # Measured on a year of 3-hourly forcing, 39% of sub-step solves reached the iteration cap
    # in a *limit cycle* — 2- and 3-cycles across the latent-heat and emissivity switches at
    # 273.15 K, not divergence. The iteration is damped for that: the step weight halves
    # whenever a step fails to shrink. These testsets pin the two properties that makes it
    # safe to have — converging solves are untouched, and cycling solves land on their best
    # iterate rather than on whichever phase the loop stopped at.
    # See `GEMB.THERMAL_IMPLICIT_DAMPING_FLOOR`.

    # --- the damping floor is a guard, not a tuned value ---------------------------------
    # It must be strictly positive (zero would freeze the iterate at the initial guess) and
    # small enough to be out of reach of the iteration cap, so it can never bind on a real
    # solve: reaching it takes 10 halvings against a 20-iteration budget also spent iterating.
    @test GEMB.THERMAL_IMPLICIT_DAMPING_FLOOR > 0.0
    @test GEMB.THERMAL_IMPLICIT_DAMPING_FLOOR < 1.0
    @test 2.0^(-GEMB.THERMAL_IMPLICIT_MAX_ITERATIONS) <
          GEMB.THERMAL_IMPLICIT_DAMPING_FLOOR

    # --- the damped iteration still solves the problem it is given ----------------------
    # Forced hard across the melting point, which is where the cycles were observed. The
    # invariants are the solver's own: the Dirichlet cell is bit-unchanged, nothing is NaN,
    # and `verbose=true` closes the energy budget internally.
    n = 40
    for (T_fill, T_air, sw) in ((GEMB.CtoK - 1e-9, 276.0, 700.0),   # sitting on the switch
                                (GEMB.CtoK, 274.0, 600.0),         # exactly on it
                                (260.0, 285.0, 800.0))             # driven hard through it
        dz, density, temperature, grain_radius = _implicit_column(n)
        fill!(temperature, T_fill)
        mp = GEMB.ModelParameters(
            thermal_solver=GEMB.ImplicitThermal(),
            emissivity_method=:grain_radius_w_threshold,
            emissivity_grain_radius_threshold=10.0,
            thermal_conductivity_method=:Sturm,
            dt_divisors=Float64[],
        )
        cfs = GEMB.ClimateForcingStep(
            10800.0, T_air, 100000.0, 0.0, 5.0, sw, 330.0,
            500.0, T_air - 1.0, 5.0, 0.0, 2.0, 2.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0
        )
        t_in = copy(temperature)
        t_out, _, _, _, _, _ = GEMB.calculate_temperature(
            copy(temperature), dz, density, 0.5, grain_radius,
            [sw * exp(-3.0 * (i - 1)) for i in 1:n], cfs, mp, true)
        @test !any(isnan, t_out)
        @test all(isfinite, t_out)
        @test t_out[end] === t_in[end]
    end

    # --- converging solves are bit-identical to the undamped iteration -------------------
    # The licence for damping unconditionally. On the quadratic path the residual falls every
    # iteration, so the weight never leaves 1.0 and `T1 + 1.0*(T1_next - T1)` is exactly
    # `T1_next`. Checked as a *fixed point*: a smooth, well-conditioned, sub-melting-point
    # column solved twice must agree to the last bit, and the steady state below is the case
    # where any spurious under-relaxation would show up as a drift away from the exact answer.
    n = N_CELLS
    dz = fill(DZ_UNIFORM, n)
    density = fill(DENSITY_UNIFORM, n)
    temperature = fill(260.0, n)
    grain_radius = fill(GRAIN_RADIUS_UNIFORM, n)
    mp = _temperature_mp(36000000; thermal_solver=GEMB.ImplicitThermal())
    cfs = GEMB.ClimateForcingStep(
        3600.0, 260.0, 100000.0, 0.0, 5.0, 0.0, _balanced_longwave(260.0),
        100.0, 260.0, 5.0, 0.0, 2.0, 2.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    )
    t1, _, _, _, _, _ = GEMB.calculate_temperature(
        copy(temperature), dz, density, 0.0, grain_radius, zeros(n), cfs, mp, true)
    t2, _, _, _, _, _ = GEMB.calculate_temperature(
        copy(temperature), dz, density, 0.0, grain_radius, zeros(n), cfs, mp, true)
    @test t1 == t2                                  # deterministic, bit-for-bit
    # The same exact steady state the "steady state is a fixed point" testset above uses, held
    # to the same drift bound: backward Euler reproduces a steady state whatever the step size,
    # so any residual under-relaxation left by the damping would show up here as a drift.
    @test all(abs.(t1 .- 260.0) .< T_STEADY_DRIFT_MAX)
end

# ---------------------------------------------------------------------------------------
# ThermalWorkspace: reusable scratch, held outside ModelParameters
#
# Both schemes draw their working arrays from a `ThermalWorkspace` and reuse them across
# calls, so the buffers outlive a single solve. That trades allocation for the risk that a
# later call reads a value the previous one left behind. Keeping the buffers out of the
# solver — and so out of `mp` — is what lets one `mp` be read by many threads at once; these
# pin both halves of that arrangement.
# ---------------------------------------------------------------------------------------

@testset "ThermalWorkspace: reuse cannot carry state between calls" begin
    # The whole correctness argument for the shared buffers is that every entry is written
    # before it is read, so a solve cannot see the previous solve's leftovers. Testing it
    # directly: run column A, then B, then A again through *one* workspace. If any buffer
    # entry were read before being written, A-after-B would differ from A-first.
    for solver in (GEMB.ExplicitThermal(), GEMB.ImplicitThermal())
        mp = _temperature_mp(108000000; thermal_solver=solver)
        cfs = GEMB.ClimateForcingStep(
            10800.0, 265.0, 100000.0, 0.0, 5.0, 200.0, 300.0,
            200.0, 260.0, 5.0, 0.0, 2.0, 2.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0
        )
        ws = ThermalWorkspace()
        function solve(n)
            dz, density, temperature, grain_radius = _implicit_column(n)
            shortwave_flux = [200.0 * exp(-3.0 * (i - 1)) for i in 1:n]
            return GEMB.calculate_temperature(copy(temperature), dz, density, 0.0,
                grain_radius, shortwave_flux, cfs, mp, true, ws)[1]
        end

        a_first = solve(40)
        solve(25)                      # a shorter column, leaving longer buffers behind
        a_again = solve(40)
        @test a_first == a_again        # bit-for-bit, not approximate

        # Same check in the other order: a long column after a short one must not read the
        # tail the short call never touched.
        b_first = solve(60)
        solve(12)
        @test b_first == solve(60)

        # And a solve drawing on the reused workspace must match one on a fresh workspace,
        # which is what the defaulted argument gives. Same physics either way.
        dz, density, temperature, grain_radius = _implicit_column(40)
        fresh = GEMB.calculate_temperature(copy(temperature), dz, density, 0.0,
            grain_radius, [200.0 * exp(-3.0 * (i - 1)) for i in 1:40], cfs, mp, true)[1]
        @test fresh == a_first
    end
end

@testset "ThermalWorkspace: instances do not share buffers" begin
    # One workspace per concurrent run is the documented rule, which is only workable if
    # separate constructions get separate buffers — pin it, since the advice is worthless
    # otherwise.
    a, b = ThermalWorkspace(), ThermalWorkspace()
    @test a.explicit.H !== b.explicit.H
    @test a.implicit.rhs !== b.implicit.rhs
    # The two schemes' buffer sets are distinct too, so switching `thermal_solver` on a
    # reused workspace cannot alias.
    @test a.explicit.H !== a.implicit.H
    # A fresh workspace starts empty and is grown on first use, so construction is cheap and
    # no column length is baked in at construction time.
    @test isempty(ThermalWorkspace().explicit.H)
    @test isempty(ThermalWorkspace().implicit.rhs)
end

@testset "ThermalWorkspace: one ModelParameters is shared across threads" begin
    # The design property this whole arrangement exists for: `mp` is a pure description, so
    # many threads can step columns concurrently against one `mp`, each with its own
    # workspace. Run the same column serially and on all available threads and demand
    # bit-identical results — a buffer reachable through `mp` would race and diverge here
    # (and, on `nthreads() == 1`, this still checks that per-workspace runs agree).
    for solver in (GEMB.ExplicitThermal(), GEMB.ImplicitThermal())
        mp = _temperature_mp(108000000; thermal_solver=solver)   # one shared, immutable mp
        cfs = GEMB.ClimateForcingStep(
            10800.0, 265.0, 100000.0, 0.0, 5.0, 200.0, 300.0,
            200.0, 260.0, 5.0, 0.0, 2.0, 2.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0
        )
        # Differing column lengths so the threads resize their buffers to different sizes.
        lengths = repeat([40, 25, 60, 12], 8)
        run(n) = begin
            dz, density, temperature, grain_radius = _implicit_column(n)
            GEMB.calculate_temperature(copy(temperature), dz, density, 0.0, grain_radius,
                [200.0 * exp(-3.0 * (i - 1)) for i in 1:n], cfs, mp, true,
                ThermalWorkspace())[1]
        end
        expected = map(run, lengths)

        threaded = Vector{Vector{Float64}}(undef, length(lengths))
        Threads.@threads for i in eachindex(lengths)
            threaded[i] = run(lengths[i])
        end
        @test threaded == expected
    end
end

@testset "ThermalWorkspace: buffers grow but never shrink" begin
    # `_resize_workspace!` only ever grows, so a long column followed by a short one leaves
    # the buffers long. Callers index by the column length, never by `length(buffer)`, which
    # is what makes that safe; this documents the behaviour so a future `resize!` down (which
    # would churn) is a deliberate change rather than an accident.
    ws = GEMB._ExplicitWorkspace()
    GEMB._resize_workspace!(ws, 40)
    @test length(ws.H) == 40
    GEMB._resize_workspace!(ws, 10)
    @test length(ws.H) == 40        # unchanged: shrinking is not a thing it does
    GEMB._resize_workspace!(ws, 100)
    @test length(ws.H) == 100
    @test all(length(b) == 100 for b in GEMB._workspace_buffers(ws))
end

@testset "Wet surface cell holds the melting point while its water freezes" begin
    # A wet surface cell at 273.15 K loses 100 W m-2: absorbed ε·LW↓ is 100 W m-2 short of
    # εσT⁴, and the air is at the melting point and saturated over water, so the turbulent
    # fluxes vanish there. Until its water has frozen the cell stays at 273.15 K inside the
    # solve, so the upward longwave is that of a 273.15 K surface throughout; the water that
    # froze, 100·dt/LF, comes back as cold content in the returned matrix temperature, and
    # calculate_melt refreezes it and counts it.
    n = N_CELLS
    ε = 0.97
    W = 1.0
    T0 = GEMB.CtoK
    lw_down = GEMB.SB * T0^4 - 100.0 / ε
    for solver in (GEMB.ExplicitThermal(), GEMB.ImplicitThermal()), dt in (1800.0, 7200.0)
        # 60 s explicit sub-steps, so the lagged surface flux resolves the cooling.
        mp = GEMB.ModelParameters(emissivity=ε, emissivity_method=:uniform,
            surface_roughness_effective_ratio=0.1, thermal_conductivity_method=:Sturm,
            thermal_solver=solver, dt_divisors=[60.0])
        cfs = GEMB.ClimateForcingStep(
            dt, T0, 100000.0, 0.0, 5.0, 0.0, lw_down,
            GEMB.saturation_vapor_pressure_water(T0), T0, 5.0, 0.0, 2.0, 2.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        water = [W; zeros(n - 1)]
        dz = fill(DZ_UNIFORM, n)
        density = fill(DENSITY_UNIFORM, n)
        t_out, lw_up, shf, lhf, _, _, _ = GEMB.calculate_temperature(fill(T0, n), dz,
            density, W, fill(GRAIN_RADIUS_UNIFORM, n), zeros(n), cfs, mp, true; water)
        M1 = DENSITY_UNIFORM * DZ_UNIFORM
        frozen = 100.0 * dt / GEMB.LF
        if frozen < W
            @test lw_up ≈ ε * GEMB.SB * T0^4 + (1 - ε) * lw_down rtol = 1e-12
            @test shf == 0.0 && lhf == 0.0
            @test all(t_out[2:end] .== T0)
            @test M1 * (GEMB.specific_enthalpy(mp, T0) - GEMB.specific_enthalpy(mp, t_out[1])) ≈
                  100.0 * dt rtol = 1e-9
            out = GEMB.calculate_melt(t_out, dz, density, copy(water), fill(0.5, n),
                zeros(n), fill(0.5, n), zeros(n), mp, true)
            @test out[1][1] ≈ T0 atol = 1e-6
            @test out[4][1] ≈ W - frozen rtol = 1e-9
            @test out[12] ≈ frozen rtol = 1e-9
        else
            # 7200 s freezes all of it after 3345 s; the surface then cools below 273.15 K.
            @test lw_up < ε * GEMB.SB * T0^4 + (1 - ε) * lw_down
            @test t_out[1] < T0
        end
    end

    # Without `water` the same column is dry and cools at once.
    mp = GEMB.ModelParameters(emissivity=ε, emissivity_method=:uniform,
        surface_roughness_effective_ratio=0.1, thermal_conductivity_method=:Sturm,
        dt_divisors=[60.0])
    cfs = GEMB.ClimateForcingStep(
        1800.0, T0, 100000.0, 0.0, 5.0, 0.0, lw_down,
        GEMB.saturation_vapor_pressure_water(T0), T0, 5.0, 0.0, 2.0, 2.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    _, lw_dry, = GEMB.calculate_temperature(fill(T0, n), fill(DZ_UNIFORM, n),
        fill(DENSITY_UNIFORM, n), 0.0, fill(GRAIN_RADIUS_UNIFORM, n), zeros(n), cfs, mp, true)
    @test lw_dry < ε * GEMB.SB * T0^4 + (1 - ε) * lw_down

    @test_throws "must equal water_surface" GEMB.calculate_temperature(fill(T0, n),
        fill(DZ_UNIFORM, n), fill(DENSITY_UNIFORM, n), 0.0, fill(GRAIN_RADIUS_UNIFORM, n),
        zeros(n), cfs, mp, false; water=[1.0; zeros(n - 1)])
end

@testset "Wet-cell temperature map" begin
    for hcm in (:CuffeyPaterson, :constant)
        mp = GEMB.ModelParameters(heat_capacity_method=hcm)
        h(T) = GEMB.specific_enthalpy(mp, T)
        M, W = 40.0, 2.0
        H_lo = (M + W) * h(GEMB.CtoK)
        H_hi = H_lo + W * GEMB.LF
        @test GEMB._wet_cell_temperature(mp, H_hi, M, W) ≈ GEMB.CtoK rtol = 1e-15
        @test GEMB._wet_cell_temperature(mp, 0.5 * (H_lo + H_hi), M, W) == GEMB.CtoK
        @test GEMB._wet_cell_temperature(mp, H_lo, M, W) == GEMB.CtoK
        # One kelvin above the plateau warms the matrix alone; one below cools matrix and
        # frozen water together.
        @test GEMB._wet_cell_temperature(mp, H_hi + M * (h(GEMB.CtoK + 1.0) - h(GEMB.CtoK)),
            M, W) ≈ GEMB.CtoK + 1.0 rtol = 1e-12
        @test GEMB._wet_cell_temperature(mp,
            H_lo - (M + W) * (h(GEMB.CtoK) - h(GEMB.CtoK - 1.0)), M, W) ≈
              GEMB.CtoK - 1.0 rtol = 1e-12
    end
end

@testset "Explicit surface flux converges in the sub-step" begin
    # A thin 300 kg m-3 surface cell under 15 m s-1 of air 15 K colder than it. The mean
    # turbulent flux over the hour at the default safety factor agrees within 2% with the
    # flux at a 16x smaller factor.
    n = 40
    dz = [0.02 * 1.08^(i - 1) for i in 1:n]
    density = fill(300.0, n)
    T0 = [255.0 + 10.0 * (i - 1) / (n - 1) for i in 1:n]
    cfs = GEMB.ClimateForcingStep(3600.0, 240.0, 80000.0, 0.0, 15.0, 0.0, 180.0, 30.0,
        250.0, 10.0, 0.0, 2.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    flux = map((GEMB.THERMAL_EXPLICIT_SAFETY_FACTOR, 0.05)) do sf
        mp = GEMB.ModelParameters(thermal_conductivity_method=:Sturm,
            thermal_explicit_safety_factor=sf,
            dt_divisors=Float64.(GEMB.fast_divisors(36000000)) ./ 10000)
        _, _, shf, lhf, = GEMB.calculate_temperature(copy(T0), dz, density, 0.0,
            fill(0.3, n), zeros(n), cfs, mp, true)
        shf + lhf
    end
    @test flux[1] ≈ flux[2] rtol = 0.02

    # The linearized surface flux uses the analytic slope of the surface balance.
    mp = GEMB.ModelParameters(thermal_conductivity_method=:Sturm,
        dt_divisors=Float64.(GEMB.fast_divisors(36000000)) ./ 10000)
    sfc = GEMB._ThermalSurface(GEMB.air_density(80000.0, 240.0), mp.emissivity, false,
        0.16 * 15.0, GEMB.GRAVITY * 10.0 / 15.0^2,
        GEMB._stability_table(10.0, 2.0, GEMB.Z0_SNOW_DRY, 0.1 * GEMB.Z0_SNOW_DRY))
    Q(T) = sum(GEMB._surface_energy_balance(T, mp.emissivity, sfc, cfs)[1:3])
    Λ = GEMB._surface_energy_balance_and_slope(T0[1], mp.emissivity, sfc, cfs)[6]
    @test Λ ≈ (Q(T0[1] + 1e-3) - Q(T0[1] - 1e-3)) / 2e-3 rtol = 1e-5
    @test Λ < 0.0
    parts = GEMB._surface_fluxes_and_slopes(T0[1], mp.emissivity, sfc, cfs)
    @test parts[6] + parts[7] + parts[8] == Λ
    @test all(iszero, GEMB._surface_fluxes_and_slopes(GEMB.CtoK, mp.emissivity, sfc, cfs)[6:8])
end

@testset "An all-low-density column stays stable under strong cooling" begin
    # At 100 kg m-3 the conductivity is so small that the diffusive limit exceeds the forcing
    # step, and only the linearized surface flux keeps the explicit scheme from oscillating. The result
    # stays physical and close to the implicit solution.
    n = 20
    dz = fill(0.1, n)
    density = fill(100.0, n)
    cfs = GEMB.ClimateForcingStep(3600.0, 230.0, 80000.0, 0.0, 12.0, 0.0, 150.0, 10.0,
        250.0, 10.0, 0.0, 2.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    T1 = map((GEMB.ExplicitThermal(), GEMB.ImplicitThermal())) do solver
        mp = GEMB.ModelParameters(thermal_conductivity_method=:Sturm, thermal_solver=solver,
            dt_divisors=Float64.(GEMB.fast_divisors(36000000)) ./ 10000)
        t, = GEMB.calculate_temperature(fill(260.0, n), dz, density, 0.0, fill(0.3, n),
            zeros(n), cfs, mp, true)
        @test all(220.0 .< t .<= 260.0)
        t[1]
    end
    @test T1[1] ≈ T1[2] atol = 2.0
end

# A forcing with no turbulent or radiative exchange at surface temperature `T_s`: air at
# `T_s`, saturated over ice, calm, no longwave. Used with `emissivity = 0`.
_quiet_cfs(dt, T_s) = GEMB.ClimateForcingStep(dt, T_s, 100000.0, 0.0, 0.0, 0.0, 0.0,
    GEMB.saturation_vapor_pressure_ice(T_s), T_s, 0.0, 0.0, 2.0, 2.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

@testset "Semi-infinite solid under constant surface flux" begin
    # F = 100 W m-2 absorbed in a 1 cm surface cell of solid ice 2 m deep, constant K and c.
    # After 6 h the profile matches
    #   T(z,t) = T₀ + (2F/k)√(κt/π)·exp(−z²/4κt) − (Fz/k)·erfc(z/2√(κt))
    # to 2% of the surface rise.
    erfc_(x) = begin   # Abramowitz & Stegun 7.1.26, |error| < 1.5e-7
        t = 1 / (1 + 0.3275911 * x)
        t * (0.254829592 + t * (-0.284496736 + t * (1.421413741 +
            t * (-1.453152027 + t * 1.061405429)))) * exp(-x^2)
    end
    n = 200
    dz = fill(0.01, n)
    ρ = 917.0
    T0 = 250.0
    F = 100.0
    t_end = 6 * 3600.0
    for solver in (GEMB.ExplicitThermal(), GEMB.ImplicitThermal())
        # The closed form needs a constant heat capacity.
        mp = GEMB.ModelParameters(emissivity=0.0, emissivity_method=:uniform,
            thermal_conductivity_method=:Sturm, thermal_solver=solver,
            heat_capacity_method=:constant,
            dt_divisors=Float64.(GEMB.fast_divisors(36000000)) ./ 10000)
        k = GEMB.thermal_conductivity(fill(T0, n), fill(ρ, n), mp)[1]
        κ = k / (ρ * mp.heat_capacity_ice)
        sw = zeros(n); sw[1] = F
        T = fill(T0, n)
        ws = ThermalWorkspace()
        for _ in 1:6
            T, = GEMB.calculate_temperature(T, dz, fill(ρ, n), 0.0, fill(1.0, n), sw,
                _quiet_cfs(3600.0, T0), mp, true, ws)
        end
        z = cumsum(dz) .- dz ./ 2
        analytic = [T0 + (2F / k) * sqrt(κ * t_end / π) * exp(-zi^2 / (4κ * t_end)) -
                    (F * zi / k) * erfc_(zi / (2 * sqrt(κ * t_end))) for zi in z]
        rise = analytic[1] - T0
        @test maximum(abs.(T[1:60] .- analytic[1:60])) < 0.02 * rise
        @test T[end] == T0
    end
end

@testset "Discrete linear steady state" begin
    # Uniform K and dz, gradient F/K downward, F absorbed in cell 1, no surface exchange:
    # every face carries F, so one step leaves T unchanged and the basal flux is −F.
    n = 30
    dz = fill(0.05, n)
    ρ = 600.0
    F = 20.0
    for solver in (GEMB.ExplicitThermal(), GEMB.ImplicitThermal())
        mp = GEMB.ModelParameters(emissivity=0.0, emissivity_method=:uniform,
            thermal_conductivity_method=:Sturm, thermal_solver=solver,
            dt_divisors=Float64.(GEMB.fast_divisors(36000000)) ./ 10000)
        K = GEMB.thermal_conductivity(fill(250.0, n), fill(ρ, n), mp)[1]
        z = cumsum(dz) .- dz ./ 2
        T_in = [250.0 - (F / K) * (z[i] - z[1]) for i in 1:n]
        sw = zeros(n); sw[1] = F
        T, _, _, _, basal, = GEMB.calculate_temperature(copy(T_in), dz, fill(ρ, n), 0.0,
            fill(1.0, n), sw, _quiet_cfs(3600.0, T_in[1]), mp, true)
        @test maximum(abs.(T .- T_in)) < 1e-8
        @test basal ≈ -F atol = 1e-8
    end
end

@testset "Melting-point clamp on the surface flux" begin
    # A surface cell above 273.15 K sees the surface fluxes of a 273.15 K surface, and its
    # excess enthalpy is exactly what it received: M₁·(h(T₁) − h(273.15)) = Q·dt.
    n = 2
    dz = [0.05, 50.0]
    ρ = 400.0
    T_in = [GEMB.CtoK, GEMB.CtoK]
    lw_down = 300.0
    cfs = GEMB.ClimateForcingStep(3600.0, 278.0, 90000.0, 0.0, 6.0, 0.0, lw_down, 500.0,
        270.0, 6.0, 0.0, 2.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    sw = [400.0, 0.0]
    for solver in (GEMB.ExplicitThermal(), GEMB.ImplicitThermal())
        mp = GEMB.ModelParameters(thermal_conductivity_method=:Sturm, thermal_solver=solver,
            dt_divisors=Float64.(GEMB.fast_divisors(36000000)) ./ 10000)
        T, lw_up, shf, lhf, basal, = GEMB.calculate_temperature(copy(T_in), dz, fill(ρ, n),
            0.0, fill(0.5, n), sw, cfs, mp, true)
        @test T[1] > GEMB.CtoK
        ε = mp.emissivity
        @test lw_up ≈ ε * GEMB.SB * GEMB.CtoK^4 + (1 - ε) * lw_down rtol = 1e-12
        z0 = GEMB.surface_roughness(mp, ρ, 0.0)
        shf_0, lhf_0, _ = GEMB.turbulent_heat_flux(GEMB.CtoK,
            GEMB.air_density(cfs.pressure_air, cfs.temperature_air), z0,
            mp.surface_roughness_effective_ratio * z0, mp.surface_roughness_effective_ratio * z0,
            cfs)
        @test shf ≈ shf_0 rtol = 1e-12
        @test lhf ≈ lhf_0 rtol = 1e-12
        Q = sw[1] + lw_down - lw_up + shf + lhf + basal
        @test ρ * dz[1] * (GEMB.specific_enthalpy(mp, T[1]) - GEMB.specific_enthalpy(mp, GEMB.CtoK)) ≈
              Q * 3600.0 rtol = 1e-9
    end
end
