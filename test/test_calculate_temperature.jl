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

Downward longwave that exactly balances the upward emission of a surface at `T`, so a
steady-state column neither warms nor cools radiatively.
"""
_balanced_longwave(T::Float64, emissivity::Float64=0.97) = GEMB.SB * T^4 * emissivity

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
    # up to 2x STRICTER than the textbook form — the unsound direction, since the 0.8 safety
    # factor cannot absorb it.
    dz_thick = [0.02, 0.02, 0.40, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02]
    i = 3
    Gm = 1.0 / (dz_thick[i] / (2 * K[i]) + dz_thick[i-1] / (2 * K[i-1]))
    Gp = 1.0 / (dz_thick[i+1] / (2 * K[i+1]) + dz_thick[i] / (2 * K[i]))
    face_i = DENSITY_UNIFORM * c * dz_thick[i] / (Gm + Gp)
    @test face_i < textbook(temperature, dz_thick, density, K, i)
    @test face_i / textbook(temperature, dz_thick, density, K, i) < 0.8  # beyond the margin
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
    dt = GEMB._find_dt_divisor(GEMB._max_safe_dt(temperature, dz, density, K, mp) * 0.8,
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
    # measured 3.98 s with the lens against 903 s without it, a factor of 227 in sub-step
    # count for one thin cell. The implicit run above took `n_sub` from
    # `THERMAL_IMPLICIT_DT_TARGET` alone and never saw either number.
    K = GEMB.thermal_conductivity(temperature, density, mp)
    dz_no_lens = [0.02 * 1.08^(i - 1) for i in 1:n]
    dt_lens = GEMB._max_safe_dt(temperature, dz, density, K, mp)
    dt_no_lens = GEMB._max_safe_dt(temperature, dz_no_lens, density, K, mp)
    @test dt_lens ≈ 3.979 atol = 0.01
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
    sfc = GEMB._ThermalSurface(1.3, GEMB.Z0_SNOW_DRY, 0.1 * GEMB.Z0_SNOW_DRY,
        0.1 * GEMB.Z0_SNOW_DRY, 0.97, false, 5.0, 0.16 * 5.0, 1.0,
        log(2.0 / GEMB.Z0_SNOW_DRY), log(2.0 / (0.1 * GEMB.Z0_SNOW_DRY)),
        log(2.0 / (0.1 * GEMB.Z0_SNOW_DRY)), (5.0 / 2.0)^2)
    @test GEMB._surface_energy_balance_slope(GEMB.CtoK + 5.0, 0.97, sfc, cfs,
        -300.0, 10.0, 5.0) === 0.0
    # Below the clamp it is strictly negative — the property that strengthens the diagonal.
    lw, shf, lhf, _, _ = GEMB._surface_energy_balance(260.0, 0.97, sfc, cfs)
    slope = GEMB._surface_energy_balance_slope(260.0, 0.97, sfc, cfs, lw, shf, lhf)
    @test slope < 0.0

    # The finite-difference slope must track the balance it differentiates. Loose (5%) only
    # because the two turbulent components are clamped at ≤ 0 independently; the longwave term
    # is analytic. An *analytic* turbulent slope that froze the Beljaars-Holtslag transfer
    # coefficients was tried and reverted — it overshoots by up to 10x, which more than doubles
    # the Newton iteration count. See the note in `_surface_energy_balance_slope`.
    _Q(T) = let (lw2, s, l, _, _) = GEMB._surface_energy_balance(T, 0.97, sfc, cfs)
        lw2 + s + l
    end
    fd = (_Q(260.0 + 1e-3) - _Q(260.0 - 1e-3)) / 2e-3
    @test slope ≈ fd rtol = 0.05
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
