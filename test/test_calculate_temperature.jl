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
function _temperature_mp(dt_total_ms::Int; emissivity::Float64=0.97)
    GEMB.ModelParameters(
        density_ice=917.0,
        emissivity=emissivity,
        emissivity_method=:uniform,
        emissivity_grain_radius_large=emissivity,
        emissivity_grain_radius_threshold=10.0,
        surface_roughness_effective_ratio=0.1,
        thermal_conductivity_method=:Sturm,
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
