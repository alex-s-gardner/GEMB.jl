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
