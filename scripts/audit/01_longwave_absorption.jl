# Known issue (excluded from the numbered findings): the surface absorbed all of LW↓ instead
# of ε·LW↓.
#
# An isothermal 1 m snow column at T_eq = 260 K under LW↓ = σT_eq⁴, no shortwave, air at T_eq
# and saturated over ice is an exact radiative steady state for a gray surface of any
# emissivity. One day is integrated with each solver.
#
# Expected: main warms the surface by ~0.45 K and returns LW↑ ≈ 252.9 against LW↓ = 259.1
# W m-2 (it absorbs the (1−ε)·LW↓ ≈ 7.8 W m-2 a gray surface reflects); fix/physics-audit stays
# within 1e-4 K and returns LW↑ = LW↓.
#
# Run (see README.md): julia --project=<GEMB checkout> scripts/audit/01_longwave_absorption.jl
# with <GEMB checkout> this branch's checkout or a worktree of main at 0cda823.
include(joinpath(@__DIR__, "common.jl"))

T_eq = 260.0
n = 10
for solver in (ExplicitThermal(), ImplicitThermal())
    mp = GEMB.ModelParameters(emissivity=0.97, emissivity_method=:uniform,
        thermal_conductivity_method=:Sturm, thermal_solver=solver,
        dt_divisors=Float64.(GEMB.fast_divisors(864000000)) ./ 10000)
    cfs = step(dt=86400.0, temperature_air=T_eq, pressure_air=100000.0, wind_speed=5.0,
        longwave_downward=GEMB.SB * T_eq^4, vapor_pressure=GEMB.saturation_vapor_pressure_ice(T_eq),
        temperature_air_mean=T_eq, wind_speed_mean=5.0, temperature_observation_height=2.0,
        wind_observation_height=2.0)
    out = thermal(fill(T_eq, n), fill(0.1, n), fill(400.0, n), zeros(n), fill(0.5, n),
        zeros(n), cfs, mp)
    @printf("%-16s surface ΔT after 1 day = %+.4f K   LW↑ = %.3f W m-2 (LW↓ = %.3f)\n",
        string(solver), out[1][1] - T_eq, out[2], GEMB.SB * T_eq^4)
end
