# Finding 17: liquid water in a cell was not a latent buffer inside the thermal solve, so a wet
# surface cooled below 0 °C within a step and radiated as a colder surface.
#
# A wet surface cell (1 kg m-2 of water) at 0 °C loses 100 W m-2 for 30 minutes, which
# freezes 0.54 kg m-2: absorbed ε·LW↓ is 100 W m-2 short of εσT⁴ and the air is at 0 °C and
# saturated over water, so the turbulent fluxes vanish at 0 °C.
#
# Expected: on main the step-mean upward longwave falls below the 0 °C value; on
# fix/physics-audit it equals the 0 °C value to the printed precision in both solvers.
#
# Run (see README.md): julia --project=<GEMB checkout> scripts/audit/07_wet_cell.jl
# with <GEMB checkout> this branch's checkout or a worktree of main at 0cda823.
include(joinpath(@__DIR__, "common.jl"))

T0 = GEMB.CtoK
ε = 0.97
lw_down = GEMB.SB * T0^4 - 100.0 / ε
n = 10
for solver in (ExplicitThermal(), ImplicitThermal())
    mp = GEMB.ModelParameters(emissivity=ε, emissivity_method=:uniform,
        thermal_conductivity_method=:Sturm, thermal_solver=solver, dt_divisors=[60.0])
    cfs = GEMB.ClimateForcingStep(1800.0, T0, 100000.0, 0.0, 5.0, 0.0, lw_down,
        GEMB.saturation_vapor_pressure_water(T0), T0, 5.0, 0.0, 2.0, 2.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    out = thermal(fill(T0, n), fill(0.1, n), fill(400.0, n), [1.0; zeros(n - 1)],
        fill(0.5, n), zeros(n), cfs, mp)
    @printf("%-16s LW↑ = %.4f W m-2 (0 °C surface: %.4f)\n", string(solver), out[2],
        ε * GEMB.SB * T0^4 + (1 - ε) * lw_down)
end
