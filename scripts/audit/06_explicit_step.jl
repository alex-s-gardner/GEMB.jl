# Explicit-solver latent item: sub-steps comparable to the surface relaxation time made the
# lagged surface flux inaccurate in windy cold conditions, and an all-low-density column
# diverged.
#
# Case A: a thin 300 kg m-3 surface cell under 15 m s-1 of air 15 K colder; hour-mean
# turbulent flux at safety factor 0.8 against 0.05.
# Case B: a 2 m column of 100 kg m-3 snow at 260 K under 12 m s-1 of 230 K air, one hour,
# explicit against implicit.
#
# Expected: main differs by ~30% in case A and its explicit case-B surface falls to ~138 K
# against an implicit ~232.5 K; fix/physics-audit differs by under 2% in case A and its
# explicit case B is within 2 K of the implicit solution.
#
# Run (see README.md): julia --project=<GEMB checkout> scripts/audit/06_explicit_step.jl
# with <GEMB checkout> this branch's checkout or a worktree of main at 0cda823.
include(joinpath(@__DIR__, "common.jl"))

divs = Float64.(GEMB.fast_divisors(36000000)) ./ 10000
n = 40
dz = [0.02 * 1.08^(i - 1) for i in 1:n]
T0 = [255.0 + 10.0 * (i - 1) / (n - 1) for i in 1:n]
cfs = GEMB.ClimateForcingStep(3600.0, 240.0, 80000.0, 0.0, 15.0, 0.0, 180.0, 30.0,
    250.0, 10.0, 0.0, 2.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
flux = map((0.8, 0.05)) do sf
    mp = GEMB.ModelParameters(thermal_conductivity_method=:Sturm,
        thermal_explicit_safety_factor=sf, dt_divisors=divs)
    out = thermal(copy(T0), dz, fill(300.0, n), zeros(n), fill(0.3, n), zeros(n), cfs, mp)
    out[3] + out[4]
end
@printf("case A: H+LE at safety 0.8 = %.2f, at 0.05 = %.2f W m-2, difference %.1f%%\n",
    flux[1], flux[2], 100 * abs(flux[1] / flux[2] - 1))

cfs_b = GEMB.ClimateForcingStep(3600.0, 230.0, 80000.0, 0.0, 12.0, 0.0, 150.0, 10.0,
    250.0, 10.0, 0.0, 2.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
for solver in (ExplicitThermal(), ImplicitThermal())
    mp = GEMB.ModelParameters(thermal_conductivity_method=:Sturm, thermal_solver=solver,
        dt_divisors=divs)
    out = thermal(fill(260.0, 20), fill(0.1, 20), fill(100.0, 20), zeros(20), fill(0.3, 20),
        zeros(20), cfs_b, mp)
    @printf("case B %-16s surface temperature after 1 h = %.2f K\n", string(solver), out[1][1])
end
