# Finding 15: the default :Arthern grain growth, a firn law, barely coarsened dry surface snow.
#
# Sixty days of dry fresh snow (r_e = 0.05 mm, 300 kg m-3) at 258 and 268 K under 5 and
# 30 K m-1 gradients, with the default `grain_growth_method`.
#
# Expected: main (default :Arthern) gives r_e ≈ 0.14-0.21 mm whatever the gradient;
# fix/physics-audit (default :hybrid, Marbouty in snow) gives 0.30-0.46 mm under the 30 K m-1
# gradient.
#
# Run (see README.md): julia --project=<GEMB checkout> scripts/audit/10_grain_growth.jl
# with <GEMB checkout> this branch's checkout or a worktree of main at 0cda823.
include(joinpath(@__DIR__, "common.jl"))

mp = initialize_parameters()
println("default grain_growth_method = :", mp.grain_growth_method)
for T0 in (258.0, 268.0), grad in (5.0, 30.0)
    n = 20
    dz = fill(0.05, n)
    rr = fill(0.05, n); gd = fill(1.0, n); gs = fill(0.5, n)
    T = collect(T0 .- grad * 0.05 .* (0:n-1))
    cfs = step(dt=3600.0)
    for _ in 1:24*60
        rr, gd, gs = GEMB.calculate_grain_size(T, dz, fill(300.0, n), zeros(n), rr, gd, gs, cfs, mp)
    end
    @printf("T = %.0f K, gradient %4.0f K/m: r_e after 60 d = %.3f mm\n", T0, grad, rr[2])
end
