# Findings 4 and 16: the stable branch switched from zL = Ri/(1−5Ri) to zL = Ri at Ri = 0.2,
# Ri used the effective height z_w²/z_T, and an Exner factor multiplied the sensible flux and Ri
#.
#
# Sensible flux over dry snow (z0 = 0.12 mm, wind at 10 m, temperature at 2 m) scanned in
# 0.01 K steps of T_air − T_s from 0.2 to 15 K, at 4 m s-1 and 80 kPa, reported as the
# exchange coefficient SHF/ΔT so a step's own 0.01 K does not count as a change; then the flux
# per unit air density at 100 and 60 kPa.
#
# Expected: main collapses to near zero below ΔT ≈ 1.6 K and jumps back by a factor of ~80 in
# one step, and gives a pressure ratio of 0.94; fix/physics-audit changes by under 1% per step
# and gives a pressure ratio of 1.0000.
#
# Run (see README.md): julia --project=<GEMB checkout> scripts/audit/03_stability.jl
# with <GEMB checkout> this branch's checkout or a worktree of main at 0cda823.
include(joinpath(@__DIR__, "common.jl"))

T_s = 265.0
z0 = GEMB.Z0_SNOW_DRY
cfs(ΔT; p=80000.0, U=4.0) = step(dt=3600.0, temperature_air=T_s + ΔT, pressure_air=p,
    wind_speed=U, vapor_pressure=300.0, temperature_observation_height=2.0,
    wind_observation_height=10.0)
ΔTs = 0.2:0.01:15.0
shf = [GEMB.turbulent_heat_flux(T_s, 1.2, z0, 0.1z0, 0.1z0, cfs(ΔT))[1] for ΔT in ΔTs]
coef = shf ./ ΔTs
jumps = abs.(diff(coef)) ./ coef[1:end-1]
k = argmax(jumps)
@printf("largest adjacent change of SHF/ΔT: %.1f%% between ΔT = %.2f and %.2f K (SHF %.2f -> %.2f W m-2)\n",
    100 * jumps[k], ΔTs[k], ΔTs[k+1], shf[k], shf[k+1])
for ΔT in (1.0, 3.0, 6.0, 10.0)
    @printf("ΔT = %4.1f K: SHF = %6.2f W m-2\n", ΔT, shf[findfirst(≈(ΔT), ΔTs)])
end
r = map((100000.0, 60000.0)) do p
    ρa = GEMB.air_density(p, T_s + 1.0)
    GEMB.turbulent_heat_flux(T_s, ρa, z0, 0.1z0, 0.1z0, cfs(1.0; p))[1] / ρa
end
@printf("SHF/ρ_air at 100 kPa over 60 kPa = %.4f\n", r[1] / r[2])
