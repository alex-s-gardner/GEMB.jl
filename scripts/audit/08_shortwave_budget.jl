# Shortwave latent items: with subsurface absorption, the Brun-Lefebre branch recomputed its own
# grain-size albedos (ignoring the firn, ice and thin-snow albedo), and the density branch
# ignored the diffuse albedo, so absorbed shortwave disagreed with the reported albedo.
#
# 4 cm of snow over ice, SW↓ = 500 W m-2 of which 200 diffuse, subsurface absorption on.
#
# Expected: on main the absorbed total departs from (1−α)(SW↓−SW↓d) + (1−α_d)SW↓d by tens of
# W m-2 for :GardnerSharp and :BrunLefebre (and their albedos ignore the ice below the snow);
# on fix/physics-audit it matches for every method.
#
# Run (see README.md): julia --project=<GEMB checkout> scripts/audit/08_shortwave_budget.jl
# with <GEMB checkout> this branch's checkout or a worktree of main at 0cda823.
include(joinpath(@__DIR__, "common.jl"))

n = 400
dz = [0.04; fill(0.1, n - 1)]
density = [350.0; fill(917.0, n - 1)]
grain_radius = fill(0.5, n)
cfs = step(dt=3600.0, shortwave_downward=500.0, shortwave_downward_diffuse=200.0,
    solar_zenith_angle=40.0, cloud_optical_thickness=1.0, cloud_fraction=0.5)
for method in (:GardnerSharp, :BrunLefebre, :GreuellKonzelmann)
    mp = GEMB.ModelParameters(albedo_method=method, shortwave_subsurface_absorption=true,
        albedo_ice=0.35)
    α, α_d = albedo(dz, density, grain_radius, cfs, mp)
    sw = GEMB.calculate_shortwave_radiation(dz, density, grain_radius, α, α_d, cfs, mp)
    expected = (1 - α) * 300.0 + (1 - α_d) * 200.0
    @printf("%-18s α = %.4f α_d = %.4f   absorbed %.3f vs expected %.3f W m-2\n", method, α,
        α_d, sum(sw), expected)
end
