# Finding 1: bare-ice albedo was α_ice + (0.58 − α_ice)·exp(−M/200) with M a few kg m-2 of
# recent melt, so exposed ice sat near 0.58 and `mp.albedo_ice` was nearly inert.
#
# Twenty days of an all-ice column with no snowfall, at albedo_ice 0.30 and 0.48.
#
# Expected: main reports a mean albedo near 0.55-0.58 for both and an absorbed-shortwave
# ratio near 1; fix/physics-audit reports 0.30 and 0.48 and a ratio of 0.70/0.52 = 1.346.
#
# Run (see README.md): julia --project=<GEMB checkout> scripts/audit/02_bare_ice_albedo.jl
# with <GEMB checkout> this branch's checkout or a worktree of main at 0cda823.
include(joinpath(@__DIR__, "common.jl"))

n = 20 * 24
time = DateTime(2000, 7, 1) .+ Hour.(0:n-1)
sw = [max(0.0, 700.0 * sin(2π * (hour(t) - 6) / 24)) for t in time]
cf = initialize_forcing(time, fill(274.0, n), fill(85000.0, n), zeros(n), fill(4.0, n), sw,
    fill(300.0, n), fill(500.0, n); temperature_air_mean=268.0, wind_speed_mean=4.0,
    precipitation_mean=0.0, temperature_observation_height=2.0, wind_observation_height=10.0)
sw_net = map((0.30, 0.48)) do albedo_ice
    mp = ModelParameters(output_frequency=:daily, albedo_ice=albedo_ice)
    out = gemb(initialize_profile(mp, cf; constant_density=true, constant_temperature=true),
        cf, mp)
    @printf("albedo_ice = %.2f   mean output albedo = %.4f   mean absorbed SW = %.2f W m-2\n",
        albedo_ice, mean(parent(out[:albedo_broadband])), mean(parent(out[:shortwave_net])))
    sum(parent(out[:shortwave_net]))
end
@printf("absorbed-SW ratio (0.30 vs 0.48) = %.4f   (expected with albedo_ice in force: %.4f)\n",
    sw_net[1] / sw_net[2], 0.70 / 0.52)
