# Finding 5: snow albedo was computed for overhead sun, clear sky and no diffuse light because
# gemb replaced the forcing's solar_zenith_angle, cloud_optical_thickness and
# shortwave_downward_diffuse with the `mp` scalars (all 0) on every step.
#
# Two days on snow, once with a forcing that carries SZA = 60° and cloud optical thickness 2,
# once without.
#
# Expected: on main the two runs have identical albedo; on fix/physics-audit the run with
# geometry has a higher mean albedo.
#
# Run (see README.md): julia --project=<GEMB checkout> scripts/audit/09_radiation_geometry.jl
# with <GEMB checkout> this branch's checkout or a worktree of main at 0cda823.
include(joinpath(@__DIR__, "common.jl"))

n = 48
time = DateTime(2020, 6, 1) .+ Hour.(0:n-1)
sw = [max(0.0, 600.0 * sin(2π * (k - 6) / 24)) for k in 0:n-1]
base = (time, fill(265.0, n), fill(80000.0, n), zeros(n), fill(3.0, n), sw, fill(250.0, n),
    fill(250.0, n))
kw = (temperature_air_mean=265.0, wind_speed_mean=3.0, precipitation_mean=500.0,
    temperature_observation_height=2.0, wind_observation_height=10.0)
cf_none = initialize_forcing(base...; kw...)
# Scalars, which both versions accept; the branch also takes per-step vectors.
cf_geom = initialize_forcing(base...; kw..., solar_zenith_angle=60.0,
    cloud_optical_thickness=2.0)
mp = initialize_parameters(output_frequency=:all)
profile = initialize_profile(mp, cf_none; constant_density=true, constant_temperature=true)
k = findfirst(>=(1.0), cumsum(parent(profile[:dz])))
parent(profile[:density])[1:k] .= 350.0
parent(profile[:grain_radius])[1:k] .= 0.3
a(cf) = mean(parent(gemb(profile, cf, mp)[:albedo_broadband]))
@printf("mean albedo without geometry %.4f, with SZA 60° and τ 2: %.4f\n", a(cf_none), a(cf_geom))
