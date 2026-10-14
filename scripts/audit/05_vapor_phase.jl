# Finding 9: at a melting surface the latent flux used LV, but the vapor mass was added to or
# taken from the ice matrix, so LF per kg exchanged was created or destroyed.
#
# One gemb_core step of a wet column at 0 °C under warm humid air, which condenses onto it.
# The column's enthalpy change is checked against the fluxes it received plus the condensate
# entering as liquid water at the melting point (its latent heat having arrived through LE).
#
# Expected: main misses by ≈ −4.6e4 J m-2, the order of −LF·EC ≈ −4.9e4 (it stores the
# condensate as ice);
# fix/physics-audit closes to round-off.
#
# Run (see README.md): julia --project=<GEMB checkout> scripts/audit/05_vapor_phase.jl
# with <GEMB checkout> this branch's checkout or a worktree of main at 0cda823.
include(joinpath(@__DIR__, "common.jl"))

n = 10
mp = GEMB.ModelParameters(column_dzmin=0.05, column_dzmax=0.10, column_depth_max=0.9,
    dt_divisors=Float64.(GEMB.fast_divisors(36000000)) ./ 10000)
state = (temperature=fill(GEMB.CtoK, n), dz=fill(0.08, n), density=fill(400.0, n),
    water=[0.3; zeros(n - 1)], grain_radius=fill(0.5, n), grain_dendricity=zeros(n),
    grain_sphericity=fill(0.5, n), age=zeros(n), evaporation_condensation=0.0,
    melt_surface=0.0)
cfs = step(dt=3600.0, temperature_air=280.0, pressure_air=100000.0, wind_speed=8.0,
    longwave_downward=320.0, vapor_pressure=900.0, temperature_air_mean=260.0,
    wind_speed_mean=5.0, precipitation_mean=200.0, temperature_observation_height=2.0,
    wind_observation_height=2.0)
new, flux = GEMB.gemb_core(deepcopy(state), cfs, mp, false)
ec = new.evaporation_condensation
E(st) = GEMB.column_enthalpy(mp, st.dz .* st.density, st.temperature, st.water)
h_w = GEMB.specific_enthalpy_water(mp)
supplied = (flux.shortwave_net + cfs.longwave_downward - flux.longwave_upward +
            flux.heat_flux_sensible + flux.heat_flux_latent + flux.heat_flux_basal) * cfs.dt +
           ec * h_w + flux.E_added
residual = E(new) + flux.runoff * h_w - E(state) - supplied
@printf("condensed EC = %.4f kg m-2, latent heat flux = %.2f W m-2\n", ec, flux.heat_flux_latent)
@printf("energy residual with the condensate as liquid = %.1f J m-2 (−LF·EC = %.1f)\n",
    residual, -GEMB.LF * ec)
