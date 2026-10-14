# Finding 3 and the rain-accounting latent item: rain joined the solid matrix of cell 1 at fixed
# dz (raising its density), and melt was reported as melt − rain, so rain frozen in a cold cell
# was never counted as refreeze.
#
# Case A: 5 kg m-2 of rain at 275.15 K on a 0.05 m, 350 kg m-3 cell at 0 °C.
# Case B: 2 kg m-2 of rain at 274.15 K on a 0.05 m, 350 kg m-3 cell at −10 °C.
#
# Expected: main leaves cell 1 of case A at 450 kg m-3 and case B with refreeze 0; on
# fix/physics-audit case A keeps 350 kg m-3 with 0.126 kg m-2 of ice melted by the rain's
# sensible heat, and case B refreezes [c·M·10 − P·c_w·1]/LF = 1.0745 kg m-2.
#
# Run (see README.md): julia --project=<GEMB checkout> scripts/audit/04_rain.jl
# with <GEMB checkout> this branch's checkout or a worktree of main at 0cda823.
include(joinpath(@__DIR__, "common.jl"))

mp = GEMB.ModelParameters(new_snow_method=:Constant150, density_ice=917.0,
    rain_temperature_threshold=273.15)
n = 20
for (name, T0, P, T_air) in (("A", GEMB.CtoK, 5.0, 275.15), ("B", GEMB.CtoK - 10.0, 2.0, 274.15))
    cfs = step(dt=86400.0, temperature_air=T_air, pressure_air=80000.0, precipitation=P,
        temperature_air_mean=260.0, precipitation_mean=200.0, temperature_observation_height=2.0,
        wind_observation_height=10.0)
    acc = GEMB.calculate_accumulation(fill(T0, n), fill(0.05, n), fill(350.0, n), zeros(n),
        fill(0.5, n), zeros(n), fill(0.5, n), zeros(n), cfs, mp, true)
    out = melt(acc[1:8]..., acc[9], mp)
    @printf("case %s: cell-1 density %.1f kg m-3 after rain, %.1f after melt; melt %.4f, refreeze %.4f, runoff %.4f kg m-2\n",
        name, acc[3][1], out[3][1], out[9], out[12], out[11])
end
c = 2102.0
@printf("expected case-B refreeze: %.4f kg m-2\n", (c * 17.5 * 10 - 2 * GEMB.HEAT_CAPACITY_WATER) / GEMB.LF)
