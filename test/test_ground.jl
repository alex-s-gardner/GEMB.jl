# The `lower_boundary = :ground` soil column: properties, analytic solutions, budgets and the
# variable snow column over it.

using Dates: month, hour, DateTime, Hour
using GEMB: ClimateForcingStep, CtoK, LF, HEAT_CAPACITY_WATER, MATERIAL_SNOW, MATERIAL_GROUND,
    MATERIAL_EMPTY, ThermalWorkspace

# erf by composite Simpson on [0, x]; agrees with erf(1) to 1e-15.
function _erf_simpson(x; n=4000)
    x == 0 && return 0.0
    h = x / n
    s = 1.0 + exp(-x^2)
    for i in 1:n-1
        s += (isodd(i) ? 4 : 2) * exp(-(i * h)^2)
    end
    return 2 / sqrt(pi) * s * h / 3
end

_ground_mp(; dt=3600.0, kw...) = GEMB.initialize_parameters(; lower_boundary=:ground,
    dt_divisors=GEMB.fast_divisors(Int(dt) * 10000) ./ 10000, kw...)

_still_air_step(dt, T_air) = ClimateForcingStep(dt=dt, temperature_air=T_air, pressure_air=9e4,
    wind_speed=1.0, longwave_downward=0.0, vapor_pressure=100.0,
    temperature_observation_height=2.0, wind_observation_height=10.0)

# One thermal step of snow-free ground whose surface exchanges nothing with the air (zero
# emissivity and transfer coefficient), with `sw` [W m-2] absorbed in the top cell. Advances `H`.
function _ground_conduction_step!(H, mp, g, T_bottom, cfs; sw=0.0)
    ng = length(g.dz)
    liquid = [GEMB.ground_state(mp, g.C_solid[j], g.water_total[j], H[j])[2] for j in 1:ng]
    K = [GEMB.ground_conductivity(g, liquid[j], g.water_total[j]) for j in 1:ng]
    gn = GEMB._GroundNodes(copy(g.C_solid), copy(g.water_total), H, K, g.dz, T_bottom)
    tab = GEMB._stability_table!(ThermalWorkspace(), 10.0, 2.0, 0.01, 0.001)
    sfc = GEMB._ThermalSurface(1.2, 0.0, false, 0.0, 0.0, tab)
    swf = zeros(ng)
    swf[1] = sw
    return GEMB._thermal_solve_ground!(Float64[], Float64[], Float64[], Float64[], gn, swf, sfc,
        0.0, false, 1.0, cfs, mp)
end

@testset "ground: Johansen conductivity and enthalpy" begin
    mp = _ground_mp(ground_temperature_deep=276.0)
    g = GEMB.ground_column(mp, GEMB.ground_grid(10.0))
    @test sum(g.dz) ≈ 10.0
    @test GEMB.ground_grid(10.0)[1:3] ≈ [0.05, 0.05, 0.1]
    W = g.water_total[1]
    # ISBA's Kersten forms at Sr = θ/φ = 0.25/0.45, by hand.
    Sr = 0.25 / 0.45
    ksat_u = 3.0^0.55 * 0.57^0.45
    ksat_f = 3.0^0.55 * 2.22^0.45
    kdry = (0.135 * 0.55 * 2700 + 64.7) / (2700 - 0.947 * 0.55 * 2700)
    @test g.conductivity_dry ≈ kdry
    @test GEMB.ground_conductivity(g, W, W) ≈ (log10(Sr) + 1) * (ksat_u - kdry) + kdry
    @test GEMB.ground_conductivity(g, 0.0, W) ≈ Sr * (ksat_f - kdry) + kdry
    @test GEMB.ground_conductivity(g, 0.0, W) > GEMB.ground_conductivity(g, W, W)
    # Enthalpy and its inverse on all three branches, and the plateau width.
    Cs = g.C_solid[1]
    for (T, l) in ((250.0, 0.0), (270.0, 0.0), (CtoK, 0.0), (CtoK, 0.4W), (CtoK, W), (280.0, W))
        T2, l2 = GEMB.ground_state(mp, Cs, W, GEMB.ground_enthalpy(mp, Cs, W, T, l))
        @test T2 ≈ T atol = 1e-9
        @test l2 ≈ (T < CtoK ? 0.0 : T > CtoK ? W : l) atol = 1e-9
    end
    @test GEMB.ground_enthalpy(mp, Cs, W, CtoK, W) - GEMB.ground_enthalpy(mp, Cs, W, CtoK, 0.0) ≈
          W * LF
    # Volumetric heat capacity (1 − φ)ρ_s c_s + θρ_w c_w above the plateau.
    H1 = GEMB.ground_enthalpy(mp, Cs, W, 280.0, W)
    H2 = GEMB.ground_enthalpy(mp, Cs, W, 281.0, W)
    @test (H2 - H1) / g.dz[1] ≈ 0.55 * 2700 * 733 + 0.25 * 1000 * HEAT_CAPACITY_WATER
end

@testset "ground: parameters are validated" begin
    @test_throws "ground_temperature_deep" GEMB.initialize_parameters(lower_boundary=:ground)
    @test_throws "ground_water must be in" GEMB.initialize_parameters(lower_boundary=:ground,
        ground_temperature_deep=275.0, ground_water=0.5, ground_porosity=0.4)
    @test_throws "ExplicitThermal" GEMB.initialize_parameters(lower_boundary=:ground,
        ground_temperature_deep=275.0, thermal_solver=ImplicitThermal())
    @test_throws "lower_boundary must be" GEMB.initialize_parameters(lower_boundary=:soil)
    @test GEMB.initialize_parameters().lower_boundary === :ice
end

# Conduction in unfrozen soil heated by a sinusoidal surface flux F0·cos(ωt): the periodic
# solution T = Tm + A·exp(−z/d)·cos(ωt − z/d − π/4), d = √(2κ/ω), A = F0·d/(k√2) (Carslaw and
# Jaeger, 1959, §2.6). The column starts on the periodic solution, 10 damping depths deep.
@testset "ground: periodic conduction against the analytic solution" begin
    for (P, F0, dt) in ((86400.0, 100.0, 900.0), (5 * 86400.0, 40.0, 3600.0))
        Tm = 283.15
        mp = _ground_mp(dt=dt, ground_temperature_deep=Tm)
        ω = 2π / P
        dzc = 0.01
        g0 = GEMB.ground_column(mp, [dzc])
        W0 = g0.water_total[1]
        k = GEMB.ground_conductivity(g0, W0, W0)
        Cv = (g0.C_solid[1] + W0 * HEAT_CAPACITY_WATER) / dzc
        d = sqrt(2k / (Cv * ω))
        n = ceil(Int, 10d / dzc)
        g = GEMB.ground_column(mp, fill(dzc, n))
        A = F0 * d / (k * sqrt(2))
        T_exact(z, t) = Tm + A * exp(-z / d) * cos(ω * t - z / d - π / 4)
        zc = [(j - 0.5) * dzc for j in 1:n]
        H = [GEMB.ground_enthalpy(mp, g.C_solid[j], W0, T_exact(zc[j], 0.0), W0) for j in 1:n]
        E0 = sum(H)
        E_in = 0.0
        t = 0.0
        err = 0.0
        steps = round(Int, 6P / dt)
        for s in 1:steps
            sw = F0 * cos(ω * (t + dt / 2))
            r = _ground_conduction_step!(H, mp, g, Tm, _still_air_step(dt, Tm); sw)
            E_in += (sw + r[4]) * dt
            t += dt
            if s > steps - round(Int, P / dt)
                for j in 1:n
                    zc[j] < 4d || continue
                    T = GEMB.ground_state(mp, g.C_solid[j], W0, H[j])[1]
                    err = max(err, abs(T - T_exact(zc[j], t)))
                end
            end
        end
        # Measured: 0.0156 K of A = 7.6 K (daily), 0.0131 K of 6.8 K (5-day).
        @test err / A < 0.005
        @test abs(sum(H) - E0 - E_in) < 1e-6 * E0 * 1e-6 + 1e-3
    end
end

# Saturated soil at the melting point frozen from a face held ΔT below it: the one-phase
# Neumann solution X = 2λ√(α t), λ·exp(λ²)·erf(λ) = St/√π, St = C_f ΔT/(ρ_w θ L_F), with α and
# C_f the frozen soil's diffusivity and volumetric heat capacity (Carslaw and Jaeger, 1959,
# §11.2). The cold face is the Dirichlet bottom cell; the surface is insulated.
@testset "ground: Stefan freezing front against the Neumann solution" begin
    ΔT = 10.0
    φ = 0.4
    dt = 3600.0
    mp = _ground_mp(dt=dt, ground_temperature_deep=CtoK - ΔT, ground_porosity=φ,
        ground_water=φ, heat_capacity_method=:constant)
    dzc = 0.01
    n = 100
    g = GEMB.ground_column(mp, vcat(fill(dzc, n), 1e-4))
    W = g.water_total[1]
    kf = GEMB.ground_conductivity(g, 0.0, W)
    Cf = (g.C_solid[1] + W * mp.heat_capacity_ice) / dzc
    St = Cf * ΔT / (W / dzc * LF)
    f(λ) = λ * exp(λ^2) * _erf_simpson(λ) - St / sqrt(pi)
    lo, hi = 1e-6, 2.0
    for _ in 1:80
        mid = (lo + hi) / 2
        f(mid) > 0 ? (hi = mid) : (lo = mid)
    end
    λ = (lo + hi) / 2
    α = kf / Cf
    H = [GEMB.ground_enthalpy(mp, g.C_solid[j], g.water_total[j], CtoK, g.water_total[j])
         for j in 1:n+1]
    t = 0.0
    errs = Float64[]
    for s in 1:240
        _ground_conduction_step!(H, mp, g, CtoK - ΔT, _still_air_step(dt, CtoK))
        t += dt
        if s % 24 == 0
            frozen = sum(g.water_total[j] -
                         GEMB.ground_state(mp, g.C_solid[j], g.water_total[j], H[j])[2] for j in 1:n)
            push!(errs, (frozen / (W / dzc) - 2λ * sqrt(α * t)) / (2λ * sqrt(α * t)))
        end
    end
    # Measured: −2.6 % after 1 day (18 cells frozen), −0.48 % after 10 days.
    @test abs(errs[1]) < 0.03
    @test abs(errs[end]) < 0.006
    @test issorted(abs.(errs[2:end]); rev=true)
end

include("synthetic_maritime_forcing.jl")
include("synthetic_seasonal_forcing.jl")

@testset "ground: a year of snow on soil closes mass and energy" begin
    cf = synthetic_seasonal_forcing(years=1)
    mp = GEMB.initialize_parameters(lower_boundary=:ground, ground_temperature_deep=276.0,
        output_frequency=:daily)
    p = initialize_profile(mp, cf)
    ng = length(GEMB.ground_grid(10.0))
    # No snow: the ground starts at row 1 and every snow slot is empty below it.
    @test all(==(MATERIAL_GROUND), parent(p[:material])[1:ng])
    @test all(==(MATERIAL_EMPTY), parent(p[:material])[ng+1:end])
    @test length(parent(p[:material])) == mp.ground_snow_cells_max + ng
    @test all(==(276.0), filter(!isnan, parent(p[:temperature])))
    # `verbose` checks the per-step mass and energy budgets and the whole-run snow mass budget.
    o = gemb(p, cf, mp; verbose=true)
    mat = parent(o[:material])
    dz = parent(o[:dz])
    # Ground cells never change thickness; the bottom cell stays at the deep temperature.
    gdz = GEMB.ground_grid(10.0)
    for j in axes(mat, 2)
        @test dz[mat[:, j] .== MATERIAL_GROUND, j] == gdz
        @test parent(o[:temperature])[findlast(==(MATERIAL_GROUND), mat[:, j]), j] == 276.0
        # Empty slots hold nothing.
        e = mat[:, j] .== MATERIAL_EMPTY
        @test all(==(0.0), dz[e, j]) && all(==(0.0), parent(o[:density])[e, j])
    end
    # Snow outputs agree with the profile.
    swe = vec(sum(dz .* parent(o[:density]) .* (mat .== MATERIAL_SNOW), dims=1) .+
              sum(parent(o[:water]) .* (mat .== MATERIAL_SNOW), dims=1))
    @test parent(o[:snow_water_equivalent]) ≈ swe
end

@testset "ground: melt-out and re-accumulation through the empty slots" begin
    cf = synthetic_seasonal_forcing(years=2)
    mp = GEMB.initialize_parameters(lower_boundary=:ground, ground_temperature_deep=276.0,
        output_frequency=:daily)
    o = gemb(initialize_profile(mp, cf), cf, mp; verbose=true)
    hs = parent(o[:snow_depth])
    mat = parent(o[:material])
    ns = vec(count(==(MATERIAL_SNOW), mat, dims=1))
    # Snow rows are the cells of non-zero thickness above the ground.
    @test ns == [count(>(0.0), parent(o[:dz])[mat[:, j] .!= MATERIAL_GROUND, j])
                 for j in axes(mat, 2)]
    @test hs ≈ vec(sum(parent(o[:dz]) .* (mat .== MATERIAL_SNOW), dims=1))
    # Snow-free days exist, and snow returns after them.
    free = findall(==(0.0), hs)
    @test !isempty(free)
    @test any(j -> j > first(free) && hs[j] > 0.1, eachindex(hs))
    # Snow-free ground is at the ground albedo; the surface rises above freezing then.
    # Albedo is an interval mean, so test days snow-free from the day before through the day.
    bare = [j for j in free if j > 1 && hs[j-1] == 0.0 && sum(parent(o[:precipitation])[j:j]) == 0.0]
    @test !isempty(bare)
    @test all(≈(mp.ground_albedo; atol=1e-12), parent(o[:albedo_broadband])[bare])
    @test maximum(parent(o[:soil_temperature_005])[free]) > CtoK + 1
    # Snow-free days carry no ground-to-snow flux sign constraint, but the snow-covered winter
    # draws heat from the ground.
    covered = findall(>(0.3), hs)
    @test sum(parent(o[:heat_flux_ground])[covered]) / length(covered) > 0
end

# Thin snow over ground: every cell below its thickness band, so each merges downward and the
# bottom cell has no surviving cell above it to merge into. It must hold the whole column.
@testset "ground: a column of only sub-band snow cells merges into its bottom cell" begin
    mp = GEMB.initialize_parameters(lower_boundary=:ground, ground_temperature_deep=275.0)
    for n in (2, 3)
        T = collect(range(269.0, 271.0; length=n)); dz = fill(0.004, n)
        ρ = collect(range(300.0, 320.0; length=n)); W = zeros(n)
        r = fill(0.1, n); d = zeros(n); s = ones(n); a = ones(n)
        M0 = sum(dz .* ρ)
        E0 = GEMB.column_mass_energy(GEMB.column_state(T, dz, ρ, W, r, d, s, a), mp)[2]
        dzmin = zeros(n); GEMB.column_bands!(dzmin, zeros(n), dz, mp)
        @test all(dz .< dzmin)
        out = GEMB.manage_layer_thickness(T, dz, ρ, W, r, d, s, a, mp, true; dt=3600.0,
            n_target=n, n_max=GEMB.ground_snow_slots(mp), dirichlet_bottom=false)
        @test length(out[2]) == 1
        @test out[2][1] * out[3][1] ≈ M0
        E1 = GEMB.column_mass_energy(GEMB.column_state(out[1:8]...), mp)[2]
        @test E1 ≈ E0 rtol = 1e-14
    end
end

# Snow-free ground in calm, clear summer weather (E2 of the SNOTEL test): without free convection
# the transfer coefficient of a calm hour vanishes and dry sunlit soil has no turbulent loss.
@testset "ground: calm sunlit dry soil stays bounded (free convection)" begin
    n = 24 * 5
    t = [DateTime(2009, 6, 5) + Hour(h) for h in 0:n-1]
    hr = hour.(t)
    sw = [max(0.0, 950.0 * sin(π * (h - 6) / 14)) for h in hr]
    Ta = [286.0 + 5.0 * sin(π * (h - 9) / 12) for h in hr]
    cf = initialize_forcing(t, Ta, fill(73000.0, n), zeros(n), zeros(n), sw, fill(287.0, n),
        fill(543.0, n); temperature_air_mean=280.0, wind_speed_mean=2.0,
        precipitation_mean=500.0, temperature_observation_height=2.7,
        wind_observation_height=10.0, solar_zenith_angle=fill(30.0, n),
        cloud_optical_thickness=zeros(n), shortwave_downward_diffuse=0.1 .* sw)
    mp = GEMB.initialize_parameters(lower_boundary=:ground, ground_temperature_deep=283.0,
        ground_water=0.05, output_frequency=:all)
    o = gemb(initialize_profile(mp, cf), cf, mp; verbose=true)
    T5 = parent(o[:soil_temperature_005])
    # Measured: the 5 cm soil peaks at 296.4 K; without free convection the surface node passes
    # 332 K on the first afternoon and the run stops.
    @test maximum(T5) < 305.0
    # Free convection carries the midday heat away as sensible and latent flux.
    noon = findall(==(13), hr)
    @test all(<(-100.0), parent(o[:heat_flux_sensible])[noon] .+ parent(o[:heat_flux_latent])[noon])
end

# Melting snow over soil warmer than 0 °C under a strong warm wind (E2): snow cannot be warmer
# than the melting point, so heat crosses the snow–ground face only from the soil upward, and the
# turbulent heat the snow absorbs melts it rather than heating the soil below.
@testset "ground: snow over warm soil does not heat the soil (interface cap)" begin
    mp = GEMB.initialize_parameters(lower_boundary=:ground, ground_temperature_deep=276.0,
        dt_divisors=GEMB.fast_divisors(3600 * 10000) ./ 10000)
    g = GEMB.ground_column(mp, GEMB.ground_grid(mp.ground_depth))
    Tg = fill(276.0, length(g.dz))
    snow = GEMB.column_state([CtoK], [0.03], [700.0], [0.1], [0.16], [0.0], [0.0], [150.0])
    cfs = ClimateForcingStep(; dt=3600.0, temperature_air=287.0, pressure_air=73760.0,
        wind_speed=22.8, longwave_downward=256.0, vapor_pressure=1049.0,
        temperature_observation_height=2.75, wind_observation_height=10.0,
        precipitation_mean=320.0, temperature_air_mean=275.0, wind_speed_mean=10.0)
    s, Tg1, _, flux = GEMB.gemb_core_ground(snow, copy(Tg), copy(g.water_total), g, cfs, mp, true)
    # Measured: G = +58.7 W m-2 and the top soil cools to 274.5 K. Without the cap the 0 °C snow
    # conducted 354 W m-2 into the soil, which warmed to 285.7 K.
    @test flux.heat_flux_ground > 0.0
    @test Tg1[1] < 276.0
    @test flux.melt > 0.0
end

@testset "ground: glacier column is untouched by the ground fields" begin
    cf = synthetic_maritime_forcing(years=1)
    mp = GEMB.initialize_parameters(output_frequency=:daily)
    mpg = GEMB.initialize_parameters(output_frequency=:daily, ground_albedo=0.3,
        ground_depth=5.0, ground_temperature_deep=270.0)
    p = initialize_profile(mp, cf; constant_density=true, constant_temperature=true)
    o1 = gemb(p, cf, mp)
    o2 = gemb(p, cf, mpg)
    @test all(k -> isequal(parent(o1[k]), parent(o2[k])), keys(o1))
    @test !haskey(o1, :material)
end

# The idealized snow–soil conduction model of the SNOTEL test design (§3.4): 5 cm cells of moist
# soil (θ = 0.25, k 1.0 unfrozen and 1.6 W m-1 K-1 frozen, the pore water freezing linearly over
# −1…0 °C as an apparent heat capacity, solids 2700 kg m-3 at 733 J kg-1 K-1), its surface coupled
# to a prescribed snow-surface temperature `Ts` [K] through a snow resistance `R` [m2 K W-1], and
# its base at `MAGT` [°C] 10 m down. Returns the hourly upward flux into the snow [W m-2] and the
# snow-base temperature [°C].
function _idealized_soil(Ts, R, dt, MAGT, T0; dz=0.05, nsub=6)
    function props(Tc)
        fliq = clamp(Tc + 1.0, 0.0, 1.0)
        kc = fliq * 1.0 + (1 - fliq) * 1.6
        cc = 2700 * 733 * 0.55 + 0.25 * (fliq * 4.18e6 + (1 - fliq) * 1.9e6)
        Lc = (Tc > -1.0 && Tc < 0.0) ? 0.25 * 1000 * 3.34e5 : 0.0
        return kc, cc + Lc
    end
    T = copy(T0)
    n = length(T)
    h = dt / nsub
    G = zeros(length(Ts))
    Tb = zeros(length(Ts))
    k = similar(T)
    C = similar(T)
    for s in eachindex(Ts)
        for _ in 1:nsub
            for i in 1:n
                k[i], C[i] = props(T[i])
            end
            Gtop = (T[1] - (Ts[s] - CtoK)) / (R[s] + dz / 2 / k[1])
            Tn = copy(T)
            for i in 1:n
                Fup = i == 1 ? Gtop : (T[i] - T[i-1]) / (dz / 2 / k[i] + dz / 2 / k[i-1])
                Fdn = i == n ? (MAGT - T[i]) / (dz / 2 / k[i]) :
                      (T[i+1] - T[i]) / (dz / 2 / k[i] + dz / 2 / k[i+1])
                Tn[i] = T[i] + h * (Fdn - Fup) / (C[i] * dz)
            end
            T = Tn
            G[s] += Gtop / nsub
        end
        Tb[s] = T[1] - G[s] * dz / 2 / props(T[1])[1]
    end
    return G, Tb
end

# Seasonal snow on soil (Colorado-like: mean air −0.2 °C, ±9 K, 900 kg m-2 yr-1, MAGT +3.1 °C),
# spun up on the deep-soil temperature. GEMB's snow base must sit near 0 °C over a zero curtain in
# the top soil, and agree with the idealized model of the design driven by GEMB's own snow-surface
# temperature and snow thermal resistance, which isolates the soil physics.
@testset "ground: seasonal snow on soil against the idealized model" begin
    MAGT = 276.25
    cf = synthetic_seasonal_forcing(years=1, temperature_mean=272.95, amplitude=9.0,
        precipitation_annual=900.0, emissivity_clear=0.85)
    mp = GEMB.initialize_parameters(lower_boundary=:ground, ground_temperature_deep=MAGT)
    sp = gemb_spinup_ground(initialize_profile(mp, cf), cf, mp; simulation_years_maximum=30)
    @test GEMB.DD.metadata(sp).spinup_converged
    o = gemb(sp, cf, mp)
    g = GEMB.ground_column(mp, GEMB.ground_grid(mp.ground_depth))
    mat = parent(o[:material]); Tp = parent(o[:temperature]); dz = parent(o[:dz])
    ρ = parent(o[:density]); W = parent(o[:water]); Gg = parent(o[:heat_flux_ground])
    nt = size(mat, 2)
    R = zeros(nt)
    Tb = zeros(nt)
    for j in 1:nt
        for i in axes(mat, 1)
            mat[i, j] == MATERIAL_SNOW &&
                (R[j] += dz[i, j] / GEMB.thermal_conductivity(Tp[i, j], ρ[i, j], mp))
        end
        i = findfirst(==(MATERIAL_GROUND), mat[:, j])
        k1 = GEMB.ground_conductivity(g, W[i, j], g.water_total[1])
        Tb[j] = Tp[i, j] + Gg[j] * g.dz[1] / (2k1) - CtoK
    end
    Tg0 = parent(sp[:temperature])[parent(sp[:material]) .== MATERIAL_GROUND]
    T0 = [GEMB.ground_soil_temperature(Tg0, g, (i - 0.5) * 0.05) - CtoK for i in 1:200]
    Gi, Tbi = _idealized_soil(Tp[1, :], R, 3600.0, MAGT - CtoK, T0)
    t = collect(GEMB.DD.dims(o, Ti))
    djf = findall(x -> month(x) in (12, 1, 2), t)
    mean(x) = sum(x) / length(x)
    # Measured: snow base −1.91 °C, G 5.20 W m-2, 43 days with the 5 cm soil within 0.01 K of
    # 0 °C.
    @test all(>(0.0), R[djf])
    @test -2.5 < mean(Tb[djf]) < 0.0
    @test abs(mean(Tb[djf]) - mean(Tbi[djf])) < 0.25
    @test abs(mean(Gg[djf]) - mean(Gi[djf])) < 0.15 * mean(Gi[djf])
    @test mean(Gg[djf]) > 0
    T5 = parent(o[:soil_temperature_005]) .- CtoK
    @test count(x -> abs(x) < 0.01, T5) / 24 > 20
end
