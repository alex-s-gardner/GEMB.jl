# The `lower_boundary = :prescribed` base: snow on a fixed-temperature reservoir whose
# temperature is the `temperature_ground` forcing layer.

using Dates: dayofyear
using GEMB: ClimateForcingStep, CtoK, MATERIAL_SNOW, MATERIAL_EMPTY, ThermalWorkspace, dims

isdefined(@__MODULE__, :synthetic_seasonal_forcing) || include("synthetic_seasonal_forcing.jl")
isdefined(@__MODULE__, :synthetic_maritime_forcing) || include("synthetic_maritime_forcing.jl")

# `cf` with its base-temperature layer replaced by `Tb` [K], every other layer and the metadata
# unchanged.
_with_base(cf, Tb) = GEMB.DD.rebuild(cf;
    data=merge(getfield(cf, :data), (; temperature_ground=Float64.(collect(Tb)))))

# A seasonal base, +3 °C in late July and −3 °C in late January.
_seasonal_base(cf) = [CtoK + 3.0 * cos(2π * (dayofyear(t) - 200) / 365)
                      for t in dims(cf.temperature_air, Ti)]

@testset "prescribed: the forcing layer and parameters" begin
    cf = synthetic_seasonal_forcing(years=1)
    # Forcing built without it carries a NaN layer, read by no other boundary.
    @test all(isnan, parent(cf.temperature_ground))
    mp = GEMB.initialize_parameters(lower_boundary=:prescribed)
    @test_throws "temperature_ground at every step" gemb(initialize_profile(mp, cf), cf, mp)
    Tb = _seasonal_base(cf)
    Tb[100] = NaN
    @test_throws "it is NaN at step 100" gemb(initialize_profile(mp, cf), _with_base(cf, Tb), mp)
    @test_throws "ExplicitThermal" GEMB.initialize_parameters(lower_boundary=:prescribed,
        thermal_solver=ImplicitThermal())
    @test_throws "lower_boundary must be ice, ground or prescribed" GEMB.initialize_parameters(
        lower_boundary=:soil)
    @test_throws "has no state to spin up" gemb_spinup(initialize_profile(mp, cf),
        _with_base(cf, _seasonal_base(cf)), mp)
    # initialize_forcing passes the layer through, from vectors or a DimStack.
    t = collect(dims(cf.temperature_air, Ti))
    v(k) = collect(parent(getproperty(cf, k)))
    cf2 = initialize_forcing(t, v(:temperature_air), v(:pressure_air), v(:precipitation),
        v(:wind_speed), v(:shortwave_downward), v(:longwave_downward), v(:vapor_pressure);
        temperature_air_mean=cf.temperature_air_mean, wind_speed_mean=cf.wind_speed_mean,
        precipitation_mean=cf.precipitation_mean,
        temperature_observation_height=cf.temperature_observation_height,
        wind_observation_height=cf.wind_observation_height, temperature_ground=Tb)
    @test isequal(parent(cf2.temperature_ground), Tb)
    @test_throws "temperature_ground values unrealistic" initialize_forcing(t,
        v(:temperature_air), v(:pressure_air), v(:precipitation), v(:wind_speed),
        v(:shortwave_downward), v(:longwave_downward), v(:vapor_pressure);
        temperature_air_mean=cf.temperature_air_mean, wind_speed_mean=cf.wind_speed_mean,
        precipitation_mean=cf.precipitation_mean,
        temperature_observation_height=cf.temperature_observation_height,
        wind_observation_height=cf.wind_observation_height, temperature_ground=Tb .- 273.15)
end

@testset "prescribed: the glacier column ignores the base layer" begin
    cf = synthetic_maritime_forcing(years=1)
    mp = GEMB.initialize_parameters(output_frequency=:daily)
    p = initialize_profile(mp, cf; constant_density=true, constant_temperature=true)
    o1 = gemb(p, cf, mp)
    o2 = gemb(p, _with_base(cf, fill(250.0, length(cf.temperature_air))), mp)
    @test all(k -> isequal(parent(o1[k]), parent(o2[k])), keys(o1))
end

# Snow slab of thickness L with a sealed surface (zero emissivity and transfer coefficient) on a
# base at Tm + A·cos(ωt). Periodic solution, x the depth below the surface:
# T = Tm + A·Re[cosh(qx)/cosh(qL)·e^{iωt}], q = √(iω/κ) (Carslaw and Jaeger, 1959, §3.6).
@testset "prescribed: snow slab on a sinusoidal base against the analytic solution" begin
    for (L, n, P, dt, tol) in ((0.5, 50, 86400.0, 300.0, 0.006), (0.3, 30, 5 * 86400.0, 900.0, 0.006))
        A = 5.0; Tm = 263.15; ρ = 300.0
        mp = GEMB.initialize_parameters(lower_boundary=:prescribed, heat_capacity_method=:constant,
            dt_divisors=GEMB.fast_divisors(Int(dt) * 10000) ./ 10000)
        dzc = L / n
        κ = GEMB.thermal_conductivity(Tm, ρ, mp) / (ρ * mp.heat_capacity_ice)
        ω = 2π / P
        q = sqrt(im * ω / κ)
        T_exact(x, t) = Tm + A * real(cosh(q * x) / cosh(q * L) * exp(im * ω * t))
        xc = [(i - 0.5) * dzc for i in 1:n]
        T = [T_exact(x, 0.0) for x in xc]
        dz = fill(dzc, n); den = fill(ρ, n); W = zeros(n)
        tab = GEMB._stability_table!(ThermalWorkspace(), 10.0, 2.0, 0.001, 0.0001)
        sfc = GEMB._ThermalSurface(1.2, 0.0, false, 0.0, 0.0, tab)
        cfs = ClimateForcingStep(dt=dt, temperature_air=Tm, pressure_air=9e4, wind_speed=1.0,
            longwave_downward=0.0, vapor_pressure=100.0, temperature_observation_height=2.0,
            wind_observation_height=10.0)
        H(T) = sum(ρ .* dz .* GEMB.specific_enthalpy.(Ref(mp), T))
        E0 = H(T); E_in = 0.0; t = 0.0; err = 0.0
        steps = round(Int, 4P / dt)
        for s in 1:steps
            Tb = Tm + A * cos(ω * (t + dt / 2))
            r = GEMB._thermal_solve_prescribed!(T, dz, den, W, Tb, zeros(n + 1), sfc, 0.0,
                false, cfs, mp)
            E_in += r[4] * dt
            t += dt
            s > steps - round(Int, P / dt) &&
                (err = max(err, maximum(abs.(T .- [T_exact(x, t) for x in xc]))))
        end
        # Measured: 0.47 % of A (daily, 300 s), 0.51 % (5-day, 900 s). The error is the base
        # being held at its step-midpoint value through each step; it falls with the step.
        @test err / A < tol
        # Energy through the base is the basal flux.
        @test abs(H(T) - E0 - E_in) < 1e-6
    end
end

@testset "prescribed: two seasonal years close mass and energy, melt out and re-accumulate" begin
    cf = _with_base(synthetic_seasonal_forcing(years=2), _seasonal_base(synthetic_seasonal_forcing(years=2)))
    mp = GEMB.initialize_parameters(lower_boundary=:prescribed, output_frequency=:daily)
    p = initialize_profile(mp, cf)
    @test length(parent(p[:dz])) == mp.ground_snow_cells_max
    @test all(==(MATERIAL_EMPTY), parent(p[:material]))
    # `verbose` checks every step's mass and energy budget, the energy through the base counted
    # as heat_flux_basal, and the whole-run snow mass budget.
    o = gemb(p, cf, mp; verbose=true)
    hs = parent(o[:snow_depth])
    mat = parent(o[:material]); dz = parent(o[:dz]); ρ = parent(o[:density]); W = parent(o[:water])
    @test hs ≈ vec(sum(dz .* (mat .== MATERIAL_SNOW), dims=1))
    @test parent(o[:snow_water_equivalent]) ≈ vec(sum((dz .* ρ .+ W) .* (mat .== MATERIAL_SNOW), dims=1))
    # Empty slots hold nothing, below the snow.
    for j in axes(mat, 2)
        ns = count(==(MATERIAL_SNOW), mat[:, j])
        @test all(==(MATERIAL_SNOW), mat[1:ns, j]) && all(==(MATERIAL_EMPTY), mat[ns+1:end, j])
        @test all(==(0.0), dz[ns+1:end, j])
    end
    # Snow-free days, then snow again.
    free = findall(==(0.0), hs)
    @test !isempty(free)
    @test any(j -> j > first(free) && hs[j] > 0.3, eachindex(hs))
    bare = [j for j in free if j > 1 && hs[j-1] == 0.0 && parent(o[:precipitation])[j] == 0.0]
    @test !isempty(bare)
    @test all(≈(mp.ground_albedo; atol=1e-12), parent(o[:albedo_broadband])[bare])
    # The output base temperature is the forcing's daily mean.
    @test parent(o[:temperature_base])[end] ≈
          sum(parent(cf.temperature_ground)[end-23:end]) / 24
    # Snow cells never exceed the melting point; under deep snow the deepest cell, whose centre
    # is half a cell above the base, follows the base (measured: within 1.6 K).
    @test all(T -> isnan(T) || T <= CtoK + 1e-9, parent(o[:temperature])[mat .== MATERIAL_SNOW])
    deep = findall(>(0.4), hs)
    Tbot = [parent(o[:temperature])[count(==(MATERIAL_SNOW), mat[:, j]), j] for j in deep]
    @test maximum(abs.(Tbot .- parent(o[:temperature_base])[deep])) < 2.0
end

@testset "prescribed: snow on a warm base melts from below" begin
    mp = GEMB.initialize_parameters(lower_boundary=:prescribed,
        dt_divisors=GEMB.fast_divisors(3600 * 10000) ./ 10000)
    snow() = GEMB.column_state([265.0, 268.0, 270.0], [0.05, 0.05, 0.05], [200.0, 220.0, 240.0],
        zeros(3), fill(0.1, 3), zeros(3), fill(0.5, 3), fill(10.0, 3))
    cfs = ClimateForcingStep(; dt=3600.0, temperature_air=263.15, pressure_air=80000.0,
        wind_speed=2.0, longwave_downward=220.0, vapor_pressure=200.0,
        temperature_observation_height=2.0, wind_observation_height=10.0,
        precipitation_mean=500.0, temperature_air_mean=270.0, wind_speed_mean=3.0)
    # A base colder than the deepest cell draws heat from the snow; a warmer one supplies it.
    _, f_cold = GEMB.gemb_core_prescribed(snow(), 268.0, cfs, mp, true)
    _, f_warm = GEMB.gemb_core_prescribed(snow(), 276.15, cfs, mp, true)
    @test f_cold.heat_flux_basal < 0.0 < f_warm.heat_flux_basal
    # Under air near the melting point a base at +5 °C melts the snow from below: the base
    # cell reaches 273.15 K and the heat it keeps receiving melts it (measured: 9.5 kg m-2 in
    # 48 h, the surface cell still below the melting point).
    cfs_mild = ClimateForcingStep(; dt=3600.0, temperature_air=271.15, pressure_air=80000.0,
        wind_speed=2.0, longwave_downward=290.0, vapor_pressure=500.0,
        temperature_observation_height=2.0, wind_observation_height=10.0,
        precipitation_mean=500.0, temperature_air_mean=270.0, wind_speed_mean=3.0)
    s = snow(); melt = 0.0; runoff = 0.0
    for _ in 1:48
        s, f = GEMB.gemb_core_prescribed(s, 278.15, cfs_mild, mp, true)
        melt += f.melt; runoff += f.runoff
    end
    @test melt > 1.0 && runoff > 0.0
    @test s.temperature[1] < CtoK
    # Trace snow on a warm base melts out in the step; on a cold one it takes the base
    # temperature.
    trace() = GEMB.column_state([270.0], [0.001], [200.0], [0.0], [0.1], [0.0], [0.5], [1.0])
    s_t, f_t = GEMB.gemb_core_prescribed(trace(), 276.15, cfs, mp, true)
    @test isempty(s_t.dz) && f_t.runoff ≈ 0.2
    s_c, _ = GEMB.gemb_core_prescribed(trace(), 266.0, cfs, mp, true)
    @test s_c.temperature == [266.0]
end
