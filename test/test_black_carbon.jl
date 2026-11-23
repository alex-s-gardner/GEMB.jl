# Tests for the per-layer equivalent impurity loading `black_carbon`: how the grid primitives and
# the accumulation carry it.

_bc_cols(; n=4, black_carbon=[0.0, 0.5, 2.0, 7.0], density=[200.0, 300.0, 450.0, 600.0],
    water=[0.0, 0.4, 1.0, 0.0]) = GEMB.column_state(fill(260.0, n), fill(0.05, n),
    copy(density), copy(water), fill(0.3, n), zeros(n), fill(0.5, n), collect(1.0:n),
    copy(black_carbon))

_impurity_mass(cols) = sum(cols.black_carbon .* cols.dz .* cols.density)

@testset "merge conserves matrix-mass-weighted black_carbon; split copies" begin
    mp = GEMB.ModelParameters()
    for merge! in (GEMB.merge_pair!, GEMB.merge_snow_pair!, GEMB.absorb_snow!)
        cols = _bc_cols()
        i, j = merge! === GEMB.absorb_snow! ? (1, 2) : (2, 3)
        Mi, Mj = cols.dz[i] * cols.density[i], cols.dz[j] * cols.density[j]
        expected = (cols.black_carbon[i] * Mi + cols.black_carbon[j] * Mj) / (Mi + Mj)
        impurity = cols.black_carbon[i] * Mi + cols.black_carbon[j] * Mj
        merge!(cols, i, j, Mi, Mj, mp)
        @test cols.black_carbon[j] ≈ expected rtol = 1e-15
        # Pore water carries no impurity: the weights are matrix masses.
        @test cols.black_carbon[j] * (Mi + Mj) ≈ impurity rtol = 1e-15
    end

    # A column at a uniform loading stays at exactly that loading through any merge.
    cols = _bc_cols(black_carbon=fill(0.3, 4))
    GEMB.merge_pair!(cols, 1, 2, 10.0, 15.0, mp)
    @test cols.black_carbon[2] === 0.3

    cols = _bc_cols()
    before = copy(cols.black_carbon)
    impurity = _impurity_mass(cols)
    GEMB.split_cell!(cols, 3)
    @test cols.black_carbon == [before[1], before[2], before[3], before[3], before[4]]
    @test _impurity_mass(cols) ≈ impurity rtol = 1e-15

    @test GEMB.mix_black_carbon(1.0, 3.0, 5.0, 1.0) == 2.0
    @test GEMB.mix_black_carbon(1.0, 0.0, 5.0, 0.0) == 1.0
end

@testset "manage_layer_thickness conserves impurity mass" begin
    mp = GEMB.ModelParameters()
    n = 12
    dz = [0.004; fill(0.05, 4); 0.2; fill(0.05, n - 6)]
    T = fill(260.0, n); ρ = [150.0; fill(350.0, n - 1)]
    bc = collect(range(0.0, 3.0; length=n))
    impurity = sum(bc .* dz .* ρ)
    out = GEMB.manage_layer_thickness(T, copy(dz), copy(ρ), zeros(n), fill(0.3, n), zeros(n),
        fill(0.5, n), zeros(n), copy(bc), mp, true; dt=10800.0, n_target=n)
    @test sum(out[9] .* out[2] .* out[3]) ≈ impurity rtol = 1e-13
end

@testset "accumulation deposits black_carbon_snow; rain leaves it" begin
    mp = GEMB.ModelParameters(new_snow_method=:Constant150, column_dzmin=0.025)
    n = 5
    column() = (fill(260.0, n), fill(0.1, n), fill(400.0, n), zeros(n), fill(0.5, n),
        zeros(n), fill(0.5, n), zeros(n), fill(3.0, n))
    step(P, T) = GEMB.ClimateForcingStep(; dt=10800.0, temperature_air=T, precipitation=P,
        black_carbon_snow=0.1, temperature_air_mean=260.0, wind_speed_mean=3.0,
        precipitation_mean=300.0, temperature_observation_height=2.0,
        wind_observation_height=10.0)

    # Snow deep enough for its own cell takes the background loading; the old cells keep theirs.
    out = GEMB.calculate_accumulation(column()..., step(10.0, 260.0), mp, true)
    @test length(out[2]) == n + 1
    @test out[9][1] == 0.1
    @test all(==(3.0), out[9][2:end])

    # Light snow merged into the top cell mixes in on matrix mass.
    out = GEMB.calculate_accumulation(column()..., step(0.5, 260.0), mp, true)
    @test length(out[2]) == n
    M = 0.1 * 400.0
    @test out[9][1] ≈ (3.0 * M + 0.1 * 0.5) / (M + 0.5) rtol = 1e-14

    # Rain does not change it.
    out = GEMB.calculate_accumulation(column()..., step(2.0, 278.0), mp, true)
    @test out[9] == fill(3.0, n)
end

@testset "albedo_observed forcing layer" begin
    n = 48
    time = collect(DateTime(2001, 6, 1):Hour(1):(DateTime(2001, 6, 1) + Hour(n - 1)))
    base = (time, fill(265.0, n), fill(80000.0, n), zeros(n), fill(3.0, n), fill(200.0, n),
        fill(250.0, n), fill(300.0, n))
    kw = (temperature_air_mean=265.0, wind_speed_mean=3.0, precipitation_mean=300.0,
        temperature_observation_height=2.0, wind_observation_height=10.0)

    # Absent: a Fill-backed NaN layer, and NaN in every ClimateForcingStep form.
    cf = initialize_forcing(base...; kw...)
    @test :albedo_observed in keys(cf)
    @test all(isnan, parent(cf.albedo_observed))
    @test isnan(GEMB.ClimateForcingStep().albedo_observed)
    @test isnan(GEMB.ClimateForcingStep(ntuple(_ -> 0.0, 23)...).albedo_observed)
    @test isnan(GEMB.ClimateForcingStep(ntuple(_ -> 0.0, 19)...).albedo_observed)

    # Present: carried through, NaN where there is no observation.
    obs = fill(NaN, n); obs[12] = 0.62
    cf = initialize_forcing(base...; kw..., albedo_observed=obs)
    @test parent(cf.albedo_observed)[12] == 0.62
    @test isnan(parent(cf.albedo_observed)[1])

    # Range-checked where finite.
    bad = copy(obs); bad[3] = 1.0
    @test_throws "albedo_observed values must be in (0, 1)" initialize_forcing(base...; kw...,
        albedo_observed=bad)
    bad[3] = 0.0
    @test_throws "albedo_observed values must be in (0, 1)" initialize_forcing(base...; kw...,
        albedo_observed=bad)

    # A DimStack layer of that name is read.
    tdim = GEMB.Ti(time)
    names = (:temperature_air, :pressure_air, :precipitation, :wind_speed,
        :shortwave_downward, :longwave_downward, :vapor_pressure)
    layers = NamedTuple{names}(Tuple(DimArray(collect(Float64, v), (tdim,)) for v in base[2:end]))
    meta = Dict{String,Any}(String(k) => v for (k, v) in pairs(kw))
    stack = DimStack(merge(layers, (; albedo_observed=DimArray(obs, (tdim,)))); metadata=meta)
    @test isequal(parent(initialize_forcing(stack).albedo_observed), obs)
end

@testset "meltwater scavenging follows CLM5" begin
    # Three cells; the top one melted out, the second lost half its matrix.
    c = [0.10, 0.20, 0.05]; M_entry = [2.0, 4.0, 10.0]
    M = [0.0, 2.0, 10.0]; water = [0.0, 0.5, 1.0]
    initial() = c .* M_entry
    load = sum(initial())

    # k = 0: no impurity moves with water; the melted-out cell passes its load down.
    bc = copy(c); m = initial()
    GEMB._scavenge_impurities!(bc, m, M, water, [0.0, 2.0, 3.0, 1.0], zeros(3), 0.0)
    @test bc[2] ≈ (0.10 * 2.0 + 0.20 * 4.0) / 2.0
    @test bc[3] == 0.05
    @test sum(m) ≈ load

    # k > 0: cell i sends k · q_out · m / (M + w) down; the deepest cell's outflow leaves.
    k = 0.2; flux_dn = [0.0, 2.0, 3.0, 1.0]
    bc = copy(c); m = initial()
    GEMB._scavenge_impurities!(bc, m, M, water, flux_dn, zeros(3), k)
    m2 = 0.10 * 2.0 + 0.20 * 4.0
    out2 = k * 3.0 * m2 / 2.5
    m3 = 0.05 * 10.0 + out2
    out3 = k * 1.0 * m3 / 11.0
    @test m[2] ≈ m2 - out2 && bc[2] ≈ (m2 - out2) / 2.0
    @test m[3] ≈ m3 - out3 && bc[3] ≈ (m3 - out3) / 10.0
    @test sum(m) ≈ load - out3

    # Ponded basal water keeps its impurity in the column.
    bc = copy(c); m = initial()
    GEMB._scavenge_impurities!(bc, m, M, water, flux_dn, zeros(3), k; basal_leaves=false)
    @test sum(m) ≈ load
    @test bc[3] ≈ m3 / 10.0

    # Runoff carries its share out of the column; the rest continues down.
    bc = copy(c); m = initial()
    GEMB._scavenge_impurities!(bc, m, M, water, [0.0, 2.0, 1.0, 0.0], [0.0, 2.0, 0.0], k)
    out2 = k * 3.0 * m2 / 2.5
    @test m[3] ≈ 0.05 * 10.0 + out2 / 3.0
    @test sum(m) ≈ load - out2 * 2.0 / 3.0

    # Outflow is capped at the cell's load (k above 1 is allowed, as in CLM5).
    bc = copy(c); m = initial()
    GEMB._scavenge_impurities!(bc, m, M, water, [0.0, 2.0, 30.0, 0.0], zeros(3), 5.0)
    @test m[2] == 0.0 && bc[2] == 0.0
    @test m[3] ≈ load

    # No water leaving: unchanged.
    bc = copy(c); m = c .* M_entry
    GEMB._scavenge_impurities!(bc, m, M_entry, zeros(3), zeros(4), zeros(3), 0.2)
    @test bc ≈ c

    @test_throws "impurity_scavenging_efficiency must be ≥ 0" GEMB.initialize_parameters(;
        impurity_scavenging_efficiency = -0.1)
end
