using Test
using GEMB
using Dates

# The `age` profile variable: mass-weighted mean age, in decimal days, of all mass in a cell
# (ice/firn matrix plus pore water), measured from column initialization. There is no MATLAB
# reference — MATLAB GEMB has no age variable — so these tests pin it against its own closed
# forms (the clock, the mass-weighted merge, the dilution formula) and against structural
# invariants that any bookkeeping error in the percolation moment algebra would break.

# A forcing helper mirroring `test_gemb_spinup.jl`: daily steps, tunable warmth, sunshine and
# precipitation so the same builder serves the dry-clock case and the melt-heavy case. Each
# of those three accepts a scalar (held constant) or a per-day vector (a seasonal cycle).
_series(x, n) = x isa AbstractVector ? collect(float.(x)) : fill(float(x), n)

function _age_forcing(; n_days=365, temperature=250.0, precipitation=0.0,
    shortwave=50.0, longwave=200.0, vapor_pressure=100.0,
    start_time=DateTime(2020, 1, 1))
    time = start_time .+ Day.(0:n_days-1)
    return initialize_forcing(
        time,
        _series(temperature, n_days),     # temperature_air
        fill(85000.0, n_days),            # pressure_air
        _series(precipitation, n_days),   # precipitation
        fill(5.0, n_days),                # wind_speed
        _series(shortwave, n_days),       # shortwave_downward
        _series(longwave, n_days),        # longwave_downward
        _series(vapor_pressure, n_days),  # vapor_pressure
        temperature_observation_height=2.0,
        wind_observation_height=10.0,
    )
end

# A `ClimateForcingStep` for the direct `calculate_accumulation` calls. Defined here rather
# than reused from `test_calculate_accumulation.jl` so this file runs standalone, and built
# with the keyword constructor so a field insertion in the struct cannot silently
# mis-assign values (every field is a `Float64`, so a positional slip is not a type error).
_age_cfs(; precipitation=0.0, temperature_air=260.0, wind_speed=5.0) =
    GEMB.ClimateForcingStep(; dt=86400.0, temperature_air, precipitation, wind_speed,
        pressure_air=80000.0, shortwave_downward=100.0, longwave_downward=250.0,
        vapor_pressure=300.0, temperature_air_mean=260.0, wind_speed_mean=5.0,
        precipitation_mean=200.0, temperature_observation_height=2.0,
        wind_observation_height=10.0, cloud_fraction=0.1)

_age_cols(; n=3, temperature=fill(260.0, n), dz=fill(0.1, n), density=fill(400.0, n),
    water=zeros(n), age=zeros(n)) = GEMB.column_state(
    copy(temperature), copy(dz), copy(density), copy(water),
    fill(0.5, n), fill(0.5, n), fill(0.5, n), copy(age))

@testset "epoch is zero at initialization" begin
    mp = initialize_parameters()
    cf = _age_forcing()

    profile = initialize_profile(mp, cf)
    @test haskey(profile, :age)
    @test all(collect(profile[:age]) .== 0.0)

    # Both fidelity flags take `_uniform_ice_profile`, a separate `layers` tuple. Missing
    # `age` there is a `KeyError` at the first `gemb_core` call rather than here, so assert
    # it explicitly.
    profile_ice = initialize_profile(mp, cf; constant_density=true, constant_temperature=true)
    @test all(collect(profile_ice[:age]) .== 0.0)
end

@testset "the clock: dry column ages by exactly dt" begin
    # No precipitation and no melt, so below the surface no mass moves at all and every cell
    # must read exactly the elapsed time. `age` carries the only absolute clock in
    # `gemb_core`, so this is the calibration of every other assertion below.
    #
    # Cell 1 is excluded: vapour deposition adds age-zero mass to the surface cell even with
    # zero precipitation, so it reads slightly *younger* than elapsed. That is the intended
    # behaviour, and it is asserted as a bound rather than an equality.
    mp = initialize_parameters(output_frequency=:daily)
    n_days = 20
    cf = _age_forcing(n_days=n_days, temperature=230.0, precipitation=0.0, shortwave=0.0)
    profile = initialize_profile(mp, cf)

    out = gemb(profile, cf, mp)
    age = Array(out[:age])

    # One daily step per output, and step 0 of `gemb_core` ages the column *before* any mass
    # moves, so the first output already reads one full day.
    for k in 1:size(age, 2)
        @test all(age[2:end, k] .≈ Float64(k))
        @test 0.0 <= age[1, k] <= Float64(k) + 1e-9
    end
end

@testset "fresh snow enters at age 0" begin
    mp = initialize_parameters(new_snow_method=Symbol("150kgm2"), column_dzmin=0.05)
    @testset "new cell" begin
        # 50 kg m-2 at 150 kg m-3 is 0.33 m, well above dzmin => its own cell.
        n = 5
        age_in = fill(100.0, n)
        cfs = _age_cfs(precipitation=50.0, temperature_air=260.0)
        (_, dz_out, _, _, _, _, _, age_out, _) = GEMB.calculate_accumulation(
            fill(260.0, n), fill(0.1, n), fill(400.0, n), zeros(n),
            fill(0.5, n), fill(0.5, n), fill(0.5, n), copy(age_in), cfs, mp, false)

        @test length(dz_out) == n + 1
        # `open_slot!` duplicates the old cell 1 into the new slot, so this is the assertion
        # that the explicit `age[1] = 0.0` is present.
        @test age_out[1] == 0.0
        # Everything below is the old column, shifted down and untouched.
        @test all(age_out[2:end] .== 100.0)
    end

    @testset "merged into cell 1 dilutes toward zero" begin
        # 2 kg m-2 is 0.013 m, below dzmin => merges into cell 1.
        n = 5
        precip = 2.0
        S_surface = 0.1 * 400.0 + 3.0          # dz*density + pore water
        age_in = fill(100.0, n)
        water = zeros(n); water[1] = 3.0
        cfs = _age_cfs(precipitation=precip, temperature_air=260.0)
        (_, dz_out, _, _, _, _, _, age_out, _) = GEMB.calculate_accumulation(
            fill(260.0, n), fill(0.1, n), fill(400.0, n), water,
            fill(0.5, n), fill(0.5, n), fill(0.5, n), copy(age_in), cfs, mp, false)

        @test length(dz_out) == n
        @test age_out[1] ≈ 100.0 * S_surface / (S_surface + precip) rtol = 1e-14
        @test all(age_out[2:end] .== 100.0)
    end

    @testset "rain dilutes on the same footing" begin
        # Rain mass joins the matrix (not `water`) in `calculate_accumulation`, so it is new
        # age-zero mass in the cell exactly as snow is.
        n = 5
        precip = 10.0
        S_surface = 0.1 * 400.0
        cfs = _age_cfs(precipitation=precip, temperature_air=280.0, wind_speed=0.0)
        (_, _, _, _, _, _, _, age_out, rain) = GEMB.calculate_accumulation(
            fill(260.0, n), fill(0.1, n), fill(400.0, n), zeros(n),
            fill(0.5, n), fill(0.5, n), fill(0.5, n), fill(50.0, n), cfs, mp, false)

        @test rain == precip
        @test age_out[1] ≈ 50.0 * S_surface / (S_surface + precip) rtol = 1e-14
    end
end

@testset "dilute_age" begin
    # Removing mass is age-neutral (the mean of a proportionally reduced cell is unchanged),
    # which is what lets the sublimation and melt sites carry no code at all.
    @test GEMB.dilute_age(42.0, 100.0, 0.0) == 42.0
    @test GEMB.dilute_age(42.0, 100.0, -10.0) == 42.0
    @test GEMB.dilute_age(42.0, 100.0, 100.0) ≈ 21.0 rtol = 1e-14
    # An empty cell gaining mass reads age 0 rather than NaN.
    @test GEMB.dilute_age(42.0, 0.0, 0.0) == 42.0
    @test GEMB.dilute_age(0.0, 0.0, 5.0) == 0.0
end

@testset "merge_pair! is mass-weighted on total cell mass" begin
    mp = initialize_parameters()

    # Two cells of known mass and age. Weighting is on `dz*density + water`, matching the
    # field's definition, so the pore water must show up in the weights.
    cols = _age_cols(n=2, dz=[0.1, 0.2], density=[400.0, 500.0],
        water=[5.0, 0.0], age=[10.0, 100.0])
    S1 = 0.1 * 400.0 + 5.0
    S2 = 0.2 * 500.0 + 0.0

    GEMB.merge_pair!(cols, 1, 2, 0.1 * 400.0, 0.2 * 500.0, mp)

    @test cols.age[2] ≈ (10.0 * S1 + 100.0 * S2) / (S1 + S2) rtol = 1e-14
    # Bounded by the two inputs — an unweighted mean or a swapped weight breaks this.
    @test 10.0 < cols.age[2] < 100.0

    # Merging equal masses is the plain mean.
    cols_eq = _age_cols(n=2, dz=[0.1, 0.1], density=[400.0, 400.0], age=[0.0, 200.0])
    GEMB.merge_pair!(cols_eq, 1, 2, 40.0, 40.0, mp)
    @test cols_eq.age[2] ≈ 100.0 rtol = 1e-14

    # A massless cell merging in leaves the target's age untouched rather than NaN.
    cols_zero = _age_cols(n=2, dz=[0.0, 0.1], density=[400.0, 400.0], age=[0.0, 55.0])
    GEMB.merge_pair!(cols_zero, 1, 2, 0.0, 40.0, mp)
    @test cols_zero.age[2] ≈ 55.0 rtol = 1e-14
end

@testset "split_cell! preserves age" begin
    # Age is intensive, so both halves inherit the parent's value unchanged. `split_cell!`
    # gets this for free via `open_slot!`; the test guards against someone "fixing" it.
    cols = _age_cols(n=3, age=[7.0, 8.0, 9.0])
    GEMB.split_cell!(cols, 2)

    @test length(cols.age) == 4
    @test cols.age == [7.0, 8.0, 8.0, 9.0]
end

@testset "manage_layer_thickness keeps age bounded through merge and split" begin
    mp = GEMB.ModelParameters(column_dzmin=0.05, column_dzmax=0.10,
        column_depth_max=1.0, column_ztop=2.0, column_zy=1.1)
    n = 10

    # A monotone age profile, as a real column has. Both the merge pass (cell 1 below dzmin)
    # and the split pass (cell 3 above dzmax) fire, plus count control in the deep column.
    age_in = collect(range(0.0, 900.0, length=n))
    dz = fill(0.08, n); dz[1] = 0.01; dz[3] = 0.20

    (_, dz_out, _, _, _, _, _, age_out, _) = GEMB.manage_layer_thickness(
        fill(260.0, n), dz, fill(400.0, n), zeros(n),
        fill(0.5, n), fill(0.5, n), fill(0.5, n), copy(age_in), mp, true)

    @test length(age_out) == n
    # Every output age is a mass-weighted mean of inputs, so it lies within the input range.
    @test all(minimum(age_in) .<= age_out .<= maximum(age_in))
    # And the ordering survives: merging and splitting only ever combine neighbours.
    @test all(diff(age_out) .>= -1e-9)
end

@testset "percolation transports the source age into the refreeze cell" begin
    # A cold column with a warm, wet surface cell. The surface cell is the only mass with a
    # nonzero age, so any age appearing at depth arrived as meltwater from it — which is the
    # whole content of the "water carries source age" decision.
    n = 8
    source_age = 500.0
    age = zeros(n)
    age[1] = source_age

    temperature = fill(255.0, n)     # cold enough that inflow refreezes
    temperature[1] = 280.0           # drives melt at the surface
    water = zeros(n)
    water[1] = 20.0                  # plus a wet surface cell to push the front down
    mp = GEMB.ModelParameters(density_ice=917.0, water_irreducible_saturation=0.07)

    (_, _, _, _, _, _, _, age_out, _, _, _, freeze, depth) =
        GEMB.calculate_melt(temperature, fill(0.1, n), fill(400.0, n), water,
            fill(0.5, n), fill(0.5, n), fill(0.5, n), copy(age), 0.0, mp, false)

    @test freeze > 0.0
    @test depth > 0.0

    # No cell may fall outside [0, source_age]: every value is a mass-weighted mean of
    # {source age, 0}. A sign error or a swapped weight in the moment update breaks this.
    @test all(0.0 .<= age_out .<= source_age + 1e-9)

    # Age moved down, and the receiving cells' ages are strictly interior — they are
    # mixtures of source water and age-zero firn, not replacements.
    receiving = findall(>(0.0), age_out[2:end])
    @test !isempty(receiving)
    @test all(age_out[1 .+ receiving] .< source_age)
end

@testset "age stays bounded over a melt-heavy run" begin
    # The real guard on the percolation moment algebra: over thousands of timesteps with
    # melt, refreeze, runoff, merges, splits and basal trimming, no cell may report an age
    # outside [0, elapsed]. A bad denominator or a dropped term shows up here immediately.
    # A seasonal cycle rather than constant warmth: melt needs summer, and refreeze needs the
    # cold winter firn beneath it, so a constant-temperature column gives one or the other but
    # not both — and it is the melt-then-refreeze pairing that exercises the moment transport.
    mp = initialize_parameters(output_frequency=:daily)
    n_days = 730
    seasonal = cos.(2π .* (1:n_days) ./ 365.0 .- π)
    cf = _age_forcing(n_days=n_days,
        temperature=258.0 .+ 17.0 .* seasonal,
        shortwave=100.0 .+ 200.0 .* max.(0.0, seasonal),
        precipitation=3.0, longwave=240.0, vapor_pressure=200.0)
    profile = initialize_profile(mp, cf)

    out = gemb(profile, cf, mp)
    age = Array(out[:age])

    @test all(isfinite, age)
    @test all(age .>= 0.0)
    # `gemb_core` ages the column at the top of the timestep, so the bound is elapsed + one
    # step. Steps are daily here.
    @test maximum(age) <= n_days + 1 + 1e-9

    # Melt actually happened, so the bound above was exercised rather than trivially met.
    @test sum(Array(out[:melt])) > 0.0
    @test sum(Array(out[:refreeze])) > 0.0

    # The surface is young and the base is old: the column has a real age gradient, not a
    # uniform clock. (Strict depth monotonicity is *not* asserted — young refrozen meltwater
    # sitting beneath older firn is physically real at a melting site.)
    final = age[:, end]
    @test final[1] < final[end]
end

@testset "age accumulates across gemb_spinup cycles" begin
    # The assertion that `profile_extract.jl` carries `age` through the profile round-trip
    # between cycles. Without it, spinup either dies on a KeyError or silently resets the
    # clock each cycle and the deep column never accumulates residence time.
    mp = initialize_parameters(output_frequency=:last)
    cf = _age_forcing(n_days=365, temperature=250.0, precipitation=1.0)
    profile = initialize_profile(mp, cf)

    spun = gemb_spinup(profile, cf, mp; max_iterations=3, verbose=false)
    age = collect(spun[:age])

    @test all(isfinite, age)
    # Three 365-day cycles: the deep column, which nothing refreshes, must exceed one cycle.
    @test age[end] > 365.0
    @test age[end] <= 3 * 365.0 + 2.0
    # The surface is still young — the clock is per-cell, not a global offset.
    @test age[1] < age[end]
end

@testset "age carries CF metadata" begin
    # The generic attribute invariants (units/long_name/cell_methods present, no
    # `standard_name` on a per-cell field, `cell_methods` dropped from extracted profiles)
    # are pinned for every layer by `test_gemb_driver.jl`. Only the unit is age-specific,
    # and `cf_attributes` has no fallback, so a missing entry is a KeyError there.
    @test GEMB.cf_attributes(:age)["units"] == "d"
end
