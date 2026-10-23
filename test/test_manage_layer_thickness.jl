# Tests for manage_layer_thickness and the grid primitives it is built from.
#
# `manage_layer_thickness` now restores the column to exactly `n_target` cells (default: the incoming
# length) using merge/split operations that conserve mass and energy exactly, so it reports
# `mass_added == 0`. Total depth is pinned separately by `trim_bottom!`, called later in the
# timestep from `gemb_core`. The tests below therefore assert the fixed cell count, exact
# conservation, and the per-operation invariants rather than the old variable output lengths.
#
# `verbose=true` throughout, so the internal mass/energy conservation checks run on every
# call and a violation surfaces as an error rather than as a silently wrong number.

# Common setup helper for manage_layer_thickness tests
function _make_layer_inputs(; n=10)
    temperature = 260.0 * ones(n)
    dz = 0.1 * ones(n)
    density = 400.0 * ones(n)
    water = zeros(n)
    grain_radius = 0.5 * ones(n)
    grain_dendricity = 0.5 * ones(n)
    grain_sphericity = 0.5 * ones(n)
    age = zeros(n)
    return temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, age
end

function _layer_mp(; column_dzmin=0.05, column_dzmax=0.10, column_depth_max=1.0,
                     column_ztop=2.0, column_zy=1.1)
    return GEMB.ModelParameters(
        column_dzmin=column_dzmin,
        column_dzmax=column_dzmax,
        column_depth_max=column_depth_max,
        column_ztop=column_ztop,
        column_zy=column_zy,
    )
end

# Column mass and enthalpy, for conservation assertions. Delegates to the same helpers the
# model's own conservation checks use, so these assertions cannot drift from the definition
# of column mass/enthalpy that `grid_ops.jl` enforces. `column_mass` needs no
# `ModelParameters`; enthalpy does, since it honours `heat_capacity_method`.
_col_mass(dz, density, water) = GEMB.column_mass((; dz, density, water))
_col_energy(temperature, dz, density, water, mp=_layer_mp()) =
    last(GEMB.column_mass_energy((; dz, density, water, temperature), mp))

@testset "No action needed" begin
    temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, age = _make_layer_inputs()
    mp = _layer_mp()
    dz_in = copy(dz)
    t_in = copy(temperature)

    (t_out, dz_out, d_out, _, _, _, _, _, e_add) =
        GEMB.manage_layer_thickness(temperature, dz, density, water, grain_radius,
            grain_dendricity, grain_sphericity, age,
            mp, true; dt=10800.0)

    # Every cell is inside its band and the count already matches, so nothing moves.
    @test length(dz_out) == 10
    @test dz_out == dz_in
    @test t_out == t_in
    @test e_add == 0.0
end

@testset "Merge small layer" begin
    temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, age = _make_layer_inputs()

    dz[1] = 0.01   # below dzmin -> merges into cell 2
    dz[2] = 0.05
    mp = _layer_mp()

    m_before = _col_mass(dz, density, water)
    e_before = _col_energy(temperature, dz, density, water)
    expected_dz = dz[1] + dz[2]
    expected_d = (dz[1] * density[1] + dz[2] * density[2]) / expected_dz

    (t_out, dz_out, d_out, w_out, _, _, _, _, e_add) =
        GEMB.manage_layer_thickness(temperature, dz, density, water, grain_radius,
            grain_dendricity, grain_sphericity, age,
            mp, true; dt=10800.0)

    # Count control restores the length the merge took away (by splitting a deep cell).
    @test length(dz_out) == 10
    # The merged pair is cell 1 of the result: combined thickness, mass-weighted density.
    @test dz_out[1] ≈ expected_dz atol = 1e-10
    @test d_out[1] ≈ expected_d atol = 1e-10
    # Merge and split are both exactly conservative; the Dirichlet BC is a no-op here
    # because the bottom cell's temperature is unchanged by a deep split.
    @test _col_mass(dz_out, d_out, w_out) ≈ m_before atol = 1e-9
    @test _col_energy(t_out, dz_out, d_out, w_out) ≈ e_before + e_add atol = 1e-6
end

@testset "Split large layer" begin
    temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, age = _make_layer_inputs()

    dz[1] = 0.2  # > dzmax
    mp = _layer_mp()

    m_before = _col_mass(dz, density, water)

    (t_out, dz_out, d_out, w_out, _, _, _, _, _) =
        GEMB.manage_layer_thickness(temperature, dz, density, water, grain_radius,
            grain_dendricity, grain_sphericity, age,
            mp, true; dt=10800.0)

    # The split adds a cell; count control reclaims the slot with a deep merge.
    @test length(dz_out) == 10
    # Splitting halves thickness and duplicates the intensive properties.
    @test dz_out[1] ≈ 0.1 atol = 1e-10
    @test dz_out[2] ≈ 0.1 atol = 1e-10
    @test t_out[1] == t_out[2]
    @test d_out[1] == d_out[2]
    @test _col_mass(dz_out, d_out, w_out) ≈ m_before atol = 1e-9
end

@testset "Fixed cell count is restored from either direction" begin
    mp = _layer_mp()

    # Start short: the count controller must split to reach the target.
    temperature, dz, density, water, gr, gd, gs, ag = _make_layer_inputs(n=8)
    (_, dz_out, _, _, _, _, _, _, _) =
        GEMB.manage_layer_thickness(temperature, dz, density, water, gr, gd, gs, ag,
            mp, true; dt=10800.0, n_target=10)
    @test length(dz_out) == 10
    @test sum(dz_out) ≈ 0.8 atol = 1e-10   # splitting preserves total depth

    # Start long: the count controller must merge to reach the target.
    temperature, dz, density, water, gr, gd, gs, ag = _make_layer_inputs(n=13)
    (_, dz_out, _, _, _, _, _, _, _) =
        GEMB.manage_layer_thickness(temperature, dz, density, water, gr, gd, gs, ag,
            mp, true; dt=10800.0, n_target=10)
    @test length(dz_out) == 10
    @test sum(dz_out) ≈ 1.3 atol = 1e-10   # merging preserves total depth
end

@testset "Bottom temperature boundary condition" begin
    temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, age = _make_layer_inputs()

    t_orig_bottom = temperature[end]
    dz[1] = 0.2  # triggers split
    mp = _layer_mp()

    (t_out, _, _, _, _, _, _, _, _) =
        GEMB.manage_layer_thickness(temperature, dz, density, water, grain_radius,
            grain_dendricity, grain_sphericity, age,
            mp, true; dt=10800.0)

    @test t_out[end] == t_orig_bottom
end

@testset "Bottom merge logic" begin
    temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, age = _make_layer_inputs()

    dz[end] = 0.01      # below dzmin; the bottom cell has no cell below to merge into
    dz[end-1] = 0.05
    m_before = _col_mass(dz, density, water)
    mp = _layer_mp()

    (_, dz_out, d_out, w_out, _, _, _, _, _) =
        GEMB.manage_layer_thickness(temperature, dz, density, water, grain_radius,
            grain_dendricity, grain_sphericity, age,
            mp, true; dt=10800.0)

    # The bottom cell merges *upward* into the deepest surviving cell above it, and count
    # control then restores the length.
    @test length(dz_out) == 10
    @test _col_mass(dz_out, d_out, w_out) ≈ m_before atol = 1e-9
end

@testset "trim_bottom! (mass controller)" begin
    mp = _layer_mp()

    # Depth pinning is signed: shrink the bottom cell when the column is too deep
    # (accumulation), grow it when too shallow (basal accretion under ablation).

    @testset "shrink — accumulation regime" begin
        n = 5
        cols = GEMB.column_state(260.0 * ones(n), 0.2 * ones(n), 400.0 * ones(n),
            zeros(n), 0.5 * ones(n), 0.5 * ones(n), 0.5 * ones(n), zeros(n))
        z_target = 0.95           # column is 1.0 m -> trim 0.05 m off the bottom cell

        mass_added, e_added = GEMB.trim_bottom!(cols, z_target, mp)

        @test sum(cols.dz) ≈ z_target atol = 1e-12
        @test cols.dz[n] ≈ 0.15 atol = 1e-12
        # Mass leaves through the base, so the reported delta is negative and matches the
        # thickness change exactly.
        @test mass_added ≈ -(0.05 * 400.0) atol = 1e-12
        @test e_added ≈ -(0.05 * 400.0 * GEMB.specific_enthalpy(mp, 260.0)) atol = 1e-6
    end

    @testset "grow — ablation regime (basal accretion)" begin
        n = 5
        cols = GEMB.column_state(260.0 * ones(n), 0.2 * ones(n), 400.0 * ones(n),
            zeros(n), 0.5 * ones(n), 0.5 * ones(n), 0.5 * ones(n), zeros(n))
        z_target = 1.05           # column is 1.0 m -> accrete 0.05 m onto the bottom cell

        mass_added, e_added = GEMB.trim_bottom!(cols, z_target, mp)

        @test sum(cols.dz) ≈ z_target atol = 1e-12
        @test cols.dz[n] ≈ 0.25 atol = 1e-12
        # Mass enters through the base: same magnitude, opposite sign.
        @test mass_added ≈ +(0.05 * 400.0) atol = 1e-12
        @test e_added ≈ +(0.05 * 400.0 * GEMB.specific_enthalpy(mp, 260.0)) atol = 1e-6
    end

    @testset "pore water is carried proportionally" begin
        n = 3
        cols = GEMB.column_state(273.15 * ones(n), [0.2, 0.2, 0.2], 400.0 * ones(n),
            [0.0, 0.0, 10.0], 0.5 * ones(n), 0.5 * ones(n), 0.5 * ones(n), zeros(n))
        # Remove half of the bottom cell -> half of its water leaves with it.
        mass_added, _ = GEMB.trim_bottom!(cols, 0.5, mp)

        @test cols.dz[n] ≈ 0.1 atol = 1e-12
        @test cols.water[n] ≈ 5.0 atol = 1e-12
        @test mass_added ≈ -(0.1 * 400.0 + 5.0) atol = 1e-12
    end

    @testset "already on target is a no-op" begin
        n = 4
        dz = 0.25 * ones(n)
        cols = GEMB.column_state(260.0 * ones(n), dz, 400.0 * ones(n),
            zeros(n), 0.5 * ones(n), 0.5 * ones(n), 0.5 * ones(n), zeros(n))
        mass_added, e_added = GEMB.trim_bottom!(cols, 1.0, mp)
        @test mass_added == 0.0
        @test e_added == 0.0
        @test dz == 0.25 * ones(n)
    end

    # An adjustment thicker than the deepest cell is not an error: a step whose accumulation
    # exceeds the bottom cell is routine at high accumulation, and the surplus comes from the cells
    # above in turn. Conservation is asserted against `column_mass_energy` rather than a re-derived
    # formula, so the bookkeeping — including the reservoir temperature a promoted cell inherits —
    # is checked rather than restated.
    @testset "adjustment spanning more than the bottom cell" begin
        _spanning_cols(temperature) = GEMB.column_state(collect(Float64, temperature),
            [0.4, 0.4, 0.2], 400.0 * ones(3), zeros(3),
            0.5 * ones(3), 0.5 * ones(3), 0.5 * ones(3), zeros(3))

        @testset "consumes the bottom cell and trims the next" begin
            cols = _spanning_cols([260.0, 260.0, 260.0])
            m0, e0 = GEMB.column_mass_energy(cols, mp)
            # 1.0 m column, 0.2 m bottom cell, 0.3 m to remove.
            mass_added, e_added = GEMB.trim_bottom!(cols, 0.7, mp)
            m1, e1 = GEMB.column_mass_energy(cols, mp)

            @test length(cols.dz) == 2
            @test sum(cols.dz) ≈ 0.7 atol = 1e-12
            @test cols.dz ≈ [0.4, 0.3] atol = 1e-12
            @test mass_added ≈ -(0.3 * 400.0) atol = 1e-9
            @test m1 - m0 ≈ mass_added atol = 1e-9
            @test e1 - e0 ≈ e_added atol = 1e-3
        end

        @testset "a promoted cell inherits the reservoir temperature" begin
            cols = _spanning_cols([250.0, 260.0, 270.0])
            m0, e0 = GEMB.column_mass_energy(cols, mp)
            mass_added, e_added = GEMB.trim_bottom!(cols, 0.7, mp)
            m1, e1 = GEMB.column_mass_energy(cols, mp)

            # Cell 3 was the Dirichlet reservoir; cell 2 takes its place and its temperature.
            @test cols.temperature ≈ [250.0, 270.0] atol = 1e-12
            @test mass_added ≈ -(0.3 * 400.0) atol = 1e-9
            @test m1 - m0 ≈ mass_added atol = 1e-9
            # Warming the promoted cell to the reservoir is part of `e_added`, so the column's own
            # enthalpy change still matches it.
            @test e1 - e0 ≈ e_added atol = 1e-3
        end

        @testset "consumes several cells" begin
            n = 6
            cols = GEMB.column_state(260.0 * ones(n), 0.1 * ones(n), 400.0 * ones(n),
                zeros(n), 0.5 * ones(n), 0.5 * ones(n), 0.5 * ones(n), zeros(n))
            m0, e0 = GEMB.column_mass_energy(cols, mp)
            # 0.6 m column -> 0.25 m: three whole cells plus half of the fourth.
            mass_added, e_added = GEMB.trim_bottom!(cols, 0.25, mp)
            m1, e1 = GEMB.column_mass_energy(cols, mp)

            @test length(cols.dz) == 3
            @test sum(cols.dz) ≈ 0.25 atol = 1e-12
            @test mass_added ≈ -(0.35 * 400.0) atol = 1e-9
            @test m1 - m0 ≈ mass_added atol = 1e-9
            @test e1 - e0 ≈ e_added atol = 1e-3
        end

        @testset "pore water of a consumed cell leaves with it" begin
            n = 3
            cols = GEMB.column_state(273.15 * ones(n), [0.4, 0.4, 0.2], 400.0 * ones(n),
                [0.0, 8.0, 10.0], 0.5 * ones(n), 0.5 * ones(n), 0.5 * ones(n), zeros(n))
            m0, _ = GEMB.column_mass_energy(cols, mp)
            mass_added, _ = GEMB.trim_bottom!(cols, 0.7, mp)
            m1, _ = GEMB.column_mass_energy(cols, mp)

            # All 10 of cell 3's water, plus a quarter of cell 2's as it goes 0.4 -> 0.3 m.
            @test cols.water ≈ [0.0, 6.0] atol = 1e-12
            @test mass_added ≈ -(0.3 * 400.0 + 10.0 + 2.0) atol = 1e-9
            @test m1 - m0 ≈ mass_added atol = 1e-9
        end
    end

    @testset "an adjustment exceeding the whole column is an error" begin
        n = 3
        cols = GEMB.column_state(260.0 * ones(n), [0.4, 0.4, 0.2], 400.0 * ones(n),
            zeros(n), 0.5 * ones(n), 0.5 * ones(n), 0.5 * ones(n), zeros(n))
        # Needs to remove more than the 1.0 m the column holds.
        @test_throws "exceeds the whole column" GEMB.trim_bottom!(cols, -0.2, mp)
    end
end

# The `ColumnWorkspace` reuse contract, pinned the same way `test_calculate_temperature.jl` pins
# the thermal one: a buffer arrives holding a previous call's values and is usually longer than the
# column, so the result must not depend on either. A regression here is a wrong answer, not a slow
# one — the surplus tail is only inert because every consumer bounds its reads by the cell count.
@testset "ColumnWorkspace reuse" begin
    mp = _layer_mp()

    run_with(ws) = GEMB.manage_layer_thickness(_make_layer_inputs(n = 10)..., mp, true; dt=10800.0, workspace = ws)

    @testset "a dirty, oversized workspace gives the fresh-workspace answer" begin
        fresh = run_with(ColumnWorkspace())

        dirty = ColumnWorkspace()
        GEMB._resize_workspace!(dirty, 400)
        fill!(dirty.dzmin, -99.0)
        fill!(dirty.dzmax, -99.0)
        fill!(dirty.mass, -99.0)
        fill!(dirty.delete_cell, true)      # every cell marked deleted, from a previous column
        append!(dirty.to_delete, 1:7)
        append!(dirty.to_split, 3:9)

        @test run_with(dirty) == fresh
    end

    @testset "one workspace across a long column then a short one" begin
        ws = ColumnWorkspace()
        long_fresh = GEMB.manage_layer_thickness(_make_layer_inputs(n = 30)..., mp, true; dt=10800.0,
                                           workspace = ColumnWorkspace())
        @test GEMB.manage_layer_thickness(_make_layer_inputs(n = 30)..., mp, true; dt=10800.0,
                                     workspace = ws) == long_fresh
        # The buffers are now 30-long; a 10-cell column must not read the stale tail.
        @test run_with(ws) == run_with(ColumnWorkspace())
    end

    @testset "buffers grow and never shrink" begin
        ws = ColumnWorkspace()
        GEMB._resize_workspace!(ws, 40)
        GEMB._resize_workspace!(ws, 10)
        @test all(length(b) == 40 for b in GEMB._workspace_buffers(ws))
    end
end

# A distinct thin surface snow layer (see `manage_layer_thickness` and `snow_layers_distinct`):
# kept while distinct and holding at least the minimum mass (3-hour step here), merged with the grain state
# mixed by mass once it is not, absorbed at the firn or ice density over firn or ice.
# Cells 2… are 0.05 m of old snow (400 kg m-3, r 0.6 mm, sphericity 0.9, 200 days) over ice
# from cell 6; `under` sets cells 2–5 to firn or ice instead.
function _surface_layer_column(; dz1=0.008, ρ1=200.0, r1=0.05, d1=1.0, s1=0.5, age1=1.0,
    T1=265.0, n=12, under=400.0)
    return ([T1; fill(260.0, n - 1)], [dz1; fill(0.05, n - 1)],
        [ρ1; fill(under, 4); fill(917.0, n - 5)], zeros(n),
        [r1; fill(0.6, 4); fill(GEMB.GRAIN_RADIUS_ICE, n - 5)], [d1; zeros(n - 1)],
        [s1; fill(0.9, n - 1)], [age1; fill(200.0, n - 1)])
end
_snow_mp(; kw...) = GEMB.ModelParameters(; column_dzmin=0.025, column_dzmax=0.075,
    column_ztop=10.0, kw...)

@testset "Mixing grain states (Crocus aggregation)" begin
    # Two dendritic states: dendricity and sphericity mass-weighted, the Brun radius of the mix.
    d, s, r = GEMB.mix_grains(0.6, 0.4, 0.09, 3.0, 1.0, 0.5, 0.05, 1.0)
    @test d ≈ (3 * 0.6 + 1.0) / 4 rtol = 1e-15
    @test s ≈ (3 * 0.4 + 0.5) / 4 rtol = 1e-15
    @test r == GEMB.dendritic_grain_radius(d, s)
    # Otherwise non-dendritic, conserving the mass-weighted Gardner–Sharp α_S of eq. 7.
    for (r1, M1, r2, M2) in ((1.0, 1.0, 0.05, 1.0), (0.6, 20.0, 0.05, 1.5), (2.0, 3.0, 0.3, 7.0))
        d, s, r = GEMB.mix_grains(0.0, 0.9, r1, M1, 1.0, 0.5, r2, M2)
        @test d == 0.0
        @test s ≈ (M1 * 0.9 + M2 * 0.5) / (M1 + M2) rtol = 1e-15
        αS(x) = big"1.48" - (3 / (big"0.091" * big(x)))^big"-0.07"
        @test GEMB._albedo_gardner(r, 0.0, 0.0, 0.0) ≈
              Float64((M1 * αS(r1) + M2 * αS(r2)) / (M1 + M2)) rtol = 1e-14
    end
    # Weight lost to rounding: the first state bit for bit.
    @test GEMB.mix_grains(0.0, 0.9, 0.6, 20.0, 1.0, 0.5, 0.05, 5e-324) === (0.0, 0.9, 0.6)
    @test GEMB.mix_grains(0.7, 0.4, 0.081, 2.0, 1.0, 0.5, 0.05, 5e-324) === (0.7, 0.4, 0.081)
end

@testset "A distinct thin surface layer is kept, then merges once it is not" begin
    mp = _snow_mp(heat_capacity_method=:constant)
    col = _surface_layer_column()
    before = deepcopy(col)
    D = GEMB.snow_layer_distance(1.0, 0.5, 0.05, false, 1.0, 200.0, 0.0, 0.9, 0.6, false, 200.0,
        400.0)
    @test GEMB.snow_layers_distinct(0.008, D, mp)
    out = GEMB.manage_layer_thickness(deepcopy(col)..., mp, true; dt=10800.0)
    @test all(out[k] == before[k] for k in 1:8)

    # The layer has aged toward the snow below (r 0.55 mm, sphericity 0.88, non-dendritic):
    # it is no longer distinct, so it merges. Closed forms, with c·T under :constant.
    col = _surface_layer_column(r1=0.55, d1=0.0, s1=0.88, age1=150.0)
    D = GEMB.snow_layer_distance(0.0, 0.88, 0.55, false, 150.0, 200.0, 0.0, 0.9, 0.6, false, 200.0,
        400.0)
    @test !GEMB.snow_layers_distinct(0.008, D, mp)
    M1, M2 = 0.008 * 200.0, 0.05 * 400.0
    m0, e0 = _col_mass(col[2], col[3], col[4]), _col_energy(col[1], col[2], col[3], col[4], mp)
    out = GEMB.manage_layer_thickness(deepcopy(col)..., mp, true; dt=10800.0, n_target=11)
    @test length(out[2]) == 11
    @test out[2][1] ≈ 0.058 rtol = 1e-14
    @test out[3][1] ≈ (M1 + M2) / 0.058 rtol = 1e-14
    @test out[1][1] ≈ (M1 * 265.0 + M2 * 260.0) / (M1 + M2) rtol = 1e-14
    @test out[8][1] ≈ (M1 * 150.0 + M2 * 200.0) / (M1 + M2) rtol = 1e-14
    @test out[7][1] ≈ (M1 * 0.88 + M2 * 0.9) / (M1 + M2) rtol = 1e-14
    @test out[6][1] == 0.0
    @test out[5][1] ≈ Float64(((M1 * big(0.55)^big"0.07" + M2 * big(0.6)^big"0.07") /
                              (M1 + M2))^(1 / big"0.07")) rtol = 1e-13
    @test _col_mass(out[2], out[3], out[4]) ≈ m0 rtol = 1e-14
    @test _col_energy(out[1], out[2], out[3], out[4], mp) ≈ e0 rtol = 1e-14
    @test all(out[k][2:end] == col[k][3:end] for k in 1:8)

    # Under :CuffeyPaterson the merge conserves enthalpy too.
    mp_cp = _snow_mp()
    e0 = _col_energy(col[1], col[2], col[3], col[4], mp_cp)
    out = GEMB.manage_layer_thickness(deepcopy(col)..., mp_cp, true; dt=10800.0, n_target=11)
    @test _col_energy(out[1], out[2], out[3], out[4], mp_cp) ≈ e0 rtol = 1e-14
end

@testset "A surface layer below the minimum mass merges however distinct" begin
    mp = _snow_mp()
    z = 0.999 * GEMB.snow_layer_min_mass(10800.0) / 200.0     # cell 1 is 200 kg m-3
    # Into snow: by mass, the fresh grains count by their mass.
    col = _surface_layer_column(dz1=z)
    m0, e0 = _col_mass(col[2], col[3], col[4]), _col_energy(col[1], col[2], col[3], col[4], mp)
    out = GEMB.manage_layer_thickness(deepcopy(col)..., mp, true; dt=10800.0, n_target=11)
    @test out[2][1] ≈ z + 0.05 rtol = 1e-14
    @test out[6][1] == 0.0
    @test out[5][1] ≈ Float64(((z * 200 * big(0.05)^big"0.07" + 20 * big(0.6)^big"0.07") /
                              (z * 200 + 20))^(1 / big"0.07")) rtol = 1e-13
    @test _col_mass(out[2], out[3], out[4]) ≈ m0 rtol = 1e-14
    @test _col_energy(out[1], out[2], out[3], out[4], mp) ≈ e0 rtol = 1e-14

    # Over firn or ice: kept at any thickness above the minimum, absorbed below it at the firn
    # or ice density, whose albedo is therefore unchanged.
    for ρ_under in (870.0, 917.0)
        col = _surface_layer_column(dz1=1.5 * GEMB.snow_layer_min_mass(10800.0) / 200.0, under=ρ_under,
            r1=0.6, d1=0.0, s1=0.9, age1=200.0)
        out = GEMB.manage_layer_thickness(deepcopy(col)..., mp, true; dt=10800.0)
        @test all(out[k] == col[k] for k in 1:8)
        col = _surface_layer_column(dz1=z, under=ρ_under)
        m0, e0 = _col_mass(col[2], col[3], col[4]), _col_energy(col[1], col[2], col[3], col[4], mp)
        out = GEMB.manage_layer_thickness(deepcopy(col)..., mp, true; dt=10800.0, n_target=11)
        @test out[3][1] === ρ_under
        @test out[2][1] ≈ 0.05 + z * 200.0 / ρ_under rtol = 1e-14
        @test (out[5][1], out[6][1], out[7][1]) == (col[5][2], col[6][2], col[7][2])
        @test _col_mass(out[2], out[3], out[4]) ≈ m0 rtol = 1e-14
        @test _col_energy(out[1], out[2], out[3], out[4], mp) ≈ e0 rtol = 1e-14
        cfs = GEMB.ClimateForcingStep(; dt=3600.0, solar_zenith_angle=60.0)
        @test GEMB.calculate_albedo(out[2], out[3], out[5], cfs, mp)[1] ===
              GEMB.firn_ice_albedo(ρ_under, mp)
    end
end
