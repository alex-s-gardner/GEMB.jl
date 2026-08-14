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
    albedo = 0.8 * ones(n)
    albedo_diffuse = 0.8 * ones(n)
    return temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, albedo, albedo_diffuse
end

function _layer_mp(; column_dzmin=0.05, column_dzmax=0.10, column_depth=1.0,
                     column_ztop=2.0, column_zy=1.1)
    return GEMB.ModelParameters(
        column_dzmin=column_dzmin,
        column_dzmax=column_dzmax,
        column_depth=column_depth,
        column_ztop=column_ztop,
        column_zy=column_zy,
    )
end

# Column mass and enthalpy, for conservation assertions. Delegates to the same helper the
# model's own conservation checks use, so these assertions cannot drift from the definition
# of column mass/enthalpy that `grid_ops.jl` enforces.
_col_mass(dz, density, water) =
    first(GEMB.column_mass_energy((; dz, density, water, temperature=zero(dz))))
_col_energy(temperature, dz, density, water) =
    last(GEMB.column_mass_energy((; dz, density, water, temperature)))

@testset "No action needed" begin
    temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, albedo, albedo_diffuse = _make_layer_inputs()
    mp = _layer_mp()
    dz_in = copy(dz)
    t_in = copy(temperature)

    (t_out, dz_out, d_out, _, _, _, _, _, _, e_add) =
        GEMB.manage_layer_thickness(temperature, dz, density, water, grain_radius,
            grain_dendricity, grain_sphericity, albedo, albedo_diffuse,
            mp, true)

    # Every cell is inside its band and the count already matches, so nothing moves.
    @test length(dz_out) == 10
    @test dz_out == dz_in
    @test t_out == t_in
    @test e_add == 0.0
end

@testset "Merge small layer" begin
    temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, albedo, albedo_diffuse = _make_layer_inputs()

    dz[1] = 0.01   # below dzmin -> merges into cell 2
    dz[2] = 0.05
    mp = _layer_mp()

    m_before = _col_mass(dz, density, water)
    e_before = _col_energy(temperature, dz, density, water)
    expected_dz = dz[1] + dz[2]
    expected_d = (dz[1] * density[1] + dz[2] * density[2]) / expected_dz

    (t_out, dz_out, d_out, w_out, _, _, _, _, _, e_add) =
        GEMB.manage_layer_thickness(temperature, dz, density, water, grain_radius,
            grain_dendricity, grain_sphericity, albedo, albedo_diffuse,
            mp, true)

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
    temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, albedo, albedo_diffuse = _make_layer_inputs()

    dz[1] = 0.2  # > dzmax
    mp = _layer_mp()

    m_before = _col_mass(dz, density, water)

    (t_out, dz_out, d_out, w_out, _, _, _, _, _, _) =
        GEMB.manage_layer_thickness(temperature, dz, density, water, grain_radius,
            grain_dendricity, grain_sphericity, albedo, albedo_diffuse,
            mp, true)

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
    temperature, dz, density, water, gr, gd, gs, al, ad = _make_layer_inputs(n=8)
    (_, dz_out, _, _, _, _, _, _, _, _) =
        GEMB.manage_layer_thickness(temperature, dz, density, water, gr, gd, gs, al, ad,
            mp, true; n_target=10)
    @test length(dz_out) == 10
    @test sum(dz_out) ≈ 0.8 atol = 1e-10   # splitting preserves total depth

    # Start long: the count controller must merge to reach the target.
    temperature, dz, density, water, gr, gd, gs, al, ad = _make_layer_inputs(n=13)
    (_, dz_out, _, _, _, _, _, _, _, _) =
        GEMB.manage_layer_thickness(temperature, dz, density, water, gr, gd, gs, al, ad,
            mp, true; n_target=10)
    @test length(dz_out) == 10
    @test sum(dz_out) ≈ 1.3 atol = 1e-10   # merging preserves total depth
end

@testset "Bottom temperature boundary condition" begin
    temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, albedo, albedo_diffuse = _make_layer_inputs()

    t_orig_bottom = temperature[end]
    dz[1] = 0.2  # triggers split
    mp = _layer_mp()

    (t_out, _, _, _, _, _, _, _, _, _) =
        GEMB.manage_layer_thickness(temperature, dz, density, water, grain_radius,
            grain_dendricity, grain_sphericity, albedo, albedo_diffuse,
            mp, true)

    @test t_out[end] == t_orig_bottom
end

@testset "Bottom merge logic" begin
    temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, albedo, albedo_diffuse = _make_layer_inputs()

    dz[end] = 0.01      # below dzmin; the bottom cell has no cell below to merge into
    dz[end-1] = 0.05
    m_before = _col_mass(dz, density, water)
    mp = _layer_mp()

    (_, dz_out, d_out, w_out, _, _, _, _, _, _) =
        GEMB.manage_layer_thickness(temperature, dz, density, water, grain_radius,
            grain_dendricity, grain_sphericity, albedo, albedo_diffuse,
            mp, true)

    # The bottom cell merges *upward* into the deepest surviving cell above it, and count
    # control then restores the length.
    @test length(dz_out) == 10
    @test _col_mass(dz_out, d_out, w_out) ≈ m_before atol = 1e-9
end

@testset "trim_bottom! (mass controller)" begin
    # Depth pinning is signed: shrink the bottom cell when the column is too deep
    # (accumulation), grow it when too shallow (basal accretion under ablation).

    @testset "shrink — accumulation regime" begin
        n = 5
        cols = GEMB.column_state(260.0 * ones(n), 0.2 * ones(n), 400.0 * ones(n),
            zeros(n), 0.5 * ones(n), 0.5 * ones(n), 0.5 * ones(n),
            0.8 * ones(n), 0.8 * ones(n))
        z_target = 0.95           # column is 1.0 m -> trim 0.05 m off the bottom cell

        mass_added, e_added = GEMB.trim_bottom!(cols, z_target)

        @test sum(cols.dz) ≈ z_target atol = 1e-12
        @test cols.dz[n] ≈ 0.15 atol = 1e-12
        # Mass leaves through the base, so the reported delta is negative and matches the
        # thickness change exactly.
        @test mass_added ≈ -(0.05 * 400.0) atol = 1e-12
        @test e_added ≈ -(0.05 * 400.0 * 260.0 * GEMB.C_ICE) atol = 1e-6
    end

    @testset "grow — ablation regime (basal accretion)" begin
        n = 5
        cols = GEMB.column_state(260.0 * ones(n), 0.2 * ones(n), 400.0 * ones(n),
            zeros(n), 0.5 * ones(n), 0.5 * ones(n), 0.5 * ones(n),
            0.8 * ones(n), 0.8 * ones(n))
        z_target = 1.05           # column is 1.0 m -> accrete 0.05 m onto the bottom cell

        mass_added, e_added = GEMB.trim_bottom!(cols, z_target)

        @test sum(cols.dz) ≈ z_target atol = 1e-12
        @test cols.dz[n] ≈ 0.25 atol = 1e-12
        # Mass enters through the base: same magnitude, opposite sign.
        @test mass_added ≈ +(0.05 * 400.0) atol = 1e-12
        @test e_added ≈ +(0.05 * 400.0 * 260.0 * GEMB.C_ICE) atol = 1e-6
    end

    @testset "pore water is carried proportionally" begin
        n = 3
        cols = GEMB.column_state(273.15 * ones(n), [0.2, 0.2, 0.2], 400.0 * ones(n),
            [0.0, 0.0, 10.0], 0.5 * ones(n), 0.5 * ones(n), 0.5 * ones(n),
            0.8 * ones(n), 0.8 * ones(n))
        # Remove half of the bottom cell -> half of its water leaves with it.
        mass_added, _ = GEMB.trim_bottom!(cols, 0.5)

        @test cols.dz[n] ≈ 0.1 atol = 1e-12
        @test cols.water[n] ≈ 5.0 atol = 1e-12
        @test mass_added ≈ -(0.1 * 400.0 + 5.0) atol = 1e-12
    end

    @testset "already on target is a no-op" begin
        n = 4
        dz = 0.25 * ones(n)
        cols = GEMB.column_state(260.0 * ones(n), dz, 400.0 * ones(n),
            zeros(n), 0.5 * ones(n), 0.5 * ones(n), 0.5 * ones(n),
            0.8 * ones(n), 0.8 * ones(n))
        mass_added, e_added = GEMB.trim_bottom!(cols, 1.0)
        @test mass_added == 0.0
        @test e_added == 0.0
        @test dz == 0.25 * ones(n)
    end

    @testset "an adjustment larger than the bottom cell is an error" begin
        n = 3
        cols = GEMB.column_state(260.0 * ones(n), [0.4, 0.4, 0.2], 400.0 * ones(n),
            zeros(n), 0.5 * ones(n), 0.5 * ones(n), 0.5 * ones(n),
            0.8 * ones(n), 0.8 * ones(n))
        # Needs to remove 0.3 m from a 0.2 m bottom cell.
        @test_throws ErrorException GEMB.trim_bottom!(cols, 0.7)
    end
end
