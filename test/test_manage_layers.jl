# Tests for manage_layers - translated from MATLAB test_manage_layers.m

# Common setup helper for manage_layers tests
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

function _layer_mp(; column_dzmin=0.05, column_dzmax=0.10, column_zmax=1.0, column_zmin=0.5,
                     column_ztop=2.0, column_zy=1.1)
    return GEMB.ModelParameters(
        column_dzmin=column_dzmin,
        column_dzmax=column_dzmax,
        column_zmax=column_zmax,
        column_zmin=column_zmin,
        column_ztop=column_ztop,
        column_zy=column_zy,
    )
end

@testset "No action needed" begin
    temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, albedo, albedo_diffuse = _make_layer_inputs()
    mp = _layer_mp(column_zmax=1.0)  # matches total depth exactly
    verbose = false

    (t_out, dz_out, d_out, _, _, _, _, _, _, m_add, e_add) =
        GEMB.manage_layers(temperature, dz, density, water, grain_radius,
            grain_dendricity, grain_sphericity, albedo, albedo_diffuse,
            mp, verbose)

    @test length(dz_out) == 10
    @test dz_out == dz
    @test t_out == temperature
    @test m_add == 0.0
    @test e_add == 0.0
end

@testset "Merge small layer" begin
    temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, albedo, albedo_diffuse = _make_layer_inputs()

    dz[1] = 0.01
    dz[2] = 0.05
    mp = _layer_mp(column_zmax=sum(dz))  # prevent padding
    verbose = false

    m1 = dz[1] * density[1]
    m2 = dz[2] * density[2]
    m_total = m1 + m2
    expected_dz = dz[1] + dz[2]
    expected_d = m_total / expected_dz

    (_, dz_out, d_out, _, _, _, _, _, _, _, _) =
        GEMB.manage_layers(temperature, dz, density, water, grain_radius,
            grain_dendricity, grain_sphericity, albedo, albedo_diffuse,
            mp, verbose)

    @test length(dz_out) == 9
    @test dz_out[1] ≈ expected_dz atol = 1e-10
    @test d_out[1] ≈ expected_d atol = 1e-10
end

@testset "Split large layer" begin
    temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, albedo, albedo_diffuse = _make_layer_inputs()

    dz[1] = 0.2  # > dzmax
    mp = _layer_mp(column_zmax=sum(dz))  # total depth = 1.1
    verbose = false

    (t_out, dz_out, d_out, _, _, _, _, _, _, _, _) =
        GEMB.manage_layers(temperature, dz, density, water, grain_radius,
            grain_dendricity, grain_sphericity, albedo, albedo_diffuse,
            mp, verbose)

    @test length(dz_out) == 11
    @test dz_out[1] ≈ 0.1 atol = 1e-10
    @test dz_out[2] ≈ 0.1 atol = 1e-10
    @test t_out[1] == t_out[2]
    @test d_out[1] == d_out[2]
end

@testset "Add bottom layer" begin
    n = 5
    temperature = 260.0 * ones(n)
    dz = 0.1 * ones(n)  # total 0.5m
    density = 400.0 * ones(n)
    water = zeros(n)
    grain_radius = 0.5 * ones(n)
    grain_dendricity = 0.5 * ones(n)
    grain_sphericity = 0.5 * ones(n)
    albedo = 0.8 * ones(n)
    albedo_diffuse = 0.8 * ones(n)

    mp = _layer_mp(column_zmax=1.0)  # target 1.0m -> adds padding
    verbose = false

    (_, dz_out, _, _, _, _, _, _, _, m_add, _) =
        GEMB.manage_layers(temperature, dz, density, water, grain_radius,
            grain_dendricity, grain_sphericity, albedo, albedo_diffuse,
            mp, verbose)

    @test length(dz_out) == 6
    @test dz_out[end] == dz_out[end-1]

    m_expected = dz[end] * density[end]
    @test m_add ≈ m_expected atol = 1e-10
end

@testset "Remove bottom layer" begin
    temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, albedo, albedo_diffuse = _make_layer_inputs()

    mp = _layer_mp(column_zmax=0.95)  # total 1.0m > 0.95 -> remove
    verbose = false

    (_, dz_out, _, _, _, _, _, _, _, m_add, _) =
        GEMB.manage_layers(temperature, dz, density, water, grain_radius,
            grain_dendricity, grain_sphericity, albedo, albedo_diffuse,
            mp, verbose)

    @test length(dz_out) == 9
    m_removed_expected = -(0.1 * 400.0)
    @test m_add ≈ m_removed_expected atol = 1e-10
end

@testset "Bottom temperature boundary condition" begin
    temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, albedo, albedo_diffuse = _make_layer_inputs()

    t_orig_bottom = temperature[end]
    dz[1] = 0.2  # triggers split
    mp = _layer_mp(column_zmax=sum(dz))  # match new depth 1.1m
    verbose = false

    (t_out, _, _, _, _, _, _, _, _, _, _) =
        GEMB.manage_layers(temperature, dz, density, water, grain_radius,
            grain_dendricity, grain_sphericity, albedo, albedo_diffuse,
            mp, verbose)

    @test t_out[end] == t_orig_bottom
end

@testset "Conservation check (verbose)" begin
    temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, albedo, albedo_diffuse = _make_layer_inputs()

    dz[1] = 0.01
    dz[2] = 0.05
    mp = _layer_mp(column_zmax=sum(dz))
    verbose = true

    # Should not throw - conservation checks pass internally
    (_, _, _, _, _, _, _, _, _, _, _) =
        GEMB.manage_layers(temperature, dz, density, water, grain_radius,
            grain_dendricity, grain_sphericity, albedo, albedo_diffuse,
            mp, verbose)
    @test true  # If we reach here, conservation passed
end

@testset "Bottom merge logic" begin
    temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, albedo, albedo_diffuse = _make_layer_inputs()

    dz[end] = 0.01
    dz[end-1] = 0.05
    expected_dz_last = dz[end] + dz[end-1]  # save before function mutates dz
    mp = _layer_mp(column_zmax=sum(dz))  # prevent padding
    verbose = false

    (_, dz_out, _, _, _, _, _, _, _, _, _) =
        GEMB.manage_layers(temperature, dz, density, water, grain_radius,
            grain_dendricity, grain_sphericity, albedo, albedo_diffuse,
            mp, verbose)

    @test length(dz_out) == 9
    @test dz_out[end] ≈ expected_dz_last atol = 1e-10
end
