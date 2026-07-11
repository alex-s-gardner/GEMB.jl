# Tests for calculate_melt - translated from MATLAB test_calculate_melt.m

# Common setup helper for melt tests
function _make_melt_inputs(; n=5)
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

function _melt_mp()
    return GEMB.ModelParameters(density_ice=920.0, water_irreducible_saturation=0.07)
end

@testset "Cold dry snow - no change" begin
    temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, albedo, albedo_diffuse = _make_melt_inputs()
    mp = _melt_mp()
    rain = 0.0
    verbose = false

    (t_out, _, d_out, w_out, _, _, _, _, _, m_tot, _, r_tot, f_tot) =
        GEMB.calculate_melt(temperature, dz, density, water, grain_radius,
            grain_dendricity, grain_sphericity, albedo, albedo_diffuse,
            rain, mp, verbose)

    @test m_tot == 0.0
    @test r_tot == 0.0
    @test f_tot == 0.0
    @test t_out == temperature
    @test d_out == density
    @test w_out == water
end

@testset "Pore water refreeze" begin
    temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, albedo, albedo_diffuse = _make_melt_inputs()
    mp = _melt_mp()
    rain = 0.0
    verbose = false

    # Add liquid water to cold snow
    water[1] = 5.0
    mass_initial = density[1] * dz[1] + water[1]

    (t_out, dz_out, d_out, w_out, _, _, _, _, _, _, _, _, f_tot) =
        GEMB.calculate_melt(temperature, dz, density, water, grain_radius,
            grain_dendricity, grain_sphericity, albedo, albedo_diffuse,
            rain, mp, verbose)

    # Verify refreeze occurred
    @test f_tot > 0.0
    @test w_out[1] < 5.0

    # Verify warming from latent heat release
    @test t_out[1] > temperature[1]

    # Verify density increase (refrozen water adds mass to matrix)
    @test d_out[1] > density[1]

    # Mass conservation
    mass_final = d_out[1] * dz_out[1] + w_out[1]
    @test mass_final ≈ mass_initial atol = 1e-10
end

@testset "Surface melt" begin
    temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, albedo, albedo_diffuse = _make_melt_inputs()
    mp = _melt_mp()
    rain = 0.0
    verbose = false

    # Hot surface layer
    temperature[1] = 280.0

    (t_out, _, _, _, _, _, _, _, _, m_tot, m_surf, _, _) =
        GEMB.calculate_melt(temperature, dz, density, water, grain_radius,
            grain_dendricity, grain_sphericity, albedo, albedo_diffuse,
            rain, mp, verbose)

    @test t_out[1] ≈ GEMB.CtoK atol = 1e-10
    @test m_surf > 0.0
    @test m_tot > 0.0
end

@testset "Runoff on ice" begin
    temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, albedo, albedo_diffuse = _make_melt_inputs()
    mp = _melt_mp()
    rain = 0.0
    verbose = false

    # Melt source at top
    temperature[1] = 280.0
    # Thick impermeable ice layer (> 0.1m threshold)
    density[2] = 830.0
    density[3] = 830.0
    # Pre-saturate top layer to trigger runoff
    water[1] = 10.0

    (_, _, _, w_out, _, _, _, _, _, _, _, r_tot, _) =
        GEMB.calculate_melt(temperature, dz, density, water, grain_radius,
            grain_dendricity, grain_sphericity, albedo, albedo_diffuse,
            rain, mp, verbose)

    @test r_tot > 0.0
    # Ice layer retains some irreducible water
    @test w_out[2] > 0.0
    # No water passes through ice layer to layer below
    @test w_out[4] == 0.0
end

@testset "Excess heat distribution" begin
    temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, albedo, albedo_diffuse = _make_melt_inputs()
    mp = _melt_mp()
    rain = 0.0
    verbose = false

    # Huge excess energy at surface
    temperature[1] = GEMB.CtoK + 500.0

    (t_out, _, d_out, _, _, _, _, _, _, _, _, _, _) =
        GEMB.calculate_melt(temperature, dz, density, water, grain_radius,
            grain_dendricity, grain_sphericity, albedo, albedo_diffuse,
            rain, mp, verbose)

    if length(d_out) < 5
        # Top cell melted completely away
        @test true
    else
        # Excess heat should have warmed underlying layer
        @test t_out[2] > 260.0
    end
end

@testset "Water squeezing" begin
    temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, albedo, albedo_diffuse = _make_melt_inputs()
    mp = _melt_mp()
    rain = 0.0
    verbose = false

    # Very wet top layer, isothermal at 0C to prevent refreeze
    water[1] = 20.0
    temperature .= GEMB.CtoK

    (_, _, _, w_out, _, _, _, _, _, _, _, r_tot, _) =
        GEMB.calculate_melt(temperature, dz, density, water, grain_radius,
            grain_dendricity, grain_sphericity, albedo, albedo_diffuse,
            rain, mp, verbose)

    # Excess water should drain from top cell
    @test w_out[1] < 20.0
    # Water should move down or run off
    @test (r_tot > 0.0 || sum(w_out[2:end]) > 0.0)
end

@testset "Rain accounting" begin
    temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, albedo, albedo_diffuse = _make_melt_inputs()
    mp = _melt_mp()
    verbose = false

    # Generate melt
    temperature[1] = 280.0

    # Baseline without rain
    (_, _, _, _, _, _, _, _, _, m_tot_base, _, _, _) =
        GEMB.calculate_melt(temperature, dz, density, water, grain_radius,
            grain_dendricity, grain_sphericity, albedo, albedo_diffuse,
            0.0, mp, verbose)

    # With rain input
    rain_input = 0.5
    (_, _, _, _, _, _, _, _, _, m_tot_rain, _, _, _) =
        GEMB.calculate_melt(temperature, dz, density, water, grain_radius,
            grain_dendricity, grain_sphericity, albedo, albedo_diffuse,
            rain_input, mp, verbose)

    # M_total with rain should subtract rain input (accounting logic)
    expected = max(0.0, m_tot_base - rain_input)
    @test m_tot_rain ≈ expected atol = 1e-6
end
