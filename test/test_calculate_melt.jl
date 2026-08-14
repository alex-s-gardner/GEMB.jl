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
    return temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity
end

function _melt_mp()
    return GEMB.ModelParameters(density_ice=920.0, water_irreducible_saturation=0.07)
end

@testset "Cold dry snow - no change" begin
    temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity = _make_melt_inputs()
    mp = _melt_mp()
    rain = 0.0
    verbose = false

    (t_out, _, d_out, w_out, _, _, _, m_tot, _, r_tot, f_tot) =
        GEMB.calculate_melt(temperature, dz, density, water, grain_radius,
            grain_dendricity, grain_sphericity,
            rain, mp, verbose)

    @test m_tot == 0.0
    @test r_tot == 0.0
    @test f_tot == 0.0
    @test t_out == temperature
    @test d_out == density
    @test w_out == water
end

@testset "Pore water refreeze" begin
    temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity = _make_melt_inputs()
    mp = _melt_mp()
    rain = 0.0
    verbose = false

    # Add liquid water to cold snow
    water[1] = 5.0
    mass_initial = density[1] * dz[1] + water[1]

    (t_out, dz_out, d_out, w_out, _, _, _, _, _, _, f_tot) =
        GEMB.calculate_melt(temperature, dz, density, water, grain_radius,
            grain_dendricity, grain_sphericity,
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
    temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity = _make_melt_inputs()
    mp = _melt_mp()
    rain = 0.0
    verbose = false

    # Hot surface layer
    temperature[1] = 280.0

    (t_out, _, _, _, _, _, _, m_tot, m_surf, _, _) =
        GEMB.calculate_melt(temperature, dz, density, water, grain_radius,
            grain_dendricity, grain_sphericity,
            rain, mp, verbose)

    @test t_out[1] ≈ GEMB.CtoK atol = 1e-10
    @test m_surf > 0.0
    @test m_tot > 0.0
end

@testset "Runoff on ice" begin
    temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity = _make_melt_inputs()
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

    (_, _, _, w_out, _, _, _, _, _, r_tot, _) =
        GEMB.calculate_melt(temperature, dz, density, water, grain_radius,
            grain_dendricity, grain_sphericity,
            rain, mp, verbose)

    @test r_tot > 0.0
    # Ice layer retains some irreducible water
    @test w_out[2] > 0.0
    # No water passes through ice layer to layer below
    @test w_out[4] == 0.0
end

@testset "Excess heat distribution" begin
    temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity = _make_melt_inputs()
    mp = _melt_mp()
    rain = 0.0
    verbose = false

    # Huge excess energy at surface
    temperature[1] = GEMB.CtoK + 500.0

    (t_out, _, d_out, _, _, _, _, _, _, _, _) =
        GEMB.calculate_melt(temperature, dz, density, water, grain_radius,
            grain_dendricity, grain_sphericity,
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
    temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity = _make_melt_inputs()
    mp = _melt_mp()
    rain = 0.0
    verbose = false

    # Very wet top layer, isothermal at 0C to prevent refreeze
    water[1] = 20.0
    temperature .= GEMB.CtoK

    (_, _, _, w_out, _, _, _, _, _, r_tot, _) =
        GEMB.calculate_melt(temperature, dz, density, water, grain_radius,
            grain_dendricity, grain_sphericity,
            rain, mp, verbose)

    # Excess water should drain from top cell
    @test w_out[1] < 20.0
    # Water should move down or run off
    @test (r_tot > 0.0 || sum(w_out[2:end]) > 0.0)
end

@testset "Rain accounting" begin
    temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity = _make_melt_inputs()
    mp = _melt_mp()
    verbose = false

    # Generate melt
    temperature[1] = 280.0

    # Baseline without rain
    (_, _, _, _, _, _, _, m_tot_base, _, _, _) =
        GEMB.calculate_melt(temperature, dz, density, water, grain_radius,
            grain_dendricity, grain_sphericity,
            0.0, mp, verbose)

    # With rain input
    rain_input = 0.5
    (_, _, _, _, _, _, _, m_tot_rain, _, _, _) =
        GEMB.calculate_melt(temperature, dz, density, water, grain_radius,
            grain_dendricity, grain_sphericity,
            rain_input, mp, verbose)

    # M_total with rain should subtract rain input (accounting logic)
    expected = max(0.0, m_tot_base - rain_input)
    @test m_tot_rain ≈ expected atol = 1e-6
end

@testset "irreducible_saturation methods" begin
    # `:constant` returns the parameter regardless of density, including above pore
    # close-off — the default path is deliberately ungated so it stays bit-identical
    # to MATLAB.
    mp_c = GEMB.ModelParameters(density_ice=917.0, water_irreducible_method=:constant,
        water_irreducible_saturation=0.07)
    for d in (200.0, 500.0, 830.0, 916.0)
        @test GEMB.irreducible_saturation(mp_c, d) === 0.07
    end

    # `:ColeouLesaffre` — Coléou & Lesaffre (1998) eq. 3 via Langen et al. (2017) eq. 4.
    mp_cl = GEMB.ModelParameters(density_ice=917.0,
        water_irreducible_method=:ColeouLesaffre, water_irreducible_saturation=0.07)

    swi(ρ) = begin
        wmi = 0.057 * (917.0 - ρ) / ρ + 0.017
        wmi / (1 - wmi) * 917.0 * ρ / (1000.0 * (917.0 - ρ))
    end
    for d in (200.0, 300.0, 500.0, 800.0, 829.0)
        @test GEMB.irreducible_saturation(mp_cl, d) ≈ swi(d) rtol = 1e-14
    end

    # The point of the scheme: retention rises with density, so the constant 0.07
    # under-retains in the percolation zone. Reference values from the CFM cross-check.
    @test GEMB.irreducible_saturation(mp_cl, 300.0) ≈ 0.0691 atol = 1e-4
    @test GEMB.irreducible_saturation(mp_cl, 800.0) ≈ 0.1630 atol = 1e-4
    @test GEMB.irreducible_saturation(mp_cl, 800.0) > GEMB.irreducible_saturation(mp_cl, 300.0)

    # S_wi is not monotone across the full range: it has a shallow minimum near
    # ρ ≈ 310 kg m-3 (0.0691) and rises on both sides — toward low density because
    # porosity grows faster than the retained mass, and toward close-off because pore
    # space vanishes. Assert the minimum is where the formula puts it, then monotone
    # increase above it (which is the firn range that matters here).
    ρ_grid = 150.0:5.0:825.0
    vals = [GEMB.irreducible_saturation(mp_cl, d) for d in ρ_grid]
    @test 300.0 <= ρ_grid[argmin(vals)] <= 325.0
    rising = [GEMB.irreducible_saturation(mp_cl, d) for d in 350.0:5.0:825.0]
    @test all(diff(rising) .> 0.0)

    # Gated at close-off (no connected pore space) — also what removes the ρ → ρᵢ
    # singularity in S_wi.
    @test GEMB.irreducible_saturation(mp_cl, GEMB.DENSITY_PORE_CLOSEOFF) == 0.0
    @test GEMB.irreducible_saturation(mp_cl, 900.0) == 0.0
    @test isfinite(GEMB.irreducible_saturation(mp_cl, 917.0))

    # Validator
    @test_throws AssertionError initialize_parameters(water_irreducible_method=:Bogus)
    @test initialize_parameters(water_irreducible_method=:ColeouLesaffre).water_irreducible_method ===
        :ColeouLesaffre
end

@testset "water_irreducible_saturation is honoured during percolation" begin
    # Upstream GEMB issue #199: MATLAB (and the Julia port before this change) declared a
    # local 0.07 that overrode `mp.water_irreducible_saturation` at every percolation
    # retention site, so the documented knob only affected the pre-percolation squeeze.
    # Retention scales linearly in S_wi, so doubling it must retain strictly more.
    function _retained(S)
        temperature, dz, density, water, grain_radius, grain_dendricity,
            grain_sphericity = _make_melt_inputs()
        temperature[1] = 280.0      # drive melt at the surface
        mp = GEMB.ModelParameters(density_ice=920.0, water_irreducible_saturation=S)
        (_, _, _, w_out, _, _, _, _, _, r_tot, _) =
            GEMB.calculate_melt(temperature, dz, density, water, grain_radius,
                grain_dendricity, grain_sphericity, 0.0, mp, false)
        return sum(w_out), r_tot
    end

    w_low, _ = _retained(0.02)
    w_high, _ = _retained(0.14)
    @test w_high > w_low
end

@testset "ColeouLesaffre retains more than the constant in dense firn" begin
    # In the percolation zone (600-830 kg m-3) Coléou & Lesaffre gives S_wi well above
    # 0.07, so the column must hold more water and shed less as runoff.
    function _run(method)
        n = 8
        temperature = fill(263.0, n)
        temperature[1] = 285.0                       # melt source
        dz = fill(0.2, n)
        density = fill(700.0, n)                     # percolation-zone firn
        water = zeros(n)
        grain_radius = fill(1.0, n)
        mp = GEMB.ModelParameters(density_ice=917.0, water_irreducible_method=method)
        (_, _, _, w_out, _, _, _, _, _, r_tot, f_tot) =
            GEMB.calculate_melt(temperature, dz, density, water, grain_radius,
                fill(0.0, n), fill(0.5, n), 0.0, mp, true)
        return sum(w_out), r_tot, f_tot
    end

    w_const, r_const, _ = _run(:constant)
    w_cl, r_cl, _ = _run(:ColeouLesaffre)

    @test w_cl > w_const
    @test r_cl <= r_const
end
