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
    age = zeros(n)
    return temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, age
end

function _melt_mp()
    return GEMB.ModelParameters(density_ice=920.0, water_irreducible_saturation=0.07)
end

@testset "Cold dry snow - no change" begin
    temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, age = _make_melt_inputs()
    mp = _melt_mp()
    rain = 0.0
    verbose = false

    (t_out, _, d_out, w_out, _, _, _, _, m_tot, _, r_tot, f_tot) =
        GEMB.calculate_melt(temperature, dz, density, water, grain_radius,
            grain_dendricity, grain_sphericity, age,
            rain, mp, verbose)

    @test m_tot == 0.0
    @test r_tot == 0.0
    @test f_tot == 0.0
    @test t_out == temperature
    @test d_out == density
    @test w_out == water
end

@testset "Pore water refreeze" begin
    temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, age = _make_melt_inputs()
    mp = _melt_mp()
    rain = 0.0
    verbose = false

    # Add liquid water to cold snow
    water[1] = 5.0
    mass_initial = density[1] * dz[1] + water[1]

    (t_out, dz_out, d_out, w_out, _, _, _, _, _, _, _, f_tot) =
        GEMB.calculate_melt(temperature, dz, density, water, grain_radius,
            grain_dendricity, grain_sphericity, age,
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
    temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, age = _make_melt_inputs()
    mp = _melt_mp()
    rain = 0.0
    verbose = false

    # Hot surface layer
    temperature[1] = 280.0

    (t_out, _, _, _, _, _, _, _, m_tot, m_surf, _, _) =
        GEMB.calculate_melt(temperature, dz, density, water, grain_radius,
            grain_dendricity, grain_sphericity, age,
            rain, mp, verbose)

    @test t_out[1] ≈ GEMB.CtoK atol = 1e-10
    @test m_surf > 0.0
    @test m_tot > 0.0
end

@testset "Runoff on ice" begin
    temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, age = _make_melt_inputs()
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

    (_, _, _, w_out, _, _, _, _, _, _, r_tot, _) =
        GEMB.calculate_melt(temperature, dz, density, water, grain_radius,
            grain_dendricity, grain_sphericity, age,
            rain, mp, verbose)

    @test r_tot > 0.0
    # Ice layer retains some irreducible water
    @test w_out[2] > 0.0
    # No water passes through ice layer to layer below
    @test w_out[4] == 0.0
end

@testset "Excess heat distribution" begin
    temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, age = _make_melt_inputs()
    mp = _melt_mp()
    rain = 0.0
    verbose = false

    # Huge excess energy at surface
    temperature[1] = GEMB.CtoK + 500.0

    (t_out, _, d_out, _, _, _, _, _, _, _, _, _) =
        GEMB.calculate_melt(temperature, dz, density, water, grain_radius,
            grain_dendricity, grain_sphericity, age,
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
    temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, age = _make_melt_inputs()
    mp = _melt_mp()
    rain = 0.0
    verbose = false

    # Very wet top layer, isothermal at 0C to prevent refreeze
    water[1] = 20.0
    temperature .= GEMB.CtoK

    (_, _, _, w_out, _, _, _, _, _, _, r_tot, _) =
        GEMB.calculate_melt(temperature, dz, density, water, grain_radius,
            grain_dendricity, grain_sphericity, age,
            rain, mp, verbose)

    # Excess water should drain from top cell
    @test w_out[1] < 20.0
    # Water should move down or run off
    @test (r_tot > 0.0 || sum(w_out[2:end]) > 0.0)
end

@testset "Rain accounting" begin
    temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, age = _make_melt_inputs()
    mp = _melt_mp()
    verbose = false

    # Generate melt
    temperature[1] = 280.0

    # Baseline without rain
    (_, _, _, _, _, _, _, _, m_tot_base, _, _, _) =
        GEMB.calculate_melt(temperature, dz, density, water, grain_radius,
            grain_dendricity, grain_sphericity, age,
            0.0, mp, verbose)

    # With rain input
    rain_input = 0.5
    (_, _, _, _, _, _, _, _, m_tot_rain, _, _, _) =
        GEMB.calculate_melt(temperature, dz, density, water, grain_radius,
            grain_dendricity, grain_sphericity, age,
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
            grain_sphericity, age = _make_melt_inputs()
        temperature[1] = 280.0      # drive melt at the surface
        mp = GEMB.ModelParameters(density_ice=920.0, water_irreducible_saturation=S)
        (_, _, _, w_out, _, _, _, _, _, _, r_tot, _) =
            GEMB.calculate_melt(temperature, dz, density, water, grain_radius,
                grain_dendricity, grain_sphericity, age, 0.0, mp, false)
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
        (_, _, _, w_out, _, _, _, _, _, _, r_tot, f_tot) =
            GEMB.calculate_melt(temperature, dz, density, water, grain_radius,
                fill(0.0, n), fill(0.5, n), zeros(n), 0.0, mp, true)
        return sum(w_out), r_tot, f_tot
    end

    w_const, r_const, _ = _run(:constant)
    w_cl, r_cl, _ = _run(:ColeouLesaffre)

    @test w_cl > w_const
    @test r_cl <= r_const
end

@testset "impermeability criterion is tunable" begin
    # RetMIP (Vandecrux et al., 2020) identifies this pair as the dominant control on
    # bucket-scheme behaviour at ice-slab sites, with the density criterion spanning
    # 810 (DMIHH) to 917 (DTU) across models. Assert direction of change, not values.
    function _run(; density_threshold=GEMB.DENSITY_PORE_CLOSEOFF, thickness_threshold=0.1,
        lens_density=850.0, lens_cells=2)
        temperature, dz, density, water, grain_radius, grain_dendricity,
            grain_sphericity, age = _make_melt_inputs()
        temperature[1] = 280.0                       # melt source at the surface
        density[2:(1+lens_cells)] .= lens_density    # candidate blocking lens
        # Pre-saturate the top cell, as "Runoff on ice" above does: melt alone refreezes in
        # the cold firn beneath and never reaches the lens, so the criterion would not be
        # exercised at all.
        water[1] = 10.0
        mp = GEMB.ModelParameters(density_ice=920.0, water_irreducible_saturation=0.07,
            impermeable_density=density_threshold,
            impermeable_thickness=thickness_threshold)
        (_, _, _, w_out, _, _, _, _, _, _, r_tot, _) =
            GEMB.calculate_melt(temperature, dz, density, water, grain_radius,
                grain_dendricity, grain_sphericity, age, 0.0, mp, false)
        return r_tot, w_out
    end

    # A 0.2 m lens at 850 kg m-3 blocks under the 830 default, so water runs off rather
    # than reaching the cells below it.
    r_default, w_default = _run()
    @test r_default > 0.0

    # Raising the criterion to 917 puts the same lens below threshold: it no longer blocks,
    # so water passes into the cold firn beneath instead of running off.
    r_high, w_high = _run(density_threshold=917.0)
    @test r_high < r_default
    @test sum(w_high[3:end]) >= sum(w_default[3:end])

    # Lowering it to 810 cannot make a 850 lens *less* blocking than 830 does.
    r_low, _ = _run(density_threshold=810.0)
    @test r_low >= r_default

    # `impermeable_thickness = 0` makes a single dense cell blocking where a 0.1 m
    # threshold lets a lone 0.1 m cell pass.
    r_thin_gate, _ = _run(thickness_threshold=0.0, lens_cells=1)
    r_thick_gate, _ = _run(thickness_threshold=0.1, lens_cells=1)
    @test r_thin_gate > r_thick_gate

    # The nonzero trickle retained inside an impermeable lens under `:constant` (the
    # behaviour pinned by "Runoff on ice" above) survives the parameterization.
    @test w_default[2] > 0.0

    # Validator
    @test_throws AssertionError initialize_parameters(impermeable_density=600.0)
    @test_throws AssertionError initialize_parameters(impermeable_density=1000.0)
    @test_throws AssertionError initialize_parameters(impermeable_thickness=-0.1)
    @test initialize_parameters(impermeable_density=810.0).impermeable_density == 810.0
    # A criterion above `density_ice` is inert, not invalid — it must not be rejected,
    # because lowering `density_ice` is a legitimate sensitivity test.
    @test initialize_parameters(density_ice=800.0).impermeable_density ==
          GEMB.DENSITY_PORE_CLOSEOFF
end

@testset "percolation_depth diagnostic" begin
    function _percolation(T_surface)
        temperature, dz, density, water, grain_radius, grain_dendricity,
            grain_sphericity, age = _make_melt_inputs(n=10)
        temperature[1] = T_surface
        mp = _melt_mp()
        out = GEMB.calculate_melt(temperature, dz, density, water, grain_radius,
            grain_dendricity, grain_sphericity, age, 0.0, mp, false)
        return out[13], sum(out[2])   # percolation_depth, Σdz
    end

    # A cold dry column: the melt block never runs, so the wetting front has zero depth.
    # This is the case that needs the diagnostic initialized outside the conditional.
    depth_dry, z_dry = _percolation(260.0)
    @test depth_dry == 0.0

    # Excess pore water in the *deepest* cell alone, with no surface melt. The percolation
    # loop still runs (there is excess water somewhere) and cannot break until it passes
    # that cell, so the loop index is a poor proxy for how far water travelled: every cell
    # above is walked without water entering it. The front depth must be 0, not the whole
    # column. The driver reduces this layer with `max` over the interval, so one such
    # timestep would otherwise pin the reported depth at the column base.
    let n = 20
        dz = fill(0.1, n)
        water = zeros(n)
        water[n] = 30.0
        mp = GEMB.ModelParameters(density_ice=917.0)
        out = GEMB.calculate_melt(fill(GEMB.CtoK, n), dz, fill(400.0, n), water,
            fill(0.5, n), fill(0.5, n), fill(0.5, n), zeros(n), 0.0, mp, false)
        @test sum(out[4][1:(n-1)]) == 0.0    # confirm no water entered the cells above
        @test out[13] == 0.0                 # so the wetting front went nowhere
    end

    # More surface energy drives the front deeper, and it can never outrun the column.
    depth_warm, z_warm = _percolation(275.0)
    depth_hot, z_hot = _percolation(295.0)
    @test depth_warm > 0.0
    @test depth_hot >= depth_warm
    @test depth_warm <= z_warm
    @test depth_hot <= z_hot
end

@testset "ice_slab_diagnostics" begin
    mp = GEMB.ModelParameters(density_ice=917.0)   # defaults: 830 kg m-3, 0.1 m

    # A 0.5 m slab of 900 kg m-3 starting at 2 m depth. Cells are 0.1 m, so cells 21-25.
    dz = fill(0.1, 40)
    density = fill(400.0, 40)
    density[21:25] .= 900.0
    thickness, depth = GEMB.ice_slab_diagnostics(dz, density, mp)
    @test thickness ≈ 0.5 atol = 1e-12
    @test depth ≈ 2.0 atol = 1e-12

    # Uniform low-density firn: no dense cells at all.
    thickness, depth = GEMB.ice_slab_diagnostics(dz, fill(400.0, 40), mp)
    @test thickness == 0.0
    @test isnan(depth)

    # A single 0.1 m dense cell counts toward `thickness` (a firn core would measure it)
    # but does not exceed the 0.1 m blocking threshold, so no depth is reported. This is
    # exactly the thin-lens-is-laterally-discontinuous assumption of the melt scheme.
    d_thin = fill(400.0, 40)
    d_thin[10] = 900.0
    thickness, depth = GEMB.ice_slab_diagnostics(dz, d_thin, mp)
    @test thickness ≈ 0.1 atol = 1e-12
    @test isnan(depth)

    # `impermeable_thickness = 0` promotes that same lens to blocking.
    mp0 = GEMB.ModelParameters(density_ice=917.0, impermeable_thickness=0.0)
    _, depth0 = GEMB.ice_slab_diagnostics(dz, d_thin, mp0)
    @test depth0 ≈ 0.9 atol = 1e-12

    # The shallowest *qualifying* run wins, and non-qualifying cells above it still count
    # toward the total thickness.
    d_two = fill(400.0, 40)
    d_two[5] = 900.0            # thin, non-blocking
    d_two[15:19] .= 900.0       # 0.5 m, blocking
    d_two[30:34] .= 900.0       # deeper, blocking — must not overwrite the first
    thickness, depth = GEMB.ice_slab_diagnostics(dz, d_two, mp)
    @test thickness ≈ 1.1 atol = 1e-12
    @test depth ≈ 1.4 atol = 1e-12

    # A slab at the surface reports depth 0.0, which is why "no slab" is NaN rather than 0.
    d_surface = fill(400.0, 40)
    d_surface[1:5] .= 900.0
    _, depth_surface = GEMB.ice_slab_diagnostics(dz, d_surface, mp)
    @test depth_surface == 0.0

    # Solid ice always counts, even when `impermeable_density` is set above `density_ice`.
    # That configuration is legal (the validator allows up to 917 so lowering `density_ice`
    # does not reject the default 830), and there `calculate_melt` blocks flow through its
    # unconditional `density_ice` clause — so reporting "no slab" would contradict the
    # physics. Verified against the runoff the same column actually produces below.
    mp_low_ice = GEMB.ModelParameters(density_ice=800.0)   # impermeable_density = 830 > 800
    d_ice = fill(400.0, 20)
    d_ice[5:10] .= 800.0                                    # solid ice under this mp
    t_ice, z_ice = GEMB.ice_slab_diagnostics(fill(0.1, 20), d_ice, mp_low_ice)
    @test t_ice ≈ 0.6 atol = 1e-12
    @test z_ice ≈ 0.4 atol = 1e-12

    # ...and that same column really does route water to runoff, so the diagnostic and the
    # percolation scheme agree rather than merely both being self-consistent.
    let n = 20
        water = zeros(n)
        water[1] = 200.0                    # overwhelm irreducible retention above the ice
        out = GEMB.calculate_melt(fill(GEMB.CtoK, n), fill(0.1, n), copy(d_ice), water,
            fill(0.5, n), fill(0.5, n), fill(0.5, n), zeros(n), 0.0, mp_low_ice, false)
        @test out[11] > 0.0                 # runoff at the slab
        @test sum(out[4][11:n]) == 0.0      # nothing got past it
        # ...and the diagnostic, run on the column `calculate_melt` returned (which is where
        # `gemb_core` runs it), sees the same slab that did the blocking.
        t_after, z_after = GEMB.ice_slab_diagnostics(out[2], out[3], mp_low_ice)
        @test t_after ≈ 0.6 atol = 1e-12
        @test z_after ≈ 0.4 atol = 1e-12
    end

    # A single cell at `density_ice` blocks unconditionally, however thin — below the 0.1 m
    # run threshold that a merely-dense cell must exceed.
    d_thin_ice = fill(400.0, 20)
    d_thin_ice[7] = 917.0
    _, z_thin_ice = GEMB.ice_slab_diagnostics(fill(0.05, 20), d_thin_ice,
        GEMB.ModelParameters(density_ice=917.0))
    @test z_thin_ice ≈ 0.3 atol = 1e-12

    # The threshold is honoured: 850 qualifies under the 830 default, not under 917.
    d_850 = fill(400.0, 40)
    d_850[21:25] .= 850.0
    t_default, z_default = GEMB.ice_slab_diagnostics(dz, d_850, mp)
    t_high, z_high = GEMB.ice_slab_diagnostics(dz, d_850,
        GEMB.ModelParameters(density_ice=917.0, impermeable_density=917.0))
    @test t_default ≈ 0.5 atol = 1e-12
    @test z_default ≈ 2.0 atol = 1e-12
    @test t_high == 0.0
    @test isnan(z_high)
end

@testset "firn_air_content depth limits" begin
    ρi = 917.0

    # Uniform column, hand-computed integral: 20 m of 400 kg m-3 firn holds
    # 20 * (1 - 400/917) m of air, and the 10 m value is exactly half of it.
    dz = fill(0.5, 40)
    density = fill(400.0, 40)
    @test firn_air_content(dz, density, ρi) ≈ 20 * (1 - 400 / ρi) rtol = 1e-14
    @test firn_air_content(dz, density, ρi, 10.0) ≈ 10 * (1 - 400 / ρi) rtol = 1e-14

    # Monotone in depth, and the ordering the driver's three layers must satisfy.
    density_grading = [400.0 + 12.0 * (i - 1) for i in 1:40]   # 400 -> 868
    fac_all = firn_air_content(dz, density_grading, ρi)
    fac_20 = firn_air_content(dz, density_grading, ρi, 20.0)
    fac_10 = firn_air_content(dz, density_grading, ρi, 10.0)
    @test fac_10 <= fac_20 <= fac_all
    @test fac_10 > 0.0

    # `z_max` beyond the column bottom is the whole column — no partial-cell artefact.
    @test firn_air_content(dz, density_grading, ρi, 1000.0) ≈ fac_all rtol = 1e-14

    # A cutoff mid-cell contributes in proportion, not as a whole cell or not at all.
    # Cutting at 10.25 m takes half of cell 21 beyond the 10 m mark.
    half_cell = 0.25 * (1 - density_grading[21] / ρi)
    @test firn_air_content(dz, density_grading, ρi, 10.25) - fac_10 ≈ half_cell rtol = 1e-12

    # Densities above ice contribute zero air rather than negative air.
    @test firn_air_content([1.0], [950.0], ρi) == 0.0
    @test firn_air_content([1.0, 1.0], [950.0, 400.0], ρi) ≈ (1 - 400 / ρi) rtol = 1e-14

    # Degenerate inputs: an empty column and a zero-thickness cell.
    @test firn_air_content(Float64[], Float64[], ρi) == 0.0
    @test firn_air_content([0.0, 1.0], [400.0, 400.0], ρi) ≈ (1 - 400 / ρi) rtol = 1e-14
end
