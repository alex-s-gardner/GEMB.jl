# Tests for calculate_melt

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

# `:constant` retention, which is what most of the assertions below are written against:
# a flat 0.07 saturation at every density, so hand-computed expectations stay tractable.
# The `:ColeouLesaffre` default is exercised by its own testsets.
function _melt_mp()
    return GEMB.ModelParameters(density_ice=920.0, water_irreducible_saturation=0.07,
        water_irreducible_method=:constant)
end

# Enthalpy of the ice matrix above `T1` at `T2` [J kg-1] under the default
# `heat_capacity_method = :CuffeyPaterson`: Cuffey & Paterson (2010) eq. 9.1, c = 152.5 + 7.122·T,
# integrated, which is (T2 − T1) times c at the midpoint.
_cp_dh(T1, T2) = (T2 - T1) * (152.5 + 7.122 * (T1 + T2) / 2)

@testset "Cold dry snow - no change" begin
    temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, age = _make_melt_inputs()
    mp = _melt_mp()
    rain = 0.0
    verbose = false

    (t_out, _, d_out, w_out, _, _, _, _, m_tot, _, r_tot, f_tot) =
        GEMB.calculate_melt(temperature, dz, density, water, grain_radius,
            grain_dendricity, grain_sphericity, age, mp, verbose)

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
            grain_dendricity, grain_sphericity, age, mp, verbose)

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
            grain_dendricity, grain_sphericity, age, mp, verbose)

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
            grain_dendricity, grain_sphericity, age, mp, verbose)

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

    # A surface cell 200 K above the melting point holds 5.62e5 J kg-1 of excess enthalpy,
    # more than LF, so it melts out and passes the surplus down; the surplus is less than the
    # 40 kg cell below can absorb by warming to 0 °C and melting.
    temperature[1] = GEMB.CtoK + 200.0
    M1 = dz[1] * density[1]

    (t_out, _, d_out, _, _, _, _, _, m_tot, _, _, _) =
        GEMB.calculate_melt(temperature, dz, density, water, grain_radius,
            grain_dendricity, grain_sphericity, age, mp, verbose)

    @test length(d_out) == 4
    @test t_out[1] > 260.0
    @test m_tot > M1
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
            grain_dendricity, grain_sphericity, age, mp, verbose)

    # Excess water should drain from top cell
    @test w_out[1] < 20.0
    # Water should move down or run off
    @test (r_tot > 0.0 || sum(w_out[2:end]) > 0.0)
end

@testset "Melt total does not depend on rain" begin
    # Rain arrives as pore water, so melt counts only ice melted by the excess enthalpy:
    # a cell at 280 K melts M·(h(280) − h(273.15))/LF, whatever water it holds.
    temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, age = _make_melt_inputs()
    mp = _melt_mp()
    temperature[1] = 280.0
    M1 = dz[1] * density[1]
    for w in (0.0, 0.5)
        water_in = copy(water); water_in[1] = w
        (_, _, _, _, _, _, _, _, m_tot, _, _, _) =
            GEMB.calculate_melt(copy(temperature), copy(dz), copy(density), water_in,
                copy(grain_radius), copy(grain_dendricity), copy(grain_sphericity), copy(age),
                mp, true)
        @test m_tot ≈ M1 * _cp_dh(GEMB.CtoK, 280.0) / GEMB.LF rtol = 1e-12
    end
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
        # `:constant` explicitly: this knob is the `:constant` method's only input, and
        # `:ColeouLesaffre` (the default) ignores it by design.
        mp = GEMB.ModelParameters(density_ice=920.0, water_irreducible_saturation=S,
            water_irreducible_method=:constant)
        (_, _, _, w_out, _, _, _, _, _, _, r_tot, _) =
            GEMB.calculate_melt(temperature, dz, density, water, grain_radius,
                grain_dendricity, grain_sphericity, age, mp, false)
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
                fill(0.0, n), fill(0.5, n), zeros(n), mp, true)
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
            water_irreducible_method=:constant,
            impermeable_density=density_threshold,
            impermeable_thickness=thickness_threshold)
        (_, _, _, w_out, _, _, _, _, _, _, r_tot, _) =
            GEMB.calculate_melt(temperature, dz, density, water, grain_radius,
                grain_dendricity, grain_sphericity, age, mp, false)
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
            grain_dendricity, grain_sphericity, age, mp, false)
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
            fill(0.5, n), fill(0.5, n), fill(0.5, n), zeros(n), mp, false)
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
            fill(0.5, n), fill(0.5, n), fill(0.5, n), zeros(n), mp_low_ice, false)
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

    # The `gemb` loop's fused one-pass form is bit-identical to the three separate calls,
    # including with a cell straddling each cutoff.
    dz_uneven = [0.3, 0.45, 0.6, 1.1, 2.3, 4.7, 6.2, 8.9]
    density_uneven = [350.0, 420.0, 505.0, 610.0, 700.0, 790.0, 880.0, 925.0]
    for (z, ρ) in ((dz, density_grading), (dz_uneven, density_uneven))
        @test GEMB._firn_air_content_whole_10_20(z, ρ, ρi) ===
              (firn_air_content(z, ρ, ρi), firn_air_content(z, ρ, ρi, 10.0),
               firn_air_content(z, ρ, ρi, 20.0))
    end

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

@testset "melt_geometry" begin
    # Fourteau et al. (2026) Sect. 2.3 argues melting should lower a cell's density at
    # constant thickness (SNOWPACK) rather than shrink its thickness at constant density
    # (Crocus, and GEMB's `:thickness` default). Assertions are on the invariants the two
    # settings are defined by, plus the mass conservation both must satisfy.
    function _run(geometry; T_surface=280.0, n=5)
        temperature, dz, density, water, grain_radius, grain_dendricity,
            grain_sphericity, age = _make_melt_inputs(; n=n)
        temperature[1] = T_surface
        mp = GEMB.ModelParameters(density_ice=920.0, water_irreducible_saturation=0.07,
            water_irreducible_method=:constant, melt_geometry=geometry)
        M_initial = sum(dz .* density) + sum(water)
        out = GEMB.calculate_melt(copy(temperature), copy(dz), copy(density), copy(water),
            copy(grain_radius), copy(grain_dendricity), copy(grain_sphericity), copy(age), mp, true)
        return out, M_initial, dz, density
    end

    # --- the defining invariant of each setting --------------------------------------
    # Cell 1 is the only one that melts here (the melt refreezes below it), so it is the
    # cell whose geometry the setting governs.
    (out_t, M0_t, dz0, rho0) = _run(:thickness)
    (out_d, M0_d, _, _) = _run(:density)
    dz_t, rho_t = out_t[2], out_t[3]
    dz_d, rho_d = out_d[2], out_d[3]

    # `:thickness` — density held, thickness shrinks.
    @test rho_t[1] ≈ rho0[1] rtol = 1e-12
    @test dz_t[1] < dz0[1]

    # `:density` — thickness held, density falls.
    @test dz_d[1] ≈ dz0[1] rtol = 1e-12
    @test rho_d[1] < rho0[1]

    # Both remove the same ice mass from the cell: only the (dz, ρ) split differs.
    @test dz_t[1] * rho_t[1] ≈ dz_d[1] * rho_d[1] rtol = 1e-12

    # And the same melt total, since melt is set by the energy excess, not the geometry.
    @test out_t[9] ≈ out_d[9] rtol = 1e-12

    # --- mass conservation under both (verbose=true already asserts it internally) ----
    for (out, M0) in ((out_t, M0_t), (out_d, M0_d))
        M_final = sum(out[2] .* out[3]) + sum(out[4]) + out[11]
        @test M_final ≈ M0 atol = 1e-9
    end

    # --- the asymmetry `:density` removes -------------------------------------------
    # A melt-refreeze cycle that restores a cell's mass should restore its geometry.
    # Refreezing is at constant thickness under both settings, so only `:density` closes
    # the loop: melt at constant density then refreeze at constant thickness leaves the
    # cell thinner and denser than it started.
    function _cycle(geometry)
        mp = GEMB.ModelParameters(density_ice=920.0, water_irreducible_saturation=0.07,
            water_irreducible_method=:constant, melt_geometry=geometry)
        M = [40.0]
        density = [400.0]
        dz = [0.1]
        # Melt 4 kg out of the cell, then put it straight back as refrozen pore water.
        dz_cell = GEMB.remove_melt!(mp, M, density, dz, 1, 4.0)
        M[1] += 4.0
        density[1] = M[1] / dz_cell
        return dz_cell, density[1]
    end
    dz_cycle_t, rho_cycle_t = _cycle(:thickness)
    dz_cycle_d, rho_cycle_d = _cycle(:density)
    @test dz_cycle_d ≈ 0.1 rtol = 1e-14
    @test rho_cycle_d ≈ 400.0 rtol = 1e-14
    @test dz_cycle_t < 0.1
    @test rho_cycle_t > 400.0

    # --- `remove_melt!` unit behaviour ----------------------------------------------
    # Returns the post-melt thickness the pore-space expressions are written against, and
    # both settings agree on it when nothing melts.
    for geometry in (:thickness, :density)
        mp = GEMB.ModelParameters(density_ice=920.0, melt_geometry=geometry)
        M = [40.0]; density = [400.0]; dz = [0.1]
        @test GEMB.remove_melt!(mp, M, density, dz, 1, 0.0) ≈ 0.1 rtol = 1e-14
        @test density[1] ≈ 400.0 rtol = 1e-14
        @test M[1] ≈ 40.0 rtol = 1e-14
    end

    # A cell that melts out entirely reports zero pore space and keeps a usable density,
    # rather than dividing by a zeroed one. `irreducible_saturation` is evaluated on the
    # density this leaves, so a zero there would give NaN retention.
    for geometry in (:thickness, :density)
        mp = GEMB.ModelParameters(density_ice=920.0, melt_geometry=geometry)
        M = [40.0]; density = [400.0]; dz = [0.1]
        dz_cell = GEMB.remove_melt!(mp, M, density, dz, 1, 40.0)
        @test dz_cell ≈ 0.0 atol = 1e-14
        @test density[1] > 0.0
        @test isfinite(GEMB.irreducible_saturation(mp, density[1]))
    end

    # --- the default is unchanged, and the option is validated -----------------------
    @test GEMB.ModelParameters().melt_geometry === :thickness
    @test initialize_parameters(melt_geometry=:density).melt_geometry === :density
    @test_throws AssertionError initialize_parameters(melt_geometry=:bogus)
end

@testset "Liquid budget over a run with rain, melt and refreeze" begin
    # Over any run, melt + rain − refreeze − runoff equals the change in column liquid water.
    # The air holds exactly the saturation vapor pressure over water at 0 °C, so a melting
    # surface exchanges no vapor with its liquid water and the budget has no fifth term.
    n = 365
    time = DateTime(2001, 1, 1) .+ Day.(0:n-1)
    seasonal = cos.(2π .* (1:n) ./ 365.0 .- π)
    cf = initialize_forcing(time, 268.0 .+ 9.0 .* seasonal, fill(85000.0, n), fill(4.0, n),
        fill(4.0, n), 120.0 .+ 180.0 .* max.(0.0, seasonal), fill(300.0, n),
        fill(GEMB.saturation_vapor_pressure_water(GEMB.CtoK), n);
        temperature_air_mean=268.0, wind_speed_mean=4.0, precipitation_mean=4.0 * 365.25,
        temperature_observation_height=2.0, wind_observation_height=10.0)
    mp = initialize_parameters(output_frequency=:all)
    profile = initialize_profile(mp, cf)
    out = gemb(profile, cf, mp; verbose=true)
    total(k) = sum(parent(out[k]))
    @test total(:melt) > 0.0 && total(:rain) > 0.0 && total(:refreeze) > 0.0
    W0 = sum(parent(profile[:water]))
    W1 = sum(parent(out[:water])[:, end])
    @test total(:melt) + total(:rain) - total(:refreeze) - total(:runoff) ≈ W1 - W0 atol = 1e-6
end

@testset "Melt from the excess enthalpy" begin
    # A cell of mass M raised ΔT above 0 °C melts M·Δh(ΔT)/LF, with Δh the enthalpy above the
    # melting point. Raised by x_full + δ, where Δh(x_full) = LF, it melts out entirely and
    # the cell below receives M·(Δh(x_full + δ) − LF).
    mp = _melt_mp()
    Δh(x) = _cp_dh(GEMB.CtoK, GEMB.CtoK + x)
    # Positive root of (b/2)x² + (a + b·CtoK)x − LF = 0 for eq. 9.1, in BigFloat.
    a, b = big"152.5", big"7.122"
    A = a + b * big(GEMB.CtoK)
    x_full = Float64(2 * big(GEMB.LF) / (A + sqrt(A^2 + 2 * b * big(GEMB.LF))))
    @test Δh(x_full) ≈ GEMB.LF rtol = 1e-14
    n = 6
    T = fill(GEMB.CtoK, n)
    dz = fill(0.1, n)
    ρ = fill(400.0, n)
    M = 0.1 * 400.0
    ΔT = 2.5
    T[1] = GEMB.CtoK + ΔT
    out = GEMB.calculate_melt(copy(T), copy(dz), copy(ρ), zeros(n), fill(0.5, n), zeros(n),
        fill(0.5, n), zeros(n), mp, true)
    @test out[9] ≈ M * Δh(ΔT) / GEMB.LF rtol = 1e-12

    δ = 3.0
    T2 = copy(T); T2[1] = GEMB.CtoK + x_full + δ
    out = GEMB.calculate_melt(T2, copy(dz), copy(ρ), zeros(n), fill(0.5, n), zeros(n),
        fill(0.5, n), zeros(n), mp, true)
    @test length(out[2]) == n - 1                     # cell 1 melted out and was removed
    # Cell 2, at 0 °C, received M·(Δh(x_full + δ) − LF) and melted that over LF of itself.
    @test out[9] ≈ M + M * (Δh(x_full + δ) - GEMB.LF) / GEMB.LF rtol = 1e-12
end

@testset "Percolation into cold firn" begin
    # Cell 1 at 0 °C holds W above its irreducible content; cell 2 is cold. The excess
    # percolates and refreezes in cell 2 up to the smallest of the water arriving, the cell's
    # cold content, and its pore space. Nothing runs off.
    mp = _melt_mp()
    n = 8
    dz = fill(0.1, n)
    ρ = fill(400.0, n)
    T = fill(GEMB.CtoK - 5.0, n)
    T[1] = GEMB.CtoK
    W = 4.0
    water = zeros(n); water[1] = W
    W_irr = (mp.density_ice - 400.0) * 0.07 * 0.1
    cold2 = 400.0 * 0.1 * _cp_dh(GEMB.CtoK - 5.0, GEMB.CtoK) / GEMB.LF
    out = GEMB.calculate_melt(copy(T), copy(dz), copy(ρ), water, fill(0.5, n), zeros(n),
        fill(0.5, n), zeros(n), mp, true)
    @test out[12] ≈ min(W - W_irr, cold2, (mp.density_ice - 400.0) * 0.1) rtol = 1e-12
    @test out[11] == 0.0
    @test out[4][1] ≈ W_irr rtol = 1e-12
end

@testset "Melt on bare ice runs off" begin
    mp = _melt_mp()
    n = 6
    T = fill(GEMB.CtoK, n); T[1] = GEMB.CtoK + 1.5
    out = GEMB.calculate_melt(T, fill(0.1, n), fill(mp.density_ice, n), zeros(n),
        fill(2.5, n), zeros(n), zeros(n), zeros(n), mp, true)
    @test out[11] == out[9]
    @test out[9] > 0.0
end

@testset "Refreeze of expelled pore water leaves the through-flow" begin
    # A cold cell above its irreducible content refreezes all the water arriving (freeze1),
    # expels pore water toward its irreducible content, and refreezes part of that with its
    # remaining cold content (freeze2). Mass: what enters is refrozen, retained, or passed on.
    mp = _melt_mp()
    dz = 0.1
    density = 500.0
    M = density * dz
    water = 5.0
    melt_input = 1.0
    freeze_max = 3.0
    f1, f2, dW, out, M_new, ρ_new = GEMB._refreeze_cell(mp, melt_input, water, freeze_max, M,
        density, dz)
    @test f1 == melt_input
    @test dW < 0.0
    @test f2 > 0.0
    @test out ≈ melt_input - f1 - dW - f2 rtol = 1e-12
    @test melt_input ≈ f1 + f2 + dW + out rtol = 1e-12
    @test M_new ≈ M + f1 + f2 rtol = 1e-12
    @test ρ_new ≈ M_new / dz rtol = 1e-12
end

@testset "SNICAR refrozen-water radius" begin
    # Under :FlannerZender2006 water that freezes in a cell enters its radius at
    # SNICAR_RADIUS_REFROZEN, weighted by the refrozen fraction of the cell's mass (CLM5.0); the
    # mass and energy budget is the same as under the default.
    function melt_run(method; hot)
        temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, age =
            _make_melt_inputs()
        density .= 300.0
        grain_dendricity .= 0.0
        grain_radius .= 0.1
        if hot
            # ~8.8 kg m-2 of melt, above the top cell's 4.3 kg m-2 of retention, percolates
            # into the cold cell below and refreezes there
            temperature[1] = 320.0
        else
            water[1] = 5.0                    # pore water in a cold cell refreezes in place
        end
        mp = GEMB.ModelParameters(density_ice=920.0, water_irreducible_saturation=0.07,
            water_irreducible_method=:constant, grain_growth_method=method)
        return GEMB.calculate_melt(temperature, dz, density, water, grain_radius,
            grain_dendricity, grain_sphericity, age, mp, false)
    end
    for hot in (false, true)
        out_h = melt_run(:hybrid; hot)
        out_f = melt_run(:FlannerZender2006; hot)
        @test out_f[12] > 0                                   # something refroze
        @test all(out_h[5] .== 0.1)                           # default: radius untouched
        for k in (1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13)
            @test out_f[k] == out_h[k]
        end
        if !hot
            # Only cell 1 froze: f = refrozen/(M + W) after freezing, 35 kg m-2 of 30 + 5.
            f = out_f[12] / 35.0
            @test out_f[5][1] ≈ (100 + f * (GEMB.SNICAR_RADIUS_REFROZEN - 100)) / 1000 rtol = 1e-12
            @test all(out_f[5][2:end] .== 0.1)
        else
            @test any(out_f[5] .> 0.1) && all(out_f[5] .<= GEMB.SNICAR_RADIUS_REFROZEN / 1000)
        end
    end
end
