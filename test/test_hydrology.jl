# Tests for the lateral-drainage / ponding runoff schemes (src/hydrology.jl and the
# ponding path in src/calculate_melt.jl). No MATLAB counterpart exists: MATLAB has only the
# `:instantaneous` behaviour, so the reference test here is that `:instantaneous` is unchanged.

function _drain_cols(; n=5, density=400.0, water=0.0, grain_radius=0.5, dz=0.1)
    return GEMB.column_state(
        260.0 * ones(n),          # temperature
        dz * ones(n),             # dz
        density * ones(n),        # density
        water * ones(n),          # water
        grain_radius * ones(n),   # grain_radius
        0.5 * ones(n),            # grain_dendricity
        0.5 * ones(n),            # grain_sphericity
        zeros(n),                 # age
    )
end

# Irreducible retention of a cell, the same expression `calculate_melt` and
# `apply_lateral_drainage!` use.
_irr(mp, density, dz) =
    (mp.density_ice - density) * GEMB.irreducible_saturation(mp, density) * dz

# A `ModelParameters` that runs through `gemb_core`: the thermal sub-stepping needs
# `dt_divisors`, which only the driver normally fills in.
_hydro_mp(; dt=3600.0, kwargs...) = GEMB.ModelParameters(;
    density_ice=917.0,
    water_irreducible_saturation=0.07,
    thermal_conductivity_method=:Sturm,
    column_dzmin=0.05,
    column_dzmax=0.10,
    column_depth_max=0.9,
    densification_method=:HerronLangway,
    dt_divisors=GEMB.fast_divisors(round(Int, dt * 10000)) ./ 10000,
    kwargs...)

# Single-timestep forcing, mild enough that melt comes from the prescribed warm surface
# rather than from the atmosphere.
_hydro_cfs(; dt=3600.0) = GEMB.ClimateForcingStep(
    dt, 265.0, 100000.0, 0.0, 5.0, 200.0, 300.0, 400.0,
    260.0, 5.0, 200.0, 2.0, 2.0, 0.0, 0.0, 0.0, 60.0, 50.0, 0.0)

_hydro_state(; n=10) = (
    temperature=260.0 * ones(n),
    dz=0.08 * ones(n),
    density=400.0 * ones(n),
    water=zeros(n),
    grain_radius=0.5 * ones(n),
    grain_dendricity=0.5 * ones(n),
    grain_sphericity=0.5 * ones(n),
    age=zeros(n),
    evaporation_condensation=0.0,
    melt_surface=0.0,
)

@testset "runoff_timescale (Zuo and Oerlemans 1996 eq. 22)" begin
    c1, c2, c3 = 0.33 * 86400.0, 25.0 * 86400.0, 140.0

    for slope in (0.0, 0.005, 0.01, 0.1)
        mp = GEMB.ModelParameters(runoff_method=:ZuoOerlemans, surface_slope=slope)
        @test GEMB.runoff_timescale(mp) ≈ c1 + c2 * exp(-c3 * slope) rtol = 1e-14
    end

    # Flat terrain saturates at c1 + c2 (~25 d) rather than diverging, so it still drains
    mp_flat = GEMB.ModelParameters(runoff_method=:ZuoOerlemans, surface_slope=0.0)
    @test GEMB.runoff_timescale(mp_flat) ≈ c1 + c2
    @test GEMB.runoff_timescale(mp_flat) / 86400 ≈ 25.33 atol = 0.01

    # Steeper drains faster, asymptotically to c1
    τ = [GEMB.runoff_timescale(GEMB.ModelParameters(runoff_method=:ZuoOerlemans,
        surface_slope=s)) for s in 0.0:0.005:0.2]
    @test all(diff(τ) .< 0)
    @test GEMB.runoff_timescale(GEMB.ModelParameters(runoff_method=:ZuoOerlemans,
        surface_slope=1.0)) ≈ c1 rtol = 1e-6

    # The RetMIP site slopes (0-0.6 degrees) all sit in the slow, c2-dominated regime
    for slope in (0.0, 0.0035, 0.0087, 0.0105)   # Summit, Dye-2, KAN_U, FA
        @test GEMB.runoff_timescale(GEMB.ModelParameters(runoff_method=:ZuoOerlemans,
            surface_slope=slope)) / 86400 > 6.0
    end
end

@testset "hydraulic_conductivity_saturated (Calonne et al. 2012 eq. 6)" begin
    # Hand-computed in metres, which is the unit trap: GEMB carries grain_radius in mm and a
    # missed conversion is a factor of 1e6 (the law goes as r^2).
    for (r_mm, ρ) in ((0.5, 400.0), (1.0, 600.0), (0.05, 300.0))
        r_m = r_mm * 1e-3
        expected = 3 * r_m^2 * 1000.0 * 9.81 / 1.792e-3 * exp(-0.013 * ρ)
        @test GEMB.hydraulic_conductivity_saturated(r_mm, ρ) ≈ expected rtol = 1e-14
    end

    # Magnitude sanity: snow lands in the 1e-4 to 1e-1 m s-1 band the permeability
    # literature reports. This is the unit guard — a missed mm/m conversion would put the
    # 0.5 mm case at ~2e4 m s-1, and a doubly-converted one at ~2e-8.
    for (r_mm, ρ) in ((0.05, 300.0), (0.5, 400.0), (1.0, 600.0))
        @test 1e-8 < GEMB.hydraulic_conductivity_saturated(r_mm, ρ) < 1e-1
    end
    @test 1e-4 < GEMB.hydraulic_conductivity_saturated(0.5, 400.0) < 1e-1

    # Denser firn is less permeable; coarser grains more so
    @test GEMB.hydraulic_conductivity_saturated(0.5, 800.0) <
          GEMB.hydraulic_conductivity_saturated(0.5, 400.0)
    @test GEMB.hydraulic_conductivity_saturated(1.0, 400.0) >
          GEMB.hydraulic_conductivity_saturated(0.5, 400.0)
end

@testset "relative_permeability (van Genuchten / Yamaguchi)" begin
    # Hand-computed from Yamaguchi et al. (2012) eqs. 6-7 and Hirashima et al. (2010) eq. 10
    function k_rel_reference(r_mm, ρ, θ_e)
        r = r_mm * 1e-3
        n = 1 + 2.7e-3 * (ρ / (2 * r))^0.61
        m = 1 - 1 / n
        return sqrt(θ_e) * (1 - (1 - θ_e^(1 / m))^m)^2
    end

    for r in (0.5, 1.0), ρ in (400.0, 700.0), θ in (0.01, 0.3, 0.9)
        @test GEMB.relative_permeability(r, ρ, θ) ≈ k_rel_reference(r, ρ, θ) rtol = 1e-14
    end

    # Bounded in [0, 1] and monotone in saturation: a barely-wet cell drains far more slowly
    # than a saturated one at the same conductivity
    ks = [GEMB.relative_permeability(0.5, 400.0, θ) for θ in 1e-9:0.01:(1 - 1e-9)]
    @test all(0 .<= ks .<= 1)
    @test all(diff(ks) .> 0)
    @test GEMB.relative_permeability(0.5, 400.0, 1.0 - 1e-12) ≈ 1.0 atol = 1e-4
end

@testset "apply_lateral_drainage! is a no-op under :instantaneous" begin
    mp = GEMB.ModelParameters()   # default runoff_method
    @test mp.runoff_method === :instantaneous

    cols = _drain_cols(water=5.0)
    water_before = copy(cols.water)
    @test GEMB.apply_lateral_drainage!(cols, 3600.0, mp) == 0.0
    @test cols.water == water_before
end

@testset "apply_lateral_drainage! :ZuoOerlemans" begin
    dt = 86400.0
    slope = 0.01
    mp = GEMB.ModelParameters(density_ice=917.0, runoff_method=:ZuoOerlemans,
        surface_slope=slope)
    irr = _irr(mp, 400.0, 0.1)
    @test irr > 0

    # A cell holding exactly its irreducible water does not drain at all
    cols = _drain_cols(water=irr)
    @test GEMB.apply_lateral_drainage!(cols, dt, mp) == 0.0
    @test cols.water ≈ fill(irr, 5)

    # Excess drains by exactly excess * dt/tau
    excess = 3.0
    cols = _drain_cols(water=irr + excess)
    drained = GEMB.apply_lateral_drainage!(cols, dt, mp)
    fraction = dt / GEMB.runoff_timescale(mp)
    @test 0 < fraction < 1
    @test drained ≈ 5 * excess * fraction rtol = 1e-12
    @test all(cols.water .≈ irr + excess * (1 - fraction))
    # Never drained below irreducible
    @test all(cols.water .> irr)

    # A step longer than the timescale drains exactly the excess and no more
    cols = _drain_cols(water=irr + excess)
    drained_long = GEMB.apply_lateral_drainage!(cols, 100 * GEMB.runoff_timescale(mp), mp)
    @test drained_long ≈ 5 * excess rtol = 1e-12
    @test all(cols.water .≈ irr)

    # Steeper slope drains more in the same step
    mp_steep = GEMB.ModelParameters(density_ice=917.0, runoff_method=:ZuoOerlemans,
        surface_slope=0.05)
    cols_a = _drain_cols(water=irr + excess)
    cols_b = _drain_cols(water=irr + excess)
    @test GEMB.apply_lateral_drainage!(cols_b, dt, mp_steep) >
          GEMB.apply_lateral_drainage!(cols_a, dt, mp)

    # Mass leaving the column equals mass lost from the cells
    cols = _drain_cols(water=irr + excess)
    before = sum(cols.water)
    drained = GEMB.apply_lateral_drainage!(cols, dt, mp)
    @test before - sum(cols.water) ≈ drained rtol = 1e-14
end

@testset "apply_lateral_drainage! :Darcy" begin
    dt = 3600.0
    slope = 0.01
    mp = GEMB.ModelParameters(density_ice=917.0, runoff_method=:Darcy, surface_slope=slope)
    irr = _irr(mp, 400.0, 0.1)

    # Nothing at or below irreducible drains
    cols = _drain_cols(water=irr)
    @test GEMB.apply_lateral_drainage!(cols, dt, mp) == 0.0

    # The flux matches Darcy's law computed independently
    excess = 0.5
    cols = _drain_cols(water=irr + excess)
    drained = GEMB.apply_lateral_drainage!(cols, dt, mp)

    capacity = mp.pore_saturation_max * (mp.density_ice - 400.0) * 0.1
    θ_e = clamp(excess / (capacity - irr), 1e-9, 1 - 1e-9)
    flux = 1000.0 * dt * slope *
           GEMB.hydraulic_conductivity_saturated(0.5, 400.0) *
           GEMB.relative_permeability(0.5, 400.0, θ_e)
    @test drained ≈ 5 * min(excess, flux) rtol = 1e-12
    @test all(cols.water .>= irr - 1e-13)

    # Capped at the excess, so a cell is never drained below irreducible even with a huge step
    cols = _drain_cols(water=irr + excess)
    @test GEMB.apply_lateral_drainage!(cols, 1e9, mp) ≈ 5 * excess rtol = 1e-12
    @test all(cols.water .≈ irr)

    # Coarse grains drain faster than fine at equal saturation, unlike the timescale law
    cols_fine = _drain_cols(water=irr + excess, grain_radius=0.1)
    cols_coarse = _drain_cols(water=irr + excess, grain_radius=2.0)
    @test GEMB.apply_lateral_drainage!(cols_coarse, dt, mp) >
          GEMB.apply_lateral_drainage!(cols_fine, dt, mp)

    # Slope scaling
    mp_flat = GEMB.ModelParameters(density_ice=917.0, runoff_method=:Darcy,
        surface_slope=0.001)
    cols_a = _drain_cols(water=irr + excess)
    cols_b = _drain_cols(water=irr + excess)
    @test GEMB.apply_lateral_drainage!(cols_b, dt, mp) >
          GEMB.apply_lateral_drainage!(cols_a, dt, mp_flat)
end

@testset "pore_capacity" begin
    mp = GEMB.ModelParameters(density_ice=917.0)
    @test GEMB.pore_capacity(mp, 400.0, 0.1) ≈ (917.0 - 400.0) * 0.1
    # No pore space at ice density, and never negative above it
    @test GEMB.pore_capacity(mp, 917.0, 0.1) == 0.0
    @test GEMB.pore_capacity(mp, 950.0, 0.1) == 0.0
    # Scales with the cap
    mp_half = GEMB.ModelParameters(density_ice=917.0, pore_saturation_max=0.5)
    @test GEMB.pore_capacity(mp_half, 400.0, 0.1) ≈ 0.5 * (917.0 - 400.0) * 0.1
    # A cell at irreducible saturation is a matching fraction of capacity, which is the
    # consistency the shared (rho_i - rho)*S*dz form exists to guarantee
    mp_c = GEMB.ModelParameters(density_ice=917.0, water_irreducible_saturation=0.07,
        water_irreducible_method=:constant)
    @test _irr(mp_c, 400.0, 0.1) / GEMB.pore_capacity(mp_c, 400.0, 0.1) ≈ 0.07
end

@testset "ponding above a barrier in calculate_melt" begin
    # A melting surface over a thick ice slab: with instantaneous runoff the blocked water
    # leaves at once; with a drainage timescale it backs up into the cells above.
    function _blocked_column(mp)
        n = 6
        temperature = 273.15 * ones(n)
        temperature[1] = 290.0            # strong surface melt
        dz = 0.1 * ones(n)
        density = 500.0 * ones(n)
        density[4] = 900.0                # barrier: thick run at/above impermeable_density
        density[5] = 900.0
        # Start the cells above the barrier at capillary saturation, so the melt this
        # timestep produces is water they cannot hold and the barrier genuinely blocks.
        water = zeros(n)
        for i in 1:3
            water[i] = _irr(mp, density[i], dz[i])
        end
        return (temperature, dz, density, water, 0.5 * ones(n), 0.5 * ones(n),
            0.5 * ones(n), zeros(n))
    end

    mp_inst = GEMB.ModelParameters(density_ice=917.0, water_irreducible_saturation=0.07,
        runoff_method=:instantaneous)
    mp_pond = GEMB.ModelParameters(density_ice=917.0, water_irreducible_saturation=0.07,
        runoff_method=:ZuoOerlemans, surface_slope=0.001)

    out_inst = GEMB.calculate_melt(_blocked_column(mp_inst)..., 0.0, mp_inst, true)
    out_pond = GEMB.calculate_melt(_blocked_column(mp_pond)..., 0.0, mp_pond, true)

    w_inst, r_inst = out_inst[4], out_inst[11]
    w_pond, r_pond = out_pond[4], out_pond[11]
    dz_inst, d_inst = out_inst[2], out_inst[3]
    dz_pond, d_pond = out_pond[2], out_pond[3]

    # Water was actually blocked, so there is something to compare
    @test r_inst > 0

    # Ponding keeps water in the column that instantaneous runoff removes
    @test r_pond < r_inst
    @test sum(w_pond) > sum(w_inst)

    # ... and it is standing water: at least one cell above the barrier now exceeds
    # irreducible, which is impossible under :instantaneous
    # Saturation above irreducible, the same measure `aquifer_diagnostics` thresholds on
    _sat(mp, w, d, dz) = (w - _irr(mp, d, dz)) / ((mp.density_ice - d) * dz)
    above_pond = [_sat(mp_pond, w_pond[i], d_pond[i], dz_pond[i]) for i in 1:3]
    above_inst = [_sat(mp_inst, w_inst[i], d_inst[i], dz_inst[i]) for i in 1:3]
    @test any(above_pond .> GEMB.AQUIFER_TOLERANCE)
    @test all(above_inst .<= GEMB.AQUIFER_TOLERANCE)

    # Filled bottom-up: the cell directly above the barrier holds more than the surface cell
    @test above_pond[3] >= above_pond[1]

    # No cell exceeds its pore capacity
    for i in eachindex(w_pond)
        @test w_pond[i] <= GEMB.pore_capacity(mp_pond, d_pond[i], dz_pond[i]) + 1e-9
    end

    # Total mass is conserved either way (verbose=true above already enforces this inside
    # calculate_melt; this checks the two runs against each other)
    mass_inst = sum(d_inst .* dz_inst) + sum(w_inst) + r_inst
    mass_pond = sum(d_pond .* dz_pond) + sum(w_pond) + r_pond
    @test mass_inst ≈ mass_pond rtol = 1e-9
end

@testset "aquifer_diagnostics" begin
    mp = GEMB.ModelParameters(density_ice=917.0, water_irreducible_saturation=0.07)
    dz = 0.1 * ones(5)
    density = 500.0 * ones(5)
    irr = _irr(mp, 500.0, 0.1)

    # A dry column has no aquifer, and NaN keeps that distinct from a surface water table
    thickness, depth = GEMB.aquifer_diagnostics(dz, density, zeros(5), mp)
    @test thickness == 0.0
    @test isnan(depth)

    # Water exactly at irreducible is capillary-held, not standing
    thickness, depth = GEMB.aquifer_diagnostics(dz, density, fill(irr, 5), mp)
    @test thickness == 0.0
    @test isnan(depth)

    # ... and so is water a hair over it. The threshold is a loose diagnostic one, in units
    # of pore space, because compaction within the timestep leaves a retaining cell marginally
    # over-saturated and that residue must not read as an aquifer.
    pore = (mp.density_ice - 500.0) * 0.1
    thickness, depth = GEMB.aquifer_diagnostics(dz, density,
        fill(irr + 0.1 * GEMB.AQUIFER_TOLERANCE * pore, 5), mp)
    @test thickness == 0.0
    @test isnan(depth)

    # A cell over the threshold does count, so the loosening has not disabled detection
    thickness, depth = GEMB.aquifer_diagnostics(dz, density,
        fill(irr + 10 * GEMB.AQUIFER_TOLERANCE * pore, 5), mp)
    @test thickness ≈ 0.5
    @test depth == 0.0

    # Two saturated cells at depth: thickness sums them, depth points at the shallower top
    water = fill(irr, 5)
    water[3] = irr + 1.0
    water[4] = irr + 1.0
    thickness, depth = GEMB.aquifer_diagnostics(dz, density, water, mp)
    @test thickness ≈ 0.2
    @test depth ≈ 0.2      # top of cell 3

    # A water table at the surface reports depth 0, not NaN
    water = fill(irr, 5)
    water[1] = irr + 1.0
    thickness, depth = GEMB.aquifer_diagnostics(dz, density, water, mp)
    @test thickness ≈ 0.1
    @test depth == 0.0
end

@testset "gemb_core: drainage joins the runoff budget conservatively" begin
    # verbose=true makes gemb_core check mass and energy conservation every timestep, so a
    # drainage term missing from either budget errors here rather than being silently wrong.
    for method in (:ZuoOerlemans, :Darcy)
        # An isothermal column at the melting point: the standing water below stays liquid
        # for the drainage to act on instead of refreezing, which would also warm the bottom
        # cell and trip gemb_core's Dirichlet check for reasons unrelated to drainage.
        state = _hydro_state()
        state.temperature .= 273.15
        state.temperature[1] = 280.0        # drive melt
        state.water .= 2.0                  # standing water to drain
        cfs = _hydro_cfs()
        mp = _hydro_mp(runoff_method=method, surface_slope=0.01)

        new_state, flux = GEMB.gemb_core(state, cfs, mp, true)

        @test flux.runoff >= 0
        @test all(new_state.water .>= 0)
        @test isfinite(flux.aquifer_thickness)
        @test flux.aquifer_thickness >= 0
    end
end

@testset "aquifer grows above an ice slab over many timesteps" begin
    # The integration test the ponding path exists for: sustained melt over a blocking slab
    # should build a water table under a delayed-runoff method, and cannot under
    # :instantaneous. Mass is audited by gemb_core itself (verbose=true) every step.
    function _run(method; slope=0.001, n_steps=48)
        n = 12
        state = (
            temperature=fill(273.15, n),
            dz=0.08 * ones(n),
            density=vcat(fill(500.0, 5), fill(900.0, 3), fill(500.0, 4)),
            water=zeros(n),
            grain_radius=0.5 * ones(n),
            grain_dendricity=0.5 * ones(n),
            grain_sphericity=0.5 * ones(n),
            age=zeros(n),
            evaporation_condensation=0.0,
            melt_surface=0.0,
        )
        # The slab must survive the run for the barrier to keep blocking, so hold the grid
        # fixed: no accumulation, and cells wide enough that nothing is merged or split.
        mp = _hydro_mp(runoff_method=method, surface_slope=slope,
            column_dzmin=0.01, column_dzmax=0.5, column_depth_max=sum(state.dz),
            impermeable_density=830.0)
        # Melting conditions: warm, wet, sunny, and no snowfall, so the surface energy
        # balance sustains melt for the whole run instead of one prescribed warm cell.
        cfs = GEMB.ClimateForcingStep(
            3600.0, 280.0, 100000.0, 0.0, 5.0, 600.0, 350.0, 900.0,
            260.0, 5.0, 200.0, 2.0, 2.0, 0.0, 0.0, 0.0, 60.0, 50.0, 0.0)

        thickness = Float64[]
        runoff = 0.0
        for _ in 1:n_steps
            state, flux = GEMB.gemb_core(state, cfs, mp, true)
            push!(thickness, flux.aquifer_thickness)
            runoff += flux.runoff
        end
        return state, mp, thickness, runoff
    end

    state_pond, mp_pond, thickness_pond, runoff_pond = _run(:ZuoOerlemans)
    state_inst, mp_inst, thickness_inst, runoff_inst = _run(:instantaneous)

    # A water table forms and thickens, where instantaneous runoff leaves none at all
    @test thickness_pond[end] > 0
    @test thickness_pond[end] >= thickness_pond[1]
    @test all(thickness_inst .== 0.0)

    # The water it holds is above irreducible in at least one cell over the slab, and in
    # none at all under :instantaneous
    _sat(mp, s, i) = (s.water[i] - _irr(mp, s.density[i], s.dz[i])) /
                     ((mp.density_ice - s.density[i]) * s.dz[i])
    excess_pond = [_sat(mp_pond, state_pond, i) for i in 1:5]
    excess_inst = [_sat(mp_inst, state_inst, i) for i in 1:5]
    @test any(excess_pond .> GEMB.AQUIFER_TOLERANCE)
    @test all(excess_inst .<= GEMB.AQUIFER_TOLERANCE)

    # Delaying runoff cannot manufacture runoff: over the same melt history the ponding run
    # sheds less, because the difference is still sitting in the column
    @test runoff_pond < runoff_inst
    @test sum(state_pond.water) > sum(state_inst.water)
end

@testset "Marchenko2019 sub-stepping stays stable" begin
    # Marchenko's conductivity is higher than Sturm's in firn, which tightens the Von Neumann
    # limit in calculate_temperature. Confirm the divisor set still resolves a stable step
    # rather than erroring or silently under-resolving.
    for method in (:Sturm, :Calonne, :Calonne2019, :Marchenko2019)
        state = _hydro_state()
        state.density .= 700.0               # firn, where Marchenko departs most from Sturm
        mp = _hydro_mp(thermal_conductivity_method=method, column_depth_max=sum(state.dz))
        new_state, flux = GEMB.gemb_core(state, _hydro_cfs(), mp, true)

        @test all(isfinite, new_state.temperature)
        # No sub-step instability: an under-resolved diffusion step oscillates outside the
        # bracketing temperatures rather than staying between them
        @test all(200.0 .< new_state.temperature .<= 273.16)
        @test isfinite(flux.shortwave_net)
    end
end

@testset "validate_parameters: runoff options" begin
    @test GEMB.validate_parameters(GEMB.ModelParameters()) === nothing
    for m in (:instantaneous, :ZuoOerlemans)
        @test GEMB.validate_parameters(GEMB.ModelParameters(runoff_method=m)) === nothing
    end
    @test GEMB.validate_parameters(
        GEMB.ModelParameters(runoff_method=:Darcy, surface_slope=0.01)) === nothing

    @test_throws AssertionError GEMB.validate_parameters(
        GEMB.ModelParameters(runoff_method=:bucket))
    # :Darcy drains nothing at zero slope while still permitting ponding, so it is rejected
    @test_throws AssertionError GEMB.validate_parameters(
        GEMB.ModelParameters(runoff_method=:Darcy))
    # Slope is a gradient, not degrees
    @test_throws AssertionError GEMB.validate_parameters(
        GEMB.ModelParameters(surface_slope=1.5))
    @test_throws AssertionError GEMB.validate_parameters(
        GEMB.ModelParameters(surface_slope=-0.1))
    @test_throws AssertionError GEMB.validate_parameters(
        GEMB.ModelParameters(pore_saturation_max=0.0))
    @test_throws AssertionError GEMB.validate_parameters(
        GEMB.ModelParameters(pore_saturation_max=1.5))
end
