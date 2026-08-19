using GEMB
using GEMB_ClimateForcing
using Test
using Dates
using GEMB: DimArray, DimStack, Ti, dims, metadata

# Blowing snow: the Crocus `SNOWDRIFT` port (`apply_blowing_snow!`) and the IMAU-FDM-style
# prescribed drift sink (`apply_prescribed_drift!`).
#
# There is no MATLAB reference for either — MATLAB GEMB has no wind-rework term — and the
# Community Firn Model has no blowing-snow physics to compare against, so these are pinned
# against, in order of strength: the closed-form solution of the compaction ODE; expressions
# re-derived here independently of the implementation, in the ordering the SURFEX source uses
# (`snowcro.F90:4685-4950`); the invariants each path claims in its docstring (exact mass
# conservation on the compaction path, `DRIFT_DENSITY_MAX` as a ceiling, monotone decay with
# depth, the loop exit); and the whole-model mass and energy budgets under `verbose=true`,
# which is what actually validates the wiring.

# Built through `column_state` for the same reason `test_horizontal_strain.jl` does: a new
# per-cell field is picked up here by changing one constructor call.
_drift_cols(; dz=[0.1, 0.2, 0.3], density=[150.0, 200.0, 250.0],
    temperature=[260.0, 262.0, 265.0], water=zeros(3),
    grain_radius=[0.15, 0.3, 0.5], dendricity=zeros(3), sphericity=fill(0.5, 3),
    age=zeros(3)) = GEMB.column_state(
    copy(temperature), copy(dz), copy(density), copy(water), copy(grain_radius),
    copy(dendricity), copy(sphericity), copy(age))

# A windy, cold, dry step. At `wind_speed = 15` the 5 m gust speed is ~17 m s-1, well above the
# ~12.4 m s-1 at which every layer with non-negative mobility drifts, so the scheme is
# unambiguously active. `vapor_pressure = 50 Pa` against a ~103 Pa ice saturation gives a
# saturation deficit of ~0.5, so the sublimation term is active too.
_drift_step(; dt=10800.0, wind_speed=15.0, temperature_air=253.0, vapor_pressure=50.0,
    pressure_air=90000.0, snow_drift=0.0) = GEMB.ClimateForcingStep(;
    dt, wind_speed, temperature_air, vapor_pressure, pressure_air, snow_drift,
    temperature_air_mean=253.0, wind_speed_mean=wind_speed,
    precipitation_mean=200.0, temperature_observation_height=2.0,
    wind_observation_height=10.0)

@testset "drift_wind_speed" begin
    cfs = _drift_step(wind_speed=10.0)

    # The log-profile extrapolation and the gust factor, written out independently.
    z0 = 1e-3
    @test GEMB.drift_wind_speed(cfs, z0) ≈
          1.25 * 10.0 * log(5.0 / z0) / log(10.0 / z0) atol = 1e-14

    # The three-way floor on the roughness length. At `z0 = 3 m` the 2.5 m half-reference-height
    # bound is the tightest of the three, so the answer must equal the one at `z0 = 2.5 m`
    # exactly.
    @test GEMB.drift_wind_speed(cfs, 3.0) == GEMB.drift_wind_speed(cfs, 2.5)

    # Rougher surfaces put more of the shear into the lowest few metres, so extrapolating the
    # same 10 m wind down to 5 m gives a *lower* speed the rougher the surface is.
    @test GEMB.drift_wind_speed(cfs, 1e-4) > GEMB.drift_wind_speed(cfs, 1e-2)

    # Linear in wind speed, and zero at zero wind.
    @test GEMB.drift_wind_speed(_drift_step(wind_speed=20.0), z0) ≈
          2 * GEMB.drift_wind_speed(_drift_step(wind_speed=10.0), z0) atol = 1e-13
    @test GEMB.drift_wind_speed(_drift_step(wind_speed=0.0), z0) == 0.0
end

@testset "drift_density_factor" begin
    # Flat at and below `DRIFT_DENSITY_MIN`: the `max(ρ_min, ρ)` of eq. 60.
    @test GEMB.drift_density_factor(50.0) == GEMB.drift_density_factor(10.0)
    @test GEMB.drift_density_factor(50.0) ≈ 1.25 atol = 1e-14

    # The source's slope, `1.25/(1000·XVMOB1)` = 0.0042373 per kg m-3, not the 0.004 the paper
    # rounds it to. Written out here independently of the implementation.
    @test GEMB.drift_density_factor(150.0) ≈ 1.25 - 1.25 * 100.0 / 1000.0 / 0.295 atol = 1e-14

    # Monotone decreasing, crossing zero at ρ_min + 1000·XVMOB1 = 345 kg m-3 — which is what
    # stops a wind slab from being re-eroded indefinitely.
    @test GEMB.drift_density_factor(400.0) < GEMB.drift_density_factor(300.0) <
          GEMB.drift_density_factor(200.0)
    @test GEMB.drift_density_factor(345.0) ≈ 0.0 atol = 1e-14
    @test GEMB.drift_density_factor(400.0) < 0.0
end

@testset "mobility_index" begin
    # Fresh dendritic snow is the most mobile; large faceted grains in dense snow the least.
    @test GEMB.mobility_index(0.15, 0.9, 0.5, 100.0, false) >
          GEMB.mobility_index(1.0, 0.0, 0.2, 350.0, false)

    # Dendritic branch, eq. 59 written out in GEMB's units (fractions on [0,1], not Crocus's
    # coded integers).
    @test GEMB.mobility_index(0.15, 0.9, 0.5, 150.0, false) ≈
          0.34 * (0.5 + 0.75 * 0.9 - 0.5 * 0.5) +
          0.66 * GEMB.drift_density_factor(150.0) atol = 1e-14

    # Non-dendritic branch. `grain_radius` is in mm, so `2·r` is the diameter in mm and the
    # `XVMOB3` term takes it directly (the source's `PSNOWGRAN2·1000` with `PSNOWGRAN2` in m).
    @test GEMB.mobility_index(0.4, 0.0, 0.6, 200.0, false) ≈
          0.34 * (0.833 - 0.833 * 0.6 - 0.583 * 0.8) +
          0.66 * GEMB.drift_density_factor(200.0) atol = 1e-14

    # Within a branch: sphericity enters negatively, and larger grains are less mobile.
    @test GEMB.mobility_index(0.3, 0.0, 0.2, 200.0, false) >
          GEMB.mobility_index(0.3, 0.0, 0.9, 200.0, false)
    @test GEMB.mobility_index(0.2, 0.0, 0.5, 200.0, false) >
          GEMB.mobility_index(0.8, 0.0, 0.5, 200.0, false)

    # The wet cap, eq. 61: a layer holding liquid water is bonded, so its mobility is pinned to
    # `DRIFT_MOB4` no matter how fresh its grains are.
    @test GEMB.mobility_index(0.15, 0.9, 0.5, 100.0, true) == GEMB.DRIFT_MOB4
    @test GEMB.DRIFT_MOB4 < 0.0

    # It is a cap, not an assignment: snow already less mobile than the cap keeps its value.
    mob_dense = GEMB.mobility_index(1.0, 0.0, 0.5, 500.0, false)
    @test mob_dense < GEMB.DRIFT_MOB4
    @test GEMB.mobility_index(1.0, 0.0, 0.5, 500.0, true) == mob_dense
end

@testset "driftability_index and its wind threshold" begin
    # The wind above which every layer with non-negative mobility drifts, from
    # `2.868·exp(-0.085·U) = 1`.
    u_all = log(GEMB.DRIFT_D1) / GEMB.DRIFT_D2
    @test u_all ≈ 12.395 atol = 1e-3
    @test GEMB.driftability_index(0.0, u_all) ≈ 0.0 atol = 1e-14

    # Clamped at zero, and monotone increasing in both arguments above the threshold.
    @test GEMB.driftability_index(0.0, 1.0) == 0.0
    @test GEMB.driftability_index(0.5, 15.0) > GEMB.driftability_index(0.2, 15.0)
    @test GEMB.driftability_index(0.5, 20.0) > GEMB.driftability_index(0.5, 15.0)

    # `drift_threshold_wind_speed` inverts `driftability_index`: at the threshold wind the
    # driftability is exactly zero. This consistency is what the sublimation rate rests on —
    # it divides by `U_t` and raises `U₅/U_t` to a power.
    for mob in (-0.5, 0.0, 0.3, 0.8)
        u_t = GEMB.drift_threshold_wind_speed(mob)
        @test GEMB.driftability_index(mob, u_t) ≈ 0.0 atol = 1e-13
        @test GEMB.driftability_index(mob, u_t + 0.1) > 0.0
    end

    # More mobile snow needs less wind.
    @test GEMB.drift_threshold_wind_speed(0.8) < GEMB.drift_threshold_wind_speed(0.2)
end

@testset "blowing_snow_sublimation_rate" begin
    mp = initialize_parameters()
    cfs = _drift_step(wind_speed=15.0)
    wind_5m = GEMB.drift_wind_speed(cfs, GEMB.surface_roughness(mp, 150.0, 0.0))
    mob = 0.3

    # Eq. 67 written out independently, in the source's exponent ordering (`**4` on the
    # temperature ratio, `**3.6` on the wind ratio) — see `DRIFT_SUBLIM_EXP_TEMPERATURE` for
    # why the source is followed rather than the paper's Table K3, which has them swapped.
    u_t = GEMB.drift_threshold_wind_speed(mob)
    rho_a = GEMB.air_density(cfs.pressure_air, cfs.temperature_air)
    q_sat = GEMB.specific_humidity(
        GEMB.saturation_vapor_pressure_ice(cfs.temperature_air), cfs.pressure_air)
    q_air = GEMB.specific_humidity(cfs.vapor_pressure, cfs.pressure_air)
    expected = 1.8e-3 * (GEMB.CtoK / cfs.temperature_air)^4 * u_t * rho_a * q_sat *
               (1.0 - q_air / q_sat) * (wind_5m / u_t)^3.6
    @test GEMB.blowing_snow_sublimation_rate(mob, wind_5m, cfs) ≈ expected rtol = 1e-14

    # Positive in dry air, and rising faster than the wind ratio itself (exponent > 1).
    @test GEMB.blowing_snow_sublimation_rate(mob, wind_5m, cfs) > 0.0
    @test GEMB.blowing_snow_sublimation_rate(mob, 2 * wind_5m, cfs) >
          4 * GEMB.blowing_snow_sublimation_rate(mob, wind_5m, cfs)

    # Zero at saturation with respect to ice, and negative above it — the caller clamps.
    cfs_sat = _drift_step(vapor_pressure=GEMB.saturation_vapor_pressure_ice(253.0))
    @test GEMB.blowing_snow_sublimation_rate(mob, wind_5m, cfs_sat) == 0.0
    cfs_super = _drift_step(vapor_pressure=1.2 * GEMB.saturation_vapor_pressure_ice(253.0))
    @test GEMB.blowing_snow_sublimation_rate(mob, wind_5m, cfs_super) < 0.0

    # Guards against a non-positive threshold wind rather than returning a nonsense rate.
    @test GEMB.drift_threshold_wind_speed(1.9) < 0.0
    @test GEMB.blowing_snow_sublimation_rate(1.9, wind_5m, cfs) == 0.0
end

@testset "apply_blowing_snow!" begin
    cfs = _drift_step()

    @testset ":none is an exact no-op" begin
        mp = initialize_parameters(blowing_snow_method=:none)
        cols = _drift_cols()
        before = deepcopy(cols)
        m, e = GEMB.apply_blowing_snow!(cols, cfs, mp)
        @test (m, e) === (0.0, 0.0)
        for k in keys(cols)
            @test cols[k] == before[k]          # bit-for-bit
        end
    end

    @testset "compaction conserves column mass exactly" begin
        mp = initialize_parameters(blowing_snow_method=:Crocus)
        cols = _drift_cols()
        mass_before = GEMB.column_mass(cols)
        m, e = GEMB.apply_blowing_snow!(cols, cfs, mp)
        # Nothing crosses the surface without `blowing_snow_sublimation`, so the scheme
        # contributes nothing to either budget.
        @test (m, e) === (0.0, 0.0)
        @test GEMB.column_mass(cols) ≈ mass_before rtol = 1e-15
        @test cols.density[1] > 150.0           # ...and it did do something
        @test cols.dz[1] < 0.1                  # at constant mass, so `dz` shrank
        @test cols.dz[1] * cols.density[1] ≈ 0.1 * 150.0 rtol = 1e-15
    end

    @testset "closed form: ρ(t) → ρ_max - (ρ_max - ρ0)·exp(-DEFF·t/τ)" begin
        # Eq. 64 with the microstructure and wind held fixed is `dρ/dt = (DEFF/τ)(ρ_max - ρ)`,
        # whose solution is the exponential above. The implementation takes the forward-Euler
        # step of that ODE, so this checks two things: that the step is exactly the Euler step
        # of the right ODE (to 1e-12, with every factor re-derived here), and that over a step
        # short enough for the Euler error to be small it agrees with the exponential.
        mp = initialize_parameters(blowing_snow_method=:Crocus)
        rho_max = GEMB.DRIFT_DENSITY_MAX
        rho0 = 150.0
        dt = 60.0

        cols = _drift_cols(dz=[0.05], density=[rho0], temperature=[260.0], water=[0.0],
            grain_radius=[0.3], dendricity=[0.0], sphericity=[0.5], age=[0.0])
        wind_5m = GEMB.drift_wind_speed(_drift_step(dt=dt),
            GEMB.surface_roughness(mp, rho0, 0.0))
        mob = GEMB.mobility_index(0.3, 0.0, 0.5, rho0, false)
        drift = GEMB.driftability_index(mob, wind_5m)
        @test drift > 0.0
        # `ZPROFEQU` for the surface cell: half of its own thickness.
        decay = 0.5 * 0.05 * 0.1 * (GEMB.DRIFT_D3 - drift)
        deff = drift * exp(-100.0 * decay)
        rate = deff * GEMB.DRIFT_WIND_EFFECT * dt / GEMB.DRIFT_GUST_COEF / GEMB.DRIFT_TAU

        GEMB.apply_blowing_snow!(cols, _drift_step(dt=dt), mp)
        @test cols.density[1] ≈ rho0 + rate * (rho_max - rho0) rtol = 1e-12
        # The Euler step and the exact exponential differ at O(rate²)·(ρ_max - ρ₀).
        @test cols.density[1] ≈ rho_max - (rho_max - rho0) * exp(-rate) rtol = 2 * rate
        @test rate < 1e-3                       # ...so the bound above is a real one
    end

    @testset "DRIFT_DENSITY_MAX is a hard ceiling" begin
        mp = initialize_parameters(blowing_snow_method=:Crocus)
        cols = _drift_cols(dz=[0.05], density=[100.0], temperature=[260.0], water=[0.0],
            grain_radius=[0.2], dendricity=[0.0], sphericity=[0.5], age=[0.0])
        for _ in 1:500
            GEMB.apply_blowing_snow!(cols, _drift_step(dt=86400.0), mp)
            @test cols.density[1] <= GEMB.DRIFT_DENSITY_MAX
        end
        # It genuinely gets there rather than stalling short of the wind-slab density.
        @test cols.density[1] > 340.0

        # Snow already denser than the ceiling is left *exactly* alone — the source's
        # `IF (PSNOWRHO < XVROMAX)` gate, which is why it is a gate rather than a `max(·, 0)`:
        # multiplying `dz` by an exact 1.0 would still be a write.
        dense = _drift_cols(dz=[0.05], density=[400.0], temperature=[260.0], water=[0.0],
            grain_radius=[0.2], dendricity=[0.0], sphericity=[0.5], age=[0.0])
        GEMB.apply_blowing_snow!(dense, cfs, mp)
        @test dense.density[1] == 400.0
        @test dense.dz[1] == 0.05
    end

    @testset "the effect decays monotonically with depth" begin
        mp = initialize_parameters(blowing_snow_method=:Crocus)
        # A uniform column, so the only thing distinguishing the cells is their depth.
        n = 6
        cols = _drift_cols(dz=fill(0.5, n), density=fill(150.0, n), temperature=fill(260.0, n),
            water=zeros(n), grain_radius=fill(0.3, n), dendricity=zeros(n),
            sphericity=fill(0.5, n), age=zeros(n))
        GEMB.apply_blowing_snow!(cols, cfs, mp)
        gains = cols.density .- 150.0
        @test gains[1] > 0.0                            # the surface drifts
        @test issorted(gains, rev=true)                 # and the effect only decays downward
        @test gains[2] < 0.1 * gains[1]                 # steeply: `exp(-100·ZPROFEQU)`
        @test gains[end] < 1e-6 * gains[1]
    end

    @testset "the loop stops at the first non-drifting layer" begin
        mp = initialize_parameters(blowing_snow_method=:Crocus)
        # A wet second cell is capped at `DRIFT_MOB4 = -0.0583`, which at a 5 m gust of
        # ~9.3 m s-1 (threshold term ~0.31) is below the drift threshold, while the dry cell 1
        # at mobility ~0.57 is above it. Cell 3 is identical to cell 1, so it *would* drift if
        # it were reached — it must not be. This is `snowcro.F90:4844`'s `EXIT`; the paper's
        # eq. 63 does not show it, and a literal port of the equation would rework cell 3.
        cols = _drift_cols(dz=[0.1, 0.1, 0.1], density=[150.0, 150.0, 150.0],
            temperature=[260.0, 273.15, 260.0], water=[0.0, 1.0, 0.0],
            grain_radius=[0.3, 0.3, 0.3], dendricity=zeros(3), sphericity=fill(0.5, 3),
            age=zeros(3))
        GEMB.apply_blowing_snow!(cols, _drift_step(wind_speed=8.0), mp)
        @test cols.density[1] > 150.0        # cell 1 drifted
        @test cols.density[2] == 150.0       # cell 2 is wet, so it stopped the loop
        @test cols.density[3] == 150.0       # and cell 3 was never reached
        @test cols.grain_radius[3] == 0.3
        @test cols.grain_sphericity[3] == 0.5
    end

    @testset "no wind is a no-op even when the scheme is on" begin
        mp = initialize_parameters(blowing_snow_method=:Crocus)
        cols = _drift_cols()
        before = deepcopy(cols)
        GEMB.apply_blowing_snow!(cols, _drift_step(wind_speed=0.0), mp)
        # At zero wind the threshold term is 1.868, which no physical mobility reaches, so the
        # loop breaks at cell 1 before writing anything.
        for k in keys(cols)
            @test cols[k] == before[k]
        end
    end

    @testset "fragmentation: dendricity falls, sphericity rises, size falls" begin
        mp = initialize_parameters(blowing_snow_method=:Crocus)

        # Dendritic branch (eq. 65). Grain radius is diagnostic while dendricity is nonzero, so
        # what is pinned is the direction of the two shape variables plus the re-derivation.
        dend = _drift_cols(dz=[0.05], density=[100.0], temperature=[260.0], water=[0.0],
            grain_radius=[0.15], dendricity=[0.9], sphericity=[0.3], age=[0.0])
        GEMB.apply_blowing_snow!(dend, _drift_step(dt=86400.0), mp)
        @test 0.0 < dend.grain_dendricity[1] < 0.9      # falls, never driven negative
        @test 0.3 < dend.grain_sphericity[1] <= 1.0     # rises, never past 1
        @test dend.grain_radius[1] ≈ GEMB.dendritic_grain_radius(
            dend.grain_dendricity[1], dend.grain_sphericity[1]) atol = 1e-15

        # Non-dendritic branch (eq. 66): grains round off toward S = 1 and break down toward
        # the `DRIFT_GRAIN_SIZE_MIN` floor, which is on the *diameter*, hence the factor 0.5.
        nond = _drift_cols(dz=[0.05], density=[100.0], temperature=[260.0], water=[0.0],
            grain_radius=[1.0], dendricity=[0.0], sphericity=[0.3], age=[0.0])
        for _ in 1:500
            GEMB.apply_blowing_snow!(nond, _drift_step(dt=86400.0), mp)
            @test nond.grain_sphericity[1] <= 1.0
            @test nond.grain_radius[1] >= 0.5 * GEMB.DRIFT_GRAIN_SIZE_MIN * 1e3
        end
        @test nond.grain_sphericity[1] ≈ 1.0 atol = 1e-6
        @test nond.grain_radius[1] ≈ 0.5 * GEMB.DRIFT_GRAIN_SIZE_MIN * 1e3 atol = 1e-12
    end

    @testset "sublimation removes mass, capped at half the surface cell" begin
        mp = initialize_parameters(blowing_snow_method=:Crocus,
            blowing_snow_sublimation=true)
        cols = _drift_cols()
        mass_before = GEMB.column_mass(cols)
        m, e = GEMB.apply_blowing_snow!(cols, cfs, mp)

        # Signed as a flux into the column, so a loss is negative.
        @test m < 0.0
        @test e < 0.0
        # The mass reported is the mass the column lost: compaction is exactly conservative, so
        # the whole difference is sublimation.
        @test GEMB.column_mass(cols) - mass_before ≈ m rtol = 1e-13

        # No latent heat: the energy is exactly the departing snow's own enthalpy at the
        # surface temperature. `LS` is deliberately *not* charged — see the docstring note.
        @test e ≈ m * GEMB.specific_enthalpy(mp, 260.0) rtol = 1e-14
        # Which is a materially different number from charging it, so the test has teeth.
        @test !isapprox(e, m * (GEMB.specific_enthalpy(mp, 260.0) + GEMB.LS); rtol=1e-3)

        # The half-cell cap. A long step at high wind would otherwise empty the cell.
        capped = _drift_cols()
        GEMB.apply_blowing_snow!(capped, _drift_step(dt=1e7, wind_speed=40.0), mp)
        # Compaction then rescales `dz` and `density` together, so the invariant to check is
        # the surface cell's *mass*, not its thickness: at most half of it may leave.
        @test capped.dz[1] * capped.density[1] ≈ 0.5 * 0.1 * 150.0 rtol = 1e-13
    end

    @testset "sublimation off is bit-identical to a zero sublimation rate" begin
        # The `ZQS_EFFECT` feedback term is proportional to the rate, so a saturated atmosphere
        # (rate exactly 0) with sublimation *on* must reproduce the compaction-only path
        # bit-for-bit. If it did not, `blowing_snow_sublimation` would be changing the
        # compaction physics rather than only adding a sink.
        cols_off = _drift_cols()
        GEMB.apply_blowing_snow!(cols_off, cfs,
            initialize_parameters(blowing_snow_method=:Crocus,
                blowing_snow_sublimation=false))

        cols_on = _drift_cols()
        m, e = GEMB.apply_blowing_snow!(cols_on,
            _drift_step(vapor_pressure=GEMB.saturation_vapor_pressure_ice(253.0)),
            initialize_parameters(blowing_snow_method=:Crocus,
                blowing_snow_sublimation=true))
        @test (m, e) === (0.0, 0.0)
        @test cols_on.density == cols_off.density
        @test cols_on.dz == cols_off.dz
        @test cols_on.grain_radius == cols_off.grain_radius
    end
end

@testset "apply_prescribed_drift!" begin
    mp = initialize_parameters()
    yr = GEMB.SECONDS_PER_YEAR

    @testset "zero rate is an exact no-op" begin
        cols = _drift_cols()
        before = deepcopy(cols)
        m, e = GEMB.apply_prescribed_drift!(cols, _drift_step(snow_drift=0.0), mp)
        @test (m, e) === (0.0, 0.0)
        for k in keys(cols)
            @test cols[k] == before[k]
        end
    end

    @testset "erosion removes mass at the surface cell's own density" begin
        # Below the cap, the mass removed is exactly the prescribed flux, and the cell's
        # density, temperature and grains are untouched — only its thickness changes.
        cols = _drift_cols()
        flux = 100.0 * 10800.0 / yr
        m, _ = GEMB.apply_prescribed_drift!(cols, _drift_step(snow_drift=100.0), mp)
        @test m ≈ -flux rtol = 1e-14
        @test cols.dz[1] ≈ 0.1 - flux / 150.0 rtol = 1e-14
        @test cols.density[1] == 150.0
        @test cols.temperature[1] == 260.0
        @test cols.grain_radius[1] == 0.15
        @test cols.age[1] == 0.0                 # age-neutral: a fraction of the cell left

        # The half-cell cap. 100 kg m-2 yr-1 over a full year would take 0.667 m from a cell
        # holding 0.1 m, so only half of it (7.5 kg m-2) may go.
        big = _drift_cols()
        m_big, _ = GEMB.apply_prescribed_drift!(big, _drift_step(dt=yr, snow_drift=100.0), mp)
        @test m_big ≈ -7.5 rtol = 1e-14
        @test big.dz[1] ≈ 0.05 rtol = 1e-15

        # Pore water leaves in proportion to the thickness removed, and is in the reported mass.
        wet = _drift_cols(water=[2.0, 0.0, 0.0])
        m_wet, _ = GEMB.apply_prescribed_drift!(wet, _drift_step(dt=yr, snow_drift=100.0), mp)
        @test wet.water[1] ≈ 1.0 rtol = 1e-14    # half the cell left, so half the water
        @test m_wet ≈ -(7.5 + 1.0) rtol = 1e-14
    end

    @testset "deposition adds mass at the fresh-snow density" begin
        cols = _drift_cols()
        cfs = _drift_step(snow_drift=-100.0)
        flux = 100.0 * 10800.0 / yr
        mass_before = GEMB.column_mass(cols)
        m, e = GEMB.apply_prescribed_drift!(cols, cfs, mp)

        @test m ≈ flux rtol = 1e-14
        @test GEMB.column_mass(cols) - mass_before ≈ m rtol = 1e-12

        # It arrives at `fresh_snow_density`, IMAU-FDM's `rho0` — not at the cell's own
        # density, which is the asymmetry that makes the sign matter.
        rho_new = clamp(GEMB.fresh_snow_density(mp, cfs.temperature_air_mean,
                cfs.precipitation_mean, cfs.wind_speed_mean, cfs.temperature_air,
                cfs.wind_speed), 1.0, mp.density_ice)
        @test cols.dz[1] ≈ 0.1 + flux / rho_new rtol = 1e-14
        # And it mixes, so the cell's density lands between the old and the new.
        @test min(150.0, rho_new) <= cols.density[1] <= max(150.0, rho_new)

        # It arrives at the air temperature, which is colder here, so the cell cools toward it
        # without overshooting.
        @test cfs.temperature_air < cols.temperature[1] < 260.0

        # Age-zero mass dilutes the cell's age.
        aged = _drift_cols(age=[100.0, 100.0, 100.0])
        GEMB.apply_prescribed_drift!(aged, cfs, mp)
        @test 0.0 < aged.age[1] < 100.0

        # No latent heat here either, but for a different reason than sublimation: this is
        # transport, so the mass arrives as snow carrying only its own enthalpy.
        @test e ≈ m * GEMB.specific_enthalpy(mp, cfs.temperature_air) rtol = 1e-14
    end

    @testset "reported mass and energy match the column change" begin
        for drift in (100.0, -100.0)
            cols = _drift_cols()
            mass_before, energy_before = GEMB.column_mass_energy(cols, mp)
            m, e = GEMB.apply_prescribed_drift!(cols, _drift_step(snow_drift=drift), mp)
            mass_after, energy_after = GEMB.column_mass_energy(cols, mp)
            @test m ≈ mass_after - mass_before atol = 1e-12
            @test e ≈ energy_after - energy_before atol =
                GEMB.energy_tolerance(energy_before)
        end
    end

    @testset "erosion and deposition are opposite in sign" begin
        eroded = _drift_cols()
        m_e, _ = GEMB.apply_prescribed_drift!(eroded, _drift_step(snow_drift=50.0), mp)
        deposited = _drift_cols()
        m_d, _ = GEMB.apply_prescribed_drift!(deposited, _drift_step(snow_drift=-50.0), mp)
        @test m_e < 0.0 < m_d
        @test eroded.dz[1] < 0.1 < deposited.dz[1]
    end
end

@testset "validator bounds" begin
    @test_throws AssertionError initialize_parameters(blowing_snow_method=:Lenaerts)
    @test_throws AssertionError initialize_parameters(drift_rate=5000.0)
    @test_throws AssertionError initialize_parameters(drift_rate=-5000.0)
    # The two flags and the rate are independent, so every combination must construct.
    @test initialize_parameters(blowing_snow_method=:Crocus,
        blowing_snow_sublimation=true, drift_rate=100.0) isa GEMB.ModelParameters
end

@testset "snow_drift forcing layer" begin
    ds = GEMB_ClimateForcing.simulate_climate_forcing("test_1", 3)
    tdim = dims(ds, Ti)
    n = length(tdim)
    _with_drift(v) = DimStack((; (k => ds[k] for k in keys(ds))...,
            snow_drift=DimArray(fill(v, n), (tdim,))); metadata=metadata(ds))

    # Absent from the producer's DimStack, so the layer is present and zero — which is what
    # makes "the forcing carries no drift" and "the forcing carries zero drift" the same state,
    # and is why the driver can resolve the `mp.drift_rate` fallback by scanning for a nonzero.
    cf = initialize_forcing(ds)
    @test haskey(cf, :snow_drift)
    @test all(iszero, cf.snow_drift)

    # A conforming DimStack that *does* carry the layer is read through.
    cf_drift = initialize_forcing(_with_drift(25.0))
    @test all(==(25.0), cf_drift.snow_drift)

    # It survives climatological averaging, which a spinup depends on.
    @test all(≈(25.0), forcing_climatology(cf_drift).snow_drift)

    # Out-of-range values are caught at construction rather than silently applied.
    # 1e5 kg m-2 yr-1 is orders of magnitude past any observed drift divergence — the
    # signature of a per-timestep mass passed where an annual rate belongs.
    @test_throws AssertionError initialize_forcing(_with_drift(1e5))
end

@testset "blowing snow through gemb_core and gemb" begin
    ds = GEMB_ClimateForcing.simulate_climate_forcing("test_1", 3)
    tdim = dims(ds, Ti)
    n = length(tdim)
    # A six-year slice of the 32-year record: this testset makes eight full runs, so it uses a
    # shorter record than `test_horizontal_strain.jl` does to keep the suite's runtime in
    # proportion. Sliced by index rather than by date so the count is exact regardless of how
    # the synthetic producer lays out leap years.
    span = Ti(1:(8*Dates.value(Date(2000, 1, 1) - Date(1994, 1, 1))))
    cf = initialize_forcing(ds)[span]

    mp0 = initialize_parameters(output_frequency=:monthly)
    profile = initialize_profile(mp0, cf)
    z_target = sum(profile[:dz])
    n_cells = length(profile[:dz])

    # Off by default: the output layer is present and identically zero. `verbose=true`
    # exercises the per-step and whole-run mass budgets, both of which now carry the term.
    out0 = gemb(profile, cf, mp0; verbose=true)
    @test haskey(out0, :blowing_snow)
    @test all(iszero, out0[:blowing_snow])

    # The guarantee the three defaults rest on. Setting them explicitly to their defaults must
    # reproduce the default run bit-for-bit — any leak surfaces here rather than as an
    # unexplained move in `test_synthetic_regression.jl`.
    out_explicit = gemb(profile, cf,
        initialize_parameters(output_frequency=:monthly, blowing_snow_method=:none,
            blowing_snow_sublimation=false, drift_rate=0.0))
    for k in (:temperature, :density, :dz, :melt, :runoff, :refreeze, :blowing_snow)
        @test parent(out_explicit[k]) == parent(out0[k])
    end

    @testset "computed scheme: densifies the near surface, conserves mass" begin
        out = gemb(profile, cf,
            initialize_parameters(output_frequency=:monthly,
                blowing_snow_method=:Crocus); verbose=true)

        # Mass-conserving, so the reported surface flux stays exactly zero even though the
        # column has been reworked. This is the property that lets the compaction path stay
        # out of the budgets entirely.
        @test all(iszero, out[:blowing_snow])
        # The near-surface density rises toward the wind-slab value.
        @test sum(view(parent(out[:density]), 1, :)) >
              sum(view(parent(out0[:density]), 1, :))
        # And both grid invariants hold, as they must for any step that changes `dz`.
        @test size(out[:dz], 1) == n_cells
        for i in axes(out[:dz], 2)
            @test sum(view(out[:dz], :, i)) ≈ z_target atol = 1e-9
        end
    end

    @testset "sublimation: removes mass, budgets close" begin
        # `verbose=true` here is the real test of the budget wiring: omitting the term from
        # either the per-step or the whole-run check makes `gemb` error rather than fail an
        # assertion below.
        out = gemb(profile, cf,
            initialize_parameters(output_frequency=:monthly,
                blowing_snow_method=:Crocus, blowing_snow_sublimation=true); verbose=true)
        @test all(<=(0.0), out[:blowing_snow])      # never a gain
        @test sum(out[:blowing_snow]) < 0.0         # and a net loss over the run
        for i in axes(out[:dz], 2)
            @test sum(view(out[:dz], :, i)) ≈ z_target atol = 1e-9
        end
    end

    @testset "prescribed drift: both signs, budgets close" begin
        out_ero = gemb(profile, cf,
            initialize_parameters(output_frequency=:monthly, drift_rate=100.0);
            verbose=true)
        @test all(<=(0.0), out_ero[:blowing_snow])
        @test sum(out_ero[:blowing_snow]) < 0.0

        out_dep = gemb(profile, cf,
            initialize_parameters(output_frequency=:monthly, drift_rate=-100.0);
            verbose=true)
        @test all(>=(0.0), out_dep[:blowing_snow])
        @test sum(out_dep[:blowing_snow]) > 0.0

        for out in (out_ero, out_dep), i in axes(out[:dz], 2)
            @test sum(view(out[:dz], :, i)) ≈ z_target atol = 1e-9
        end
    end

    @testset "the forcing layer takes precedence over mp.drift_rate" begin
        cf_drift = initialize_forcing(DimStack(
            (; (k => ds[k] for k in keys(ds))...,
                snow_drift=DimArray(fill(100.0, n), (tdim,)));
            metadata=metadata(ds)))[span]

        # A conflicting scalar is ignored rather than added to the layer, so a real drift
        # product can never be perturbed by a stray parameter.
        out_layer = gemb(profile, cf_drift,
            initialize_parameters(output_frequency=:monthly); verbose=true)
        out_both = gemb(profile, cf_drift,
            initialize_parameters(output_frequency=:monthly, drift_rate=-500.0);
            verbose=true)
        @test parent(out_both[:blowing_snow]) == parent(out_layer[:blowing_snow])

        # And a constant layer is equivalent to the scalar fallback at the same rate, which is
        # what makes `mp.drift_rate` a genuine shorthand rather than a second code path.
        out_scalar = gemb(profile, cf,
            initialize_parameters(output_frequency=:monthly, drift_rate=100.0))
        @test parent(out_layer[:blowing_snow]) ≈ parent(out_scalar[:blowing_snow]) rtol = 1e-12
    end

    @testset "both paths together" begin
        out = gemb(profile, cf,
            initialize_parameters(output_frequency=:monthly,
                blowing_snow_method=:Crocus, blowing_snow_sublimation=true,
                drift_rate=100.0); verbose=true)
        @test sum(out[:blowing_snow]) < 0.0
        for i in axes(out[:dz], 2)
            @test sum(view(out[:dz], :, i)) ≈ z_target atol = 1e-9
        end
    end

    @testset "CF attributes are declared for the new layer" begin
        attrs = cf_attributes(:blowing_snow)
        @test attrs["units"] == "kg m-2"
        @test attrs["cell_methods"] == "time: sum"
        @test !isempty(attrs["long_name"])
        # No CF standard name fits this quantity, so none must be claimed.
        @test !haskey(attrs, "standard_name")
    end
end
