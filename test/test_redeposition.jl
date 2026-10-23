using GEMB
using Test

# `blowing_snow_method = :Keenan2021`: Keenan et al. (2021) eqs. 3-4, pp. 1067-1068, with the
# erosion threshold from the Crocus mobility index (Lafaysse et al., 2026, eqs. 59-62).
#
# Expected values are computed here in BigFloat from the equations as printed:
#   eq. 3  Φ = 0.0014·ρ_a·u*·(u* − u*_t)·(u* + 7.6·u*_t + 205)/L
#   eq. 4  ρ_rd = 361·log10(U₁₀) + 33 for U₁₀ > 1 m s-1, else 33
#   u* = κ·U₁₀/ln(10/z₀), κ = 0.4, z₀ = 2 mm (p. 1067, neutral)
#   mobility (dendritic) = 0.34·(0.5 + 0.75·d − 0.5·s) + 0.66·F(ρ),
#   F(ρ) = 1.25 − 1.25·(max(ρ, 50) − 50)/(1000·0.295),  U₅,t = −ln((MO + 1)/2.868)/0.085

_rd_big(U) = big(U) > 1 ? big"361" * log10(big(U)) + big"33" : big"33"
_ustar_big(U) = big"0.4" * big(U) / log(big"10" / big"0.002")
_flux_big(ρa, us, ust, L) = big"0.0014" * big(ρa) * big(us) * (big(us) - big(ust)) *
                            (big(us) + big"7.6" * big(ust) + big"205") / big(L)
_F_big(ρ) = big"1.25" - big"1.25" * (max(big(ρ), big"50") - 50) / (1000 * big"0.295")
_mob_dendritic_big(d, s, ρ) = big"0.34" * (big"0.5" + big"0.75" * big(d) - big"0.5" * big(s)) +
                              big"0.66" * _F_big(ρ)
_ut5_big(mob) = -log((mob + 1) / big"2.868") / big"0.085"

_rd_step(; wind_speed=10.0, dt=3600.0, temperature_air=253.0, pressure_air=80000.0,
    wind_observation_height=10.0) = GEMB.ClimateForcingStep(; dt, temperature_air,
    pressure_air, precipitation=0.0, wind_speed, vapor_pressure=50.0,
    temperature_air_mean=temperature_air, wind_speed_mean=wind_speed, precipitation_mean=200.0,
    temperature_observation_height=2.0, wind_observation_height)

_rd_cols(; dz=[0.05, 0.05, 0.1, 0.2], density=[150.0, 180.0, 250.0, 500.0],
    water=zeros(4), grain_radius=[0.05, 0.05, 0.2, 0.4], dendricity=[1.0, 1.0, 0.0, 0.0],
    sphericity=fill(0.5, 4)) = GEMB.column_state(fill(253.0, length(dz)), copy(dz),
    copy(density), copy(water), copy(grain_radius), copy(dendricity), copy(sphericity),
    zeros(length(dz)))

const _MP_RD = GEMB.ModelParameters(blowing_snow_method=:Keenan2021, grain_growth_method=:hybrid)

@testset "Keenan 2021 eq. 4: redeposited density" begin
    setprecision(BigFloat, 256) do
        for U in (0.5, 1.0, 5.0, 10.0, 15.0, 20.0)
            @test GEMB.redeposition_density(U) ≈ Float64(_rd_big(U)) rtol = 1e-14
        end
        @test GEMB.redeposition_density(0.5) == 33.0
        @test GEMB.redeposition_density(1.0) == 33.0
        @test Float64(_rd_big(5)) ≈ 285.3281715653028 rtol = 1e-14
        @test Float64(_rd_big(10)) ≈ 394.0 rtol = 1e-15
        @test Float64(_rd_big(20)) ≈ 502.6718284346972 rtol = 1e-14
    end
    # Continuous at 1 m s-1 and non-decreasing.
    @test GEMB.redeposition_density(1.0 + 1e-9) ≈ 33.0 atol = 1e-6
    ρ = [GEMB.redeposition_density(U) for U in 0.0:0.1:30.0]
    @test all(diff(ρ) .>= 0)
    @test_throws "wind speed must be finite and non-negative" GEMB.redeposition_density(-1.0)
end

@testset "Keenan 2021 eq. 3: friction velocity and saltation flux" begin
    setprecision(BigFloat, 256) do
        @test GEMB.saltation_friction_velocity(10.0) ≈ Float64(_ustar_big(10)) rtol = 1e-14
        @test Float64(_ustar_big(10)) ≈ 0.46963828459723833 rtol = 1e-14
        @test GEMB.saltation_mass_flux(1.0, 0.5, 0.3, 10.0) ≈
              Float64(_flux_big(1, big"0.5", big"0.3", 10)) rtol = 1e-14
        @test Float64(_flux_big(1, big"0.5", big"0.3", 10)) ≈ 2.90892e-3 rtol = 1e-12
        @test GEMB.saltation_mass_flux(0.9, 0.62, 0.21, 3.0) ≈
              Float64(_flux_big(big"0.9", big"0.62", big"0.21", 3)) rtol = 1e-14
    end
    # Zero at and below the threshold; inversely proportional to the fetch length.
    @test GEMB.saltation_mass_flux(1.0, 0.3, 0.3, 10.0) == 0.0
    @test GEMB.saltation_mass_flux(1.0, 0.2, 0.3, 10.0) == 0.0
    @test GEMB.saltation_mass_flux(1.0, 0.5, 0.3, 5.0) ≈
          2 * GEMB.saltation_mass_flux(1.0, 0.5, 0.3, 10.0) rtol = 1e-15
    @test_throws "fetch length must be positive" GEMB.saltation_mass_flux(1.0, 0.5, 0.3, 0.0)
    @test_throws "friction velocities must be non-negative" GEMB.saltation_mass_flux(1.0, -0.5, 0.3, 10.0)
end

@testset "erosion threshold from the Crocus mobility index" begin
    setprecision(BigFloat, 256) do
        # Fresh dendritic snow, d = 1, s = 0.5, at 180 kg m-3.
        mob_big = _mob_dendritic_big(1, big"0.5", 180)
        mob = GEMB.mobility_index(0.05, 1.0, 0.5, 180.0, false)
        @test mob ≈ Float64(mob_big) rtol = 1e-14
        @test Float64(mob_big) ≈ 0.8014406779661016 rtol = 1e-14
        @test GEMB.redeposition_threshold_wind_speed(mob) ≈ Float64(_ut5_big(mob_big)) rtol = 1e-14
        @test Float64(_ut5_big(mob_big)) ≈ 5.470920015505936 rtol = 1e-13
    end
    # Dense snow never drifts; the formula's undefined range fails.
    @test GEMB.redeposition_threshold_wind_speed(-1.0) == Inf
    @test GEMB.redeposition_threshold_wind_speed(GEMB.mobility_index(0.4, 0.0, 0.5, 900.0, false)) == Inf
    @test_throws "outside the range of the drift threshold" GEMB.redeposition_threshold_wind_speed(1.9)
end

@testset "one-cell update, closed form" begin
    setprecision(BigFloat, 256) do
        # The cell update of eq. 4: 10 kg m-2 at 150 kg m-3 with 4 kg m-2 reworked at
        # ρ_rd(10 m s-1) = 394: dz = 6/150 + 4/394, ρ = 10/dz.
        dz = 6 / big"150" + 4 / _rd_big(10)
        @test Float64(dz) ≈ 0.05015228426395939 rtol = 1e-14
        @test Float64(10 / dz) ≈ 199.39271255060729 rtol = 1e-13

        # Through `apply_redeposition!`: one thick fresh-snow cell over firn that does not
        # drift, so all of Φ·dt comes from cell 1. Every quantity written out in BigFloat.
        U = big"9"
        T = big"253"
        p = big"80000"
        dt = big"3600"
        z0 = big"0.00012"                            # dry snow, Z0_SNOW_DRY
        U5 = big"1.25" * U * log(5 / z0) / log(10 / z0)
        ut5 = _ut5_big(_mob_dendritic_big(1, big"0.5", 150))
        @test U5 > ut5
        us = _ustar_big(U)
        ust = us * ut5 / U5
        ρa = big"0.029" * p / (big"8.314" * T)
        M = _flux_big(ρa, us, ust, 10) * dt
        mass = big"150" * big"0.5"
        @test M < mass
        dz_e = (mass - M) / 150 + M / _rd_big(U)
        ρ_e = mass / dz_e

        cols = _rd_cols(; dz=[0.5, 0.2], density=[150.0, 600.0], water=zeros(2),
            grain_radius=[0.05, 0.4], dendricity=[1.0, 0.0], sphericity=[0.5, 0.5])
        out = GEMB.apply_redeposition!(cols, _rd_step(; wind_speed=9.0), _MP_RD)
        @test out == (0.0, 0.0)
        @test cols.dz[1] ≈ Float64(dz_e) rtol = 1e-12
        @test cols.density[1] ≈ Float64(ρ_e) rtol = 1e-12
        @test cols.density[1] * cols.dz[1] ≈ 75.0 rtol = 1e-14
        # Microstructure of the reworked fraction (Wever et al., 2023, Table 1).
        f = Float64(M / mass)
        @test cols.grain_dendricity[1] ≈ (1 - f) * 1.0 + f * 0.325 rtol = 1e-12
        @test cols.grain_sphericity[1] ≈ (1 - f) * 0.5 + f * 0.875 rtol = 1e-12
        @test cols.grain_radius[1] ≈
              GEMB.dendritic_grain_radius(cols.grain_dendricity[1], cols.grain_sphericity[1])
        # The firn below does not drift and is untouched.
        @test cols.dz[2] == 0.2 && cols.density[2] == 600.0
    end
end

@testset "redeposition: limits and invariants" begin
    # No wind: no packing.
    c0 = _rd_cols()
    c = _rd_cols()
    GEMB.apply_redeposition!(c, _rd_step(; wind_speed=0.0), _MP_RD)
    @test c == c0
    # Below the fresh-snow threshold (U₁₀ ≈ 4.7 m s-1 over dry snow): no packing.
    c = _rd_cols()
    GEMB.apply_redeposition!(c, _rd_step(; wind_speed=4.0), _MP_RD)
    @test c == c0
    # A wet surface cell does not erode.
    c = _rd_cols(; water=[0.1, 0.0, 0.0, 0.0])
    cw = _rd_cols(; water=[0.1, 0.0, 0.0, 0.0])
    GEMB.apply_redeposition!(c, _rd_step(; wind_speed=15.0), _MP_RD)
    @test c == cw

    # Windy: mass, temperature, water and age are unchanged in every cell, so column mass and
    # enthalpy are conserved; density never falls; the reworked mass is at most Φ·dt.
    for U in (6.0, 10.0, 15.0, 25.0)
        c = _rd_cols()
        cfs = _rd_step(; wind_speed=U)
        GEMB.apply_redeposition!(c, cfs, _MP_RD)
        @test c.density .* c.dz ≈ c0.density .* c0.dz rtol = 1e-14
        @test c.temperature == c0.temperature && c.water == c0.water && c.age == c0.age
        @test all(c.density .>= c0.density)
        @test all(c.density .<= max.(c0.density, GEMB.redeposition_density(U)) .+ 1e-9)
        # Reworked mass of cell i is m·(1/ρ0 − 1/ρ1)/(1/ρ0 − 1/ρ_rd) wherever it densified.
        ρrd = GEMB.redeposition_density(U)
        m = c0.density .* c0.dz
        reworked = sum(m[i] * (1 / c0.density[i] - 1 / c.density[i]) /
                       (1 / c0.density[i] - 1 / ρrd)
                       for i in eachindex(m) if c.density[i] > c0.density[i]; init=0.0)
        z0 = GEMB.Z0_SNOW_DRY
        U5 = GEMB.drift_wind_speed(cfs, z0)
        ut5 = GEMB.redeposition_threshold_wind_speed(GEMB.mobility_index(0.05, 1.0, 0.5, 150.0, false))
        us = GEMB.saltation_friction_velocity(U)
        Φ = GEMB.saltation_mass_flux(GEMB.air_density(80000.0, 253.0), us, us * ut5 / U5, 10.0)
        @test reworked <= Φ * cfs.dt * (1 + 1e-12)
        # Erosion stops at the first cell whose own threshold the gust speed does not exceed.
        # The 500 kg m-3 firn cell (non-dendritic, 0.4 mm) has a threshold of 19.4 m s-1.
        ut5_firn = GEMB.redeposition_threshold_wind_speed(
            GEMB.mobility_index(0.4, 0.0, 0.5, 500.0, false))
        @test ut5_firn ≈ 19.4 atol = 0.05
        U5 > ut5_firn || @test c.density[4] == 500.0
    end

    # Density non-decreasing in wind speed at a fixed state.
    ρ1 = map(0.0:0.5:30.0) do U
        c = _rd_cols()
        GEMB.apply_redeposition!(c, _rd_step(; wind_speed=U), _MP_RD)
        c.density[1]
    end
    @test all(diff(ρ1) .>= 0)
    @test ρ1[end] > 150.0

    # A cell denser than ρ_rd is reworked (its grains change) but keeps its density. A
    # dendritic cell at 400 kg m-3 drifts above U₅ = 10.4 m s-1; at U₁₀ = 9.5 m s-1 the gust
    # speed is 11.1 and ρ_rd = 386.
    c = _rd_cols(; density=[400.0, 180.0, 250.0, 500.0])
    GEMB.apply_redeposition!(c, _rd_step(; wind_speed=9.5), _MP_RD)
    @test GEMB.redeposition_density(9.5) < 400.0
    @test c.density[1] == 400.0 && c.dz[1] == 0.05
    @test c.grain_dendricity[1] < 1.0

    # A non-10 m observation height goes through the 10 m profile.
    c_a = _rd_cols()
    c_b = _rd_cols()
    GEMB.apply_redeposition!(c_a, _rd_step(; wind_speed=10.0), _MP_RD)
    U2 = 10.0 * log(2.0 / GEMB.Z0_SNOW_DRY) / log(10.0 / GEMB.Z0_SNOW_DRY)
    GEMB.apply_redeposition!(c_b, _rd_step(; wind_speed=U2, wind_observation_height=2.0), _MP_RD)
    @test c_b.density ≈ c_a.density rtol = 1e-12

    @test_throws "wind speed must be finite and non-negative" GEMB.apply_redeposition!(
        _rd_cols(), _rd_step(; wind_speed=-1.0), _MP_RD)
end

@testset "Keenan2021 selection and validation" begin
    # `apply_blowing_snow!` dispatches to it and reports no surface flux.
    c_a = _rd_cols()
    c_b = _rd_cols()
    @test GEMB.apply_blowing_snow!(c_a, _rd_step(; wind_speed=12.0), _MP_RD) == (0.0, 0.0)
    GEMB.apply_redeposition!(c_b, _rd_step(; wind_speed=12.0), _MP_RD)
    @test c_a == c_b
    # It replaces the Crocus compaction rather than adding to it.
    c_c = _rd_cols()
    GEMB.apply_blowing_snow!(c_c, _rd_step(; wind_speed=12.0),
        GEMB.ModelParameters(blowing_snow_method=:Crocus))
    @test c_c != c_a

    @test GEMB.validate_parameters(_MP_RD) === nothing
    @test_throws "cannot be combined with :Keenan2021" GEMB.validate_parameters(
        GEMB.ModelParameters(blowing_snow_method=:Keenan2021, blowing_snow_sublimation=true))
    @test_throws "drift_fetch_length must be in" GEMB.validate_parameters(
        GEMB.ModelParameters(blowing_snow_method=:Keenan2021, drift_fetch_length=0.0))
    @test_throws "unrecognized blowing_snow_method" GEMB.apply_blowing_snow!(_rd_cols(),
        _rd_step(), GEMB.ModelParameters(blowing_snow_method=:Bogus))
    # The fetch length scales the reworked mass.
    c_L3 = _rd_cols()
    GEMB.apply_redeposition!(c_L3, _rd_step(; wind_speed=12.0),
        GEMB.ModelParameters(blowing_snow_method=:Keenan2021, drift_fetch_length=3.0))
    @test sum(c_L3.dz) < sum(c_b.dz)
end
