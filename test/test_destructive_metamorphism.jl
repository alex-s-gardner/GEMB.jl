using GEMB
using Test

# `destructive_metamorphism = true`: Anderson (1976) settling as van Kampenhout et al. (2017)
# eq. 5, p. 2587 (c3 = 2.777e-6 s-1, c4 = 0.04 K-1, c2 = 2 with liquid water, ρ_dm = 175)
# with the taper of Schneider et al. (2022) eq. 3, p. 3 (c5 = 0.046 m3 kg-1), summed with
# the `:Crocus`/`:CrocusPure` viscous strain rate.
#
# Expected values are computed here in BigFloat from those equations, written out
# independently of the implementation.

function _dm_rate_big(ρ, T; wet=false, ρdm=big"175")
    c3 = big"2.777e-6"
    c4 = big"0.04"
    c5 = big"0.046"
    Tf = big"273.15"
    ρ = big(ρ)
    c1 = ρ < ρdm ? big"1" : exp(-c5 * (ρ - ρdm))
    c2 = wet ? big"2" : big"1"
    return c3 * c2 * c1 * exp(-c4 * (Tf - min(big(T), Tf)))
end

_dm(ρ, T; wet=false, ρdm=175.0) = GEMB.destructive_metamorphism_rate(ρ, T, wet, ρdm)

@testset "destructive metamorphism rate at stated states" begin
    setprecision(BigFloat, 256) do
        # Dry, -10 °C, below ρ_dm: c3·exp(-0.4).
        @test _dm(150.0, 263.15) ≈ Float64(_dm_rate_big(150, big"263.15")) rtol = 1e-14
        @test Float64(_dm_rate_big(150, big"263.15")) ≈ 1.86147877e-6 rtol = 1e-8
        # Dry, -50 °C: c3·exp(-2).
        @test _dm(150.0, 223.15) ≈ Float64(_dm_rate_big(150, big"223.15")) rtol = 1e-14
        # Above ρ_dm, 275 kg m-3: an extra factor exp(-0.046·100) = exp(-4.6).
        @test _dm(275.0, 263.15) ≈ Float64(_dm_rate_big(275, big"263.15")) rtol = 1e-14
        @test Float64(_dm_rate_big(275, big"263.15") / _dm_rate_big(150, big"263.15")) ≈
              Float64(exp(big"-4.6")) rtol = 1e-14
        # Liquid water doubles it: wet at the melting point is 2·c3.
        @test _dm(150.0, 273.15; wet=true) ≈ Float64(2 * big"2.777e-6") rtol = 1e-14
        @test _dm(300.0, 268.15; wet=true) ≈
              Float64(_dm_rate_big(300, big"268.15"; wet=true)) rtol = 1e-14
        @test _dm(300.0, 268.15; wet=true) == 2 * _dm(300.0, 268.15)
        # A configured ρ_dm moves the taper.
        @test _dm(250.0, 263.15; ρdm=300.0) == _dm(150.0, 263.15)
        @test _dm(350.0, 263.15; ρdm=300.0) ≈
              Float64(_dm_rate_big(350, big"263.15"; ρdm=big"300")) rtol = 1e-14
    end
end

@testset "destructive metamorphism rate: limits and shape" begin
    # Continuous at ρ_dm: c1 = exp(0) = 1 on the taper branch.
    @test _dm(175.0, 263.15) == _dm(174.999, 263.15)
    @test _dm(175.0001, 263.15) ≈ _dm(175.0, 263.15) rtol = 1e-5
    # Temperature is capped at T_frz.
    @test _dm(150.0, 280.0) == _dm(150.0, 273.15)
    # Non-increasing in density, non-decreasing in temperature up to T_frz.
    r_ρ = [_dm(ρ, 263.15) for ρ in 50.0:5.0:600.0]
    @test all(diff(r_ρ) .<= 0)
    r_T = [_dm(150.0, T) for T in 200.0:0.5:273.15]
    @test all(diff(r_T) .>= 0)
    # Fail fast on states the formula is not defined for.
    @test_throws "density must be positive" _dm(0.0, 263.15)
    @test_throws "temperature must be finite and positive" _dm(150.0, NaN)
end

# A one-cell column. Fresh dendritic grains (radius 0.05 mm) close eq. 9's gate, so f2 = 1,
# and a dry cell has f1 = 1, so the Vionnet et al. (2012) eq. 7 viscosity is
# η0·(ρ/c_η)·exp(a_η·(T_frz − T)
# + b_η·ρ) with η0 = 7.62237e6 kg s-1 m-1, a_η = 0.1 K-1, b_η = 0.023 m3 kg-1, c_η = 250
# kg m-3, and the overburden is half the cell's own weight, σ = ½·ρ·dz·g with g = 9.81.
function _one_cell(mp; ρ=150.0, dz=0.08, T=263.15, water=0.0, dt=3600.0)
    cfs = GEMB.ClimateForcingStep(; dt, temperature_air=T, pressure_air=80000.0,
        precipitation=0.0, wind_speed=0.0, vapor_pressure=100.0, temperature_air_mean=T,
        wind_speed_mean=0.0, precipitation_mean=300.0, temperature_observation_height=2.0,
        wind_observation_height=10.0)
    return GEMB.calculate_density([T], [dz], [ρ], [0.05], [water], cfs, mp;
        grain_dendricity=[1.0], grain_sphericity=[0.5])
end

function _one_cell_big(; ρ, dz, T, dt, dm::Bool)
    ρ, dz, T, dt = big(ρ), big(dz), big(T), big(dt)
    η = big"7.62237e6" * (ρ / 250) * exp(big"0.1" * (big"273.15" - T) + big"0.023" * ρ)
    σ = big"0.5" * ρ * dz * big"9.81"
    rate = σ / η + (dm ? _dm_rate_big(ρ, T) : big"0")
    dz_new = dz * exp(-rate * dt)
    return dz_new, ρ * dz / dz_new
end

@testset "one-step closed form, single layer" begin
    mp_on = GEMB.ModelParameters(densification_method=:CrocusPure, destructive_metamorphism=true,
        destructive_metamorphism_density=175.0)
    mp_off = GEMB.ModelParameters(densification_method=:CrocusPure, destructive_metamorphism=false)
    setprecision(BigFloat, 256) do
        for (ρ, T, dt) in ((150.0, 263.15, 3600.0), (100.0, 243.15, 10800.0),
                           (300.0, 268.15, 86400.0))
            dz_e, ρ_e = _one_cell_big(; ρ, dz=0.08, T, dt, dm=true)
            dz_o, ρ_o = _one_cell(mp_on; ρ, T, dt)
            @test dz_o[1] ≈ Float64(dz_e) rtol = 1e-12
            @test ρ_o[1] ≈ Float64(ρ_e) rtol = 1e-12
            # Mass conserved exactly enough to be invisible: ρ·dz is the input mass.
            @test ρ_o[1] * dz_o[1] ≈ ρ * 0.08 rtol = 1e-14

            dz_e0, ρ_e0 = _one_cell_big(; ρ, dz=0.08, T, dt, dm=false)
            dz_f, ρ_f = _one_cell(mp_off; ρ, T, dt)
            @test ρ_f[1] ≈ Float64(ρ_e0) rtol = 1e-12
            @test ρ_o[1] > ρ_f[1]
        end
        # The destructive term alone, isolated as the ratio of the two steps: exp(r·dt).
        _, ρ_on = _one_cell(mp_on)
        _, ρ_off = _one_cell(mp_off)
        @test ρ_on[1] / ρ_off[1] ≈ Float64(exp(_dm_rate_big(150, big"263.15") * 3600)) rtol = 1e-12
        # One hour at -10 °C below ρ_dm: dz × exp(-1.8615e-6·3600) = 0.99332.
        @test Float64(exp(-_dm_rate_big(150, big"263.15") * 3600)) ≈ 0.9933211 rtol = 1e-7
    end
    # No change over a zero-length step.
    dz0, ρ0 = _one_cell(mp_on; dt=0.0)
    @test dz0[1] == 0.08
    @test ρ0[1] == 150.0
end

@testset "under :Crocus, cells handed to the firn law are untouched" begin
    mp_on = GEMB.ModelParameters(densification_method=:Crocus, destructive_metamorphism=true,
        destructive_metamorphism_density=175.0)
    mp_off = GEMB.ModelParameters(densification_method=:Crocus, destructive_metamorphism=false)
    for ρ in (450.0, 600.0)
        @test _one_cell(mp_on; ρ) == _one_cell(mp_off; ρ)
    end
    # Below the handover it acts.
    @test _one_cell(mp_on; ρ=200.0)[2][1] > _one_cell(mp_off; ρ=200.0)[2][1]
    # Liquid water selects c2 = 2 (and f1 < 1 in the viscosity, so compare against the same
    # wet cell without the term).
    on_wet = _one_cell(mp_on; ρ=150.0, water=0.5, T=273.15)[2][1]
    off_wet = _one_cell(mp_off; ρ=150.0, water=0.5, T=273.15)[2][1]
    @test on_wet / off_wet ≈ exp(2 * 2.777e-6 * 3600.0) rtol = 1e-12
end

@testset "destructive metamorphism: switch is inert when off and guarded" begin
    # Off, the density threshold is unread: bit-identical outputs.
    mp_a = GEMB.ModelParameters(densification_method=:Crocus, destructive_metamorphism=false)
    mp_b = GEMB.ModelParameters(densification_method=:Crocus, destructive_metamorphism=false,
        destructive_metamorphism_density=300.0)
    @test _one_cell(mp_a) == _one_cell(mp_b)

    # Only the two viscous schemes accept it.
    @test GEMB.validate_parameters(GEMB.ModelParameters(densification_method=:Crocus,
        destructive_metamorphism=true)) === nothing
    @test_throws "requires densification_method :Crocus or :CrocusPure" GEMB.validate_parameters(
        GEMB.ModelParameters(densification_method=:Arthern, destructive_metamorphism=true))
    @test_throws "destructive_metamorphism_density must be in" GEMB.validate_parameters(
        GEMB.ModelParameters(destructive_metamorphism_density=20.0))
    # A ModelParameters built without validation fails at the first step instead of
    # silently dropping the term.
    @test_throws "destructive_metamorphism requires densification_method" _one_cell(
        GEMB.ModelParameters(densification_method=:Arthern, destructive_metamorphism=true))
end
