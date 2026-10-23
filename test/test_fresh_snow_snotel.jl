using GEMB
using Test

# `new_snow_method = :vanKampenhout2017SNOTEL`: the van Kampenhout et al. (2017) fresh-snow density with its
# temperature term replaced by 56 + 0.85·(clamp(ΔT, -8, 2) + 15)^1.5. Expected values are evaluated here in
# BigFloat from that form, independently of the implementation.
_snotel_rho_T(T_kelvin) = (ΔT = clamp(big(T_kelvin) - big"273.15", big"-8", big"2");
                           big"56" + big"0.85" * (ΔT + 15)^(big"3" / 2))

@testset "vanKampenhout2017SNOTEL: temperature term" begin
    for T in (243.15, 263.15, 265.15, 265.65, 268.15, 271.15, 273.15, 275.15, 290.0)
        @test GEMB.snotel_density_temperature(T) ≈ Float64(_snotel_rho_T(T)) rtol = 1e-14
    end
    # Level at and below -8 °C, held above +2 °C, continuous at both ends.
    floor_value = GEMB.snotel_density_temperature(265.15)
    @test floor_value ≈ 71.7422 rtol = 1e-5
    @test GEMB.snotel_density_temperature(213.15) === floor_value
    @test GEMB.snotel_density_temperature(265.15 + 1e-9) ≈ floor_value rtol = 1e-9
    @test GEMB.snotel_density_temperature(300.0) === GEMB.snotel_density_temperature(275.15)
    @test GEMB.snotel_density_temperature(275.15 - 1e-9) ≈ GEMB.snotel_density_temperature(275.15) rtol = 1e-9
    # Rises with temperature between the two limits.
    ρ = [GEMB.snotel_density_temperature(T) for T in 265.15:0.5:275.15]
    @test all(diff(ρ) .> 0)
    @test_throws "air temperature must be finite and positive" GEMB.snotel_density_temperature(NaN)
    @test_throws "air temperature must be finite and positive" GEMB.snotel_density_temperature(0.0)
end

@testset "vanKampenhout2017SNOTEL: fresh_snow_density" begin
    mp = GEMB.ModelParameters(new_snow_method=:vanKampenhout2017SNOTEL)
    @test_nowarn GEMB.validate_parameters(mp)
    # Wind term carried to 10 m as for :vanKampenhout2017.
    for (T, U, z) in ((263.15, 0.0, 10.0), (268.15, 6.0, 2.0), (272.15, 12.0, 3.5))
        ρW = GEMB.vankampenhout_density_wind(GEMB.wind_speed_10m(U, z, GEMB.Z0_SNOW_DRY))
        @test GEMB.fresh_snow_density(mp, 250.0, 200.0, 3.0, T, U; wind_observation_height=z) ≈
              GEMB.snotel_density_temperature(T) + ρW rtol = 1e-14
    end
end

@testset "vanKampenhout2017SNOTEL: calculate_accumulation opens a cell at the formula density" begin
    mp = GEMB.ModelParameters(new_snow_method=:vanKampenhout2017SNOTEL, densification_method=:Crocus,
        column_dzmin=0.05)
    n = 5
    cfs = GEMB.ClimateForcingStep(; dt=3600.0, temperature_air=262.0, pressure_air=80000.0,
        precipitation=30.0, wind_speed=2.0, vapor_pressure=200.0,
        temperature_air_mean=255.0, wind_speed_mean=4.0, precipitation_mean=500.0,
        temperature_observation_height=2.0, wind_observation_height=2.0)
    out = GEMB.calculate_accumulation(fill(260.0, n), fill(0.1, n), fill(300.0, n), zeros(n),
        fill(0.5, n), zeros(n), fill(0.5, n), zeros(n), cfs, mp, true)
    expected = GEMB.snotel_density_temperature(262.0) +
               GEMB.vankampenhout_density_wind(GEMB.wind_speed_10m(2.0, 2.0, GEMB.Z0_SNOW_DRY))
    @test out[3][1] ≈ expected rtol = 1e-14
    @test out[2][1] ≈ 30.0 / expected rtol = 1e-14
    @test out[6][1] == GEMB.GDN_NEW_SNOW
    @test out[7][1] == GEMB.GSP_NEW_SNOW
end

@testset "default snow set" begin
    mp = GEMB.ModelParameters()
    @test mp.new_snow_method === :vanKampenhout2017SNOTEL
    @test mp.densification_method === :Crocus
    @test mp.destructive_metamorphism === true
    @test mp.destructive_metamorphism_density === 100.0
    @test mp.grain_growth_method === :FlannerZender2006
    # Destructive metamorphism follows the densification method.
    @test GEMB.ModelParameters(densification_method=:Arthern).destructive_metamorphism === false
end
