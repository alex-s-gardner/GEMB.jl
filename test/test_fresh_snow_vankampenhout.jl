using GEMB
using Test

# `new_snow_method = :vanKampenhout2017`: van Kampenhout et al. (2017) eqs. 2-4, p. 2587.
#
# Expected values are computed here in BigFloat from the equations as printed, written out
# independently of the implementation:
#   eq. 2  ρ_fs = ρ_T + ρ_W
#   eq. 3  ρ_T  = 50 + 1.7·17^(3/2)                  if T > T_frz + 2
#                 50 + 1.7·(T − T_frz + 15)^(3/2)    if T_frz − 15 < T ≤ T_frz + 2
#                 −3.8328·(T − T_frz) − 0.0333·(T − T_frz)²   if T ≤ T_frz − 15
#   eq. 4  ρ_W  = 266.861·(½(1 + tanh(U/5)))^8.8, U the 10 m wind speed

const _VK_TFRZ = big"273.15"

function _vk_rho_T(T_kelvin)
    ΔT = big(T_kelvin) - _VK_TFRZ
    if ΔT > 2
        return big"50" + big"1.7" * big"17"^(big"3" / 2)
    elseif ΔT > -15
        return big"50" + big"1.7" * (ΔT + 15)^(big"3" / 2)
    else
        return big"-3.8328" * ΔT - big"0.0333" * ΔT^2
    end
end

_vk_rho_W(U) = big"266.861" * (big"0.5" * (1 + tanh(big(U) / 5)))^big"8.8"

@testset "vanKampenhout2017: temperature term, each branch" begin
    setprecision(BigFloat, 256) do
        # Middle branch, -10 °C: 50 + 1.7·5^1.5.
        @test GEMB.vankampenhout_density_temperature(263.15) ≈
              Float64(_vk_rho_T(big"263.15")) rtol = 1e-14
        @test Float64(_vk_rho_T(big"263.15")) ≈ 69.00657780874821 rtol = 1e-14
        # Middle branch, 0 °C: 50 + 1.7·15^1.5.
        @test GEMB.vankampenhout_density_temperature(273.15) ≈
              Float64(_vk_rho_T(big"273.15")) rtol = 1e-14
        # Warm branch, +5 °C: 50 + 1.7·17^1.5, and the same at +30 °C (T ≫ 0 cap).
        warm = Float64(big"50" + big"1.7" * big"17"^(big"3" / 2))
        @test GEMB.vankampenhout_density_temperature(278.15) ≈ warm rtol = 1e-14
        @test GEMB.vankampenhout_density_temperature(303.15) ==
              GEMB.vankampenhout_density_temperature(278.15)
        # Cold branch, -40 °C: 153.312 − 53.28.
        @test GEMB.vankampenhout_density_temperature(233.15) ≈
              Float64(_vk_rho_T(big"233.15")) rtol = 1e-13
        @test Float64(_vk_rho_T(big"233.15")) ≈ 100.032 rtol = 1e-14
        # Cold branch, -60 °C.
        @test GEMB.vankampenhout_density_temperature(213.15) ≈
              Float64(_vk_rho_T(big"213.15")) rtol = 1e-13
    end
end

@testset "vanKampenhout2017: temperature term at the branch joins" begin
    setprecision(BigFloat, 256) do
        # +2 °C: continuous. The middle branch at exactly +2 equals the warm constant.
        at2 = GEMB.vankampenhout_density_temperature(275.15)
        @test at2 ≈ Float64(big"50" + big"1.7" * big"17"^(big"3" / 2)) rtol = 1e-14
        @test GEMB.vankampenhout_density_temperature(275.15 + 1e-9) ≈ at2 rtol = 1e-12
        @test GEMB.vankampenhout_density_temperature(275.15 - 1e-9) ≈ at2 rtol = 1e-10

        # -15 °C: the cold branch applies at exactly -15 and gives 49.9995; the middle branch
        # tends to 50 from above. The published step is 5e-4 kg m-3.
        cold = _vk_rho_T(big"258.15")
        @test Float64(cold) ≈ 49.9995 rtol = 1e-14
        @test GEMB.vankampenhout_density_temperature(258.15) ≈ Float64(cold) rtol = 1e-12
        middle_limit = big"50"                     # 50 + 1.7·0^1.5
        step = middle_limit - cold
        @test Float64(step) ≈ 5.0e-4 rtol = 1e-10
        @test GEMB.vankampenhout_density_temperature(258.15 + 1e-6) -
              GEMB.vankampenhout_density_temperature(258.15) ≈ Float64(step) atol = 1e-6
    end
end

@testset "vanKampenhout2017: wind term" begin
    setprecision(BigFloat, 256) do
        for U in (0, 5, 10, 20)
            @test GEMB.vankampenhout_density_wind(Float64(U)) ≈ Float64(_vk_rho_W(U)) rtol = 1e-13
        end
        # Two reference values, 266.861·0.5^8.8 and 266.861·(½(1 + tanh 1))^8.8.
        @test Float64(_vk_rho_W(0)) ≈ 0.5987163900641868 rtol = 1e-13
        @test Float64(_vk_rho_W(5)) ≈ 87.33588449248248 rtol = 1e-13
        # Limit: approaches 266.861 from below.
        @test GEMB.vankampenhout_density_wind(60.0) < 266.861
        @test GEMB.vankampenhout_density_wind(60.0) ≈ 266.861 rtol = 1e-9
        # Non-decreasing in U.
        ρ = [GEMB.vankampenhout_density_wind(U) for U in 0.0:0.25:40.0]
        @test all(diff(ρ) .>= 0)
    end
end

@testset "vanKampenhout2017: fresh_snow_density and the 10 m wind" begin
    mp = GEMB.ModelParameters(new_snow_method=:vanKampenhout2017, densification_method=:Crocus)
    @test GEMB.validate_parameters(mp) === nothing

    setprecision(BigFloat, 256) do
        # Eq. 2 at -10 °C and 5 m s-1 (10 m wind; the default height is 10 m).
        expected = Float64(_vk_rho_T(big"263.15") + _vk_rho_W(5))
        @test GEMB.fresh_snow_density(mp, 250.0, 200.0, 3.0, 263.15, 5.0) ≈ expected rtol = 1e-13
        @test expected ≈ 156.34246230123068 rtol = 1e-13
    end

    # Reads the instantaneous arguments, not the climatological means.
    @test GEMB.fresh_snow_density(mp, 250.0, 200.0, 3.0, 263.15, 5.0) ==
          GEMB.fresh_snow_density(mp, 230.0, 900.0, 9.0, 263.15, 5.0)
    @test GEMB.fresh_snow_density(mp, 250.0, 200.0, 3.0, 263.15, 8.0) >
          GEMB.fresh_snow_density(mp, 250.0, 200.0, 3.0, 263.15, 5.0)

    # The 10 m conversion is the identity at a 10 m observation height, whatever z0.
    for z0 in (1.2e-4, 1.3e-3, 3.2e-3)
        @test GEMB.wind_speed_10m(7.3, 10.0, z0) == 7.3
    end
    # Below 10 m it is the neutral log profile U·ln(10/z0)/ln(z/z0), computed in BigFloat.
    setprecision(BigFloat, 256) do
        z0 = big"1.2e-4"
        U2 = big"4"
        U10 = U2 * log(10 / z0) / log(2 / z0)
        @test GEMB.wind_speed_10m(4.0, 2.0, 1.2e-4) ≈ Float64(U10) rtol = 1e-14
        @test GEMB.fresh_snow_density(mp, 250.0, 200.0, 3.0, 263.15, 4.0;
                  wind_observation_height=2.0, z0=1.2e-4) ≈
              Float64(_vk_rho_T(big"263.15") + _vk_rho_W(U10)) rtol = 1e-13
    end

    # Bounds over the forcing range: 49.9995 at -15 °C calm, below 436.1 everywhere.
    ρ = [GEMB.fresh_snow_density(mp, 250.0, 200.0, 3.0, T, U)
         for T in 193.15:1.0:283.15, U in 0.0:1.0:50.0]
    @test minimum(ρ) ≥ 49.99
    @test maximum(ρ) < 436.1

    # No wind: the deposition density is the temperature term plus 266.861·0.5^8.8.
    @test GEMB.fresh_snow_density(mp, 250.0, 200.0, 3.0, 263.15, 0.0) ≈
          GEMB.vankampenhout_density_temperature(263.15) + 266.861 * 0.5^8.8 rtol = 1e-14
end

@testset "vanKampenhout2017: fail fast" begin
    mp = GEMB.ModelParameters(new_snow_method=:vanKampenhout2017)
    @test_throws "wind speed must be finite and non-negative" GEMB.fresh_snow_density(
        mp, 250.0, 200.0, 3.0, 263.15, -0.1)
    @test_throws "wind speed must be finite and non-negative" GEMB.vankampenhout_density_wind(NaN)
    @test_throws "air temperature must be finite and positive" GEMB.vankampenhout_density_temperature(-1.0)
    @test_throws "air temperature must be finite and positive" GEMB.vankampenhout_density_temperature(Inf)
    @test_throws "roughness length must be positive" GEMB.wind_speed_10m(5.0, 10.0, 0.0)
    @test_throws "must exceed the roughness length" GEMB.wind_speed_10m(5.0, 1e-4, 1e-3)
    # An unrecognized method is an error, not a zero density.
    @test_throws "unrecognized new_snow_method" GEMB.fresh_snow_density(
        GEMB.ModelParameters(new_snow_method=:NotAMethod), 250.0, 200.0, 3.0)
end

@testset "vanKampenhout2017: calculate_accumulation opens a cell at the formula density" begin
    mp = GEMB.ModelParameters(new_snow_method=:vanKampenhout2017, densification_method=:Crocus,
        column_dzmin=0.05)
    n = 5
    # 2 m wind over a dry snow surface: calculate_accumulation passes the forcing height and
    # the roughness of the surface cell (dry, below ice density -> Z0_SNOW_DRY).
    cfs = GEMB.ClimateForcingStep(; dt=3600.0, temperature_air=265.0, pressure_air=80000.0,
        precipitation=30.0, wind_speed=6.0, vapor_pressure=200.0,
        temperature_air_mean=255.0, wind_speed_mean=4.0, precipitation_mean=500.0,
        temperature_observation_height=2.0, wind_observation_height=2.0)
    out = GEMB.calculate_accumulation(fill(260.0, n), fill(0.1, n), fill(300.0, n), zeros(n),
        fill(0.5, n), zeros(n), fill(0.5, n), zeros(n), cfs, mp, true)
    expected = GEMB.vankampenhout_density_temperature(265.0) +
               GEMB.vankampenhout_density_wind(GEMB.wind_speed_10m(6.0, 2.0, GEMB.Z0_SNOW_DRY))
    @test out[3][1] ≈ expected rtol = 1e-14
    @test out[2][1] ≈ 30.0 / expected rtol = 1e-14
    # Default fresh-grain properties.
    @test out[6][1] == GEMB.GDN_NEW_SNOW
    @test out[7][1] == GEMB.GSP_NEW_SNOW
end
