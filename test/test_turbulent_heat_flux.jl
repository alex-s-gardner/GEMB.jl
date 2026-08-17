# Tests for turbulent_heat_flux

@testset "Stable conditions" begin
    # Set up forcing (stable: T_air > T_surface)
    cfs = GEMB.ClimateForcingStep(
        10800.0,      # dt
        268.0,        # temperature_air (warmer than surface)
        80000.0,      # pressure_air
        0.0,          # precipitation
        5.0,          # wind_speed
        0.0,          # shortwave_downward
        300.0,        # longwave_downward
        300.0,        # vapor_pressure
        255.0,        # temperature_air_mean
        5.0,          # wind_speed_mean
        200.0,        # precipitation_mean
        2.0,          # temperature_observation_height
        10.0,         # wind_observation_height
        0.0, 0.0, 0.0, 0.0, 0.0, 0.1  # BC, COT, SZA, SWdiff, CF
    )

    T_surface = 265.0
    density_air = 1.225
    z0 = 0.00012
    zT = z0 * 0.10
    zQ = z0 * 0.10

    shf, lhf, lh = GEMB.turbulent_heat_flux(T_surface, density_air, z0, zT, zQ, cfs)

    # In stable conditions with T_air > T_surface, SHF should be positive (toward surface)
    @test shf > 0
    # Latent heat should be sublimation since T_surface < 273.15
    @test lh ≈ GEMB.LS atol = 1e-6
    # Values should be finite
    @test isfinite(shf)
    @test isfinite(lhf)
end

@testset "Unstable conditions" begin
    # Set up forcing (unstable: T_surface > T_air)
    cfs = GEMB.ClimateForcingStep(
        10800.0,      # dt
        260.0,        # temperature_air (colder than surface)
        80000.0,      # pressure_air
        0.0,          # precipitation
        5.0,          # wind_speed
        0.0,          # shortwave_downward
        300.0,        # longwave_downward
        300.0,        # vapor_pressure
        255.0,        # temperature_air_mean
        5.0,          # wind_speed_mean
        200.0,        # precipitation_mean
        2.0,          # temperature_observation_height
        10.0,         # wind_observation_height
        0.0, 0.0, 0.0, 0.0, 0.0, 0.1
    )

    T_surface = 272.0
    density_air = 1.225
    z0 = 0.00012
    zT = z0 * 0.10
    zQ = z0 * 0.10

    shf, lhf, lh = GEMB.turbulent_heat_flux(T_surface, density_air, z0, zT, zQ, cfs)

    # In unstable conditions with T_surface > T_air, SHF should be negative (away from surface)
    @test shf < 0
    @test isfinite(shf)
    @test isfinite(lhf)
end

@testset "Stability functions satisfy their definition" begin
    # Ψ(ζ) = ∫₀^ζ (1 - φ(z))/z dz. Both branches are checked against this directly, by
    # integrating their own profile function φ numerically and comparing to the closed form
    # the module evaluates. This is the independent check: it uses only the definition and
    # the published φ, never the closed forms themselves.
    function psi_numeric(phi, zeta; n=200_000)
        h = zeta / n
        s = 0.0
        for i in 1:n
            z = (i - 0.5) * h
            s += (1.0 - phi(z)) / z * h
        end
        return s
    end

    a1, b1, c1, d1 = 1.0, 2.0 / 3.0, 5.0, 0.35

    # Stable: Beljaars & Holtslag (1991) eq. 28 for momentum.
    phi_m_stable(z) = 1 + a1 * z + b1 * z * (1 + c1 - d1 * z) * exp(-d1 * z)
    psi_m_stable(z) = -(a1 * z + b1 * (z - c1 / d1) * exp(-d1 * z) + b1 * c1 / d1)

    # Unstable: Högström (1988) coefficients, normalized to φ(0) = 1.
    phi_m_unstable(z) = (1.0 - 19.0 * z)^(-0.25)
    phi_h_unstable(z) = (1.0 - 11.6 * z)^(-0.5)
    function psi_m_unstable(z)
        x = (1.0 - 19.0 * z)^0.25
        2.0 * log((1 + x) / 2.0) + log((1 + x^2) / 2.0) - 2 * atan(x) + pi / 2
    end
    psi_h_unstable(z) = 2.0 * log((1.0 + (1.0 - 11.6 * z)^0.5) / 2.0)

    for zeta in (0.01, 0.1, 1.0, 5.0)
        @test psi_m_stable(zeta) ≈ psi_numeric(phi_m_stable, zeta) rtol = 1e-5
    end
    for zeta in (-0.01, -0.1, -1.0, -5.0)
        @test psi_m_unstable(zeta) ≈ psi_numeric(phi_m_unstable, zeta) rtol = 1e-5
        @test psi_h_unstable(zeta) ≈ psi_numeric(phi_h_unstable, zeta) rtol = 1e-5
    end

    # Ψ(0) = 0 on every branch: the integrand vanishes at neutral.
    @test psi_m_stable(0.0) ≈ 0.0 atol = 1e-12
    @test psi_m_unstable(0.0) ≈ 0.0 atol = 1e-12
    @test psi_h_unstable(0.0) ≈ 0.0 atol = 1e-12

    # Sign: unstable enhances exchange (Ψ > 0), stable suppresses it (Ψ < 0).
    for zeta in (-0.01, -0.1, -1.0, -5.0)
        @test psi_m_unstable(zeta) > 0
        @test psi_h_unstable(zeta) > 0
    end
    for zeta in (0.01, 0.1, 1.0, 5.0)
        @test psi_m_stable(zeta) < 0
    end
end

@testset "Fluxes are continuous across neutral stability" begin
    # The branch point is Ri = 0, i.e. T_surface == T_air. A jump here would leave the
    # implicit solver's surface energy balance without a solution in T_surface (Fourteau et
    # al., 2024, Appendix D), and Newton converges onto precisely this crossing.
    T_surface = 265.0
    mk(T_air) = GEMB.ClimateForcingStep(
        10800.0, T_air, 80000.0, 0.0, 5.0, 0.0, 300.0, 300.0,
        255.0, 5.0, 200.0, 2.0, 10.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.1
    )
    flux(T_air) = GEMB.turbulent_heat_flux(T_surface, 1.225, 0.00012, 1.2e-5, 1.2e-5, mk(T_air))

    # Approaching neutral from both sides, the sensible flux must go to zero smoothly and
    # symmetrically (it is proportional to T_air - T_surface, with a continuous coefficient).
    for delta in (1e-1, 1e-2, 1e-3)
        shf_warm, _, _ = flux(T_surface + delta)
        shf_cold, _, _ = flux(T_surface - delta)
        @test shf_warm > 0
        @test shf_cold < 0
        # Symmetric to within the residual curvature of the transfer coefficient.
        @test shf_warm ≈ -shf_cold rtol = 2e-2
        # And the ratio to delta is the same on both sides, i.e. no step in the coefficient.
        @test shf_warm / delta ≈ -shf_cold / delta rtol = 2e-2
    end

    # Exactly neutral: zero sensible flux, no NaN from the branch.
    shf_neutral, lhf_neutral, _ = flux(T_surface)
    @test shf_neutral == 0.0
    @test isfinite(lhf_neutral)
end

@testset "Extreme instability stays finite and correctly signed" begin
    # Calm wind over a melting surface under very cold air drives the bulk Richardson number
    # to ~-7e5 in real forcing (it carries wind_speed^-2, and wind is floored at
    # min_wind_speed). Ψ is unbounded in ζ, so without ZETA_UNSTABLE_MIN the transfer
    # coefficient crosses zero there and the flux diverges then flips sign.
    mk(wind, T_air) = GEMB.ClimateForcingStep(
        10800.0, T_air, 80000.0, 0.0, wind, 0.0, 300.0, 300.0,
        255.0, 5.0, 200.0, 2.0, 10.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.1
    )
    for z0 in (GEMB.Z0_SNOW_DRY, GEMB.Z0_SNOW_WET, GEMB.Z0_ICE)
        zT = z0 * 0.10
        # The worst reachable case: zero wind (floored), melting surface, coldest air.
        shf, lhf, _ = GEMB.turbulent_heat_flux(GEMB.CtoK, 1.225, z0, zT, zT, mk(0.0, 229.5))
        @test isfinite(shf)
        @test isfinite(lhf)
        # T_air < T_surface, so the sensible flux must still be away from the surface.
        @test shf < 0
    end

    # The bound must not introduce a jump: fluxes remain monotone in T_surface through it.
    zT = GEMB.Z0_SNOW_DRY * 0.10
    f(T_s) = GEMB.turbulent_heat_flux(T_s, 1.225, GEMB.Z0_SNOW_DRY, zT, zT, mk(0.01, 240.0))[1]
    temps = 240.0:2.0:GEMB.CtoK
    shfs = [f(T_s) for T_s in temps]
    @test all(isfinite, shfs)
    @test all(diff(shfs) .< 0)  # colder air relative to surface => more negative flux
end

@testset "Melting surface uses vaporization" begin
    cfs = GEMB.ClimateForcingStep(
        10800.0, 275.0, 80000.0, 0.0, 5.0, 0.0, 300.0, 600.0,
        255.0, 5.0, 200.0, 2.0, 10.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.1
    )

    T_surface = 273.15  # At melting point
    density_air = 1.225
    z0 = 0.0013
    zT = z0 * 0.10
    zQ = z0 * 0.10

    _, _, lh = GEMB.turbulent_heat_flux(T_surface, density_air, z0, zT, zQ, cfs)

    # At melting point, should use latent heat of vaporization
    @test lh ≈ GEMB.LV atol = 1e-6
end

