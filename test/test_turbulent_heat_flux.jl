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
    # Ψ(ζ) = ∫₀^ζ (1 - φ(z))/z dz. GEMB's closed forms are checked against this definition by
    # integrating the published profile functions φ numerically. The φ below are retyped from
    # the papers, not taken from GEMB, so this is the independent check; GEMB's own φ are then
    # checked against them.
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

    # Stable: Beljaars & Holtslag (1991) eqs. 28 (momentum) and 32 (heat).
    phi_m_stable(z) = 1 + a1 * z + b1 * z * (1 + c1 - d1 * z) * exp(-d1 * z)
    phi_h_stable(z) = 1 + z * (a1 * sqrt(1 + 2 * a1 * z / 3) + b1 * (1 + c1 - d1 * z) * exp(-d1 * z))

    # Unstable: Högström (1988) coefficients, normalized to φ(0) = 1.
    phi_m_unstable(z) = (1.0 - 19.0 * z)^(-0.25)
    phi_h_unstable(z) = (1.0 - 11.6 * z)^(-0.5)

    for zeta in (0.01, 0.1, 1.0, 5.0)
        @test GEMB.stability_psi_m(zeta) ≈ psi_numeric(phi_m_stable, zeta) rtol = 1e-5
        @test GEMB.stability_psi_h(zeta) ≈ psi_numeric(phi_h_stable, zeta) rtol = 1e-5
        @test GEMB.stability_phi_m(zeta) ≈ phi_m_stable(zeta) rtol = 1e-14
        @test GEMB.stability_phi_h(zeta) ≈ phi_h_stable(zeta) rtol = 1e-14
    end
    for zeta in (-0.01, -0.1, -1.0, -5.0)
        @test GEMB.stability_psi_m(zeta) ≈ psi_numeric(phi_m_unstable, zeta) rtol = 1e-5
        @test GEMB.stability_psi_h(zeta) ≈ psi_numeric(phi_h_unstable, zeta) rtol = 1e-5
        @test GEMB.stability_phi_m(zeta) ≈ phi_m_unstable(zeta) rtol = 1e-14
        @test GEMB.stability_phi_h(zeta) ≈ phi_h_unstable(zeta) rtol = 1e-14
    end

    # Ψ(0) = 0 on every branch: the integrand vanishes at neutral.
    for f in (GEMB.stability_psi_m, GEMB.stability_psi_h)
        @test f(0.0) ≈ 0.0 atol = 1e-12
        @test f(1e-12) ≈ 0.0 atol = 1e-10
        @test f(-1e-12) ≈ 0.0 atol = 1e-10
    end

    # Sign: unstable enhances exchange (Ψ > 0), stable suppresses it (Ψ < 0).
    for zeta in (-0.01, -0.1, -1.0, -5.0)
        @test GEMB.stability_psi_m(zeta) > 0
        @test GEMB.stability_psi_h(zeta) > 0
    end
    for zeta in (0.01, 0.1, 1.0, 5.0)
        @test GEMB.stability_psi_m(zeta) < 0
        @test GEMB.stability_psi_h(zeta) < 0
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


# ---------------------------------------------------------------------------------------
# Analytical checks: neutral limits, pressure invariance, continuity, and agreement with an
# iterated Monin-Obukhov solution.
# ---------------------------------------------------------------------------------------

_thf_cfs(; T_air, pressure=80000.0, wind=5.0, vapor=300.0, z_T=2.0, z_w=10.0) =
    GEMB.ClimateForcingStep(; dt=3600.0, temperature_air=T_air, pressure_air=pressure,
        wind_speed=wind, vapor_pressure=vapor, temperature_observation_height=z_T,
        wind_observation_height=z_w)

@testset "Neutral sensible coefficient" begin
    T_s, ΔT, U, p = 260.0, 1e-4, 5.0, 60000.0
    z0 = GEMB.Z0_SNOW_DRY
    zT = 0.1 * z0
    cfs = _thf_cfs(T_air=T_s + ΔT, pressure=p, wind=U)
    ρa = GEMB.air_density(p, T_s + ΔT)
    shf, _, _ = GEMB.turbulent_heat_flux(T_s, ρa, z0, zT, zT, cfs)
    expected = ρa * GEMB.HEAT_CAPACITY_AIR * 0.4^2 * U / (log(10.0 / z0) * log(2.0 / zT))
    @test shf / ΔT ≈ expected rtol = 1e-3
end

@testset "Neutral latent flux" begin
    for (T, L, es) in ((260.0, GEMB.LS, GEMB.saturation_vapor_pressure_ice),
                       (GEMB.CtoK, GEMB.LV, GEMB.saturation_vapor_pressure_water))
        e_a = 150.0
        z0 = GEMB.Z0_SNOW_WET
        zQ = 0.1 * z0
        cfs = _thf_cfs(T_air=T, wind=4.0, vapor=e_a)
        _, lhf, lh = GEMB.turbulent_heat_flux(T, 1.2, z0, zQ, zQ, cfs)
        @test lh == L
        @test lhf ≈ 0.4^2 * 4.0 * L * (e_a - es(T)) /
                    (461.9 * T * log(10.0 / z0) * log(2.0 / zQ)) rtol = 1e-4
    end
end

@testset "Sensible flux per unit air density is pressure-invariant" begin
    z0 = GEMB.Z0_SNOW_DRY
    zT = 0.1 * z0
    for ΔT in (-3.0, 0.0, 1e-4, 2.0)
        r = map((100000.0, 60000.0)) do p
            ρa = GEMB.air_density(p, 262.0)
            GEMB.turbulent_heat_flux(260.0, ρa, z0, zT, zT,
                _thf_cfs(T_air=260.0 + ΔT, pressure=p))[1] / ρa
        end
        @test r[1] ≈ r[2] rtol = 1e-6 atol = 1e-12
    end
end

@testset "Stable sensible flux is continuous in ΔT" begin
    # A scan in 0.01 K steps from neutral to 15 K: adjacent fluxes differ by less than 5%
    # once ΔT ≥ 0.2 K (below that a 0.01 K step is itself more than 5% of ΔT).
    T_s = 265.0
    z0 = GEMB.Z0_SNOW_DRY
    zT = 0.1 * z0
    for U in (2.0, 4.0, 6.0, 8.0), z_w in (2.0, 10.0)
        ΔTs = 0.0:0.01:15.0
        shf = [GEMB.turbulent_heat_flux(T_s, 1.2, z0, zT, zT,
            _thf_cfs(T_air=T_s + ΔT, wind=U, z_w=z_w))[1] for ΔT in ΔTs]
        @test all(isfinite, shf)
        @test shf[1] == 0.0
        @test all(>(0.0), shf[2:end])
        for k in eachindex(ΔTs)[2:end]
            ΔTs[k-1] >= 0.2 || continue
            @test abs(shf[k] - shf[k-1]) / shf[k-1] < 0.05
        end
    end
end

@testset "Agreement with an iterated Monin-Obukhov solution" begin
    # Independent reference: the Beljaars-Holtslag / Högström functions as published, the
    # Obukhov length found by bisection on z_w/L = Ri·F_m²/F_h, and the fluxes from u* and θ*.
    a, b, c, d = 1.0, 2.0 / 3.0, 5.0, 0.35
    ψm(ζ) = ζ > 0 ? -(a * ζ + b * (ζ - c / d) * exp(-d * ζ) + b * c / d) :
            (x = (1 - 19ζ)^0.25; 2log((1 + x) / 2) + log((1 + x^2) / 2) - 2atan(x) + π / 2)
    ψh(ζ) = ζ > 0 ? -((1 + 2a * ζ / 3)^1.5 + b * (ζ - c / d) * exp(-d * ζ) + b * c / d - 1) :
            2log((1 + sqrt(1 - 11.6ζ)) / 2)
    function reference_shf(T_s, T_a, U, ρa, z_w, z_T, z0, zT)
        Fm(ζ) = log(z_w / z0) - ψm(ζ) + ψm(ζ * z0 / z_w)
        Fh(ζ) = log(z_T / zT) - ψh(ζ * z_T / z_w) + ψh(ζ * zT / z_w)
        Ri = 9.81 * (T_a - T_s) * z_w / (0.5 * (T_a + T_s) * U^2)
        g(ζ) = ζ * Fh(ζ) - Ri * Fm(ζ)^2           # increasing in ζ, zero at the solution
        lo, hi = Ri >= 0 ? (0.0, 1.0) : (-1.0, 0.0)
        Ri >= 0 ? (while g(hi) < 0; hi *= 2; end) : (while g(lo) > 0; lo *= 2; end)
        for _ in 1:200
            m = 0.5 * (lo + hi)
            g(m) > 0 ? (hi = m) : (lo = m)
        end
        ζ = 0.5 * (lo + hi)
        u_star = 0.4 * U / Fm(ζ)
        θ_star = 0.4 * (T_a - T_s) / Fh(ζ)
        return ρa * GEMB.HEAT_CAPACITY_AIR * u_star * θ_star, ζ
    end

    n_checked = 0
    for z0 in (GEMB.Z0_SNOW_DRY, GEMB.Z0_SNOW_WET, GEMB.Z0_ICE), (z_T, z_w) in ((2.0, 10.0),
            (2.0, 2.0), (10.0, 2.0)), U in (1.0, 3.0, 8.0), ΔT in (-8.0, -2.0, -0.3, 0.3, 2.0, 8.0)
        zT = 0.1 * z0
        T_s = 263.0
        cfs = _thf_cfs(T_air=T_s + ΔT, wind=U, z_T=z_T, z_w=z_w)
        shf, _, _ = GEMB.turbulent_heat_flux(T_s, 1.2, z0, zT, zT, cfs)
        ref, ζ = reference_shf(T_s, T_s + ΔT, U, 1.2, z_w, z_T, z0, zT)
        (-100 < ζ < 1) || continue
        n_checked += 1
        @test shf ≈ ref rtol = 0.05
        # The fixed-step closure is in fact far tighter than the 5% requirement.
        @test shf ≈ ref rtol = 1e-8
    end
    @test n_checked > 100
end

@testset "Saturation vapor pressure over ice" begin
    # Murphy & Koop (2005) eq. 7 in closed form, 611.657 Pa at the triple point.
    T = 273.16
    @test GEMB.saturation_vapor_pressure_ice(T) ≈
          exp(9.550426 - 5723.265 / T + 3.53068 * log(T) - 0.00728332 * T) rtol = 1e-15
    @test GEMB.saturation_vapor_pressure_ice(273.16) ≈ 611.657 rtol = 1e-4
    # Against the IAPWS (2011) sublimation curve (Wagner et al., 2011, J. Phys. Chem. Ref.
    # Data 40, 043103), which Murphy & Koop's eq. 7 reproduces to better than 0.1%.
    function iapws(T)
        θ = T / 273.16
        a = (-0.212144006e2, 0.273203819e2, -0.610598130e1)
        b = (0.333333333e-2, 0.120666667e1, 0.170333333e1)
        return 611.657 * exp(sum(a[i] * θ^b[i] for i in 1:3) / θ)
    end
    for T in 200.0:5.0:273.15
        @test GEMB.saturation_vapor_pressure_ice(T) ≈ iapws(T) rtol = 1e-3
    end
end


@testset "Stability table against the converged Monin-Obukhov solution" begin
    # Independent reference: bisection on z_w/L = Ri·F_m²/F_h with the published functions.
    a, b, c, d = 1.0, 2.0 / 3.0, 5.0, 0.35
    ψm(ζ) = ζ > 0 ? -(a * ζ + b * (ζ - c / d) * exp(-d * ζ) + b * c / d) :
            (x = (1 - 19ζ)^0.25; 2log((1 + x) / 2) + log((1 + x^2) / 2) - 2atan(x) + π / 2)
    ψh(ζ) = ζ > 0 ? -((1 + 2a * ζ / 3)^1.5 + b * (ζ - c / d) * exp(-d * ζ) + b * c / d - 1) :
            2log((1 + sqrt(1 - 11.6ζ)) / 2)
    function reference(Ri, z_w, z_T, z0, zT)
        Fm(ζ) = log(z_w / z0) - ψm(ζ) + ψm(ζ * z0 / z_w)
        Fh(ζ) = log(z_T / zT) - ψh(ζ * z_T / z_w) + ψh(ζ * zT / z_w)
        g(ζ) = ζ * Fh(ζ) - Ri * Fm(ζ)^2
        lo, hi = Ri >= 0 ? (0.0, 1.0) : (-1.0, 0.0)
        Ri >= 0 ? (while g(hi) < 0; hi *= 2; end) : (while g(lo) > 0 && lo > -1e3; lo *= 2; end)
        for _ in 1:300
            m = 0.5 * (lo + hi)
            g(m) > 0 ? (hi = m) : (lo = m)
        end
        ζ = max(0.5 * (lo + hi), GEMB.ZETA_UNSTABLE_MIN)
        return ζ, Fm(ζ) * Fh(ζ)
    end
    worst_ζ = 0.0
    worst_F = 0.0
    worst_dG = 0.0
    for (z_w, z_T) in ((10.0, 2.0), (2.0, 2.0)), z0 in (GEMB.Z0_SNOW_DRY, GEMB.Z0_ICE)
        zT = 0.1 * z0
        tab = GEMB._stability_table(z_w, z_T, z0, zT)
        # A grid that falls between the table's nodes, spanning both signs and every regime.
        for sgn in (1.0, -1.0), q in range(log(1e-12), log(1e6); length=1777)
            Ri = sgn * exp(q)
            ζ_ref, FF_ref = reference(Ri, z_w, z_T, z0, zT)
            ζ = GEMB._stability_zeta(tab, Ri)
            G, dG = GEMB._stability_lookup(tab, Ri)
            worst_ζ = max(worst_ζ, abs(ζ / ζ_ref - 1))
            worst_F = max(worst_F, abs(exp(-G) * FF_ref - 1))
            # The tabulated derivative against a central difference of the reference.
            if abs(Ri) > 1e-6 && ζ_ref > GEMB.ZETA_UNSTABLE_MIN + 1
                δ = 1e-5 * abs(Ri)
                Gp = log(reference(Ri + δ, z_w, z_T, z0, zT)[2])
                Gm = log(reference(Ri - δ, z_w, z_T, z0, zT)[2])
                ref = (Gp - Gm) / (2δ) * Ri          # dG/d ln|Ri|
                worst_dG = max(worst_dG, abs(dG * Ri - ref) / max(1.0, abs(ref)))
            end
        end
    end
    @test worst_ζ < 1e-6
    @test worst_F < 1e-6
    @test worst_dG < 1e-6

    # The tabulated fluxes agree with the exact path of `turbulent_heat_flux`, and their
    # T_surface derivatives with a central difference.
    tab = GEMB._stability_table(10.0, 2.0, GEMB.Z0_SNOW_WET, 0.1 * GEMB.Z0_SNOW_WET)
    for (T_s, T_a, U) in ((260.0, 268.0, 3.0), (265.0, 255.0, 6.0), (GEMB.CtoK, 278.0, 0.5))
        cfs = _thf_cfs(T_air=T_a, wind=U)
        C = 0.16 * U
        ri_scale = GEMB.GRAVITY * 10.0 / U^2
        shf, lhf, L, dshf, dlhf = GEMB._table_turbulent_heat_flux(tab, T_s, 1.2, C, ri_scale, cfs)
        shf_e, lhf_e, L_e = GEMB.turbulent_heat_flux(T_s, 1.2, GEMB.Z0_SNOW_WET,
            0.1 * GEMB.Z0_SNOW_WET, 0.1 * GEMB.Z0_SNOW_WET, cfs)
        @test L == L_e
        @test shf ≈ shf_e rtol = 1e-6
        @test lhf ≈ lhf_e rtol = 1e-6
        T_d = min(T_s, GEMB.CtoK - 1e-3)
        f(T) = GEMB._table_turbulent_heat_flux(tab, T, 1.2, C, ri_scale, cfs)
        h = 1e-4
        @test f(T_d)[4] ≈ (f(T_d + h)[1] - f(T_d - h)[1]) / 2h rtol = 1e-5
        @test f(T_d)[5] ≈ (f(T_d + h)[2] - f(T_d - h)[2]) / 2h rtol = 1e-5
    end
end
