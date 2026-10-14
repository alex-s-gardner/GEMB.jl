"""
    turbulent_heat_flux(T_surface, density_air, z0, zT, zQ, cfs::ClimateForcingStep; min_wind_speed=0.01)

Compute sensible and latent heat fluxes using Monin-Obukhov similarity theory.

`z0`, `zT` and `zQ` are the roughness lengths for momentum, heat and moisture [m]. Wind is
observed at `cfs.wind_observation_height` and temperature and humidity at
`cfs.temperature_observation_height`; the two may differ (see [`obukhov_stability`](@ref)).

`min_wind_speed` floors `cfs.wind_speed` (which would otherwise drive the bulk
Richardson number to infinity at zero wind); clamping here avoids rebuilding the
`ClimateForcingStep` in the caller.

Returns (heat_flux_sensible, heat_flux_latent, latent_heat) [W m-2, W m-2, J kg-1].
"""
function turbulent_heat_flux(T_surface::Float64, density_air::Float64,
    z0::Float64, zT::Float64, zQ::Float64, cfs::ClimateForcingStep;
    min_wind_speed::Float64=0.01)

    wind_speed = max(cfs.wind_speed, min_wind_speed)

    # Neutral bulk-transfer factor κ²U
    C = VON_KARMAN^2 * wind_speed

    # Neutral roughness log terms (invariant in T_surface).
    logM = log(cfs.wind_observation_height / z0)
    logHT = log(cfs.temperature_observation_height / zT)
    logHQ = log(cfs.temperature_observation_height / zQ)

    # g·z_w/U², the part of the bulk Richardson number that does not depend on T_surface.
    ri_scale = GRAVITY * cfs.wind_observation_height / wind_speed^2

    return _turbulent_heat_flux(T_surface, density_air, z0, zT, zQ, cfs,
        wind_speed, C, logM, logHT, logHQ, ri_scale)
end

"""
    _turbulent_heat_flux(T_surface, density_air, z0, zT, zQ, cfs, wind_speed, C, logM, logHT, logHQ, ri_scale)

Core of [`turbulent_heat_flux`](@ref) with the sub-timestep-invariant quantities
(`wind_speed`, the neutral factor `C = κ²U`, the neutral roughness logs `logM`/`logHT`/`logHQ`,
and `ri_scale = g·z_w/U²`) hoisted to the caller. `calculate_temperature` computes these once
per timestep and reuses them across all thermal sub-steps, where only `T_surface` changes.

## Fluxes

With `ζ = z_w/L` from [`obukhov_stability`](@ref), the integrated transfer functions are

    F_m = ln(z_w/z0) − Ψ_m(ζ) + Ψ_m(ζ·z0/z_w)
    F_h = ln(z_T/zT) − Ψ_h(ζ·z_T/z_w) + Ψ_h(ζ·zT/z_w)
    F_q = ln(z_T/zQ) − Ψ_h(ζ·z_T/z_w) + Ψ_h(ζ·zQ/z_w)

and the fluxes are

    H  = ρ_a·c_p·κ²U·(T_a − T_s)/(F_m·F_h)
    LE = κ²U·L·(e_a − e_s(T_s))/(R_v·T̄·F_m·F_q),    T̄ = (T_a + T_s)/2

The potential temperature is referenced to the surface pressure, so `θ_a − θ_s = T_a − T_s`
and no pressure factor enters either flux or the Richardson number. The dry-adiabatic
difference over the observation height, `g·z_T/c_p` (0.02 K at 2 m), is neglected.

`L` is `LS` with saturation over ice below the melting point and `LV` with saturation over
water at it.

# References
- Monin, A. S. and Obukhov, A. M. (1954). Basic laws of turbulent mixing in the surface
  layer of the atmosphere. *Tr. Akad. Nauk SSSR Geofiz. Inst.* 24, 163-187.
"""
@inline function _turbulent_heat_flux(T_surface::Float64, density_air::Float64,
    z0::Float64, zT::Float64, zQ::Float64, cfs::ClimateForcingStep,
    wind_speed::Float64, C::Float64,
    logM::Float64, logHT::Float64, logHQ::Float64, ri_scale::Float64)

    z_w = cfs.wind_observation_height
    z_T = cfs.temperature_observation_height
    T_mean = 0.5 * (cfs.temperature_air + T_surface)

    # Bulk Richardson number with the wind height as reference: Ri_b = g·Δθ·z_w/(θ̄·U²).
    Ri = ri_scale * (cfs.temperature_air - T_surface) / T_mean

    ζ = obukhov_stability(Ri, z_w, z_T, z0, zT, logM, logHT)

    psi_h_T = stability_psi_h(ζ * z_T / z_w)
    coefM = logM - stability_psi_m(ζ) + stability_psi_m(ζ * z0 / z_w)
    coefHT = logHT - psi_h_T + stability_psi_h(ζ * zT / z_w)
    # GEMB's scalar roughness lengths are equal, which makes the two coefficients identical.
    coefHQ = zQ == zT ? coefHT : logHQ - psi_h_T + stability_psi_h(ζ * zQ / z_w)

    # Sensible Heat Flux [W m-2]
    heat_flux_sensible = density_air * C * HEAT_CAPACITY_AIR *
                         (cfs.temperature_air - T_surface) / (coefM * coefHT)

    # Latent Heat Flux [W m-2]
    if T_surface >= CtoK - T_TOLERANCE
        # Liquid water surface
        latent_heat = LV
        # Saturation Vapor Pressure (Murray 1967)
        eS = saturation_vapor_pressure_water(T_surface)
    else
        # Ice surface
        latent_heat = LS
        # Saturation Vapor Pressure (Murray 1967)
        eS = saturation_vapor_pressure_ice(T_surface)
    end

    # 461.9 is the specific gas constant for water vapor [J kg-1 K-1]
    heat_flux_latent = C * latent_heat * (cfs.vapor_pressure - eS) /
                       (461.9 * T_mean) / (coefM * coefHQ)

    return heat_flux_sensible, heat_flux_latent, latent_heat
end

"""
    obukhov_stability(Ri, z_w, z_T, z0, zT, logM=log(z_w/z0), logHT=log(z_T/zT)) -> ζ

Stability parameter `ζ = z_w/L` at the wind observation height `z_w`, from the bulk Richardson
number `Ri = g·Δθ·z_w/(θ̄·U²)` (Δθ between the temperature height `z_T` and the surface, `U`
at `z_w`).

Monin-Obukhov similarity with `u* = κU/F_m(ζ)` and `θ* = κΔθ/F_h(ζ)` (see
[`_turbulent_heat_flux`](@ref) for `F_m` and `F_h`, which evaluate `Ψ_m` at `z_w/L` and `Ψ_h`
at `z_T/L`) gives `z_w/L = Ri·F_m(ζ)²/F_h(ζ)`. This function solves that relation for `ζ` by
`MO_NEWTON_STEPS` Newton steps on `ln|ζ|`,

    R(s) = s + ln F_h − 2 ln F_m − ln|Ri| = 0,    s = ln|ζ|,
    dR/ds = 1 + (φ_h(ζ z_T/z_w) − φ_h(ζ zT/z_w))/F_h − 2(φ_m(ζ) − φ_m(ζ z0/z_w))/F_m,

starting from the neutral estimate `ζ₀ = Ri·ln(z_w/z0)²/ln(z_T/zT)`. The number of steps is
fixed and there is no convergence test, so `ζ` is a smooth function of `Ri` with the sign of
`Ri`, zero at `Ri = 0`, on both sides of neutral. Against the converged solution the transfer
coefficient `1/(F_m F_h)` agrees to 1e-12 for `-100 < ζ < 1` and to 5e-5 for any `ζ` reached
with `|Ri| ≤ 10^4`, over `z0` from 1e-5 to 1e-2 m, `zT/z0` from 0.01 to 1, and wind and
temperature heights of 2 and 10 m in any combination.

The unstable result is bounded below by `ZETA_UNSTABLE_MIN`.

The stable functions are Beljaars & Holtslag (1991); the unstable ones Paulson (1970) with
Högström's (1988) coefficients (see [`stability_psi_m`](@ref)).
"""
@inline function obukhov_stability(Ri::Float64, z_w::Float64, z_T::Float64, z0::Float64,
    zT::Float64, logM::Float64=log(z_w / z0), logHT::Float64=log(z_T / zT))

    isfinite(Ri) || error("obukhov_stability: bulk Richardson number is $(Ri)")
    Ri == 0.0 && return 0.0

    sgn = sign(Ri)
    log_abs_Ri = log(abs(Ri))
    s = log_abs_Ri + 2 * log(logM) - log(logHT)
    r0 = z0 / z_w
    rT = z_T / z_w
    rTs = zT / z_w
    for _ in 1:MO_NEWTON_STEPS
        ζ = sgn * exp(s)
        ψm, φm = _stability_m(ζ)
        ψm0, φm0 = _stability_m(ζ * r0)
        ψhT, φhT = _stability_h(ζ * rT)
        ψh0, φh0 = _stability_h(ζ * rTs)
        Fm = logM - ψm + ψm0
        Fh = logHT - ψhT + ψh0
        R = s + log(Fh / (Fm * Fm)) - log_abs_Ri
        dR = 1 + (φhT - φh0) / Fh - 2 * (φm - φm0) / Fm
        s -= R / dR
    end
    return max(sgn * exp(s), ZETA_UNSTABLE_MIN)
end

# Beljaars & Holtslag (1991) stable-function constants.
const _BH_A = 1.0
const _BH_B = 2.0 / 3.0
const _BH_C = 5.0
const _BH_D = 0.35

"""
    stability_psi_m(ζ), stability_psi_h(ζ)

Integrated stability functions `Ψ(ζ) = ∫₀^ζ (1 − φ(x))/x dx` for momentum and heat, which
enter the transfer functions as `ln(z/z₀) − Ψ(z/L) + Ψ(z₀/L)`.

- **Stable (`ζ > 0`)**: Beljaars & Holtslag (1991) eqs. 28 and 32 in closed form, with
  `a = 1, b = 2/3, c = 5, d = 0.35`:
  `Ψ_m = −[aζ + b(ζ − c/d)e^{−dζ} + bc/d]` and
  `Ψ_h = −[(1 + 2aζ/3)^{3/2} + b(ζ − c/d)e^{−dζ} + bc/d − 1]`.
- **Unstable (`ζ ≤ 0`)**: Paulson (1970) with Högström's (1988) coefficients, in terms of the
  inverse profile functions `x_m = (1 − 19ζ)^{1/4}` and `x_h = (1 − 11.6ζ)^{1/2}`:
  `Ψ_m = 2ln((1 + x_m)/2) + ln((1 + x_m²)/2) − 2atan(x_m) + π/2` and
  `Ψ_h = 2ln((1 + x_h)/2)`.

`Ψ(0) = 0` on both branches, so the transfer functions are continuous through neutral;
`Ψ > 0` when unstable and `Ψ < 0` when stable. Högström's `φ_h` carries a `0.95` prefactor (his
`κ_H/κ_M` ratio) that makes `φ_h(0) = 0.95`; it is dropped here so that `Ψ_h(0) = 0`.

# References
- Paulson, C. A. (1970). The mathematical representation of wind speed and temperature profiles
  in the unstable atmospheric surface layer. *J. Appl. Meteorol.* 9, 857-861.
- Högström, U. (1988). Non-dimensional wind and temperature profiles in the atmospheric surface
  layer: a re-evaluation. *Boundary-Layer Meteorol.* 42, 55-78.
- Beljaars, A. C. M. & Holtslag, A. A. M. (1991). Flux parameterization over land surfaces for
  atmospheric models. *J. Appl. Meteorol.* 30, 327-341.
"""
@inline stability_psi_m(ζ::Float64) = _stability_m(ζ)[1]
@inline stability_psi_h(ζ::Float64) = _stability_h(ζ)[1]

"""
    stability_phi_m(ζ), stability_phi_h(ζ)

Profile functions `φ = 1 − ζ·dΨ/dζ` matching [`stability_psi_m`](@ref) and
`stability_psi_h`: Beljaars & Holtslag (1991) when stable,
`φ_m = 1 + ζ[a + b e^{−dζ}(1 + c − dζ)]` and
`φ_h = 1 + ζ[a(1 + 2aζ/3)^{1/2} + b e^{−dζ}(1 + c − dζ)]`; Högström (1988) when unstable,
`φ_m = (1 − 19ζ)^{−1/4}` and `φ_h = (1 − 11.6ζ)^{−1/2}`.
"""
@inline stability_phi_m(ζ::Float64) = _stability_m(ζ)[2]
@inline stability_phi_h(ζ::Float64) = _stability_h(ζ)[2]

# `(Ψ, φ)` together: the Newton step in `obukhov_stability` needs both at each argument, and
# they share their one transcendental (`exp(−dζ)` stable, the square roots unstable).
@inline function _stability_m(ζ::Float64)
    if ζ > 0.0
        e = exp(-_BH_D * ζ)
        ψ = -(_BH_A * ζ + _BH_B * (ζ - _BH_C / _BH_D) * e + _BH_B * _BH_C / _BH_D)
        φ = 1 + ζ * (_BH_A + _BH_B * e * (1 + _BH_C - _BH_D * ζ))
        return ψ, φ
    end
    x2 = sqrt(1.0 - 19.0 * ζ)          # x_m²
    x = sqrt(x2)                       # x_m = (1 − 19ζ)^(1/4)
    ψ = 2.0 * log((1 + x) / 2.0) + log((1 + x2) / 2.0) - 2 * atan(x) + pi / 2
    return ψ, 1.0 / x
end

@inline function _stability_h(ζ::Float64)
    if ζ > 0.0
        e = exp(-_BH_D * ζ)
        y = 1 + 2 * _BH_A * ζ / 3
        sy = sqrt(y)
        ψ = -(y * sy + _BH_B * (ζ - _BH_C / _BH_D) * e + _BH_B * _BH_C / _BH_D - 1.0)
        φ = 1 + ζ * (_BH_A * sy + _BH_B * e * (1 + _BH_C - _BH_D * ζ))
        return ψ, φ
    end
    x = sqrt(1.0 - 11.6 * ζ)           # x_h
    return 2.0 * log((1.0 + x) / 2.0), 1.0 / x
end
