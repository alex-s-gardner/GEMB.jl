"""
    turbulent_heat_flux(T_surface, density_air, z0, zT, zQ, cfs::ClimateForcingStep; min_wind_speed=0.01)

Compute sensible and latent heat fluxes using Monin-Obukhov similarity theory.

`min_wind_speed` floors `cfs.wind_speed` (which would otherwise drive the bulk
Richardson number to infinity at zero wind); clamping here avoids rebuilding the
`ClimateForcingStep` in the caller.

Returns (heat_flux_sensible, heat_flux_latent, latent_heat) [W m-2, W m-2, J kg-1].
"""
function turbulent_heat_flux(T_surface::Float64, density_air::Float64,
    z0::Float64, zT::Float64, zQ::Float64, cfs::ClimateForcingStep;
    min_wind_speed::Float64=0.01)

    wind_speed = max(cfs.wind_speed, min_wind_speed)

    # Bulk-transfer coefficient (Neutral)
    An = VON_KARMAN^2  # 0.4^2 = 0.16
    C = An * wind_speed

    # Exner-like pressure factor (100000/p)^0.286; reused below for the
    # sensible-heat flux, so compute the pow only once.
    pressure_factor = (100000 / cfs.pressure_air)^0.286

    # Neutral roughness log terms (invariant in T_surface).
    logM = log(cfs.wind_observation_height / z0)
    logHT = log(cfs.temperature_observation_height / zT)
    logHQ = log(cfs.temperature_observation_height / zQ)

    # Squared wind-to-height ratio in the bulk Richardson number, also invariant in
    # T_surface. Integer exponent: `^2` and `^2.0` agree bit-for-bit but the float
    # exponent compiles to a `pow` call rather than a multiply.
    wind_ratio_sq = (wind_speed / cfs.wind_observation_height)^2

    return _turbulent_heat_flux(T_surface, density_air, z0, zT, zQ, cfs,
        wind_speed, C, pressure_factor, logM, logHT, logHQ, wind_ratio_sq)
end

"""
    _turbulent_heat_flux(T_surface, density_air, z0, zT, zQ, cfs, wind_speed, C, pressure_factor, logM, logHT, logHQ, wind_ratio_sq)

Core of [`turbulent_heat_flux`](@ref) with the sub-timestep-invariant quantities
(`wind_speed`, bulk coefficient `C`, Exner `pressure_factor`, the neutral
roughness logs `logM`/`logHT`/`logHQ`, and the squared wind-to-height ratio
`wind_ratio_sq`) hoisted to the caller. `calculate_temperature`
computes these once per timestep and reuses them across all thermal sub-steps,
where only `T_surface` changes. Numerically identical to the inline form.

## Integrated stability functions

Both branches return `Ψ(ζ) = ∫₀^ζ (1 - φ(z))/z dz`, the integrated correction that enters the
transfer coefficients as `coef = log(z/z₀) - Ψ`. Two properties follow from that definition and
are what the branches are checked against in `test_turbulent_heat_flux.jl`:

- **`Ψ(0) = 0`.** The integrand vanishes at neutral, so both branches must meet at `Ri = 0`.
  Without this the fluxes step discontinuously as `T_surface` crosses `T_air` — and since the
  implicit solver's Newton iteration converges on exactly that crossing, a jump there can leave
  the surface energy balance with no solution in `T_surface` (Fourteau et al., 2024, Appendix D,
  who move their own branch point to `Ri_b = 0` for this reason).
- **`Ψ > 0` when unstable, `Ψ < 0` when stable.** Convection enhances exchange, so subtracting a
  positive `Ψ` shrinks `coef` and raises the flux; stable stratification does the reverse.

The stable branch is Beljaars & Holtslag (1991) eqs. 28 and 32, integrated in closed form. The
unstable branch is Paulson (1970) with Högström's (1988) coefficients. Paulson's argument is the
*inverse* profile function `x = φ⁻¹`, so both exponents are positive:
`x_m = (1 - 19ζ)^{1/4}` inverts `φ_m = (1 - 19ζ)^{-1/4}`, and `x_h = (1 - 11.6ζ)^{1/2}` inverts
`φ_h = (1 - 11.6ζ)^{-1/2}`. `Ψ_h = 2·ln((1 + x_h)/2)` takes `x_h` to the first power, the square
root already being carried by `x_h` itself.

Both closed forms were verified by numerically integrating their own `φ` back through the
definition above, agreeing to 6 decimal places at `ζ = ±0.01, ±0.1, ±1, ±5`.

Högström's `φ_h` carries a `0.95` prefactor (his `κ_H/κ_M` ratio) which makes `φ_h(0) = 0.95`
rather than 1. It is dropped here: `Ψ` is defined against a profile function that is 1 at
neutral, and retaining it would put `Ψ_h(0) = -0.0999` instead of 0 — reintroducing the very
discontinuity described above. The ratio belongs to the neutral transfer coefficient, not
inside the integral.

`ζ` is bounded below by `ZETA_UNSTABLE_MIN` on the unstable branch. This is a numerical
guard, not physics: `Ri` carries `wind_speed^-2`, so at the `min_wind_speed` floor it reaches
values far outside the range Monin-Obukhov theory was ever fit over, where `Ψ_h` would grow past
`log(z_T/z_Q)` and drive the transfer coefficient through zero. See that constant's comment.

# References
- Paulson, C. A. (1970). The mathematical representation of wind speed and temperature profiles
  in the unstable atmospheric surface layer. *J. Appl. Meteorol.* 9, 857-861.
- Högström, U. (1988). Non-dimensional wind and temperature profiles in the atmospheric surface
  layer: a re-evaluation. *Boundary-Layer Meteorol.* 42, 55-78.
- Beljaars, A. C. M. & Holtslag, A. A. M. (1991). Flux parameterization over land surfaces for
  atmospheric models. *J. Appl. Meteorol.* 30, 327-341.
- Fourteau, K., Brondex, J., Brun, F. & Dumont, M. (2024). A novel numerical implementation for
  the surface energy budget of melting snowpacks and glaciers. *Geosci. Model Dev.* 17,
  1903-1929. https://doi.org/10.5194/gmd-17-1903-2024
"""
@inline function _turbulent_heat_flux(T_surface::Float64, density_air::Float64,
    z0::Float64, zT::Float64, zQ::Float64, cfs::ClimateForcingStep,
    wind_speed::Float64, C::Float64, pressure_factor::Float64,
    logM::Float64, logHT::Float64, logHQ::Float64, wind_ratio_sq::Float64)

    # Bulk Richardson Number
    Ri = pressure_factor *
         (2.0 * GRAVITY * (cfs.temperature_air - T_surface)) /
         (cfs.temperature_observation_height * (cfs.temperature_air + T_surface) *
          wind_ratio_sq)

    # Constants for Beljaars and Holtslag (1991)
    a1 = 1.0
    b1 = 2.0 / 3.0
    c1 = 5.0
    d1 = 0.35
    PhiMz0 = 0.0
    PhiHzT = 0.0
    PhiHzQ = 0.0

    if Ri > 0.0 + T_TOLERANCE  # STABLE
        if Ri < 0.2 - T_TOLERANCE
            zL = Ri / (1.0 - 5.0 * Ri)
        else
            zL = Ri
        end

        zLM = max(zL / cfs.wind_observation_height * z0, 1e-3)
        zLT = max(zL / cfs.temperature_observation_height * zT, 1e-3)

        # Integrated Stability Functions (Psi)
        PhiMz = -(a1 * zL + b1 * (zL - c1 / d1) * exp(-d1 * zL) + b1 * c1 / d1)
        PhiHz = -((1 + 2 * a1 * zL / 3)^1.5 + b1 * (zL - c1 / d1) * exp(-d1 * zL) + b1 * c1 / d1 - 1.0)

        PhiMz0 = -(a1 * zLM + b1 * (zLM - c1 / d1) * exp(-d1 * zLM) + b1 * c1 / d1)
        PhiHzT = -((1 + 2 * a1 * zLT / 3)^1.5 + b1 * (zLT - c1 / d1) * exp(-d1 * zLT) + b1 * c1 / d1 - 1.0)

        PhiHzQ = PhiHzT
    else  # UNSTABLE
        # `Ri` carries `wind_speed^-2`, so at the `min_wind_speed` floor it reaches
        # magnitudes that are an artifact of that floor rather than a stability regime.
        # See `ZETA_UNSTABLE_MIN` for why it is bounded and why -100 is safe.
        zL = max(Ri / 1.5, ZETA_UNSTABLE_MIN)
        # Paulson (1970) integrated forms. Their argument is the *inverse* profile
        # function `x = φ⁻¹`, not `φ` itself: `x_m = (1 - 19ζ)^(1/4)` is `φ_m⁻¹` for
        # Högström's (1988) `φ_m = (1 - 19ζ)^(-1/4)`, and likewise for heat. Both
        # exponents are therefore positive here, and `x_h` enters `Ψ_h` to the first
        # power because it already carries the square root. See the docstring.
        xm = (1.0 - 19.0 * zL)^0.25
        PhiMz = 2.0 * log((1 + xm) / 2.0) + log((1 + xm^2) / 2.0) - 2 * atan(xm) + pi / 2
        xh = (1.0 - 11.6 * zL)^0.5
        PhiHz = 2.0 * log((1.0 + xh) / 2.0)
    end

    # Final Transfer Coefficients
    coefM = logM - PhiMz + PhiMz0
    coefHT = logHT - PhiHz + PhiHzT
    coefHQ = logHQ - PhiHz + PhiHzQ

    # Sensible Heat Flux [W m-2]
    heat_flux_sensible = density_air * C * HEAT_CAPACITY_AIR * (cfs.temperature_air - T_surface) * pressure_factor
    heat_flux_sensible = heat_flux_sensible / (coefM * coefHT)

    # Latent Heat Flux [W m-2]
    if T_surface >= CtoK - T_TOLERANCE
        # Liquid water surface
        latent_heat = LV
        # Saturation Vapor Pressure (Murray 1967)
        eS = saturation_vapor_pressure_water(T_surface)
    else
        # Ice surface
        latent_heat = LS
        # Saturation Vapor Pressure (Bolton 1980)
        eS = saturation_vapor_pressure_ice(T_surface)
    end

    # 461.9 is the specific gas constant for water vapor [J kg-1 K-1]
    heat_flux_latent = C * latent_heat * (cfs.vapor_pressure - eS) /
                       (461.9 * (cfs.temperature_air + T_surface) / 2.0)
    heat_flux_latent = heat_flux_latent / (coefM * coefHQ)

    return heat_flux_sensible, heat_flux_latent, latent_heat
end
