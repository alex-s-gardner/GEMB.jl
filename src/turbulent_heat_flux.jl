"""
    turbulent_heat_flux(T_surface, density_air, z0, zT, zQ, cfs::ClimateForcingStep)

Compute sensible and latent heat fluxes using Monin-Obukhov similarity theory.
Matches MATLAB's `turbulent_heat_flux.m`.

Returns (heat_flux_sensible, heat_flux_latent, latent_heat) [W m-2, W m-2, J kg-1].
"""
function turbulent_heat_flux(T_surface::Float64, density_air::Float64,
    z0::Float64, zT::Float64, zQ::Float64, cfs::ClimateForcingStep)

    T_tolerance = 1e-10

    # Bulk-transfer coefficient (Neutral)
    An = VON_KARMAN^2  # 0.4^2 = 0.16
    C = An * cfs.wind_speed

    # Bulk Richardson Number
    Ri = ((100000 / cfs.pressure_air)^0.286) *
         (2.0 * GRAVITY * (cfs.temperature_air - T_surface)) /
         (cfs.temperature_observation_height * (cfs.temperature_air + T_surface) *
          ((cfs.wind_speed / cfs.wind_observation_height)^2.0))

    # Constants for Beljaars and Holtslag (1991)
    a1 = 1.0
    b1 = 2.0 / 3.0
    c1 = 5.0
    d1 = 0.35
    PhiMz0 = 0.0
    PhiHzT = 0.0
    PhiHzQ = 0.0

    if Ri > 0.0 + T_tolerance  # STABLE
        if Ri < 0.2 - T_tolerance
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
        zL = Ri / 1.5
        xm = (1.0 - 19.0 * zL)^(-0.25)
        PhiMz = 2.0 * log((1 + xm) / 2.0) + log((1 + xm^2) / 2.0) - 2 * atan(xm) + pi / 2
        xh = 0.95 * (1.0 - 11.6 * zL)^(-0.5)
        PhiHz = 2.0 * log((1.0 + xh^2) / 2.0)
    end

    # Final Transfer Coefficients
    coefM = log(cfs.wind_observation_height / z0) - PhiMz + PhiMz0
    coefHT = log(cfs.temperature_observation_height / zT) - PhiHz + PhiHzT
    coefHQ = log(cfs.temperature_observation_height / zQ) - PhiHz + PhiHzQ

    # Sensible Heat Flux [W m-2]
    heat_flux_sensible = density_air * C * C_AIR * (cfs.temperature_air - T_surface) * (100000 / cfs.pressure_air)^0.286
    heat_flux_sensible = heat_flux_sensible / (coefM * coefHT)

    # Latent Heat Flux [W m-2]
    if T_surface >= CtoK - T_tolerance
        # Liquid water surface
        latent_heat = LV
        # Saturation Vapor Pressure (Murray 1967)
        eS = 610.78 * exp(17.2693882 * (T_surface - CtoK - 0.01) / (T_surface - 35.86))
    else
        # Ice surface
        latent_heat = LS
        # Saturation Vapor Pressure (Bolton 1980)
        eS = 610.78 * exp(21.8745584 * (T_surface - CtoK - 0.01) / (T_surface - 7.66))
    end

    # 461.9 is the specific gas constant for water vapor [J kg-1 K-1]
    heat_flux_latent = C * latent_heat * (cfs.vapor_pressure - eS) /
                       (461.9 * (cfs.temperature_air + T_surface) / 2.0)
    heat_flux_latent = heat_flux_latent / (coefM * coefHQ)

    return heat_flux_sensible, heat_flux_latent, latent_heat
end
