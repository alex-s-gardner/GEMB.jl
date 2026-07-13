"""
    calculate_shortwave_radiation(dz, density, grain_radius, albedo_surface, albedo_diffuse_surface, cfs::ClimateForcingStep, mp::ModelParameters)

Distribute absorbed shortwave radiation vertically within snow/ice.

Depending on model configuration:
1. Surface Absorption: All net shortwave energy is absorbed by the top grid cell
   (`shortwave_subsurface_absorption = false`).
2. Subsurface Penetration: Shortwave energy penetrates and is absorbed by deeper layers
   (`shortwave_subsurface_absorption = true`), using either:
   - Density-dependent extinction (Bassford, 2002)
   - Spectral-dependent extinction for "BrunLefebre" (Lefebre et al., 2003)

Returns `shortwave_flux` vector [W m-2] of absorbed shortwave radiation per grid cell.

# References
- Lefebre, F., et al. (2003). J. Geophys. Res., 108, 4231.
- Greuell, W. and Konzelmann, T. (1994). Global Planet. Change, 9, 91-114.
"""
function calculate_shortwave_radiation(dz::Vector{Float64}, density::Vector{Float64},
    grain_radius::Vector{Float64}, albedo_surface::Float64,
    albedo_diffuse_surface::Float64,
    cfs::ClimateForcingStep, mp::ModelParameters)

    d_tolerance = 1e-11

    # Initialize variables
    m = length(density)
    shortwave_flux = zeros(m)

    if (!mp.shortwave_subsurface_absorption) ||
       ((mp.density_ice - density[1]) < d_tolerance)
        # all sw radiation is absorbed by the top grid cell

        if mp.albedo_method == :GardnerSharp
            shortwave_flux[1] = (1.0 - albedo_surface) * max(0.0, (cfs.shortwave_downward - cfs.shortwave_downward_diffuse)) +
                (1.0 - albedo_diffuse_surface) * cfs.shortwave_downward_diffuse
        else
            shortwave_flux[1] = (1 - albedo_surface) * cfs.shortwave_downward
        end

    else  # sw radiation is absorbed at depth within the glacier

        if mp.albedo_method == :BrunLefebre
            # convert effective radius [mm] to grain size [m]
            gsz = (grain_radius .* 2) ./ 1000

            # Spectral fractions [0.3-0.8um 0.8-1.5um 1.5-2.8um]
            sF = [0.606, 0.301, 0.093]

            # initialize variables
            B1_cum = ones(m + 1)
            B2_cum = ones(m + 1)

            # spectral albedos:
            a1 = min(0.98, 0.95 - 1.58 * gsz[1]^0.5)
            a2 = max(0.0, 0.95 - 15.4 * gsz[1]^0.5)
            a3 = max(0.127, 0.88 + 346.3 * gsz[1] - 32.31 * gsz[1]^0.5)

            # separate net shortwave radiative flux into spectral ranges
            swfS = (sF .* cfs.shortwave_downward) .* (1 .- [a1, a2, a3])

            # absorption coefficient for spectral range
            h = density ./ (gsz .^ 0.5)
            B1 = 0.0192 .* h                 # 0.3 - 0.8um
            B2 = 0.1098 .* h                 # 0.8 - 1.5um

            # cumulative extinction factors
            B1_cum[2:end] = cumprod(exp.(-B1 .* dz))
            B2_cum[2:end] = cumprod(exp.(-B2 .* dz))

            # flux across grid cell boundaries
            Qs1 = swfS[1] .* B1_cum
            Qs2 = swfS[2] .* B2_cum

            # net energy flux to each grid cell
            shortwave_flux = (Qs1[1:m] .- Qs1[2:m+1]) .+ (Qs2[1:m] .- Qs2[2:m+1])

            # add flux absorbed at surface
            shortwave_flux[1] = shortwave_flux[1] + swfS[3]

        else  # function of grid cell density
            # fraction of sw radiation absorbed in top grid cell (wavelength > 0.8um)
            SWs = 0.36

            # calculate surface shortwave radiation fluxes [W m-2]
            swf_s = SWs * (1 - albedo_surface) * cfs.shortwave_downward

            # calculate subsurface shortwave radiation fluxes [W m-2]
            swf_ss = (1 - SWs) * (1 - albedo_surface) * cfs.shortwave_downward

            # SW extinction coefficients
            Bs = 10.0    # snow SW extinction coefficient [m-1] (Bassford, 2006)
            Bi = 1.3     # ice SW extinction coefficient [m-1] (Bassford, 2006)

            # calculate extinction coefficient B [m-1] vector
            B = Bs .+ (300 .- density) .* ((Bs - Bi) / (mp.density_ice - 300))

            # cumulative extinction factor
            B_cum = vcat([1.0], cumprod(exp.(-B .* dz)))

            # flux across grid cell boundaries
            Qs = swf_ss .* B_cum

            # net energy flux to each grid cell
            shortwave_flux = Qs[1:m] .- Qs[2:m+1]

            # add flux absorbed at surface
            shortwave_flux[1] = shortwave_flux[1] + swf_s
        end
    end

    return shortwave_flux
end
