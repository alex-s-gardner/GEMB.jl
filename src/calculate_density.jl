"""
    calculate_density(temperature, dz, density, grain_radius, cfs::ClimateForcingStep, mp::ModelParameters)

Compute the densification of snow/firn using one of several models:
- "HerronLangway": Herron and Langway (1980)
- "Arthern": semi-empirical model of Arthern et al. (2010) [default]
- "ArthernB": physical model from Appendix B of Arthern et al. (2010)
- "LiZwally": empirical model of Li and Zwally (2004)
- "Helsen": modified empirical model by Helsen et al. (2008)
- "Ligtenberg": semi-empirical model of Ligtenberg et al. (2011)

Returns `(dz, density)` as new vectors.

# References
- Arthern, R. J., et al. (2010). J. Geophys. Res., 115, F03011.
- Herron, M. and Langway, C. (1980). J. Glaciol., 25, 373-385.
- Li, J. and Zwally, H. (2004). Ann. Glaciol., 38, 309-313.
- Helsen, M. M., et al. (2008). Science, 320, 1626-1629.
- Ligtenberg, S. R. M., et al. (2011). The Cryosphere, 5, 809-819.
"""
function calculate_density(temperature::Vector{Float64}, dz::Vector{Float64},
    density::Vector{Float64}, grain_radius::Vector{Float64},
    cfs::ClimateForcingStep, mp::ModelParameters)

    # Note: density is modified in-place. dz is recomputed from mass_init / density.

    d_tolerance = 1e-11

    # specify constants
    dt = cfs.dt / 86400.0   # convert from [s] to [d]
    R = 8.314               # gas constant [mol-1 K-1]

    # initial mass
    mass_init = density .* dz

    # calculate new snow/firn density for:
    #   snow with densities <= 550 [kg m-3]
    #   snow with densities > 550 [kg m-3]
    idx = density .<= (550.0 + d_tolerance)

    local c0::Vector{Float64}
    local c1::Vector{Float64}

    if mp.densification_method == :HerronLangway
        c0 = (11 .* exp.(-10160 ./ (temperature[idx] .* R))) .* cfs.precipitation_mean / 1000
        c1 = (575 .* exp.(-21400 ./ (temperature[.!idx] .* R))) .* (cfs.precipitation_mean / 1000)^0.5

    elseif mp.densification_method == :Arthern
        H = exp.((-60000.0 ./ (temperature .* R)) .+ (42400.0 ./ (cfs.temperature_air_mean .* R))) .*
            (cfs.precipitation_mean * 9.81)

        c0 = 0.07 .* H[idx]
        c1 = 0.03 .* H[.!idx]

    elseif mp.densification_method == :ArthernB
        # calculate overburden pressure
        obp = vcat([0.0], cumsum(dz[1:end-1]) .* density[1:end-1])

        H = exp.((-60000.0 ./ (temperature .* R))) .* obp ./ (grain_radius ./ 1000) .^ 2
        c0 = 9.2e-9 .* H[idx]
        c1 = 3.7e-9 .* H[.!idx]

    elseif mp.densification_method == :LiZwally
        c_all = (cfs.precipitation_mean ./ mp.density_ice) .*
            max.(139.21 .- 0.542 .* cfs.temperature_air_mean, 1.0) .* 8.36 .* max.(CtoK .- temperature, 1.0) .^ (-2.061)
        c0 = c_all[idx]
        c1 = c_all[.!idx]

    elseif mp.densification_method == :Helsen
        c_all = (cfs.precipitation_mean ./ mp.density_ice) .*
            max.(76.138 .- 0.28965 .* cfs.temperature_air_mean, 1.0) .* 8.36 .* max.(CtoK .- temperature, 1.0) .^ (-2.061)
        c0 = c_all[idx]
        c1 = c_all[.!idx]

    elseif mp.densification_method == :Ligtenberg
        H = exp.((-60000.0 ./ (temperature .* R)) .+ (42400.0 ./ (cfs.temperature_air_mean .* R))) .*
            (cfs.precipitation_mean .* 9.81)

        c0arth = 0.07 .* H
        c1arth = 0.03 .* H

        M01 = densification_lookup_M01(mp.densification_coeffs_M01)

        if length(M01) == 4
            M0 = max(M01[1] - (M01[2] * log(cfs.precipitation_mean)), 0.25)
            M1 = max(M01[3] - (M01[4] * log(cfs.precipitation_mean)), 0.25)
        else
            if abs(mp.density_ice - 820.0) < d_tolerance
                M0 = max(M01[1, 1] - (M01[1, 2] * log(cfs.precipitation_mean)), 0.25)
                M1 = max(M01[1, 3] - (M01[1, 4] * log(cfs.precipitation_mean)), 0.25)
            else
                M0 = max(M01[2, 1] - (M01[2, 2] * log(cfs.precipitation_mean)), 0.25)
                M1 = max(M01[2, 3] - (M01[2, 4] * log(cfs.precipitation_mean)), 0.25)
            end
        end

        c0 = M0 .* c0arth[idx]
        c1 = M1 .* c1arth[.!idx]
    else
        error("unrecognized densification method")
    end

    # new snow density
    density[idx] = density[idx] .+ (c0 .* (mp.density_ice .- density[idx]) ./ 365 .* dt)
    density[.!idx] = density[.!idx] .+ (c1 .* (mp.density_ice .- density[.!idx]) ./ 365 .* dt)

    # do not allow densities to exceed the density of ice
    density[density .> (mp.density_ice - d_tolerance)] .= mp.density_ice

    # calculate new grid cell length
    dz = mass_init ./ density

    return dz, density
end
