"""
    calculate_grain_size(temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, cfs::ClimateForcingStep, mp::ModelParameters)

Model the evolution of effective snow grain size, dendricity, and sphericity.

Accounts for different physical processes depending on snow state:
- Dendritic Snow (fresh): Evolves based on temperature gradients and liquid water content.
- Nondendritic Dry Snow: Temperature gradient metamorphism using Marbouty (1980).
- Wet Snow: Rapid grain growth due to liquid water using Brun (1989).

Only executes if `mp.albedo_method` is "GardnerSharp" or "BrunLefebre".

Returns `(grain_radius, grain_dendricity, grain_sphericity)` as new vectors.

# References
- Brun, E., et al. (1992). J. Glaciol., 38, 13-22.
- Marbouty, D. (1980). J. Glaciol., 26, 303-312.
- Brun, E. (1989). Ann. Glaciol., 13, 22-26.
"""
function calculate_grain_size(temperature::Vector{Float64}, dz::Vector{Float64},
    density::Vector{Float64}, water::Vector{Float64},
    grain_radius::Vector{Float64}, grain_dendricity::Vector{Float64},
    grain_sphericity::Vector{Float64},
    cfs::ClimateForcingStep, mp::ModelParameters)

    # Note: grain_dendricity and grain_sphericity are modified in-place.
    # grain_radius is returned as a new array (derived from gsz).

    T_tolerance = 1e-10
    gdn_tolerance = 1e-10
    water_tolerance = 1e-13

    # Only run grain growth for these albedo methods
    if !(mp.albedo_method == :GardnerSharp || mp.albedo_method == :BrunLefebre)
        return grain_radius, grain_dendricity, grain_sphericity
    end

    gsz = grain_radius .* 2

    # Convert dt from seconds to days
    dt_days = cfs.dt / 86400.0

    # Determine liquid-water content in terms of mass fraction [%]
    lwc = water ./ (density .* dz) .* 100

    # Set maximum water content by mass to 9 percent (Brun, 1980)
    lwc[lwc .> (9 + water_tolerance)] .= 9

    # Calculate temperature gradient across grid cells
    dT = zeros(length(temperature))

    # depth of grid point center from surface
    z_center = cumsum(dz) .- dz ./ 2

    # Take forward differences on left and right edges
    m = length(z_center)
    if m > 2
        dT[1] = (temperature[3] - temperature[1]) / (z_center[3] - z_center[1])
        dT[m] = (temperature[m] - temperature[m-2]) / (z_center[m] - z_center[m-2])
    elseif m > 1
        dT[1] = (temperature[2] - temperature[1]) / (z_center[2] - z_center[1])
        dT[m] = (temperature[m] - temperature[m-1]) / (z_center[m] - z_center[m-1])
    end

    # Take centered differences on interior points
    if m > 2
        z_center_diff = z_center[3:end] .- z_center[1:end-2]
        dT[2:end-1] = (temperature[3:end] .- temperature[1:end-2]) ./ z_center_diff
    end

    # take absolute value of temperature gradient
    dT = abs.(dT)

    # index for dendricity > 0 & == 0
    G = grain_dendricity .> (0 + gdn_tolerance)
    J = .!G

    ## DENDRITIC SNOW METAMORPHISM
    # FOR SNOW DENTRICITY > 0

    if sum(G) != 0
        # index for dentricity > 0 and temperature gradients < 5 degC m-1 and >= 5 degC m-1
        H = (abs.(dT) .<= 5 + T_tolerance) .& G .& (water .<= 0 + water_tolerance)
        I = (abs.(dT) .> 5 + T_tolerance) .& G .& (water .<= 0 + water_tolerance)

        # determine coefficients
        A = -2e8 .* exp.(-6e3 ./ temperature[H]) .* dt_days
        B = 1e9 .* exp.(-6e3 ./ temperature[H]) .* dt_days
        C = (-2e8 .* exp.(-6e3 ./ temperature[I]) .* dt_days) .* (abs.(dT[I]) .^ 0.4)

        # new dendricity and sphericity for dT < 5 degC m-1
        grain_dendricity[H] = grain_dendricity[H] .+ A
        grain_sphericity[H] = grain_sphericity[H] .+ B

        # new dendricity and sphericity for dT >= 5 degC m-1
        grain_dendricity[I] = grain_dendricity[I] .+ C
        grain_sphericity[I] = grain_sphericity[I] .+ C

        # WET SNOW METAMORPHISM
        # index for dendritic wet snow
        L = (water .> (0 + water_tolerance)) .& G

        # check if snowpack is wet
        if sum(L) != 0
            # determine coefficient
            D = (1 / 16) .* (lwc[L] .^ 3) .* dt_days

            # new dendricity and sphericity for wet snow
            grain_dendricity[L] = grain_dendricity[L] .- D
            grain_sphericity[L] = grain_sphericity[L] .+ D
        end

        # dendricity and sphericity can not be > 1 or < 0
        grain_dendricity[grain_dendricity .<= 0 + gdn_tolerance] .= 0
        grain_sphericity[grain_sphericity .<= 0 + gdn_tolerance] .= 0
        grain_dendricity[grain_dendricity .>= 1 - gdn_tolerance] .= 1
        grain_sphericity[grain_sphericity .>= 1 - gdn_tolerance] .= 1

        # determine new grain size (mm)
        gsz[G] = max.(0.1 .* (grain_dendricity[G] ./ 0.99 .+ (1.0 .- 1.0 .* grain_dendricity[G] ./ 0.99) .* (grain_sphericity[G] ./ 0.99 .* 3.0 .+ (1.0 .- grain_sphericity[G] ./ 0.99) .* 4.0)), gdn_tolerance * 2)
    end

    # if there is snow dendricity == 0
    if sum(J) != 0
        # When wet-snow grains (class 6) are submitted to a temperature gradient
        # higher than 5 degC m-1, their sphericity decreases according to Equations (4).
        P1 = J .& (grain_sphericity .> gdn_tolerance) .& (grain_sphericity .< 1 - gdn_tolerance) .& (abs.(dT) .> 5 + T_tolerance)
        P2 = J .& (grain_sphericity .> gdn_tolerance) .& (grain_sphericity .< 1 - gdn_tolerance) .& ((abs.(dT) .<= 5 + T_tolerance) .& (water .> 0 + water_tolerance))
        P3 = J .& (grain_sphericity .> gdn_tolerance) .& (grain_sphericity .< 1 - gdn_tolerance) .& .!P1 .& .!P2

        F1 = (-2e8 .* exp.(-6e3 ./ temperature[P1]) .* dt_days) .* abs.(dT[P1]) .^ 0.4
        F2 = (1.0 / 16.0) .* lwc[P2] .^ 3.0 .* dt_days
        F3 = 1e9 .* exp.(-6e3 ./ temperature[P3]) .* dt_days

        grain_sphericity[P1] = grain_sphericity[P1] .+ F1
        grain_sphericity[P2] = grain_sphericity[P2] .+ F2
        grain_sphericity[P3] = grain_sphericity[P3] .+ F3

        # sphericity can not be > 1 or < 0
        grain_sphericity[grain_sphericity .<= 0 + gdn_tolerance] .= 0
        grain_sphericity[grain_sphericity .>= 1 - gdn_tolerance] .= 1

        # DRY SNOW METAMORPHISM (Marbouty, 1980)
        P = J .& ((water .<= 0 + water_tolerance) .| ((grain_sphericity .<= 0 + gdn_tolerance) .& (abs.(dT) .> 5 + T_tolerance)))
        Q = _Marbouty(temperature[P], density[P], dT[P])

        # Calculate grain growth
        gsz[P] = gsz[P] .+ Q .* dt_days

        # WET SNOW METAMORPHISM (Brun, 1989)
        # Index for nondendritic wet snow
        K = J .& .!((water .<= 0 + water_tolerance) .| ((grain_sphericity .<= 0 + gdn_tolerance) .& (abs.(dT) .> 5 + T_tolerance)))

        # check if snowpack is wet
        if sum(K) != 0
            # wet rate of change coefficient
            E = (1.28e-8 .+ 4.22e-10 .* (lwc[K] .^ 3)) .* (dt_days * 86400)   # [mm^3 s^-1]

            # calculate change in grain volume and convert to grain size
            gsz[K] = 2 .* (3 / (pi * 4) .* ((4 / 3) .* pi .* (gsz[K] ./ 2) .^ 3 .+ E)) .^ (1 / 3)
        end

        # grains with sphericity == 1 can not have grain sizes > 2 mm (Brun, 1992)
        gsz[(abs.(grain_sphericity .- 1) .< water_tolerance) .& (gsz .> 2 - water_tolerance)] .= 2

        # grains with sphericity == 0 can not have grain sizes > 5 mm (Brun, 1992)
        gsz[(abs.(grain_sphericity .- 1) .>= water_tolerance) .& (gsz .> 5 - water_tolerance)] .= 5
    end

    # Convert grain size back to effective grain radius
    grain_radius = gsz ./ 2

    return grain_radius, grain_dendricity, grain_sphericity
end

"""
    _Marbouty(temperature, density, dT)

Calculate grain growth according to Fig. 9 of Marbouty (1980).
No grain growth for density > 400 kg m-3 (H is set to zero).
"""
function _Marbouty(temperature::Vector{Float64}, density::Vector{Float64}, dT::Vector{Float64})
    T_tolerance = 1e-10
    d_tolerance = 1e-11

    F = zeros(length(temperature))
    H = zeros(length(temperature))
    G = zeros(length(temperature))

    E = 0.09       # model time growth constant [mm d-1]
    temperature = temperature .- 273.15  # converts temperature from K to C
    dT = dT ./ 100.0   # convert dT from degC/m to degC/cm

    ## Temperature coefficient F
    I = temperature .> -6 + T_tolerance
    F[I] = 0.7 .+ ((temperature[I] ./ -6) .* 0.3)

    I = (temperature .<= -6 + T_tolerance) .& (temperature .> -22 + T_tolerance)
    F[I] = 1 .- ((temperature[I] .+ 6) ./ -16 .* 0.8)

    I = (temperature .<= -22 + T_tolerance) .& (temperature .> -40 + T_tolerance)
    F[I] = 0.2 .- ((temperature[I] .+ 22) ./ -18 .* 0.2)

    ## Density coefficient H
    H[density .< 150 - d_tolerance] .= 1

    I = (density .>= 150 - d_tolerance) .& (density .< 400 - d_tolerance)
    H[I] = 1 .- ((density[I] .- 150) ./ 250)

    ## Temperature gradient coefficient G
    I = (dT .>= 0.16 - T_tolerance) .& (dT .< 0.25 - T_tolerance)
    G[I] = ((dT[I] .- 0.16) ./ 0.09) .* 0.1

    I = (dT .>= 0.25 - T_tolerance) .& (dT .< 0.40 - T_tolerance)
    G[I] = 0.10 .+ (((dT[I] .- 0.25) ./ 0.15) .* 0.57)

    I = (dT .>= 0.40 - T_tolerance) .& (dT .< 0.50 - T_tolerance)
    G[I] = 0.67 .+ (((dT[I] .- 0.40) ./ 0.10) .* 0.23)

    I = (dT .>= 0.50 - T_tolerance) .& (dT .< 0.70 - T_tolerance)
    G[I] = 0.90 .+ (((dT[I] .- 0.50) ./ 0.20) .* 0.1)

    G[dT .>= 0.7 - T_tolerance] .= 1

    ## Grouped coefficient Q
    Q = F .* H .* G .* E

    return Q
end
