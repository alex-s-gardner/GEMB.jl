"""
    calculate_melt(temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, albedo, albedo_diffuse, rain, mp::ModelParameters, verbose::Bool)

Compute meltwater generation, percolation, refreezing, and runoff using a tipping bucket approach.

Processes:
1. Initial Refreeze: Existing pore water in cold layers is refrozen.
2. Melt Generation: Excess energy above 0 C is converted to liquid meltwater.
3. Percolation: Liquid water percolates downward, refreezing in cold layers,
   being retained as pore water, or running off at impermeable ice lenses.

Returns `(temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, albedo, albedo_diffuse, melt_total, melt_surface, runoff_total, freeze_total)`.
Arrays may shrink (cells deleted when mass=0).
"""
function calculate_melt(temperature::Vector{Float64}, dz::Vector{Float64},
    density::Vector{Float64}, water::Vector{Float64},
    grain_radius::Vector{Float64}, grain_dendricity::Vector{Float64},
    grain_sphericity::Vector{Float64}, albedo::Vector{Float64},
    albedo_diffuse::Vector{Float64}, rain::Float64,
    mp::ModelParameters, verbose::Bool)

    # Note: arrays are modified in-place. May shrink via deleteat! when cells lose all mass.

    T_tolerance = 1e-10
    d_tolerance = 1e-11
    water_tolerance = 1e-13

    # Specify constants
    d_phc = 830.0           # pore hole close off density [kg m-3]
    ice_layer_dzmin = 0.1   # minimum ice layer thickness for runoff [m]

    m = length(temperature)
    water_delta = zeros(m)

    # store initial mass [kg]
    M = dz .* density

    if verbose
        M_total_initial = sum(water) + sum(M)
        E_total_initial = sum(M .* temperature .* C_ICE) +
            sum(water .* (LF + CtoK * C_ICE))
    end

    # initialize melt and runoff scalars
    runoff_total = 0.0
    melt_total = 0.0
    melt_surface = 0.0

    # calculate temperature excess above 0 degC
    T_excess = max.(0.0, temperature .- CtoK)

    # new grid point center temperature
    temperature = min.(temperature, CtoK)

    # specify irreducible water content saturation [fraction]
    water_irreducible_saturation = 0.07

    ## REFREEZE PORE WATER

    if sum(water) > water_tolerance
        # calculate maximum freeze amount [kg]
        freeze_max = max.(0.0, -((temperature .- CtoK) .* M .* C_ICE) ./ LF)

        # freeze pore water and change snow/ice properties
        water_delta = min.(freeze_max, water)
        water = water .- water_delta
        M = M .+ water_delta
        density = M ./ dz
        temperature = temperature .+ Float64.(M .> water_tolerance) .* (water_delta .* (LF .+ (CtoK .- temperature) .* C_ICE) ./ (M .* C_ICE))

        # if pore water froze in ice then adjust density and dz thickness
        density[density .> mp.density_ice - d_tolerance] .= mp.density_ice
        dz = M ./ density
    end

    # squeeze water from snow pack
    water_irreducible = (mp.density_ice .- density) .* mp.water_irreducible_saturation .* (M ./ density)
    water_excess = max.(0.0, water .- water_irreducible)

    ## MELT, PERCOLATION AND REFREEZE

    freeze = zeros(m)

    # Add previous freeze to freeze and reset water_delta
    freeze = freeze .+ water_delta
    water_delta .= 0.0

    # run melt algorithm if there is melt water or excess pore water
    if (sum(T_excess) > T_tolerance) || (sum(water_excess) > water_tolerance)

        # Check to see if thermal energy exceeds energy to melt entire cell
        T_surplus = max.(0.0, T_excess .- LF / C_ICE)

        if sum(T_surplus) > T_tolerance
            # calculate surplus energy
            E_surplus = T_surplus .* C_ICE .* M
            i = 1

            while (sum(E_surplus) > T_tolerance) && (i < (m + 1))
                if i < m
                    # use surplus energy to increase the temperature of lower cell
                    temperature[i+1] = E_surplus[i] / M[i+1] / C_ICE + temperature[i+1]

                    T_excess[i+1] = max(0.0, temperature[i+1] - CtoK) + T_excess[i+1]
                    temperature[i+1] = min(CtoK, temperature[i+1])

                    T_surplus[i+1] = max(0.0, T_excess[i+1] - LF / C_ICE)
                    E_surplus[i+1] = T_surplus[i+1] * C_ICE * M[i+1]
                else
                    error("surplus energy reached the base of gemb column (i.e. entire column melted out in a single time step)")
                end

                # adjust current cell properties
                T_excess[i] = LF / C_ICE
                E_surplus[i] = 0.0
                i += 1
            end
        end

        # convert temperature excess to melt [kg]
        melt_maximum = T_excess .* density .* dz .* C_ICE ./ LF
        melt = min.(melt_maximum, M)
        melt_surface = melt[1]
        melt_total = max(0.0, sum(melt) - rain)

        # calculate maximum refreeze amount [kg]
        freeze_max = max.(0.0, -((temperature .- CtoK) .* density .* dz .* C_ICE) ./ LF)

        # initialize refreeze, runoff, flux_dn and water_delta vectors
        runoff = zeros(m)
        flux_dn = zeros(m + 1)

        # determine the deepest grid cell where melt/pore water is generated
        X_idx = findlast((melt .> water_tolerance) .| (water_excess .> water_tolerance))
        X = isnothing(X_idx) ? 1 : X_idx

        Xi = 1
        m = length(temperature)

        # meltwater percolation
        for i in 1:m
            # calculate total melt water entering cell
            melt_input = melt[i] + flux_dn[i]

            ice_depth = 0.0
            # If this grid cell's density exceeds the pore closeoff density:
            if density[i] >= d_phc - d_tolerance
                for l in i:m
                    if density[l] >= d_phc - d_tolerance
                        ice_depth += dz[l]
                        if ice_depth > ice_layer_dzmin + d_tolerance
                            break
                        end
                    else
                        break
                    end
                end
            end

            # break loop if there is no meltwater and if depth is > mw_depth
            if abs(melt_input) < water_tolerance && i > X
                break

            # if reaches impermeable ice layer all liquid water runs off
            elseif (density[i] >= (mp.density_ice - d_tolerance)) ||
                   ((density[i] >= d_phc - d_tolerance) && (ice_depth > ice_layer_dzmin + d_tolerance))

                M[i] = M[i] - melt[i]
                water_irr = (mp.density_ice - density[i]) * water_irreducible_saturation * (M[i] / density[i])
                water_delta[i] = max(min(melt_input, water_irr - water[i]), -water[i])
                runoff[i] = max(0.0, melt_input - water_delta[i])

            # check if no energy to refreeze meltwater
            elseif abs(freeze_max[i]) < d_tolerance

                M[i] = M[i] - melt[i]
                water_irr = (mp.density_ice - density[i]) * water_irreducible_saturation * (M[i] / density[i])
                water_delta[i] = max(min(melt_input, water_irr - water[i]), -1 * water[i])
                flux_dn[i+1] = max(0.0, melt_input - water_delta[i])
                runoff[i] = 0.0

            # some or all meltwater refreezes
            else
                M[i] = M[i] - melt[i]
                dz_0 = M[i] / density[i]
                d_max = (mp.density_ice - density[i]) * dz_0
                freeze1 = min(min(melt_input, d_max), freeze_max[i])
                M[i] = M[i] + freeze1
                density[i] = M[i] / dz_0

                # pore water
                water_irr = (mp.density_ice - density[i]) * water_irreducible_saturation * dz_0
                water_delta[i] = max(min(melt_input - freeze1, water_irr - water[i]), -1 * water[i])
                freeze2 = 0.0

                if water_delta[i] < 0.0 - water_tolerance
                    d_max = (mp.density_ice - density[i]) * dz_0
                    freeze2_max = min(d_max, freeze_max[i] - freeze1)
                    freeze2 = min(-1.0 * water_delta[i], freeze2_max)
                    M[i] = M[i] + freeze2
                    density[i] = M[i] / dz_0
                end

                freeze[i] = freeze[i] + freeze1 + freeze2
                flux_dn[i+1] = max(0.0, melt_input - freeze1 - water_delta[i])

                if M[i] > water_tolerance
                    temperature[i] = temperature[i] +
                        ((freeze1 + freeze2) * (LF + (CtoK - temperature[i]) * C_ICE) / (M[i] * C_ICE))
                end

                # check if an ice layer forms
                if abs(density[i] - mp.density_ice) < d_tolerance
                    runoff[i] = flux_dn[i+1]
                    flux_dn[i+1] = 0.0
                end
            end

            Xi = Xi + 1
        end

        # Check for negative pore water
        if verbose
            if any(water .< 0.0 - water_tolerance)
                error("Negative pore water generated in melt equations.")
            end
        end

        # adjust pore water
        water = water .+ water_delta

        # calculate runoff_total
        runoff_total = sum(runoff) + flux_dn[Xi]

        # delete all cells with zero mass
        to_delete = findall(M .<= water_tolerance)
        if !isempty(to_delete)
            deleteat!(M, to_delete)
            deleteat!(water, to_delete)
            deleteat!(density, to_delete)
            deleteat!(temperature, to_delete)
            deleteat!(albedo, to_delete)
            deleteat!(grain_radius, to_delete)
            deleteat!(grain_dendricity, to_delete)
            deleteat!(grain_sphericity, to_delete)
            deleteat!(albedo_diffuse, to_delete)
        end

        # calculate new grid lengths
        dz = M ./ density
    end

    freeze_total = sum(freeze)

    ## CHECK FOR MASS AND ENERGY CONSERVATION
    if verbose
        E_total_runoff = runoff_total * (LF + CtoK * C_ICE)

        M_total_final = sum(water) + sum(M) + runoff_total
        E_total_final = sum(M .* temperature .* C_ICE) + sum(water .* (LF + CtoK * C_ICE))

        M_delta = M_total_initial - M_total_final
        E_delta = E_total_initial - E_total_final - E_total_runoff

        if (abs(M_delta) > 1e-3) || (abs(E_delta) > 1e-3)
            error("Mass and/or energy are not conserved in melt equations:\n M_delta: $(M_delta) E_delta: $(E_delta)\n")
        end

        if any(water .< 0.0 - water_tolerance)
            error("Negative pore water generated in melt equations.")
        end
    end

    return temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, albedo, albedo_diffuse, melt_total, melt_surface, runoff_total, freeze_total
end
