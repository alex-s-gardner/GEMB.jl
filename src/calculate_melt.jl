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

    # new grid point center temperature. Rebind to a fresh array (do not mutate
    # the caller's temperature vector) to preserve the previous behavior.
    temperature = min.(temperature, CtoK)

    # specify irreducible water content saturation [fraction]
    water_irreducible_saturation = 0.07

    ## REFREEZE PORE WATER

    if sum(water) > water_tolerance
        # Fused per-cell refreeze: freeze pore water and update snow/ice
        # properties. Numerically identical (same per-element ops and
        # left-to-right ordering) to the previous broadcast chain, but computes
        # everything in a single pass. density, dz, and water are rebound to
        # fresh arrays (matching the previous `density = M ./ dz` etc.) so the
        # caller's inputs are not mutated; M and temperature are mutated in
        # place since they are already function-local fresh arrays.
        dz_orig = dz
        density = similar(density)
        dz = similar(dz_orig)
        water = copy(water)
        @inbounds for i in 1:m
            # maximum freeze amount [kg]
            freeze_max_i = max(0.0, -((temperature[i] - CtoK) * M[i] * C_ICE) / LF)

            # freeze pore water and change snow/ice properties
            wd = min(freeze_max_i, water[i])
            water_delta[i] = wd
            water[i] = water[i] - wd
            M[i] = M[i] + wd
            density[i] = M[i] / dz_orig[i]
            mask = M[i] > water_tolerance ? 1.0 : 0.0
            temperature[i] = temperature[i] + mask *
                (wd * (LF + (CtoK - temperature[i]) * C_ICE) / (M[i] * C_ICE))

            # if pore water froze in ice then adjust density and dz thickness
            if density[i] > mp.density_ice - d_tolerance
                density[i] = mp.density_ice
            end
            dz[i] = M[i] / density[i]
        end
    end

    # squeeze water from snow pack (compute water_excess without materializing
    # the water_irreducible temporary)
    water_excess = Vector{Float64}(undef, m)
    @inbounds for i in 1:m
        water_irreducible = (mp.density_ice - density[i]) * mp.water_irreducible_saturation * (M[i] / density[i])
        water_excess[i] = max(0.0, water[i] - water_irreducible)
    end

    ## MELT, PERCOLATION AND REFREEZE

    # Seed freeze with the pore-water refreeze accumulated above, then reset
    # water_delta for reuse in the percolation loop.
    freeze = copy(water_delta)
    fill!(water_delta, 0.0)

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

        # convert temperature excess to melt [kg] and compute the max refreeze
        # amount in a single fused pass. melt[i] = min(melt_maximum[i], M[i]);
        # freeze_max[i] = max(0, -((T-CtoK)*density*dz*C_ICE)/LF). Also track the
        # running melt sum and the deepest cell with melt/excess pore water,
        # avoiding the melt_maximum, freeze_max, and findlast BitVector
        # temporaries. Numerically identical to the previous broadcasts.
        melt = Vector{Float64}(undef, m)
        freeze_max = Vector{Float64}(undef, m)
        melt_sum = 0.0
        X = 1
        @inbounds for i in 1:m
            melt_max_i = T_excess[i] * density[i] * dz[i] * C_ICE / LF
            mi = min(melt_max_i, M[i])
            melt[i] = mi
            melt_sum += mi
            freeze_max[i] = max(0.0, -((temperature[i] - CtoK) * density[i] * dz[i] * C_ICE) / LF)
            if mi > water_tolerance || water_excess[i] > water_tolerance
                X = i
            end
        end
        melt_surface = melt[1]
        melt_total = max(0.0, melt_sum - rain)

        # initialize refreeze, runoff, flux_dn and water_delta vectors
        runoff = zeros(m)
        flux_dn = zeros(m + 1)

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
