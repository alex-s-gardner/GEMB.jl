"""
    manage_layers(temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, albedo, albedo_diffuse, mp::ModelParameters, verbose::Bool)

Adjust the depth and number of vertical layers to maintain proper grid discretization.

Performs three main operations:
1. Merging: Cells thinner than `mp.column_dzmin` are merged with neighbors.
2. Splitting: Cells thicker than `mp.column_dzmax` are split in half.
3. Depth Adjustment: Ensures total column depth stays within limits (`mp.column_zmax`).
   Adds or removes layers at the bottom boundary as necessary.

The Dirichlet temperature boundary condition is enforced at the bottom.

Returns `(temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, albedo, albedo_diffuse, mass_added, E_added)`.
Arrays may grow/shrink (split/merge/depth adjustment).
"""
function manage_layers(temperature::Vector{Float64}, dz::Vector{Float64},
    density::Vector{Float64}, water::Vector{Float64},
    grain_radius::Vector{Float64}, grain_dendricity::Vector{Float64},
    grain_sphericity::Vector{Float64}, albedo::Vector{Float64},
    albedo_diffuse::Vector{Float64},
    mp::ModelParameters, verbose::Bool)

    # Copy inputs to avoid mutation
    temperature = copy(temperature)
    dz = copy(dz)
    density = copy(density)
    water = copy(water)
    grain_radius = copy(grain_radius)
    grain_dendricity = copy(grain_dendricity)
    grain_sphericity = copy(grain_sphericity)
    albedo = copy(albedo)
    albedo_diffuse = copy(albedo_diffuse)

    d_tolerance = 1e-11

    # store initial mass [kg] and energy [J]
    M = dz .* density
    M_total_initial = sum(water) + sum(M)
    E_total_initial = sum(M .* temperature .* C_ICE) +
        sum(water .* (LF + CtoK * C_ICE))

    T_bottom = temperature[end]
    m = length(temperature)

    z_cumulative = cumsum(dz)

    # A logical mask that indicates which cells are in the top layers
    top_layers = z_cumulative .<= (mp.column_ztop + d_tolerance)

    # Define column_dzmin2 array
    column_dzmin2 = mp.column_dzmin .* ones(m)

    # Overwrite the bottom layers with stretched values
    n_bottom = sum(.!top_layers)
    if n_bottom > 0
        column_dzmin2[.!top_layers] = cumprod(mp.column_zy .* ones(n_bottom)) .* mp.column_dzmin
    end

    # Define column_dzmax2 array
    column_dzmax2 = mp.column_dzmax .* ones(m)

    if n_bottom > 0
        column_dzmax2[.!top_layers] = cumprod(mp.column_zy .* ones(n_bottom)) .* mp.column_dzmax
    end

    # Preallocate a logical array for cells to be deleted
    delete_cell = falses(m)

    # Check to see if any cells are too small and need to be merged
    for i in 1:m
        if dz[i] < (column_dzmin2[i] - d_tolerance)
            # dz has not met minimum thickness requirements, delete and merge
            delete_cell[i] = true

            # Determine the target location for the cell contents
            if i == m
                i_target = findlast(.!delete_cell)
            else
                i_target = i + 1
            end

            # Move quantities to target (linearly weighted by mass)
            m_new = M[i] + M[i_target]
            temperature[i_target] = (temperature[i] * M[i] + temperature[i_target] * M[i_target]) / m_new
            albedo[i_target] = (albedo[i] * M[i] + albedo[i_target] * M[i_target]) / m_new
            albedo_diffuse[i_target] = (albedo_diffuse[i] * M[i] + albedo_diffuse[i_target] * M[i_target]) / m_new

            # Use grain properties from lower cell
            grain_radius[i_target] = grain_radius[i]
            grain_dendricity[i_target] = grain_dendricity[i]
            grain_sphericity[i_target] = grain_sphericity[i]

            # Merge with underlying grid cell and delete old cell
            dz[i_target] = dz[i] + dz[i_target]
            density[i_target] = m_new / dz[i_target]
            water[i_target] = water[i] + water[i_target]
            M[i_target] = m_new
        end
    end

    # Delete combined cells
    keep = .!delete_cell
    water = water[keep]
    dz = dz[keep]
    density = density[keep]
    temperature = temperature[keep]
    albedo = albedo[keep]
    grain_radius = grain_radius[keep]
    grain_dendricity = grain_dendricity[keep]
    grain_sphericity = grain_sphericity[keep]
    albedo_diffuse = albedo_diffuse[keep]
    column_dzmax2 = column_dzmax2[keep]

    # Calculate new length of cells
    m = length(temperature)

    ## Split cells
    # Find the cells that exceed tolerances
    f = findall(dz .> (column_dzmax2 .+ d_tolerance))

    # Conserve quantities among the cells that will be split
    dz[f] = dz[f] ./ 2
    water[f] = water[f] ./ 2

    # Sort the indices of all the cells including the ones that will be duplicated
    fs = sort(vcat(collect(1:m), f))

    # Recreate the variables with split cells
    dz = dz[fs]
    water = water[fs]
    temperature = temperature[fs]
    density = density[fs]
    albedo = albedo[fs]
    albedo_diffuse = albedo_diffuse[fs]
    grain_radius = grain_radius[fs]
    grain_dendricity = grain_dendricity[fs]
    grain_sphericity = grain_sphericity[fs]

    ## CORRECT FOR TOTAL MODEL DEPTH

    # Calculate total model depth
    z_total = sum(dz)

    local mass_added::Float64
    local E_added::Float64

    if z_total < (mp.column_zmax - d_tolerance)
        # Mass and energy to be added
        mass_added = (dz[end] * density[end]) + water[end]
        E_added = temperature[end] * (dz[end] * density[end]) * C_ICE + water[end] * (LF + CtoK * C_ICE)

        # Add a grid cell of the same size and temperature to the bottom
        # Optimized: use push! instead of vcat for single elements
        push!(dz, dz[end])
        push!(temperature, temperature[end])
        push!(water, water[end])
        push!(density, density[end])
        push!(albedo, albedo[end])
        push!(albedo_diffuse, albedo_diffuse[end])
        push!(grain_radius, grain_radius[end])
        push!(grain_dendricity, grain_dendricity[end])
        push!(grain_sphericity, grain_sphericity[end])

    elseif z_total > mp.column_zmax + d_tolerance
        # Mass and energy loss
        mass_added = -((dz[end] * density[end]) + water[end])
        E_added = -(temperature[end] * (dz[end] * density[end]) * C_ICE) - water[end] * (LF + CtoK * C_ICE)

        # Remove a grid cell from the bottom
        dz = dz[1:end-1]
        temperature = temperature[1:end-1]
        water = water[1:end-1]
        density = density[1:end-1]
        albedo = albedo[1:end-1]
        grain_radius = grain_radius[1:end-1]
        grain_dendricity = grain_dendricity[1:end-1]
        grain_sphericity = grain_sphericity[1:end-1]
        albedo_diffuse = albedo_diffuse[1:end-1]
    else
        # No mass or energy is added or removed
        mass_added = 0.0
        E_added = 0.0
    end

    # Enforce Dirichlet boundary condition at bottom
    E_added = E_added + ((T_bottom - temperature[end]) * (dz[end] * density[end]) * C_ICE)
    temperature[end] = T_bottom

    ## CHECK FOR MASS AND ENERGY CONSERVATION
    if verbose
        M = dz .* density
        M_total_final = sum(water) + sum(M)
        E_total_final = sum(M .* temperature .* C_ICE) +
            sum(water .* (LF + CtoK * C_ICE))

        M_delta = M_total_initial - M_total_final + mass_added
        E_delta = E_total_initial - E_total_final + E_added

        if (abs(M_delta) > 1e-3) || (abs(E_delta) > 1e-3)
            error("Mass and/or energy are not conserved in manage_layers:\n M_delta: $(M_delta) E_delta: $(E_delta)\n")
        end
    end

    return temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, albedo, albedo_diffuse, mass_added, E_added
end
