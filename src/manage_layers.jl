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

    # Note: arrays are modified in-place. May grow/shrink via layer management.

    d_tolerance = 1e-11

    # store initial mass [kg] and energy [J]
    M_total_initial = 0.0
    E_total_initial = 0.0
    @inbounds for i in eachindex(dz)
        mi = dz[i] * density[i]
        M_total_initial += mi + water[i]
        E_total_initial += mi * temperature[i] * C_ICE + water[i] * (LF + CtoK * C_ICE)
    end

    T_bottom = temperature[end]
    m = length(temperature)

    # Compute per-cell dzmin/dzmax thresholds (replaces cumsum + ones + cumprod allocations)
    column_dzmin2 = Vector{Float64}(undef, m)
    column_dzmax2 = Vector{Float64}(undef, m)
    z_cum = 0.0
    zy_power = 1.0
    @inbounds for i in 1:m
        z_cum += dz[i]
        if z_cum <= mp.column_ztop + d_tolerance
            column_dzmin2[i] = mp.column_dzmin
            column_dzmax2[i] = mp.column_dzmax
        else
            zy_power *= mp.column_zy
            column_dzmin2[i] = mp.column_dzmin * zy_power
            column_dzmax2[i] = mp.column_dzmax * zy_power
        end
    end

    # Compute mass for merge calculations
    M = dz .* density

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

    # Delete merged cells in-place using deleteat!
    to_delete = findall(delete_cell)
    if !isempty(to_delete)
        deleteat!(temperature, to_delete)
        deleteat!(dz, to_delete)
        deleteat!(density, to_delete)
        deleteat!(water, to_delete)
        deleteat!(grain_radius, to_delete)
        deleteat!(grain_dendricity, to_delete)
        deleteat!(grain_sphericity, to_delete)
        deleteat!(albedo, to_delete)
        deleteat!(albedo_diffuse, to_delete)
        deleteat!(column_dzmax2, to_delete)
    end

    # Calculate new length of cells
    m = length(temperature)

    ## Split cells
    # Find the cells that exceed tolerances
    f = findall(dz .> (column_dzmax2 .+ d_tolerance))

    if !isempty(f)
        # Halve dz and water at split positions
        @inbounds for idx in f
            dz[idx] /= 2
            water[idx] /= 2
        end

        # Insert duplicates back-to-front to preserve indices
        @inbounds for k in length(f):-1:1
            idx = f[k]
            insert!(dz, idx, dz[idx])
            insert!(temperature, idx, temperature[idx])
            insert!(density, idx, density[idx])
            insert!(water, idx, water[idx])
            insert!(grain_radius, idx, grain_radius[idx])
            insert!(grain_dendricity, idx, grain_dendricity[idx])
            insert!(grain_sphericity, idx, grain_sphericity[idx])
            insert!(albedo, idx, albedo[idx])
            insert!(albedo_diffuse, idx, albedo_diffuse[idx])
        end
    end

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
        pop!(dz)
        pop!(temperature)
        pop!(water)
        pop!(density)
        pop!(albedo)
        pop!(grain_radius)
        pop!(grain_dendricity)
        pop!(grain_sphericity)
        pop!(albedo_diffuse)
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
        M_total_final = 0.0
        E_total_final = 0.0
        @inbounds for i in eachindex(dz)
            mi = dz[i] * density[i]
            M_total_final += mi + water[i]
            E_total_final += mi * temperature[i] * C_ICE + water[i] * (LF + CtoK * C_ICE)
        end

        M_delta = M_total_initial - M_total_final + mass_added
        E_delta = E_total_initial - E_total_final + E_added

        if (abs(M_delta) > 1e-3) || (abs(E_delta) > 1e-3)
            error("Mass and/or energy are not conserved in manage_layers:\n M_delta: $(M_delta) E_delta: $(E_delta)\n")
        end
    end

    return temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, albedo, albedo_diffuse, mass_added, E_added
end
